from formatizer import f
from .exceptions import PipelineError
from . import utils

import os
import gzip
import shutil


class Mapping:

    def __init__(self, analysis, fastq_dir):
        self.analysis = analysis
        self.fastq_dir = fastq_dir

        self.sample_base = os.path.join(self.fastq_dir, self.analysis.sample)
        self.sample_base_out = os.path.join(
            self.fastq_dir,
            "REPORTS",
            self.analysis.sample)
        self.output_basename = os.path.join("REPORTS", self.analysis.basename)

        if not os.path.exists(os.path.join(self.analysis.bam_dir, "REPORTS")):
            os.makedirs(os.path.join(self.analysis.bam_dir, "REPORTS"))

        if not os.path.exists(os.path.join(self.fastq_dir, "REPORTS")):
            os.makedirs(os.path.join(self.fastq_dir, "REPORTS"))

        self.gatk_threads = self.analysis.parameters["gatk_threads"]
        self.max_records_str = utils.get_picard_max_records_string(
            self.analysis.parameters["picard_max_records"])

    def chdir(self):
        os.chdir(self.analysis.bam_dir)

    def cutadapt(self):
        self.analysis.logger.info("Cutting adapters")
        self.chdir()

        input_files = utils.find_fastqs_by_organism(
            self.analysis.sample,
            self.fastq_dir)
        output_files = {}
        for organism, filenames_and_indexes in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_" + organism
            out_filenames = sorted(
                [
                    "%s%s.clipped.R%d.fastq" %
                    (self.analysis.sample,
                     organism_str,
                     read_index) for _,
                    read_index in filenames_and_indexes])
            in_filenames = sorted([
                "\"%s\"" % os.path.join(self.fastq_dir, filename)
                for filename, _ in filenames_and_indexes])
            in_filenames_str = " ".join(in_filenames)

            retval = utils.run_and_log(f(
                'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAG '
                '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG '
                '-m 20 -o "{out_filenames[0]}" -p '
                '"{out_filenames[1]}" {in_filenames_str} '
                '> "{self.sample_base_out}.cutadapt.txt"'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Cutadapt exited with status %d" %
                    retval)
                raise PipelineError("cutadapt error")

            output_files[organism] = out_filenames

        self.analysis.last_operation_filenames = output_files
        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self):
        self.analysis.logger.info("Running fastqc")
        self.chdir()
        config = self.analysis.config

        qcdir = os.path.join(self.fastq_dir, "FASTQC")
        if not os.path.exists(qcdir):
            os.mkdir(qcdir)

        for input_filename in utils.get_sample_filenames(
                self.analysis.last_operation_filenames):
            retval = utils.run_and_log(f(
                '{config.fastqc} "{input_filename}"'
                ' --outdir {qcdir}'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "fastqc exited with status %d" %
                    retval)
                raise PipelineError("fastqc error")

        self.analysis.logger.info("Finished fastqc")

    def trim(self):
        trim_end = False
        if self.analysis.parameters["run_xenome"]:
            trim_end = True

        if trim_end:
            self.analysis.logger.info("Trimming first 5 bp")
            trim_end_cmd = ""
        else:
            self.analysis.logger.info("Trimming first 5 bp and last 10 bp")
            trim_end_cmd = "-e 10 "

        self.chdir()
        config = self.analysis.config

        out_filenames = []
        for filename in utils.get_sample_filenames(
                self.analysis.last_operation_filenames):
            organism, read_index, _ = utils.get_params_from_filename(
                filename, self.analysis.sample)
            if organism is None or organism == "":
                organism_str = ""
                organism = "hg19"
            else:
                organism_str = "_" + organism

            out_filename = os.path.join(
                self.fastq_dir, "%s%s.trimmed.R%d.fastq" %
                (self.analysis.sample, organism_str, read_index))

            retval = utils.run_and_log(f(
                '{config.seqtk} trimfq -b 5 '
                '{trim_end_cmd}'
                '"{filename}" '
                '> "{out_filename}"'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Trimming with seqtk exited with status %d" %
                    retval)
                raise PipelineError("trimming error")

            out_filenames.append(out_filename)
            os.unlink(filename)

        self.analysis.last_operation_filenames = out_filenames

        self.analysis.logger.info("Finished trimming")

    def align(self):
        self.analysis.logger.info("Running alignment")
        self.chdir()
        config = self.analysis.config
        samples_by_organism = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)

        output_files = {}
        for organism, filenames in samples_by_organism.items():
            input_files = " ".join(
                ["\"%s\"" % filename
                 for filename in filenames])

            if len(samples_by_organism) == 1:
                output_file = self.analysis.basename + ".sam"
            else:
                output_file = "%s_%s.sam" % (self.analysis.basename, organism)

            genome_index = \
                utils.get_genome_ref_index_by_organism(config, organism)[1]

            retval = utils.run_and_log(f(
                '{config.novoalign} -oSAM "@RG\tID:{self.analysis.basename}\t'
                'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
                '-d {genome_index} '
                '-i PE {config.mean_len_library},{config.sd_len_library} '
                '-t 90 -f {input_files}> {output_file}'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Novoalign exited with status %d" %
                    retval)
                raise PipelineError("novoalign error")

            output_files[organism] = os.path.join(os.getcwd(), output_file)

        self.analysis.last_operation_filenames = output_files

        self.analysis.logger.info("Alignment SAM -> BAM")
        output_files = {}
        for organism, filename in self.analysis.last_operation_filenames.items():
            if len(self.analysis.last_operation_filenames) == 1:
                output_file = self.analysis.basename + ".bam"
            else:
                output_file = "%s_%s.bam" % (self.analysis.basename, organism)

            retval = utils.run_and_log(f(
                '{config.picard} SamFormatConverter '
                'I={filename} '
                'O={output_file}'
                '{self.max_records_str}'
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Picard SamFormatConverter exited with status %d" %
                    retval)
                raise PipelineError("picard SamFormatConverter error")

            os.unlink(filename)
            output_files[organism] = os.path.join(os.getcwd(), output_file)

        self.analysis.last_operation_filenames = output_files

        output_files = {}
        for organism, filename in self.analysis.last_operation_filenames.items():
            if len(self.analysis.last_operation_filenames) == 1:
                output_file = self.analysis.basename + ".rg.bam"
            else:
                output_file = "%s_%s.rg.bam" % (
                    self.analysis.basename, organism)

            retval = utils.run_and_log(f(
                '{config.picard} AddOrReplaceReadGroups '
                'I={filename} '
                'O={output_file} RGID={self.analysis.basename} '
                'RGLB=lib1 RGPL=ILLUMINA RGPU={config.kit} '
                'RGSM={self.analysis.basename}'
                '{self.max_records_str}',
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Picard AddOrReplaceReadGroup exited with status %d" %
                    retval)
                raise PipelineError("picard AddOrReplaceReadGroup error")

            output_files[organism] = os.path.join(os.getcwd(), output_file)

        self.analysis.last_operation_filenames = output_files

        output_files = {}
        for organism, filename in self.analysis.last_operation_filenames.items():
            if len(self.analysis.last_operation_filenames) == 1:
                output_file = self.analysis.basename + ".bam"
            else:
                output_file = "%s_%s.bam" % (self.analysis.basename, organism)
            os.rename(filename, output_file)

            if not os.path.exists(output_file):
                self.analysis.logger.error(
                    "BAM file still does not exist, when it should have "
                    "been created")
                raise PipelineError("picard failed")

            output_files[organism] = os.path.join(os.getcwd(), output_file)

        self.analysis.last_operation_filenames = output_files

        self.analysis.bamfiles = self.analysis.last_operation_filenames
        self.analysis.logger.info("Finished alignment")

    def sort_bam(self):
        self.analysis.logger.info("Sorting BAM(s)")
        self.chdir()
        config = self.analysis.config

        input_files = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)

        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                output_file = self.analysis.basename + ".srt.bam"
            else:
                output_file = "%s_%s.srt.bam" % (
                    self.analysis.basename, organism)

            retval = utils.run_and_log(f(
                '{config.picard} SortSam '
                'I={filename} '
                'O={output_file} SO=coordinate'
                '{self.max_records_str}'),
                self.analysis.logger
            )
            output_files[organism] = os.path.join(os.getcwd(), output_file)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard SortSam exited with status %d" %
                    retval)
                raise PipelineError("picard SortSam error")

            os.unlink(filename)

        self.analysis.last_operation_filenames = output_files

        input_files = self.analysis.last_operation_filenames
        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                output_file = self.analysis.basename + ".srt.reorder.bam"
            else:
                output_file = "%s_%s.srt.reorder.bam" % (
                    self.analysis.basename, organism)

            output_file = self.analysis.basename + ".srt.reorder.bam"
            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            retval = utils.run_and_log(f(
                '{config.picard} ReorderSam '
                'I={filename} '
                'O={output_file} R={genome_ref} '
                'CREATE_INDEX=true'
                '{self.max_records_str}',
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Picard SortSam exited with status %d" %
                    retval)
                raise PipelineError("picard SortSam")

            if not os.path.exists(output_file):
                self.analysis.logger.error(
                    "BAM file still does not exist, when it should have "
                    "been created")
                raise PipelineError("picard failed")

            output_files[organism] = os.path.join(os.getcwd(), output_file)
            os.unlink(filename)

        self.analysis.last_operation_filenames = output_files
        self.analysis.bamfiles = output_files
        self.analysis.logger.info("Finished sorting")

    def mark_duplicates(self):
        self.analysis.logger.info("Marking duplicates")
        self.chdir()
        config = self.analysis.config

        input_files = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)

        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            output_file = self.output_basename + \
                organism_str + ".srt.marked.dup.bam"
            retval = utils.run_and_log(f(
                '{config.picard} MarkDuplicates '
                'I={filename} '
                'O={output_file} '
                'M={self.output_basename}{organism_str}.marked_dup_metrics.txt '
                'CREATE_INDEX=true '
                'SO=coordinate'
                '{self.max_records_str}'),
                self.analysis.logger
            )
            output_files[organism] = os.path.join(os.getcwd(), output_file)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard MarkDuplicates exited with status %d" %
                    retval)
                raise PipelineError("picard MarkDuplicates error")

            os.unlink(filename)
            index_file = filename[:-4] + ".bai"
            if os.path.exists(index_file):
                os.unlink(index_file)

        self.analysis.last_operation_filenames = output_files
        self.analysis.bamfiles = output_files
        self.analysis.logger.info("Finished marking duplicates")

    def indel_realign(self):
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config

        input_files = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)

        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            output_file = self.output_basename + \
                organism_str + ".realignment.intervals"
            retval = utils.run_and_log(f(
                '{config.gatk} -T RealignerTargetCreator -R {genome_ref} '
                '-I {filename} -nt {self.gatk_threads} '
                '-known {config.indel_1} -known {config.indel_2} '
                '-L {config.target_list} '
                '-ip 50 -o {output_file}'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk RealignerTargetCreator exited with status %d" %
                    retval)
                raise PipelineError("gatk RealignerTargetCreator error")

        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            output_file = self.analysis.basename + \
                organism_str + ".srt.realigned.bam"
            retval = utils.run_and_log(
                f('{config.gatk} -T IndelRealigner -R {genome_ref} '
                  '-I {filename} '
                  '-known {config.indel_1} -known {config.indel_2} '
                  '-targetIntervals {self.output_basename}{organism_str}.realignment.intervals '
                  '-o {output_file}'), self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk IndelRealigner exited with status %d" %
                    retval)
                raise PipelineError("gatk IndelRealigner failed")

            if not os.path.exists(output_file):
                self.analysis.logger.error(
                    "BAM file still does not exist, when it should have "
                    "been created")
                raise PipelineError("gatk failed")

            output_files[organism] = os.path.join(
                os.getcwd(),
                output_file)

            os.unlink(filename)
            os.unlink(filename[:-4] + ".bai")

        self.analysis.last_operation_filenames = output_files
        self.analysis.bamfiles = self.analysis.last_operation_filenames
        self.analysis.logger.info("Finished indel realignment")

    def recalibration(self):
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config

        input_files = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)
        ignored_files = {}
        for key in input_files:
            if not key.startswith("hg"):
                ignored_files[key] = input_files[key]
                del input_files[key]

        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            dbsnp = utils.get_dbsnp_by_organism(config, organism)
            retval = utils.run_and_log(f(
                '{config.gatk} -T BaseRecalibrator -R {genome_ref} -I '
                '{filename} -nct {self.gatk_threads} '
                '-knownSites {dbsnp} '
                '-o {self.output_basename}{organism_str}.recalibration.table',
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk BaseRecalibrator exited with status %d" %
                    retval)
                raise PipelineError("gatk BaseRecalibrator error")

        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            output_file = self.analysis.basename + \
                organism_str + ".srt.realigned.recal.bam"
            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            retval = utils.run_and_log(f(
                '{config.gatk} -T PrintReads -R {genome_ref} -I '
                '{filename} -nct {self.gatk_threads} '
                '-BQSR {self.output_basename}{organism_str}.recalibration.table -o '
                '{output_file}',
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk PrintReads exited with status %d" %
                    retval)
                raise PipelineError("gatk PrintReads error")

            if not os.path.exists(output_file):
                self.analysis.logger.error(
                    "BAM file still does not exist, when it should have "
                    "been created")
                raise PipelineError("gatk failed")

            output_files[organism] = os.path.join(os.getcwd(), output_file)

            os.unlink(filename)
            os.unlink(filename[:-4] + ".bai")

        self.analysis.last_operation_filenames = output_files
        self.analysis.last_operation_filenames.update(ignored_files)

        if not self.analysis.parameters["run_post_recalibration"]:
            self.analysis.logger.info("Finished recalibration")
            return

        input_files = self.analysis.last_operation_filenames
        for key in input_files.keys():
            if not key.startswith("hg"):
                del input_files[key]

        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            dbsnp = utils.get_dbsnp_by_organism(config, organism)
            retval = utils.run_and_log(
                f('{config.gatk} -T BaseRecalibrator -R {genome_ref} '
                  '-I {filename} '
                  '-knownSites {dbsnp} '
                  '-L {config.target_list} '
                  '-ip 50 '
                  '-nct {self.gatk_threads} '
                  '-o {self.output_basename}{organism_str}'
                  '.post_realignment.table'),
                self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk BaseRecalibrator exited with status %d" %
                    retval)
                raise PipelineError("gatk BaseRecalibrator failed")

        input_files = self.analysis.last_operation_filenames
        for key in input_files.keys():
            if not key.startswith("hg"):
                del input_files[key]

        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            retval = utils.run_and_log(
                f('{config.gatk} -T AnalyzeCovariates -R {genome_ref} '
                  '-before {self.output_basename}{organism_str}.recalibration.table '
                  '-after {self.output_basename}{organism_str}.post_realignment.table '
                  '-plots {self.output_basename}{organism_str}.recalibration_plots.pdf '
                  ), self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Gatk AnalyzeCovariates exited with status %d" %
                    retval)
                raise PipelineError("gatk AnalyzeCovariates failed")

        output_files = {}
        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            output_file = self.analysis.basename + \
                organism_str + ".srt.realigned.recal.no_dup.bam"
            retval = utils.run_and_log(
                f('{config.picard} MarkDuplicates '
                  'I={filename} '
                  'O={output_file} '
                  'REMOVE_DUPLICATES=true '
                  'M={self.output_basename}{organism_str}.no_dup_metrics.txt '
                  'CREATE_INDEX=true'
                  '{self.max_records_str}'), self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard MarkDuplicates exited with status %d" %
                    retval)
                raise PipelineError("picard MarkDuplicates failed")

            if not os.path.exists(output_file):
                self.analysis.logger.error(
                    "BAM file still does not exist, when it should have "
                    "been created")
                raise PipelineError("gatk failed")

            output_files[organism] = os.path.join(
                os.getcwd(),
                output_file)

            os.unlink(filename)
            os.unlink(filename[:-4] + ".bai")

        self.analysis.last_operation_filenames = output_files
        self.analysis.bamfiles = self.analysis.last_operation_filenames
        self.analysis.bamfiles = self.analysis.last_operation_filenames
        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self):
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        input_files = utils.get_samples_by_organism(
            self.analysis.last_operation_filenames)

        for organism, filename in input_files.items():
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(
                config,
                organism)[0]
            retval = utils.run_and_log(f(
                '{config.picard} CollectHsMetrics '
                'I={filename} '
                'BI={config.bait_list} '
                'TI={config.target_list} R={genome_ref} '
                'O={self.output_basename}{organism_str}.hs_metrics.txt'
                '{self.max_records_str}'),
                self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard CollectHsMetrics exited with status %d" %
                    retval)
                raise PipelineError("picard CollectHsMetrics error")

            retval = utils.run_and_log(
                f('{config.picard} CollectTargetedPcrMetrics '
                  'I={filename} AMPLICON_INTERVALS='
                  '{config.bait_list} TARGET_INTERVALS={config.target_list} '
                  'R={genome_ref} '
                  'O={self.output_basename}{organism_str}.targeted_metrics.txt'
                  '{self.max_records_str} '
                  'PER_BASE_COVERAGE={self.output_basename}{organism_str}.coverage.txt',),
                self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard CollectTargetedPcrMetrics exited with status %d" %
                    retval)
                raise PipelineError("picard CollectTargetedPcrMetrics error")

            retval = utils.run_and_log(
                f('{config.picard} CollectGcBiasMetrics '
                  'R={genome_ref} '
                  'I={filename} '
                  'O={self.output_basename}{organism_str}.gcbias.metrics.txt '
                  'CHART={self.output_basename}{organism_str}.gcbias_metrics.pdf '
                  'S={self.output_basename}{organism_str}.gcbias_summ_metrics.txt'
                  '{self.max_records_str}'),
                self.analysis.logger)

            if retval != 0:
                self.analysis.logger.error(
                    "Picard CollectGcBiasMetrics exited with status %d" %
                    retval)
                raise PipelineError("picard CollectGcBiasMetrics error")

        self.analysis.logger.info("Finished metrics collection")

    def compress_fastq(self):
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = utils.find_fastqs_by_organism(
            self.analysis.sample,
            self.fastq_dir)
        for filenames in fastq_files.values():
            for filename, _ in filenames:
                filename = os.path.join(self.fastq_dir, filename)
                compressed_filename = filename + ".gz"
                with open(filename, "rb") as in_fd, \
                        gzip.open(compressed_filename, "wb") as out_fd:
                    shutil.copyfileobj(in_fd, out_fd)
                os.unlink(filename)
        self.analysis.logger.info("Finished compressing fastq files")

    def variant_calling(self):
        self.analysis.logger.info("Running variant calling")
        self.chdir()
        retval = utils.run_and_log(
            f('{config.gatk} -T HaplotypeCaller -R {config.genome_ref} -I '
              '{self.analysis.last_operation_filenames} -ERC GVCF '
              '-nct {self.gatk_threads} '
              '--genotyping_mode DISCOVERY -L {config.target_list} '
              '-stand_emit_conf 10 -stand_call_conf 30 --output_mode '
              'EMIT_VARIANTS_ONLY '
              '-o {self.analysis.basename}.variants.g.vcf',),
            self.analysis.logger)

        if retval != 0:
            self.analysis.logger.error(
                "Gatk HaplotypeCaller exited with status %d" %
                retval)
            raise PipelineError("gatk HaplotypeCaller error")

        retval = utils.run_and_log(f(
            '{config.gatk} -T VariantFiltration -R {config.genome_ref} '
            '-nct {self.gatk_threads} '
            '--variant {self.analysis.basename}.variants.vcf '
            '-o {self.analysis.basename}.variants.filtered.vcf '
            '--clusterWindowSize 10 '
            '--clusterSize 3 --filterExpression "DP < 10" '
            '--filterName "LowCoverage" '
            '--filterExpression "QUAL < 30.0" --filterName '
            '"VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" '
            '--filterName "LowQual" --filterExpression "QD < 1.5" '
            '--filterName "LowQD" ',
            # --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"
            # --filterName "HARD_TO_VALIDATE"
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error(
                "Gatk VariantFiltration exited with status %d" %
                retval)
            raise PipelineError("gatk VariantFiltration error")

        self.analysis.logger.info("Finished variant calling")

    def run(self):
        if self.analysis.parameters["run_cutadapt"]:
            self.cutadapt()
        self.fastqc()
        self.trim()
        self.align()
        self.sort_bam()
        self.indel_realign()
        self.recalibration()
        self.metrics_collection()
        if self.analysis.parameters["compress_fastq"]:
            self.compress_fastq()