from formatizer import f
from . import utils
from .executor import Executor

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

        input_filenames = [
            os.path.join(self.fastq_dir, filename)
            for pair in utils.find_fastqs_by_organism(
                self.analysis.sample,
                self.fastq_dir).values()
            for filename, _ in pair]
        executor = Executor(self.analysis)
        executor(
            f('cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAG '
              '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG '
              '-m 20 -o "{{output_filename[0]}}" -p '
              '"{{output_filename[1]}}" {{input_filename}} '
              '> "{self.sample_base_out}.cutadapt.txt"'),
            output_format=f(
                "{self.analysis.sample}{{organism_str}}.clipped.R%d.fastq"),
            input_filenames=input_filenames,
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            output_path=self.fastq_dir,
            output_function=lambda filename: [filename % (index + 1)
                                              for index in range(2)])

        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self):
        self.analysis.logger.info("Running fastqc")
        self.chdir()

        qcdir = os.path.join(self.fastq_dir, "FASTQC")
        if not os.path.exists(qcdir):
            os.mkdir(qcdir)

        executor = Executor(self.analysis)
        executor(f('{self.analysis.config.fastqc} '
                   '"{{input_filename}}" --outdir {qcdir}'),
                 override_last_files=False)

        self.analysis.logger.info("Finished fastqc")

    def trim(self):
        trim_end = False
        if self.analysis.parameters["run_xenome"]:
            trim_end = True

        if trim_end:
            self.analysis.logger.info("Trimming first 5 bp and last 10 bp")
            trim_end_cmd = "-e 10 "
        else:
            self.analysis.logger.info("Trimming first 5 bp")
            trim_end_cmd = ""

        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f('{config.seqtk} trimfq -b 5 '
                   '{trim_end_cmd}'
                   '"{{input_filename}}" '
                   '> "{{output_filename}}"'),
                 output_format=os.path.join(
                     self.fastq_dir,
                     "%s{organism_str}.trimmed.R{read_index}.fastq"
                     % self.analysis.sample),
                 output_path=self.fastq_dir,
                 error_string="Trimming with seqtk exited with status "
                              "{status}",
                 exception_string="trimming error",
                 unlink_inputs=True)

        self.analysis.logger.info("Finished trimming")

    def align(self):
        self.analysis.logger.info("Running alignment")
        self.chdir()

        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.novoalign} -oSAM "@RG\tID:{self.analysis.basename}\t'
            'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
            '-d {{genome_index}} '
            '-i PE {config.mean_len_library},{config.sd_len_library} '
            '-t 90 -f {{input_filename}}> {{output_filename}}'),
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            output_format=f("{self.analysis.basename}{{organism_str}}.sam"),
            split_by_organism=True,
            only_human=True
        )

        self.analysis.logger.info("Alignment SAM -> BAM")
        executor(f(
            '{config.picard} SamFormatConverter '
            'I={{input_filename}} '
            'O={{output_filename}}'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.bam"),
            error_string="Picard SamFormatConverter exited with status "
                         "{status}",
            exception_string="picard SamFormatConverter error",
            unlink_inputs=True
        )

        executor(f(
            '{config.picard} AddOrReplaceReadGroups '
            'I={{input_filename}} '
            'O={{output_filename}} RGID={self.analysis.basename} '
            'RGLB=lib1 RGPL=ILLUMINA RGPU={config.kit} '
            'RGSM={self.analysis.basename}'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.rg.bam"),
            error_string="Picard AddOrReplaceReadGroups exited with status {status}",
            exception_string="picard AddOrReplaceReadGroups error",
            unlink_inputs=True
        )

        executor(lambda **kwargs: os.rename(kwargs["input_filename"], kwargs["output_filename"]),
                 output_format=f("{self.analysis.basename}{{organism_str}}.bam"))

        self.analysis.logger.info("Finished alignment")

    def sort_bam(self):
        self.analysis.logger.info("Sorting BAM(s)")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.picard} SortSam '
            'I={{input_filename}} '
            'O={{output_filename}} SO=coordinate'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.srt.bam"),
            error_string="Picard SortSam exited with status {status}",
            exception_string="picard SortSam error",
            unlink_inputs=True
        )

        executor(f(
            '{config.picard} ReorderSam '
            'I={{input_filename}} '
            'O={{output_filename}} R={{genome_ref}} '
            'CREATE_INDEX=true'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.srt.reorder.bam"),
            error_string="Picard ReorderSam exited with status {status}",
            exception_string="picard ReorderSam error",
            unlink_inputs=True
        )

        self.analysis.logger.info("Finished sorting")

    def mark_duplicates(self):
        self.analysis.logger.info("Marking duplicates")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.picard} MarkDuplicates '
            'I={{input_filename}} '
            'O={{output_filename}} '
            'M={self.output_basename}{{organism_str}}.marked_dup_metrics.txt '
            'CREATE_INDEX=true'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.srt.marked.dup.bam"),
            error_string="Picard MarkDuplicates exited with status {status}",
            exception_string="picard MarkDuplicates error",
            unlink_inputs=True
        )

        self.analysis.logger.info("Finished marking duplicates")

    def indel_realign(self):
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.gatk} -T RealignerTargetCreator -R {{genome_ref}} '
            '-I {{input_filename}} -nt {self.gatk_threads} '
            '-known {config.indel_1} -known {config.indel_2} '
            '-L {config.target_list} '
            '-ip 50 -o {{output_filename}}'),
            output_format=f("{self.output_basename}{{organism_str}}"
                            ".realignment.intervals"),
            error_string="Gatk RalignerTargetCreator exited with status "
                         "{status}",
            exception_string="gatk RealignerTargetCreator error",
            override_last_files=False
        )

        executor(f(
            '{config.gatk} -T IndelRealigner -R {{genome_ref}} '
            '-I {{input_filename}} '
            '-known {config.indel_1} -known {config.indel_2} '
            '-targetIntervals {self.output_basename}{{organism_str}}'
            '.realignment.intervals -o {{output_filename}}'),
            output_format=f("{self.analysis.basename}{{organism_str}}"
                            ".srt.realigned.bam"),
            error_string="Gatk IndelRealigner exited with status {status}",
            exception_string="gatk IndelRealigner error",
            unlink_inputs=True
        )

        self.analysis.logger.info("Finished indel realignment")

    def _filter_non_hg(self, filename):
        organism = utils.get_params_from_filename(filename)[8]
        if organism is None or organism.lower().startswith("hg"):
            return filename
        else:
            return None

    def recalibration(self):
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.gatk} -T BaseRecalibrator -R {{genome_ref}} '
            '-I {{input_filename}} -nct {self.gatk_threads} '
            '-knownSites {{dbsnp}} '
            '-o {self.output_basename}{{organism_str}}.recalibration.table'),
            input_function=lambda filename: self._filter_non_hg(filename),
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            override_last_files=False
        )

        executor(f(
            '{config.gatk} -T PrintReads -R {{genome_ref}} '
            '-I {{input_filename}} -nct {self.gatk_threads} '
            '-BQSR {self.output_basename}{{organism_str}}.recalibration.table '
            '-o {{output_filename}}'),
            output_format=f("{self.analysis.basename}{{organism_str}}"
                            ".srt.realigned.recal.bam"),
            error_string="Gatk PrintReads exited with status {status}",
            exception_string="gatk PrintReads error",
            unlink_inputs=True
        )

        if not self.analysis.parameters["run_post_recalibration"]:
            self.analysis.logger.info("Finished recalibration")
            return

        executor(f(
            '{config.gatk} -T BaseRecalibrator -R {{genome_ref}} '
            '-I {{input_filename}} -knownSites {{dbsnp}} '
            '-L {config.target_list} -ip 50 '
            '-nct {self.gatk_threads} '
            '-o {self.output_basename}{{organism_str}}'
            '.post_realignment.table'),
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            override_last_files=False
        )

        executor(f(
            '{config.gatk} -T AnalyzeCovariates -R {{genome_ref}} '
            '-before {self.output_basename}{{organism_str}}.recalibration.table '
            '-after {self.output_basename}{{organism_str}}.post_realignment.table '
            '-plots {self.output_basename}{{organism_str}}.recalibration_plots.pdf'),
            error_string="Gatk AnalyzeCovariates exited with status {status}",
            exception_string="gatk AnalyzeCovariates error",
            override_last_files=False
        )

        executor(f(
            '{config.picard} MarkDuplicates '
            'I={{input_filename}} O={{output_filename}} '
            'REMOVE_DUPLICATES=true '
            'M={self.output_basename}{{organism_str}}.no_dup_metrics.txt '
            'CREATE_INDEX=true'
            '{self.max_records_str}'),
            output_format=f("{self.analysis.basename}{{organism_str}}.srt.realigned.recal.no_dup.bam"),
            error_string="Picard MarkDuplicates exited with status {status}",
            exception_string="picard MarkDuplicates error",
            unlink_inputs=True
        )

        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self):
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(f(
            '{config.picard} CollectHsMetrics '
            'I={{input_filename}} BI={config.bait_list} '
            'TI={config.target_list} R={{genome_ref}} '
            'O={self.output_basename}{{organism_str}}.hs_metrics.txt'
            '{self.max_records_str}'),
            error_string="Picard CollectHsMetrics exited with status {status}",
            exception_string="picard CollectHsMetrics error",
            override_last_files=False
        )

        executor(f(
            '{config.picard} CollectTargetedPcrMetrics '
            'I={{input_filename}} AMPLICON_INTERVALS='
            '{config.bait_list} TARGET_INTERVALS={config.target_list} '
            'R={{genome_ref}} '
            'O={self.output_basename}{{organism_str}}.targeted_metrics.txt'
            '{self.max_records_str} '
            'PER_BASE_COVERAGE={self.output_basename}{{organism_str}}.coverage.txt'),
            error_string="Picard CollectTargetedPcrMetrics exited with status {status}",
            exception_string="picard CollectTargetedPcrMetrics error",
            override_last_files=False
        )

        executor(f(
            '{config.picard} CollectGcBiasMetrics '
            'R={{genome_ref}} I={{input_filename}} '
            'O={self.output_basename}{{organism_str}}.gcbias.metrics.txt '
            'CHART={self.output_basename}{{organism_str}}.gcbias_metrics.pdf '
            'S={self.output_basename}{{organism_str}}.gcbias_summ_metrics.txt'
            '{self.max_records_str}'),
            error_string="Picard CollectGcBiasMetrics exited with status {status}",
            exception_string="picard CollectGcBiasMetrics error",
            override_last_files=False
        )

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

    def run(self):
        if self.analysis.parameters["run_cutadapt"]:
            self.cutadapt()
        self.fastqc()
        self.trim()
        self.align()
        self.sort_bam()
        if self.analysis.parameters["mark_duplicates"]:
            self.mark_duplicates()
        self.indel_realign()
        self.recalibration()
        self.metrics_collection()
        if self.analysis.parameters["compress_fastq"]:
            self.compress_fastq()
