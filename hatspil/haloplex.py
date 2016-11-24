import os
from formatizer import f
from .exceptions import PipelineError
from . import utils


class Haloplex:

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

    def chdir(self):
        os.chdir(self.analysis.bam_dir)

    def cutadapt(self):
        self.analysis.logger.info("Cutting adapters")
        self.chdir()
        output_filenames = [
            "%s.clipped.R%d.fastq" % (self.sample_base, index + 1)
            for index in range(2)]
        retval = utils.run_and_log(f(
            'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAG '
            '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG '
            '-m 20 -o "{output_filenames[0]}" -p '
            '"{output_filenames[1]}" "{self.sample_base}_R1.fastq" '
            '"{self.sample_base}_R2.fastq" '
            '> "{self.sample_base_out}.cutadapt.txt"'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error(
                "Cutadapt exited with status %d" %
                retval)
            raise PipelineError("cutadapt error")

        self.analysis.last_operation_filenames = output_filenames
        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self):
        self.analysis.logger.info("Running fastqc")
        self.chdir()
        config = self.analysis.config

        qcdir = os.path.join(self.fastq_dir, "FASTQC")
        if not os.path.exists(qcdir):
            os.mkdir(qcdir)

        for input_filename in self.analysis.last_operation_filenames:
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
        self.analysis.logger.info("Trimming first 5 bp")
        self.chdir()
        config = self.analysis.config

        out_filenames = []
        for filename in self.analysis.last_operation_filenames:
            read_index = utils.get_read_index(filename)
            out_filename = "%s.trimmed.R%d.fastq" % (filename, read_index)
            retval = utils.run_and_log(f(
                '{config.seqtk} trimfq -b 5 "{filename}" '
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
        input_files = " ".join(
            ["\"%s\"" % filename
             for filename in self.analysis.last_operation_filenames])

        output_file = self.analysis.basename + ".sam"
        retval = utils.run_and_log(f(
            '{config.novoalign} -oSAM "@RG\tID:{self.analysis.basename}\tSM:'
            '{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
            '-d {config.hg19_index} '
            '-i PE {config.mean_len_library},{config.sd_len_library} -t 90 '
            '-f {input_files}> {output_file}'),
            self.analysis.logger
        )
        self.analysis.last_operation_filenames = output_file

        if retval != 0:
            self.analysis.logger.error(
                "Novoalign exited with status %d" %
                retval)
            raise PipelineError("novoalign error")

        self.analysis.logger.info("Alignment SAM -> BAM")
        output_file = self.analysis.basename + ".bam"
        retval = utils.run_and_log(f(
            '{config.picard} SamFormatConverter '
            'I={self.analysis.last_operation_filenames} '
            'O={output_file}',  # MAX_RECORDS_IN_RAM=1000000
            ),
            self.analysis.logger
        )
        os.unlink(self.analysis.last_operation_filenames)
        self.analysis.last_operation_filenames = output_file

        if retval != 0:
            self.analysis.logger.error(
                "Picard SamFormatConverter exited with status %d" %
                retval)
            raise PipelineError("picard SamFormatConverter error")

        output_file = self.analysis.basename + ".rg.bam"
        retval = utils.run_and_log(f(
            '{config.picard} AddOrReplaceReadGroups '
            'I={self.analysis.last_operation_filenames} '
            'O={output_file} RGID={self.analysis.basename} '
            'RGLB=lib1 RGPL=ILLUMINA RGPU={config.kit} '
            'RGSM={self.analysis.basename}',  # MAX_RECORDS_IN_RAM=100000
            ),
            self.analysis.logger
        )
        self.analysis.last_operation_filenames = output_file

        if retval != 0:
            self.analysis.logger.error(
                "Picard AddOrReplaceReadGroup exited with status %d" %
                retval)
            raise PipelineError("picard AddOrReplaceReadGroup error")

        output_file = self.analysis.basename + ".bam"
        os.rename(self.analysis.last_operation_filenames, output_file)
        self.analysis.last_operation_filenames = output_file

        if not os.path.exists(output_file):
            self.analysis.logger.error(
                "BAM file still does not exist, when it should have "
                "been created")
            raise PipelineError("picard failed")

        self.analysis.bamfile = os.path.join(
            os.getcwd(),
            self.analysis.last_operation_filenames)
        self.analysis.logger.info("Finished alignment")

    def sort_bam(self):
        self.analysis.logger.info("Sorting BAM")
        self.chdir()
        config = self.analysis.config

        output_file = self.analysis.basename + ".srt.bam"
        retval = utils.run_and_log(f(
            '{config.picard} SortSam '
            'I={self.analysis.last_operation_filenames} '
            'O={output_file} SO=coordinate'),
            self.analysis.logger
        )
        os.unlink(self.analysis.last_operation_filenames)
        self.analysis.last_operation_filenames = output_file

        if retval != 0:
            self.analysis.logger.error(
                "Picard SortSam exited with status %d" %
                retval)
            raise PipelineError("picard SortSam error")

        output_file = self.analysis.basename + ".srt.reorder.bam"
        retval = utils.run_and_log(f(
            '{config.picard} ReorderSam '
            'I={self.analysis.last_operation_filenames} '
            'O={output_file} R={config.genome_ref} '
            'CREATE_INDEX=true ',  # MAX_RECORDS_IN_RAM=1000000
            ),
            self.analysis.logger
        )
        os.unlink(self.analysis.last_operation_filenames)
        self.analysis.last_operation_filenames = output_file

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

        self.analysis.bamfile = os.path.join(
            os.getcwd(),
            self.analysis.last_operation_filenames)
        self.analysis.logger.info("Finished sorting")

    def indel_realign(self):
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.gatk} -T RealignerTargetCreator -R {config.genome_ref} '
            '-I {self.analysis.last_operation_filenames} -nt 15 -known '
            '{config.indel_1} -known {config.indel_2} -L {config.target_list} '
            '-ip 50 -o {self.output_basename}.realignment.intervals'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error(
                "Gatk RealignerTargetCreator exited with status %d" %
                retval)
            raise PipelineError("gatk RealignerTargetCreator error")

        output_file = self.analysis.basename + ".srt.realigned.bam"
        retval = utils.run_and_log(
            f('{config.gatk} -T IndelRealigner -R {config.genome_ref} -I '
              '{self.analysis.last_operation_filenames} '
              '-known {config.indel_1} -known {config.indel_2} '
              '-targetIntervals {self.output_basename}.realignment.intervals '
              '-o {output_file}'), self.analysis.logger)
        os.unlink(self.analysis.last_operation_filenames)
        os.unlink(self.analysis.last_operation_filenames[:-4] + ".bai")
        self.analysis.last_operation_filenames = output_file

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

        self.analysis.bamfile = os.path.join(
            os.getcwd(),
            self.analysis.last_operation_filenames)
        self.analysis.logger.info("Finished indel realignment")

    def recalibration(self):
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.gatk} -T BaseRecalibrator -R {config.genome_ref} -I '
            '{self.analysis.last_operation_filenames} -nct 20 -knownSites '
            '{config.dbsnp138} -o {self.output_basename}.recalibration.table',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error(
                "Gatk BaseRecalibrator exited with status %d" %
                retval)
            raise PipelineError("gatk BaseRecalibrator error")

        output_file = self.analysis.basename + ".srt.realigned.recal.bam"
        retval = utils.run_and_log(f(
            '{config.gatk} -T PrintReads -R {config.genome_ref} -I '
            '{self.analysis.last_operation_filenames} -nct 20 -BQSR '
            '{self.output_basename}.recalibration.table -o '
            '{output_file}',
            ),
            self.analysis.logger
        )
        os.unlink(self.analysis.last_operation_filenames)
        os.unlink(self.analysis.last_operation_filenames[:-4] + ".bai")
        self.analysis.last_operation_filenames = output_file

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

        self.analysis.bamfile = os.path.join(
            os.getcwd(),
            self.analysis.last_operation_filenames)
        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self):
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(
            f('{config.picard} CollectHsMetrics '
              'I={self.analysis.last_operation_filenames} '
              'BI={config.bait_list} '
              'TI={config.target_list} R={config.genome_ref} '
              'O={self.output_basename}.hs_metrics.txt '
              'MAX_RECORDS_IN_RAM=1500000',),
            self.analysis.logger)

        if retval != 0:
            self.analysis.logger.error(
                "Picard CollectHsMetrics exited with status %d" %
                retval)
            raise PipelineError("picard CollectHsMetrics error")

        retval = utils.run_and_log(
            f('{config.picard} CollectTargetedPcrMetrics '
              'I={self.analysis.last_operation_filenames} AMPLICON_INTERVALS='
              '{config.bait_list} TARGET_INTERVALS={config.target_list} '
              'R={config.genome_ref} '
              'O={self.output_basename}.targeted_metrics.txt '
              'MAX_RECORDS_IN_RAM=1500000 '
              'PER_BASE_COVERAGE={self.output_basename}.coverage.txt',),
            self.analysis.logger)

        if retval != 0:
            self.analysis.logger.error(
                "Picard CollectTargetedPcrMetrics exited with status %d" %
                retval)
            raise PipelineError("picard CollectTargetedPcrMetrics error")

        self.analysis.logger.info("Finished metrics collection")

    def variant_calling(self):
        self.analysis.logger.info("Running variant calling")
        self.chdir()
        retval = utils.run_and_log(
            f('{config.gatk} -T HaplotypeCaller -R {config.genome_ref} -I '
              '{self.analysis.last_operation_filenames} -ERC GVCF -nct 20 '
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
        self.cutadapt()
        self.fastqc()
        self.trim()
        self.align()
        self.sort_bam()
        self.indel_realign()
        self.recalibration()
        self.metrics_collection()
