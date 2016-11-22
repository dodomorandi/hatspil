import os
from formatizer import f
from .exceptions import PipelineError
from . import utils


class Haloplex:

    def __init__(self, analysis, fastq_dir):
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.sample_base = "%s/%s" % (self.fastq_dir, self.analysis.sample)
        self.fastq_R1 = self.sample_base + "_R1.fastq"
        self.fastq_R2 = self.sample_base + "_R2.fastq"

    def chdir(self):
        os.chdir(self.analysis.bam_dir)

    def cutadapt(self):
        self.analysis.logger.info("Cutting adapters")
        self.chdir()
        retval = utils.run_and_log(f(
            'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG '
            '-m 20 -o "{self.sample_base}.clipped.R1.fastq" -p '
            '"{self.sample_base}.clipped.R2.fastq" {self.fastq_R1} {self.fastq_R2} > '
            '"{self.sample_base}.cutadapt.txt"'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Cutadapt exited with status %d" % retval)
            raise PipelineError("cutadapt error")

        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self):
        self.analysis.logger.info("Running fastqc")
        self.chdir()
        config = self.analysis.config

        qcdir = os.path.join(self.fastq_dir, "FASTQC")
        if not os.path.exists(qcdir):
            os.mkdir(qcdir)

        for read_value in ("R1", "R2"):
            retval = utils.run_and_log(f(
                '{config.fastqc} "{self.sample_base}.clipped.{read_value}.fastq"'
                ' --outdir {qcdir}'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error("fastqc exited with status %d" % retval)
                raise PipelineError("fastqc error")

        self.analysis.logger.info("Finished fastqc")

    def trim(self):
        self.analysis.logger.info("Trimming first 5 bp")
        self.chdir()
        config = self.analysis.config

        for read_value in ("R1", "R2"):
            retval = utils.run_and_log(f(
                '{config.seqtk} trimfq -b 5 "{self.sample_base}.clipped.{read_value}.fastq" '
                '> "{self.sample_base}.trimmed.{read_value}.fastq"'),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error("Trimming with seqtk exited with status %d" % retval)
                raise PipelineError("trimming error")

            os.unlink(f(
                "{self.sample_base}.clipped.{read_value}.fastq")
            )

        self.analysis.logger.info("Finished trimming")

    def align(self):
        self.analysis.logger.info("Running alignment")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.novoalign} -oSAM "@RG\tID:{self.analysis.basename}\tSM:'
            '{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" -d {config.hg19_index} '
            '-i PE {config.mean_len_library},{config.sd_len_library} -t 90 '
            '-f "{self.sample_base}.trimmed.R1.fastq" '
            '"{self.sample_base}.trimmed.R2.fastq"> {self.analysis.basename}.sam'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Novoalign exited with status %d" % retval)
            raise PipelineError("novoalign error")

        self.analysis.logger.info("Alignment SAM -> BAM")
        retval = utils.run_and_log(f(
            '{config.picard} SamFormatConverter I={self.analysis.basename}.sam '
            'O={self.analysis.basename}.bam',  # MAX_RECORDS_IN_RAM=1000000
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard SamFormatConverter exited with status %d" % retval)
            raise PipelineError("picard SamFormatConverter error")

        retval = utils.run_and_log(f(
            '{config.picard} AddOrReplaceReadGroups I={self.analysis.basename}.bam '
            'O={self.analysis.basename}.rg.bam RGID={self.analysis.basename} '
            'RGLB=lib1 RGPL=ILLUMINA RGPU={config.kit} '
            'RGSM={self.analysis.basename}',  # MAX_RECORDS_IN_RAM=100000
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard AddOrReplaceReadGroup exited with status %d" % retval)
            raise PipelineError("picard AddOrReplaceReadGroup error")

        os.rename("%s.rg.bam" % self.analysis.basename, "%s.bam" % self.analysis.basename)

        if not os.path.exists("%s.bam" % self.analysis.basename):
            self.analysis.logger.error("BAM file still does not exist, when it should have been created")
            raise PipelineError("picard failed")

        self.analysis.bamfile = os.path.join(os.getcwd(), "%s.bam" % self.analysis.basename)

        os.unlink("%s.sam" % self.analysis.basename)
        self.analysis.logger.info("Finished alignment")

    def sort_bam(self):
        self.analysis.logger.info("Sorting BAM")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.picard} SortSam I={self.analysis.basename}.bam '
            'O={self.analysis.basename}.srt.bam SO=coordinate'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard SortSam exited with status %d" % retval)
            raise PipelineError("picard SortSam error")

        retval = utils.run_and_log(f(
            '{config.picard} ReorderSam I={self.analysis.basename}.srt.bam '
            'O={self.analysis.basename}.srt.reorder.bam R={config.genome_ref} '
            'CREATE_INDEX=true ',  # MAX_RECORDS_IN_RAM=1000000
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard SortSam exited with status %d" % retval)
            raise PipelineError("picard SortSam")

        if not os.path.exists("%s.srt.reorder.bam" % self.analysis.basename):
            self.analysis.logger.error("BAM file still does not exist, when it should have been created")
            raise PipelineError("picard failed")

        self.analysis.bamfile = os.path.join(os.getcwd(), "%s.srt.reorder.bam" % self.analysis.basename)

        os.unlink("%s.bam" % self.analysis.basename)
        os.unlink("%s.srt.bam" % self.analysis.basename)

        self.analysis.logger.info("Finished sorting")

    def indel_realign(self):
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.gatk} -T RealignerTargetCreator -R {config.genome_ref} '
            '-I {self.analysis.basename}.srt.reorder.bam -nt 15 -known {config.indel_1} '
            '-known {config.indel_2} -L {config.target_list} -ip 50 -o '
            '{self.analysis.basename}.realignment.intervals'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk RealignerTargetCreator exited with status %d" % retval)
            raise PipelineError("gatk RealignerTargetCreator error")

        retval = utils.run_and_log(f(
            '{config.gatk} -T IndelRealigner -R {config.genome_ref} -I '
            '{self.analysis.basename}.srt.reorder.bam -known {config.indel_1} -known '
            '{config.indel_2} -targetIntervals '
            '{self.analysis.basename}.realignment.intervals -o '
            '{self.analysis.basename}.srt.realigned.bam'),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk IndelRealigner exited with status %d" % retval)
            raise PipelineError("gatk IndelRealigner failed")

        if not os.path.exists("%s.srt.realigned.bam" % self.analysis.basename):
            self.analysis.logger.error("BAM file still does not exist, when it should have been created")
            raise PipelineError("gatk failed")

        self.analysis.bamfile = os.path.join(os.getcwd(), "%s.srt.realigned.bam" % self.analysis.basename)

        os.unlink("%s.srt.reorder.bam" % self.analysis.basename)
        os.unlink("%s.srt.reorder.bai" % self.analysis.basename)

        self.analysis.logger.info("Finished indel realignment")

    def recalibration(self):
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.gatk} -T BaseRecalibrator -R {config.genome_ref} -I '
            '{self.analysis.basename}.srt.realigned.bam -nct 20 -knownSites '
            '{config.dbsnp138} -o {self.analysis.basename}.recalibration.table',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk BaseRecalibrator exited with status %d" % retval)
            raise PipelineError("gatk BaseRecalibrator error")

        retval = utils.run_and_log(f(
            '{config.gatk} -T PrintReads -R {config.genome_ref} -I '
            '{self.analysis.basename}.srt.realigned.bam -nct 20 -BQSR '
            '{self.analysis.basename}.recalibration.table -o '
            '{self.analysis.basename}.srt.realigned.recal.bam',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk PrintReads exited with status %d" % retval)
            raise PipelineError("gatk PrintReads error")

        if not os.path.exists("%s.srt.realigned.recal.bam" % self.analysis.basename):
            self.analysis.logger.error("BAM file still does not exist, when it should have been created")
            raise PipelineError("gatk failed")

        self.analysis.bamfile = os.path.join(os.getcwd(), "%s.srt.realigned.recal.bam" % self.analysis.basename)

        os.unlink("%s.srt.realigned.bam" % self.analysis.basename)
        os.unlink("%s.srt.realigned.bai" % self.analysis.basename)

        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self):
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.picard} CollectHsMetrics '
            'I={self.analysis.basename}.srt.realigned.recal.bam BI={config.bait_list} '
            'TI={config.target_list} R={config.genome_ref} '
            'O={self.analysis.basename}.hs_metrics.txt MAX_RECORDS_IN_RAM=1500000',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard CollectHsMetrics exited with status %d" % retval)
            raise PipelineError("picard CollectHsMetrics error")

        retval = utils.run_and_log(f(
            '{config.picard} CollectTargetedPcrMetrics '
            'I={self.analysis.basename}.srt.realigned.recal.bam AMPLICON_INTERVALS='
            '{config.bait_list} TARGET_INTERVALS={config.target_list} '
            'R={config.genome_ref} O={self.analysis.basename}.targeted_metrics.txt '
            'MAX_RECORDS_IN_RAM=1500000 '
            'PER_BASE_COVERAGE={self.analysis.basename}.coverage.txt',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Picard CollectTargetedPcrMetrics exited with status %d" % retval)
            raise PipelineError("picard CollectTargetedPcrMetrics error")

        self.analysis.logger.info("Finished metrics collection")

    def variant_calling(self):
        self.analysis.logger.info("Running variant calling")
        self.chdir()
        retval = utils.run_and_log(f(
            '{config.gatk} -T HaplotypeCaller -R {config.genome_ref} -I '
            '{self.analysis.basename}.srt.realigned.recal.bam -ERC GVCF -nct 20 '
            '--genotyping_mode DISCOVERY -L {config.target_list} '
            '-stand_emit_conf 10 -stand_call_conf 30 --output_mode '
            'EMIT_VARIANTS_ONLY -o {self.analysis.basename}.variants.g.vcf',
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk HaplotypeCaller exited with status %d" % retval)
            raise PipelineError("gatk HaplotypeCaller error")

        retval = utils.run_and_log(f(
            '{config.gatk} -T VariantFiltration -R {config.genome_ref} '
            '--variant {self.analysis.basename}.variants.vcf -o '
            '{self.analysis.basename}.variants.filtered.vcf --clusterWindowSize 10 '
            '--clusterSize 3 --filterExpression "DP < 10" --filterName '
            '"LowCoverage" --filterExpression "QUAL < 30.0" --filterName '
            '"VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" '
            '--filterName "LowQual" --filterExpression "QD < 1.5" --filterName '
            '"LowQD" ',  # --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE"
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Gatk VariantFiltration exited with status %d" % retval)
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
