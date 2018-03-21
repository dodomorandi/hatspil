import csv
import logging
import math
import os
import re
import shutil
import tempfile
from enum import Enum, auto

from . import utils
from .barcoded_filename import Analyte, BarcodedFilename
from .executor import Executor


class Aligner(Enum):
    NOVOALIGN = auto()
    BWA = auto()


class Mapping:
    def __init__(self, analysis, fastq_dir):
        self.analysis = analysis
        self.fastq_dir = fastq_dir

        self.sample_base = os.path.join(self.fastq_dir, self.analysis.sample)
        self.sample_base_out = os.path.join(self.fastq_dir, "REPORTS",
                                            self.analysis.sample)
        self.output_basename = os.path.join("REPORTS", self.analysis.basename)

        try:
            os.makedirs(
                os.path.join(self.analysis.get_bam_dir(), "REPORTS"),
                exist_ok=True)
        except OSError:
            pass

        try:
            os.makedirs(os.path.join(self.fastq_dir, "REPORTS"), exist_ok=True)
        except OSError:
            pass

        self.gatk_threads = self.analysis.parameters["gatk_threads"]
        self.max_records_str = utils.get_picard_max_records_string(
            self.analysis.parameters["picard_max_records"])

        self.sort_tempdir = os.path.join(self.analysis.get_bam_dir(),
                                         "%s_sort_tmp" % self.analysis.sample)

    def chdir(self):
        os.chdir(self.analysis.get_bam_dir())

    def cutadapt(self):
        self.analysis.logger.info("Cutting adapters")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            f'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAG '
            '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG '
            f'-m 20 -o "{{output_filename[0]}}" -p '
            f'"{{output_filename[1]}}" {{input_filename}} '
            f'> "{self.sample_base_out}.cutadapt.txt"',
            output_format=
            f"{self.analysis.sample}{{organism_str}}.clipped.R%d.fastq",
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            split_by_organism=True,
            output_path=self.fastq_dir,
            unlink_inputs=True,
            output_function=
            lambda filename: [filename % (index + 1) for index in range(2)])

        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self):
        self.analysis.logger.info("Running fastqc")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            f'{self.analysis.config.fastqc} '
            f'"{{input_filename}}" --outdir REPORTS',
            override_last_files=False)

        self.analysis.logger.info("Finished fastqc")

    def trim(self):
        trim_end = False
        if self.analysis.parameters["use_xenome"]:
            trim_end = True

        if trim_end:
            self.analysis.logger.info("Trimming first 5 bp and last 10 bp")
            trim_3 = self.analysis.parameters["trim_3"]
            if trim_3 is None:
                trim_end_cmd = "-e 10 "
            else:
                trim_end_cmd = "-e %d " % trim_3
        else:
            self.analysis.logger.info("Trimming first 5 bp")
            trim_end_cmd = ""

        self.chdir()
        config = self.analysis.config

        trim_5 = self.analysis.parameters["trim_5"]
        executor = Executor(self.analysis)
        executor(
            f'{config.seqtk} trimfq -b {trim_5} '
            f'{trim_end_cmd}'
            f'"{{input_filename}}" '
            f'> "{{output_filename}}"',
            output_format=os.path.join(
                self.fastq_dir, "%s{organism_str}.trimmed.R{read_index}.fastq"
                % self.analysis.sample),
            output_path=self.fastq_dir,
            unlink_inputs=True,
            error_string="Trimming with seqtk exited with status "
            "{status}",
            exception_string="trimming error")

        self.analysis.logger.info("Finished trimming")

    def filter_alignment(*args, **kwargs):
        """
        keep only aligned reads with maximum of N mismatches and without
        Ns, hard clipping and padding
        """
        if len(kwargs['input_filename']) != 1:
            raise "Expected a list with only one file"
        input_filename = kwargs["input_filename"][0]
        tmp_filename = input_filename + ".tmp"
        reSpaces = re.compile(R"\s+")
        reCigar = re.compile(R"N|H|P")
        with open(input_filename) as fd, open(tmp_filename, "w") as tmp_fd:
            for line in fd:
                if (line[0] == "@"):
                    tmp_fd.write(line)
                    continue

                params = reSpaces.split(line)
                if reCigar.search(params[5]):
                    continue

                read_length = len(params[9])

                params = params[11:]
                mutations = float("inf")
                for param in params:
                    if param.startswith("NM:i:"):
                        mutations = int(param[5:])
                        break

                if mutations > math.floor(read_length * 0.04):
                    continue

                tmp_fd.write(line)
        os.rename(tmp_filename, input_filename)

    def convert_alignment(self):
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        executor(
            self.filter_alignment,
            input_split_reads=False,
            split_by_organism=True,
            only_human=True,
            override_last_files=False)

        self.analysis.logger.info("Alignment SAM -> BAM")
        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'SamFormatConverter '
            f'I={{input_filename}} '
            f'O={{output_filename}}'
            f'{self.max_records_str}',
            output_format=f"{self.analysis.basename}{{organism_str}}.bam",
            error_string="Picard SamFormatConverter exited with status "
            "{status}",
            exception_string="picard SamFormatConverter error",
            unlink_inputs=True)

        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'AddOrReplaceReadGroups '
            f'I={{input_filename}} '
            f'O={{output_filename}} RGID={self.analysis.basename} '
            f'RGLB=lib1 RGPL=ILLUMINA RGPU={config.kit} '
            f'RGSM={self.analysis.basename}'
            f'{self.max_records_str}',
            output_format=f"{self.analysis.basename}{{organism_str}}.rg.bam",
            error_string=
            "Picard AddOrReplaceReadGroups exited with status {status}",
            exception_string="picard AddOrReplaceReadGroups error",
            unlink_inputs=True)

        executor(
            lambda **kwargs: os.rename(kwargs["input_filename"], kwargs["output_filename"]),
            output_format=f"{self.analysis.basename}{{organism_str}}.bam")

        self.analysis.logger.info("Finished alignment SAM->BAM")

    def align_novoalign(self):
        self.analysis.logger.info("Running alignment with NovoAlign")
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, "align.log")
            fh = logging.FileHandler(filename)
            self.analysis.logger.addHandler(fh)
            if barcoded.analyte == Analyte.WHOLE_EXOME:
                executor(
                    f'{config.novoalign} -oSAM "@RG\tID:{self.analysis.basename}\t'
                    f'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
                    f'-d {{genome_index}} '
                    f'-i PE {config.mean_len_library},{config.sd_len_library} '
                    f'-t 90 -f {{input_filename}}> {{output_filename}}',
                    input_function=lambda l: " ".join(sorted(l)),
                    input_split_reads=False,
                    output_format=
                    f"{self.analysis.basename}{{organism_str}}.sam",
                    split_by_organism=True,
                    only_human=True,
                    unlink_inputs=True)
            elif barcoded.analyte == Analyte.GENE_PANEL:
                executor(
                    f'{config.novoalign} '
                    f'-C '
                    f'-oSAM "@RG\tID:{self.analysis.basename}\t'
                    f'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
                    f'-d {{genome_index}} '
                    f'-i 50-500 -h 8 -H 20 --matchreward 3 -t 90 '
                    f'-f {{input_filename}}> {{output_filename}}',
                    input_function=lambda l: " ".join(sorted(l)),
                    input_split_reads=False,
                    output_format=
                    f"{self.analysis.basename}{{organism_str}}.sam",
                    split_by_organism=True,
                    only_human=True,
                    unlink_inputs=True)
            else:
                raise Exception("Unnhandled analyte")
            #  CSV NOVOALIGN
            with open(filename, 'r') as file_log, \
                    open(self.output_basename + "_novoalign.csv", 'w') \
                    as csv_file, \
                    open(self.output_basename + "_stat_novoalign.csv", "w") \
                    as stat_csv_file:
                writer = csv.writer(csv_file)
                writer_stat = csv.writer(stat_csv_file)
                is_csv = False
                is_stat = False
                values = []
                labels = []
                for line in file_log:
                    fields = line.split(":")
                    label = fields[0][1:].strip()

                    if is_stat is True:
                        if label == "No Mapping Found":
                            is_stat = False
                        values.append(fields[1].strip().split()[0])
                        labels.append(label)
                    elif label == "Paired Reads":
                        values.append(fields[1].strip().split()[0])
                        labels.append(label)
                        is_stat = True
                    else:
                        fields = line.split()
                        if is_csv is True:
                            if fields[1] == "Mean":
                                break
                            else:
                                writer.writerow(fields[1:4])
                        elif fields[1] == "From":
                            writer.writerow(fields[1:4])
                            is_csv = True
                writer_stat.writerow(labels)
                writer_stat.writerow(values)
            self.analysis.logger.removeHandler(fh)
            fh.close()
        self.analysis.logger.info(
            "Alignment finished. Aligner used: NovoAlign")

    def align_bwa(self):
        self.analysis.logger.info("Running alignment with BWA")
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        executor(
            f'{config.bwa} mem -t 6 -L 5,10 -v 1 {{genome_index}}.fasta '
            f'{{input_filename}}> {{output_filename}}',
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            output_format=f"{self.analysis.basename}{{organism_str}}.sam",
            split_by_organism=True,
            only_human=True,
            unlink_inputs=True)
        self.analysis.logger.info("Alignment finished. Aligner used: BWA")

    def sort_bam(self):
        self.analysis.logger.info("Sorting BAM(s)")
        self.chdir()
        config = self.analysis.config

        try:
            os.makedirs(self.sort_tempdir, exist_ok=True)
        except OSError:
            pass

        executor = Executor(self.analysis)
        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'SortSam '
            f'I={{input_filename}} '
            f'O={{output_filename}} SO=coordinate '
            f"TMP_DIR={self.sort_tempdir}"
            f'{self.max_records_str}',
            output_format=f"{self.analysis.basename}{{organism_str}}.srt.bam",
            error_string="Picard SortSam exited with status {status}",
            exception_string="picard SortSam error",
            unlink_inputs=True)

        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'ReorderSam '
            f'I={{input_filename}} '
            f'O={{output_filename}} R={{genome_ref}} '
            f'CREATE_INDEX=true'
            f'{self.max_records_str}',
            output_format=
            f"{self.analysis.basename}{{organism_str}}.srt.reorder.bam",
            error_string="Picard ReorderSam exited with status {status}",
            exception_string="picard ReorderSam error",
            unlink_inputs=True)

        if os.path.exists(self.sort_tempdir):
            shutil.rmtree(self.sort_tempdir)
        self.analysis.logger.info("Finished sorting")

    def mark_duplicates(self):
        self.analysis.logger.info("Marking duplicates")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        if barcoded.analyte == Analyte.WHOLE_EXOME:
            executor(
                f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
                f'MarkDuplicates '
                f'I={{input_filename}} '
                f'O={{output_filename}} '
                f'M={self.output_basename}{{organism_str}}.marked_dup_metrics.txt '
                f'CREATE_INDEX=true '
                f'{self.max_records_str}',
                output_format=
                f"{self.analysis.basename}{{organism_str}}.srt.marked.dup.bam",
                error_string=
                "Picard MarkDuplicates exited with status {status}",
                exception_string="picard MarkDuplicates error",
                unlink_inputs=True)
        elif barcoded.analyte == Analyte.GENE_PANEL:
            executor(
                f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
                f'FixMateInformation '
                f'I={{input_filename}} '
                f'O={{output_filename}} '
                f'ADD_MATE_CIGAR=true '
                f'IGNORE_MISSING_MATES=true '
                f'{self.max_records_str}',
                output_format=
                f"{self.analysis.basename}{{organism_str}}.srt.mc.bam",
                error_string=
                "Picard FixMateInformation exited with status {status}",
                exception_string="picard FixMateInformation error",
                unlink_inputs=True)

            executor(
                f'{config.samtools} view -H '
                f'{{input_filename}} > {{output_filename}}',
                output_format=
                f"{self.analysis.basename}{{organism_str}}.srt.mc.filtered.sam",
                error_string="samtools view exited with status {status}",
                exception_string="samtools view error",
                override_last_files=False)

            executor(
                f'{config.samtools} view '
                f'{{input_filename}} | grep "MC:" >> {{output_filename}}',
                output_format=
                f"{self.analysis.basename}{{organism_str}}.srt.mc.filtered.sam",
                error_string="samtools view exited with status {status}",
                exception_string="samtools view error",
                unlink_inputs=True)

            executor(
                f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
                f'UmiAwareMarkDuplicatesWithMateCigar '
                f'I={{input_filename}} '
                f'O={{output_filename}} '
                f'UMI_METRICS_FILE={self.output_basename}{{organism_str}}.UMI_metrics.txt '
                f'METRICS_FILE={self.output_basename}{{organism_str}}.marked_dup_metrics.txt '
                f'UMI_TAG_NAME=BX '
                f'CREATE_INDEX=true '
                f'TAGGING_POLICY=All '
                f'REMOVE_DUPLICATES=true '
                f'{self.max_records_str}',
                output_format=
                f"{self.analysis.basename}{{organism_str}}.srt.no_duplicates.bam",
                error_string=
                "Picard UmiAwareMarkDuplicatesWithMateCigar exited with status {status}",
                exception_string=
                "picard UmiAwareMarkDuplicatesWithMateCigar error",
                unlink_inputs=True)
        else:
            raise Exception("Unhandled analyte")

        self.analysis.logger.info("Finished marking duplicates")

    def indel_realign(self):
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T RealignerTargetCreator -R {{genome_ref}} '
            f'-I {{input_filename}} -nt {self.gatk_threads} '
            f'-known {config.indel_1} -known {config.indel_2} '
            f'-L {config.target_list} '
            f'-ip 50 -o {{output_filename}}',
            output_format=f"{self.output_basename}{{organism_str}}"
            ".realignment.intervals",
            error_string="Gatk RalignerTargetCreator exited with status "
            "{status}",
            exception_string="gatk RealignerTargetCreator error",
            override_last_files=False)

        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T IndelRealigner -R {{genome_ref}} '
            f'-I {{input_filename}} '
            f'-known {config.indel_1} -known {config.indel_2} '
            f'-targetIntervals {self.output_basename}{{organism_str}}'
            f'.realignment.intervals -o {{output_filename}}',
            output_format=f"{self.analysis.basename}{{organism_str}}"
            ".srt.realigned.bam",
            error_string="Gatk IndelRealigner exited with status {status}",
            exception_string="gatk IndelRealigner error",
            unlink_inputs=True)

        self.analysis.logger.info("Finished indel realignment")

    def _filter_non_hg(self, filename):
        organism = BarcodedFilename(filename).organism
        if organism is None or organism.lower().startswith("hg"):
            return filename
        else:
            return None

    def recalibration(self):
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T BaseRecalibrator -R {{genome_ref}} '
            f'-I {{input_filename}} -nct {self.gatk_threads} '
            f'-knownSites {{dbsnp}} '
            f'-o {self.output_basename}{{organism_str}}.recalibration.table',
            input_function=lambda filename: self._filter_non_hg(filename),
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            override_last_files=False)

        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T PrintReads -R {{genome_ref}} '
            f'-I {{input_filename}} -nct {self.gatk_threads} '
            f'-BQSR {self.output_basename}{{organism_str}}.recalibration.table '
            f'-o {{output_filename}}',
            output_format=f"{self.analysis.basename}{{organism_str}}"
            ".srt.realigned.recal.bam",
            error_string="Gatk PrintReads exited with status {status}",
            exception_string="gatk PrintReads error",
            unlink_inputs=True)

        if not self.analysis.parameters["run_post_recalibration"]:
            self.analysis.logger.info("Finished recalibration")
            return

        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T BaseRecalibrator -R {{genome_ref}} '
            f'-I {{input_filename}} -knownSites {{dbsnp}} '
            f'-L {config.target_list} -ip 50 '
            f'-nct {self.gatk_threads} '
            f'-o {self.output_basename}{{organism_str}}'
            f'.post_realignment.table',
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            override_last_files=False)

        executor(
            f'{config.java} {config.gatk_jvm_args} -jar {config.gatk} '
            f'-T AnalyzeCovariates -R {{genome_ref}} '
            f'-before {self.output_basename}{{organism_str}}.recalibration.table '
            f'-after {self.output_basename}{{organism_str}}.post_realignment.table '
            f'-plots {self.output_basename}{{organism_str}}.recalibration_plots.pdf',
            error_string="Gatk AnalyzeCovariates exited with status {status}",
            exception_string="gatk AnalyzeCovariates error",
            override_last_files=False)

        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'MarkDuplicates '
            f'I={{input_filename}} O={{output_filename}} '
            f'REMOVE_DUPLICATES=true '
            f'M={self.output_basename}{{organism_str}}.no_dup_metrics.txt '
            f'CREATE_INDEX=true'
            f'{self.max_records_str}',
            output_format=
            f"{self.analysis.basename}{{organism_str}}.srt.realigned.recal.no_dup.bam",
            error_string="Picard MarkDuplicates exited with status {status}",
            exception_string="picard MarkDuplicates error",
            unlink_inputs=True)

        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self):
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'CollectHsMetrics '
            f'I={{input_filename}} BI={config.bait_list} '
            f'TI={config.target_list} R={{genome_ref}} '
            f'O={self.output_basename}{{organism_str}}.hs_metrics.txt '
            f'MINIMUM_MAPPING_QUALITY=0 '
            f'MINIMUM_BASE_QUALITY=0 '
            f'COVERAGE_CAP=10000 '
            f'CLIP_OVERLAPPING_READS=false '
            f'PER_BASE_COVERAGE={self.output_basename}{{organism_str}}.coverage.txt'
            f'{self.max_records_str}',
            error_string="Picard CollectHsMetrics exited with status {status}",
            exception_string="picard CollectHsMetrics error",
            override_last_files=False)

        executor(
            f'{config.java} {config.picard_jvm_args} -jar {config.picard} '
            f'CollectGcBiasMetrics '
            f'R={{genome_ref}} I={{input_filename}} '
            f'O={self.output_basename}{{organism_str}}.gcbias.metrics.txt '
            f'CHART={self.output_basename}{{organism_str}}.gcbias_metrics.pdf '
            f'S={self.output_basename}{{organism_str}}.gcbias_summ_metrics.txt'
            f'{self.max_records_str}',
            error_string=
            "Picard CollectGcBiasMetrics exited with status {status}",
            exception_string="picard CollectGcBiasMetrics error",
            override_last_files=False)

        self.analysis.logger.info("Finished metrics collection")

    def bam2tdf(self):
        self.analysis.logger.info("Converting BAM to TDF")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f'{config.java} -jar {config.bam2tdf} -m 10 {{input_filename}}',
            error_string="Java bam2tdf exited with status {status}",
            exception_string="bam2tdf error",
            override_last_files=False)

    def compress_fastq(self):
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = utils.find_fastqs_by_organism(self.analysis.sample,
                                                    self.fastq_dir)
        for filenames in fastq_files.values():
            for filename, _ in filenames:
                utils.gzip(filename)
        self.analysis.logger.info("Finished compressing fastq files")

    def run(self):
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        if barcoded.analyte == Analyte.WHOLE_EXOME:
            if self.analysis.parameters["use_cutadapt"]:
                self.cutadapt()

        self.fastqc()
        self.trim()

        if self.analysis.parameters["aligner"] == Aligner.NOVOALIGN:
            self.align_novoalign()
        elif self.analysis.parameters["aligner"] == Aligner.BWA:
            self.align_bwa()
        else:
            raise Exception("unexpected aligner")
        self.convert_alignment()
        self.sort_bam()

        if self.analysis.parameters["mark_duplicates"]:
            self.mark_duplicates()

        if barcoded.analyte == Analyte.WHOLE_EXOME:
            self.indel_realign()
            self.recalibration()

        self.metrics_collection()
        if self.analysis.parameters["use_tdf"]:
            self.bam2tdf()
        if self.analysis.parameters["compress_fastq"]:
            self.compress_fastq()
