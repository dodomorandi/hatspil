import csv
import logging
import os
import shutil
import tempfile
from enum import Enum, auto

from . import utils
from .analysis import Analysis
from .barcoded_filename import Analyte, BarcodedFilename
from .executor import Executor


class GenericAligner(Enum):
    NOVOALIGN = auto()
    BWA = auto()


class RnaSeqAligner(Enum):
    STAR = auto()


class Aligner:
    def __init__(self, analysis: Analysis) -> None:
        self.analysis = analysis
        self.output_basename = os.path.join("REPORTS", self.analysis.basename)
        self.only_human = True

        self.max_records_str = utils.get_picard_max_records_string(
            self.analysis.parameters["picard_max_records"]
        )
        self.sort_tempdir = os.path.join(
            self.analysis.get_bam_dir(), "%s_sort_tmp" % self.analysis.sample
        )

    def chdir(self) -> None:
        os.chdir(self.analysis.get_bam_dir())

    def novoalign(self) -> None:
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
                    f"{config.novoalign} "
                    f'-oSAM "@RG\tID:{self.analysis.basename}\t'
                    f'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
                    f"-d {{genome_index}} "
                    f"-i PE {config.mean_len_library},{config.sd_len_library} "
                    f"-t 90 -f {{input_filename}}> {{output_filename}}",
                    input_function=lambda l: " ".join(sorted(l)),
                    input_split_reads=False,
                    output_format=f"{self.analysis.basename}{{organism_str}}.sam",
                    split_by_organism=True,
                    only_human=self.only_human,
                    unlink_inputs=True,
                )
            elif barcoded.analyte == Analyte.GENE_PANEL:
                executor(
                    f"{config.novoalign} "
                    f"-C "
                    f'-oSAM "@RG\tID:{self.analysis.basename}\t'
                    f'SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA" '
                    f"-d {{genome_index}} "
                    f"-i 50-500 -h 8 -H 20 --matchreward 3 -t 90 "
                    f"-f {{input_filename}}> {{output_filename}}",
                    input_function=lambda l: " ".join(sorted(l)),
                    input_split_reads=False,
                    output_format=f"{self.analysis.basename}{{organism_str}}.sam",
                    split_by_organism=True,
                    only_human=self.only_human,
                    unlink_inputs=True,
                )
            else:
                raise Exception("Unnhandled analyte")
            #  CSV NOVOALIGN
            with open(filename, "r") as file_log, open(
                self.output_basename + "_novoalign.csv", "w"
            ) as csv_file, open(
                self.output_basename + "_stat_novoalign.csv", "w"
            ) as stat_csv_file:
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
        self.analysis.logger.info("Alignment finished. Aligner used: NovoAlign")

    def bwa(self) -> None:
        self.analysis.logger.info("Running alignment with BWA")
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        executor(
            f"{config.bwa} mem -t 6 -L 5,10 -v 1 {{genome_ref}} "
            f"{{input_filename}}> {{output_filename}}",
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            output_format=f"{self.analysis.basename}{{organism_str}}.sam",
            split_by_organism=True,
            only_human=self.only_human,
            unlink_inputs=True,
        )
        self.analysis.logger.info("Alignment finished. Aligner used: BWA")

    def star(self) -> None:
        self.analysis.logger.info("Running alignment with STAR")
        config = self.analysis.config
        executor = Executor(self.analysis)

        bam_directory = self.analysis.get_bam_dir()

        def create_dirs(*args: str, **kwargs: str) -> None:
            output_path = os.path.join(bam_directory, kwargs["output_filename"])
            os.makedirs(output_path, exist_ok=True)
            os.makedirs(os.path.join(output_path, "index"), exist_ok=True)

        output_format = f"{self.analysis.basename}{{organism_str}}"
        star_dir_format = os.path.join(bam_directory, output_format + "_star")
        star_index_dir_format = os.path.join(star_dir_format, "index")
        executor(
            create_dirs,
            input_function=lambda l: " ".join(sorted(l)),
            output_format=star_dir_format,
            input_split_reads=False,
            override_last_files=False,
            split_by_organism=True,
            only_human=self.only_human,
        )
        star_index_format = "{getattr(config, 'star_index_' + organism)}"
        features_format = "{getattr(config, 'features_' + organism)}"

        self.analysis.logger.info("Step 1: Alignment 1st Pass:")

        executor(
            f"{config.star} --genomeDir {star_index_format} "
            f"--readFilesIn {{input_filename}} "
            f"--outFileNamePrefix {star_dir_format}/{output_format}. "
            f"--runThreadN 5 "
            f"--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 "
            f"--outFilterMismatchNmax 10 --alignIntronMax 500000 "
            f"--alignMatesGapMax 1000000 --sjdbScore 2 "
            f"--alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory "
            f"--outFilterMatchNminOverLread 0.33 "
            f"--outFilterScoreMinOverLread 0.33 "
            f"--sjdbOverhang 100 --outSAMstrandField intronMotif "
            f"--outSAMtype None --outSAMmode None",
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            override_last_files=False,
            split_by_organism=True,
            only_human=self.only_human,
        )
        self.analysis.logger.info("Finished step 1")

        self.analysis.logger.info("Step 2: Intermediate Index Generation:")

        executor(
            f"{config.star} "
            f"--runMode genomeGenerate "
            f"--genomeDir {star_index_dir_format} "
            f"--genomeFastaFiles {{genome_ref}} "
            f"--sjdbOverhang 100 "
            f"--runThreadN 5 "
            f"--sjdbFileChrStartEnd {star_dir_format}/"
            f"{output_format}.SJ.out.tab",
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            override_last_files=False,
            split_by_organism=True,
            only_human=self.only_human,
        )
        self.analysis.logger.info("Finished step 2")
        self.analysis.logger.info("Step 3: Alignment 2nd Pass")

        executor(
            f"{config.star} "
            f"--genomeDir {star_index_dir_format} "
            f"--readFilesIn {{input_filename}} "
            f"--outFileNamePrefix {star_dir_format}/{output_format}. "
            f"--runThreadN 5 "
            f"--outFilterMultimapScoreRange 1 "
            f"--outFilterMultimapNmax 20 "
            f"--outFilterMismatchNmax 10 "
            f"--alignIntronMax 500000 "
            f"--alignMatesGapMax 1000000 "
            f"--sjdbScore 2 "
            f"--alignSJDBoverhangMin 1 "
            f"--genomeLoad NoSharedMemory "
            f"--limitBAMsortRAM 0 "
            f"--outSAMattrRGline ID:{self.analysis.basename}\t"
            f"SM:{self.analysis.sample}\tLB:lib1\tPL:ILLUMINA "
            f"--outFilterMatchNminOverLread 0.33 "
            f"--outFilterScoreMinOverLread 0.33 "
            f"--sjdbOverhang 100 "
            f"--outSAMstrandField intronMotif "
            f"--outSAMattributes NH HI NM MD AS XS "
            f"--outSAMunmapped Within",
            input_function=lambda l: " ".join(sorted(l)),
            input_split_reads=False,
            output_format=f"{star_dir_format}/{output_format}.Aligned.out.sam",
            split_by_organism=True,
            only_human=self.only_human,
        )

        self.analysis.logger.info("Finished step 3")

        executor(
            lambda *args, **kwargs: os.rename(
                str(kwargs["input_filename"]), str(kwargs["output_filename"])
            ),
            split_by_organism=True,
            output_format=lambda *args, **kwargs: kwargs[
                "input_filename"
            ].filename.replace(".Aligned.out", ""),
            only_human=self.only_human,
        )

        executor(
            lambda *args, **kwargs: os.rename(
                str(kwargs["input_filename"]), str(kwargs["output_filename"])
            ),
            split_by_organism=True,
            output_format=os.path.join(bam_directory, f"{output_format}.sam"),
            only_human=self.only_human,
        )

        self.analysis.logger.info("Step 4: get HTseq count")

        counts_format = os.path.join(bam_directory, f"{output_format}_counts.txt")
        executor(
            f"{config.samtools} view -F 4 {{input_filename}} |"
            f"htseq-count "
            f"-m intersection-nonempty "
            f"-i gene_id "
            f"-r pos "
            f"-s no "
            f"- {features_format} "
            f"> {counts_format}",
            override_last_files=False,
            split_by_organism=True,
            only_human=self.only_human,
        )
        self.analysis.logger.info("Finished HTseq count")
        self.analysis.logger.info("Alignment finished. Aligner used: STAR")

    def sort_bam(self) -> None:
        self.analysis.logger.info("Sorting BAM(s)")
        self.chdir()
        config = self.analysis.config

        os.makedirs(self.sort_tempdir, exist_ok=True)

        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"SortSam "
            f"I={{input_filename}} "
            f"O={{output_filename}} SO=coordinate "
            f"TMP_DIR={self.sort_tempdir}"
            f"{self.max_records_str}",
            output_format=f"{self.analysis.basename}.srt{{organism_str}}.bam",
            error_string="Picard SortSam exited with status {status}",
            exception_string="picard SortSam error",
            only_human=self.only_human,
            split_by_organism=True,
            unlink_inputs=True,
        )

        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"ReorderSam "
            f"I={{input_filename}} "
            f"O={{output_filename}} R={{genome_ref}} "
            f"CREATE_INDEX=true"
            f"{self.max_records_str}",
            output_format=f"{self.analysis.basename}"
            f".srt.reorder{{organism_str}}.bam",
            error_string="Picard ReorderSam exited with status {status}",
            exception_string="picard ReorderSam error",
            only_human=self.only_human,
            split_by_organism=True,
            unlink_inputs=True,
        )

        if os.path.exists(self.sort_tempdir):
            shutil.rmtree(self.sort_tempdir)
        self.analysis.logger.info("Finished sorting")

    def run(self) -> None:
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)

        if barcoded.analyte == Analyte.RNASEQ:
            if self.analysis.parameters["rnaseq_aligner"] == RnaSeqAligner.STAR:
                self.star()
            else:
                raise Exception("unexpected aligner for this type of sample")
        else:
            if self.analysis.parameters["aligner"] == GenericAligner.NOVOALIGN:
                self.novoalign()
            elif self.analysis.parameters["aligner"] == GenericAligner.BWA:
                self.bwa()
            else:
                raise Exception("unexpected aligner for this type of sample")
