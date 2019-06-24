"""Module to perform sequences mapping.

NGS data evaluation is not just about alignment, there are many steps
that are necessary or extremely informative and that have to be handled
carefully. This module is about all these steps.
"""
import math
import os
import re
from typing import Any, List, Optional, Sequence, Union

from .aligner import Aligner, RnaSeqAligner
from .core import utils
from .core.analysis import Analysis
from .core.barcoded_filename import Analyte, BarcodedFilename
from .core.exceptions import PipelineError
from .core.executor import AnalysisFileData, Executor
from .db import Db
from .db.cutadapt import Cutadapt
from .db.picard_metrics import PicardMetrics, PicardMetricsType
from .xenograft import Xenograft, XenograftClassifier


class Mapping:
    """Handle NGS to produce meaningful results.

    This class does not strictly handle sequence alignment, but all the
    steps that should be performed to produce meaningful results. There
    are both pre-alignment and post-alignment steps, and some of them
    depends on the sample type and/or command-line parameters.
    """

    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        """Create a new instance.

        The reports directories are also created if needed.
        """
        self.analysis = analysis
        self.fastq_dir = fastq_dir

        self.sample_base = os.path.join(self.fastq_dir, self.analysis.sample)
        self.sample_base_out = os.path.join(
            self.fastq_dir, "REPORTS", self.analysis.sample
        )
        self.output_basename = os.path.join("REPORTS", self.analysis.basename)

        os.makedirs(os.path.join(self.analysis.get_bam_dir(), "REPORTS"), exist_ok=True)
        os.makedirs(os.path.join(self.fastq_dir, "REPORTS"), exist_ok=True)

        self.gatk_threads = self.analysis.parameters["gatk_threads"]
        self.max_records_str = utils.get_picard_max_records_string(
            self.analysis.parameters["picard_max_records"]
        )

    def chdir(self) -> None:
        """Change current directory to the BAM folder."""
        os.chdir(self.analysis.get_bam_dir())

    def cutadapt(self) -> None:
        """Run cutadapt software on the current sample data.

        Cutadapt software is run depending on the fact the NGS data is
        single-end or paired-end.

        The process uses the adapter information specified in
        configuration under the kit needed for the sample.

        Results are saved in the database if available.
        """
        self.analysis.logger.info("Cutting adapters")
        self.chdir()
        config = self.analysis.config

        report_filename = f"{self.sample_base_out}.cutadapt.txt"
        if isinstance(self.analysis.last_operation_filenames, list):
            is_paired_end = len(self.analysis.last_operation_filenames) == 2
        elif isinstance(self.analysis.last_operation_filenames, dict):
            first_filenames = next(
                iter(self.analysis.last_operation_filenames.values())
            )
            is_paired_end = len(first_filenames) == 2
            assert all(
                map(
                    lambda filenames: len(filenames) == len(first_filenames),
                    self.analysis.last_operation_filenames.values(),
                )
            )
        else:
            is_paired_end = False

        executor = Executor(self.analysis)
        barcoded_sample = BarcodedFilename.from_sample(self.analysis.sample)
        kit = utils.get_kit_from_barcoded(self.analysis.config, barcoded_sample)
        assert kit
        if is_paired_end:
            if not kit.adapter_r1 or not kit.adapter_r2:
                self.analysis.logger.warning(
                    "cutadapt needs both adapter to perform "
                    "a paired-end trimming. Skipping cutadapt"
                )
                return

            executor(
                f"cutadapt -a {{kit.adapter_r1}} "
                f"-A {{kit.adapter_r2}} "
                f'-m 20 -o "{{output_filename[0]}}" -p '
                f'"{{output_filename[1]}}" {{input_filename}} '
                f'> "{report_filename}"',
                output_format=f"{self.analysis.sample}.clipped{{organism_str}}"
                ".R%d.fastq",
                input_function=lambda l: " ".join(sorted(l)),
                input_split_reads=False,
                split_by_organism=True,
                output_path=self.fastq_dir,
                unlink_inputs=True,
                output_function=lambda filename: [
                    filename % (index + 1) for index in range(2)
                ],
            )
        else:
            if not kit.adapter_r1:
                self.analysis.logger.warning(
                    "cutadapt needs the R1 adapter to perform "
                    "a single end trimming. Skipping cutadapt"
                )
                return

            executor(
                f"cutadapt -a {{kit.adapter_r1}} "
                f'-m 20 -o "{{output_filename}}" '
                f" {{input_filename}} "
                f'> "{report_filename}"',
                output_format=f"{self.analysis.sample}.clipped{{organism_str}}"
                ".R1.fastq",
                input_split_reads=False,
                split_by_organism=True,
                output_path=self.fastq_dir,
                unlink_inputs=True,
            )

        db = Db(config)
        cutadapt = Cutadapt(db)
        cutadapt.store_from_file(self.analysis, report_filename)

        self.analysis.logger.info("Finished cutting adapters")

    def fastqc(self) -> None:
        """Run FastQC tool on the input files.

        Results are placed inside BAM reports directory.
        """
        self.analysis.logger.info("Running fastqc")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            f'{self.analysis.config.fastqc} "{{input_filename}}" --outdir REPORTS',
            override_last_files=False,
        )

        self.analysis.logger.info("Finished fastqc")

    def trim(self) -> None:
        """Trim the 5' and 3' ends of the reads.

        Uses SeqTK TrimFQ tool to perform a trimming operation depending
        on the command line parameters.
        """
        trim_end = False
        if self.analysis.parameters["use_xenograft_classifier"]:
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
            f"{config.seqtk} trimfq -b {trim_5} "
            f"{trim_end_cmd}"
            f'"{{input_filename}}" '
            f'> "{{output_filename}}"',
            output_format=os.path.join(
                self.fastq_dir,
                "%s.trimmed{organism_str}"
                ".R{input_filename.barcode.read_index}.fastq" % self.analysis.sample,
            ),
            output_path=self.fastq_dir,
            split_by_organism=True,
            unlink_inputs=True,
            error_string="Trimming with seqtk exited with status {status}",
            exception_string="trimming error",
        )

        self.analysis.logger.info("Finished trimming")

    def _filter_alignment(*args, **kwargs: Sequence[AnalysisFileData]) -> None:
        """Perform some SAM filtering.

        Keep only aligned reads with maximum of N mismatches and without
        Ns, hard clipping and padding.
        """
        if len(kwargs["input_filenames"]) != 1:
            raise PipelineError("Expected a list with only one file")
        input_filename: str = kwargs["input_filenames"][0].filename
        tmp_filename = input_filename + ".tmp"
        reSpaces = re.compile(r"\s+")
        reCigar = re.compile(r"N|H|P")
        with open(input_filename) as fd, open(tmp_filename, "w") as tmp_fd:
            for line in fd:
                if line[0] == "@":
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

    def convert_to_bam(self) -> None:
        """Filter the SAM file and convert the content to BAM.

        The filtering process discards all the reads with a number of
        mutations above 4% of the length of the read and in case of Ns,
        hard clipping or padding.

        The conversion to BAM is performed by Picard SamFormatConverter.
        """
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        executor(
            self._filter_alignment,
            input_split_reads=False,
            split_by_organism=True,
            only_human=True,
            override_last_files=False,
        )

        self.analysis.logger.info("Alignment SAM -> BAM")
        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"SamFormatConverter "
            f"I={{input_filename}} "
            f"O={{output_filename}} "
            f"TMP_DIR={config.temporary_dir}"
            f"{self.max_records_str}",
            output_format=f"{self.analysis.basename}{{organism_str}}.bam",
            error_string="Picard SamFormatConverter exited with status {status}",
            exception_string="picard SamFormatConverter error",
            split_by_organism=True,
            unlink_inputs=True,
        )
        self.analysis.logger.info("Finished alignment SAM->BAM")

    def create_bam_index(self) -> None:
        """Sort BAM and create BAI.

        Both processes are performed by Samtools.
        """
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        self.analysis.logger.info("Sorting BAM and creating BAI file")
        executor(
            f"{config.samtools} sort {{input_filename}} -O BAM -o {{output_filename}}",
            unlink_inputs=True,
            split_by_organism=True,
            output_format=f"{self.analysis.basename}.srt{{organism_str}}.bam",
        )

        executor(
            f"{config.samtools} index {{input_filename}}",
            split_by_organism=True,
            override_last_files=False,
        )

        def rename_files(*args: Any, **kwargs: Any) -> None:
            bam_filename = kwargs["input_filename"].filename
            assert os.path.splitext(bam_filename)[1].lower() == ".bam"
            output_bam_filename = kwargs["output_filename"]
            os.rename(bam_filename, output_bam_filename)
            os.rename(f"{bam_filename}.bai", f"{output_bam_filename}.bai")

        executor(
            rename_files,
            split_by_organism=True,
            output_format=f"{self.analysis.basename}{{organism_str}}.bam",
        )

        self.analysis.logger.info("Finished sorting and creating BAI file")

    def add_bam_groups(self) -> None:
        """Add/replace read groups to BAM.

        It uses Picard AddOrReplaceReadGroups to correctly handle the
        groups for the analysis. At the end of the process
        `Mapping.create_bam_index` is run.
        """
        self.chdir()
        config = self.analysis.config
        executor = Executor(self.analysis)
        self.analysis.logger.info("Adding/replacing read groups to BAM")
        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"AddOrReplaceReadGroups "
            f"I={{input_filename}} "
            f"O={{output_filename}} RGID={self.analysis.basename} "
            f"RGLB=lib1 RGPL=ILLUMINA RGPU={{kit.name}} "
            f"RGSM={self.analysis.basename} "
            f"TMP_DIR={config.temporary_dir}"
            f"{self.max_records_str}",
            output_format=f"{self.analysis.basename}.rg{{organism_str}}.bam",
            error_string="Picard AddOrReplaceReadGroups exited with status {status}",
            exception_string="picard AddOrReplaceReadGroups error",
            split_by_organism=True,
            unlink_inputs=True,
        )

        executor(
            lambda **kwargs: os.rename(
                kwargs["input_filename"].filename, kwargs["output_filename"]
            ),
            output_format=f"{self.analysis.basename}{{organism_str}}.bam",
            split_by_organism=True,
        )
        self.create_bam_index()

        self.analysis.logger.info("Finished adding/replacing read groups to BAM")

    def mark_duplicates(self) -> None:
        """Mark duplicates depending on the analyte.

        In case the sample is a whole exome, Picard MarkDuplicates is
        used. Otherwise, when the sample is a gene panel, there are more
        steps involved. Picard FixMateInformation is used to perform the
        mate-pair information, then the reads with the mate-pair tag are
        given to Picard UmiAwareMarkDuplicatesWithMateCigar.

        In both cases the results are stored in the database if
        possible.

        At to date only these two type of analyte are handled.
        """
        self.analysis.logger.info("Marking duplicates")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        if barcoded.analyte == Analyte.WHOLE_EXOME:
            executor(
                f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
                f"MarkDuplicates "
                f"I={{input_filename}} "
                f"O={{output_filename}} "
                f"M={self.output_basename}.marked_dup_metrics"
                f"{{organism_str}}.txt "
                f"CREATE_INDEX=true "
                f"TMP_DIR={config.temporary_dir}"
                f"{self.max_records_str}",
                output_format=f"{self.analysis.basename}.srt.marked.dup"
                f"{{organism_str}}.bam",
                error_string="Picard MarkDuplicates exited with status {status}",
                exception_string="picard MarkDuplicates error",
                split_by_organism=True,
                unlink_inputs=True,
            )

            db = Db(self.analysis.config)
            picard_metrics = PicardMetrics(db)
            executor(
                lambda *args, **kwargs: picard_metrics.store_from_file(
                    self.analysis,
                    f"{self.output_basename}.marked_dup_metrics{{}}.txt".format(
                        kwargs["organism_str"]
                    ),
                    PicardMetricsType.marked_duplicates,
                ),
                split_by_organism=True,
                override_last_files=False,
            )
        elif barcoded.analyte == Analyte.GENE_PANEL:
            executor(
                f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
                f"FixMateInformation "
                f"I={{input_filename}} "
                f"O={{output_filename}} "
                f"ADD_MATE_CIGAR=true "
                f"IGNORE_MISSING_MATES=true "
                f"TMP_DIR={config.temporary_dir}"
                f"{self.max_records_str}",
                output_format=f"{self.analysis.basename}.srt.mc"
                f"{{organism_str}}.bam",
                error_string="Picard FixMateInformation exited with status {status}",
                exception_string="picard FixMateInformation error",
                split_by_organism=True,
                unlink_inputs=True,
            )

            executor(
                f"{config.samtools} view -H "
                f"{{input_filename}} > {{output_filename}}",
                output_format=f"{self.analysis.basename}.srt.mc.filtered"
                f"{{organism_str}}.sam",
                error_string="samtools view exited with status {status}",
                exception_string="samtools view error",
                split_by_organism=True,
                override_last_files=False,
            )

            executor(
                f"{config.samtools} view "
                f'{{input_filename}} | grep "MC:" >> {{output_filename}}',
                output_format=f"{self.analysis.basename}.str.mc.filtered"
                f"{{organism_str}}.sam",
                error_string="samtools view exited with status {status}",
                exception_string="samtools view error",
                split_by_organism=True,
                unlink_inputs=True,
            )

            executor(
                f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
                f"UmiAwareMarkDuplicatesWithMateCigar "
                f"I={{input_filename}} "
                f"O={{output_filename}} "
                f"UMI_METRICS_FILE={self.output_basename}.UMI_metrics"
                f"{{organism_str}}.txt "
                f"METRICS_FILE={self.output_basename}.marked_dup_metrics"
                f"{{organism_str}}.txt"
                f"UMI_TAG_NAME=BX "
                f"CREATE_INDEX=true "
                f"TAGGING_POLICY=All "
                f"REMOVE_DUPLICATES=true "
                f"TMP_DIR={config.temporary_dir}"
                f"{self.max_records_str}",
                output_format=f"{self.analysis.basename}.srt.no_duplicates"
                f"{{organism_str}}.bam",
                error_string="Picard UmiAwareMarkDuplicatesWithMateCigar "
                "exited with status {status}",
                exception_string="picard UmiAwareMarkDuplicatesWithMateCigar error",
                split_by_organism=True,
                unlink_inputs=True,
            )

            db = Db(self.analysis.config)
            picard_metrics = PicardMetrics(db)
            executor(
                lambda *args, **kwargs: picard_metrics.store_from_file(
                    self.analysis,
                    f"{self.output_basename}.marked_dup_metrics{{}}.txt".format(
                        kwargs["organism_str"]
                    ),
                    PicardMetricsType.marked_duplicates,
                ),
                split_by_organism=True,
                override_last_files=False,
            )

            executor(
                lambda *args, **kwargs: picard_metrics.store_from_file(
                    self.analysis,
                    f"{self.output_basename}.UMI_metrics{{}}.txt".format(
                        kwargs["organism_str"]
                    ),
                    PicardMetricsType.umi,
                ),
                split_by_organism=True,
                override_last_files=False,
            )

        else:
            raise Exception("Unhandled analyte")

        self.analysis.logger.info("Finished marking duplicates")

    def indel_realign(self) -> None:
        """Perform an indel realignment.

        A realigner target is created using GATK RealignerTargetCreator
        and the realigned is peformed using GATK IndelRealigner.
        The parameters are extracted from the configuration in
        relationship of the kit and the analyte for the current sample.
        """
        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        self.analysis.logger.info("Running indel realignment")
        self.chdir()
        config = self.analysis.config
        rnaseq_parameter = ""
        if barcoded.analyte == Analyte.RNASEQ:
            rnaseq_parameter = "-U ALLOW_N_CIGAR_READS "
        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T RealignerTargetCreator -R {{genome_ref}} "
            f"-I {{input_filename}} -nt {self.gatk_threads} "
            f"-known {{indel_1}} -known {{indel_2}} "
            f"{rnaseq_parameter}"
            f"-L {{kit.target_list}} "
            f"-ip 50 -o {{output_filename}}",
            output_format=f"{self.output_basename}.realignment"
            f"{{organism_str}}.intervals",
            error_string="Gatk RalignerTargetCreator exited with status {status}",
            exception_string="gatk RealignerTargetCreator error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T IndelRealigner -R {{genome_ref}} "
            f"-I {{input_filename}} "
            f"-known {{indel_1}} -known {{indel_2}} "
            f"{rnaseq_parameter}"
            f"-targetIntervals {self.output_basename}.realignment"
            f"{{organism_str}}.intervals "
            f"-o {{output_filename}}",
            output_format=f"{self.analysis.basename}.srt.realigned"
            f"{{organism_str}}.bam",
            error_string="Gatk IndelRealigner exited with status {status}",
            exception_string="gatk IndelRealigner error",
            split_by_organism=True,
            unlink_inputs=True,
        )

        self.analysis.logger.info("Finished indel realignment")

    def _filter_non_hg(self, filename: Union[str, List[str]]) -> Optional[str]:
        assert isinstance(filename, str)
        organism = BarcodedFilename(filename).organism
        if organism is None or organism.lower().startswith("hg"):
            return filename
        else:
            return None

    def recalibration(self) -> None:
        """Perform a base recalibration step.

        Uses GATK BaseRecalibrator to create a recalibration table, then
        GATK PrintReads is used to obtain the recalibrated file.

        In case the "run_post_recalibration" is set from command-line,
        a further step is performed. The BaseRecalibrator is run again
        with a target list, GATK AnalyzeCovariates is run to generate
        useful plots and Picard MarkDuplicates is run. The metrics data
        of this last step is collected in the database, if possible.
        """
        self.analysis.logger.info("Running base recalibration")
        self.chdir()
        config = self.analysis.config
        rnaseq_parameter = ""
        if self.analysis.parameters["aligner"] in RnaSeqAligner:
            rnaseq_parameter = "-U ALLOW_N_CIGAR_READS "
        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T BaseRecalibrator -R {{genome_ref}} "
            f"-I {{input_filename}} -nct {self.gatk_threads} "
            f"-knownSites {{dbsnp}} "
            f"{rnaseq_parameter}"
            f"-o {self.output_basename}.recalibration"
            f"{{organism_str}}.table",
            input_function=lambda filename: self._filter_non_hg(filename),
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T PrintReads -R {{genome_ref}} "
            f"{rnaseq_parameter}"
            f"-I {{input_filename}} -nct {self.gatk_threads} "
            f"-BQSR {self.output_basename}.recalibration"
            f"{{organism_str}}.table "
            f"-o {{output_filename}}",
            output_format=f"{self.analysis.basename}.srt.realigned.recal"
            f"{{organism_str}}.bam",
            error_string="Gatk PrintReads exited with status {status}",
            exception_string="gatk PrintReads error",
            split_by_organism=True,
            unlink_inputs=True,
        )

        if not self.analysis.parameters["run_post_recalibration"]:
            self.analysis.logger.info("Finished recalibration")
            return

        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T BaseRecalibrator -R {{genome_ref}} "
            f"{rnaseq_parameter}"
            f"-I {{input_filename}} -knownSites {{dbsnp}} "
            f"-L {{kit.target_list}} -ip 50 "
            f"-nct {self.gatk_threads} "
            f"-o {self.output_basename}.post_realignment"
            f"{{organism_str}}.table",
            error_string="Gatk BaseRecalibrator exited with status {status}",
            exception_string="gatk BaseRecalibrator error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            f"{config.java} {config.gatk_jvm_args} -jar {config.gatk} "
            f"-T AnalyzeCovariates -R {{genome_ref}} "
            f"{rnaseq_parameter}"
            f"-before {self.output_basename}.recalibration"
            f"{{organism_str}}.table "
            f"-after {self.output_basename}.post_realignment"
            f"{{organism_str}}.table "
            f"-plots {self.output_basename}.recalibration_plots"
            f"{{organism_str}}.pdf",
            error_string="Gatk AnalyzeCovariates exited with status {status}",
            exception_string="gatk AnalyzeCovariates error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"MarkDuplicates "
            f"{rnaseq_parameter}"
            f"I={{input_filename}} O={{output_filename}} "
            f"REMOVE_DUPLICATES=true "
            f"M={self.output_basename}.no_dup_metrics{{organism_str}}.txt "
            f"CREATE_INDEX=true "
            f"TMP_DIR={config.temporary_dir}"
            f"{self.max_records_str}",
            output_format=f"{self.analysis.basename}"
            ".srt.realigned.recal.no_dup"
            f"{{organism_str}}.bam",
            error_string="Picard MarkDuplicates exited with status {status}",
            exception_string="picard MarkDuplicates error",
            split_by_organism=True,
            unlink_inputs=True,
        )

        db = Db(self.analysis.config)
        picard_metrics = PicardMetrics(db)
        executor(
            lambda *args, **kwargs: picard_metrics.store_from_file(
                self.analysis,
                f"{self.output_basename}.no_dup_metrics{{}}.txt".format(
                    kwargs["organism_str"]
                ),
                PicardMetricsType.no_duplicates,
            ),
            split_by_organism=True,
            override_last_files=False,
        )

        self.analysis.logger.info("Finished recalibration")

    def metrics_collection(self) -> None:
        """Collect some metrics from the data.

        Picard CollectHsMetrics and Picard CollectGcBiasMetrics are run
        and the output is stored into the database, if possible.
        """
        self.analysis.logger.info("Running metrics collection")
        self.chdir()
        config = self.analysis.config

        db = Db(self.analysis.config)
        picard_metrics = PicardMetrics(db)

        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"CollectHsMetrics "
            f"I={{input_filename}} BI={{kit.bait_list}} "
            f"TI={{kit.target_list}} R={{genome_ref}} "
            f"O={self.output_basename}.hs_metrics{{organism_str}}.txt "
            f"MINIMUM_MAPPING_QUALITY=0 "
            f"MINIMUM_BASE_QUALITY=0 "
            f"COVERAGE_CAP=10000 "
            f"CLIP_OVERLAPPING_READS=false "
            f"PER_BASE_COVERAGE={self.output_basename}.coverage"
            f"{{organism_str}}.txt "
            f"TMP_DIR={config.temporary_dir}"
            f"{self.max_records_str}",
            error_string="Picard CollectHsMetrics exited with status {status}",
            exception_string="picard CollectHsMetrics error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            lambda *args, **kwargs: picard_metrics.store_from_file(
                self.analysis,
                f"{self.output_basename}.hs_metrics{{}}.txt".format(
                    kwargs["organism_str"]
                ),
                PicardMetricsType.hs,
            ),
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            f"{config.java} {config.picard_jvm_args} -jar {config.picard} "
            f"CollectGcBiasMetrics "
            f"R={{genome_ref}} I={{input_filename}} "
            f"O={self.output_basename}.gcbias.metrics{{organism_str}}.txt "
            f"CHART={self.output_basename}.gcbias_metrics{{organism_str}}.pdf "
            f"S={self.output_basename}.gcbias_summ_metrics{{organism_str}}.txt "
            f"TMP_DIR={config.temporary_dir}"
            f"{self.max_records_str}",
            error_string="Picard CollectGcBiasMetrics exited with status {status}",
            exception_string="picard CollectGcBiasMetrics error",
            split_by_organism=True,
            override_last_files=False,
        )

        executor(
            lambda *args, **kwargs: picard_metrics.store_from_file(
                self.analysis,
                f"{self.output_basename}.gcbias.metrics{{}}.txt".format(
                    kwargs["organism_str"]
                ),
                PicardMetricsType.gcbias,
            ),
            split_by_organism=True,
            override_last_files=False,
        )

        self.analysis.logger.info("Finished metrics collection")

    def bam2tdf(self) -> None:
        """Convert the input BAM to TDF.

        This is pretty useful to allow manual inspection of the data
        using IGV, even if the process is time consuming.
        """
        self.analysis.logger.info("Converting BAM to TDF")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f"{config.java} -jar {config.bam2tdf} -m 10 {{input_filename}}",
            error_string="Java bam2tdf exited with status {status}",
            exception_string="bam2tdf error",
            split_by_organism=True,
            override_last_files=False,
        )

    def compress_fastq(self) -> None:
        """Compress the fastq data using gzip compression."""
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = utils.find_fastqs_by_organism(
            self.analysis.sample,
            self.fastq_dir,
            utils.get_human_annotation(self.analysis.config),
        )
        for filenames in fastq_files.values():
            for filename, _ in filenames:
                utils.gzip(filename)
        self.analysis.logger.info("Finished compressing fastq files")

    def run(self) -> None:
        """Run the mapping process.

        The steps that are performed depends on the characteristics of
        the sample and the command line parameters. The whole process
        involves reads alignment using the `Aligner` class.

        At the end of the function the current `Analysis` instance is
        has a BAM file that can be furtherly analysed with other modules
        and tools in order to obtain valuable results.
        """
        if self.analysis.parameters["skip_mapping"]:
            self.analysis.run_fake = True

        barcoded = BarcodedFilename.from_sample(self.analysis.sample)

        parsing_xenograft = (
            self.analysis.parameters["use_xenograft_classifier"]
            and barcoded.is_xenograft()
        )
        if parsing_xenograft:
            xenograft = Xenograft(self.analysis, self.fastq_dir)

            if xenograft.classifier == XenograftClassifier.XENOME:
                xenograft.run()

        if barcoded.analyte == Analyte.WHOLE_EXOME:
            if self.analysis.parameters["use_cutadapt"]:
                self.cutadapt()

        self.fastqc()
        self.trim()

        if parsing_xenograft and xenograft.classifier != XenograftClassifier.XENOME:
            xenograft.run()
            self.add_bam_groups()
        else:
            aligner = Aligner(self.analysis)
            aligner.run()

            self.convert_to_bam()
            self.add_bam_groups()
            aligner.sort_bam()

        if (
            self.analysis.parameters["mark_duplicates"]
            and barcoded.analyte != Analyte.RNASEQ
        ):
            self.mark_duplicates()

        if barcoded.analyte == Analyte.WHOLE_EXOME:
            self.indel_realign()
            self.recalibration()

        self.metrics_collection()
        if self.analysis.parameters["use_tdf"]:
            self.bam2tdf()
        if self.analysis.parameters["compress_fastq"]:
            self.compress_fastq()

        self.analysis.run_fake = False
