"""The module to handle xenograft NGS data.

Xenograft tissues result in mixed information of host and graft. There
are a few tools to handle these cases, and their approach are pretty
different.

The purpose of this module is to abstract away some parts and ease the
integration with other modules.
"""
import argparse
import itertools
import os
import re
import shutil
import sys
from copy import deepcopy
from enum import Enum, auto
from typing import Any, Dict, Iterable, List, Optional, Union, cast

from .core import utils
from .aligner import Aligner, GenericAligner, RnaSeqAligner
from .analysis import Analysis
from .barcoded_filename import Analyte, BarcodedFilename
from .config import Config
from .exceptions import PipelineError
from .executor import AnalysisFileData, Executor, SingleAnalysis


class XenograftClassifier(Enum):
    """The Xenograft classifier software."""

    XENOME = auto()
    DISAMBIGUATE = auto()


class SampleFileAvailability:
    """A class to track the status of the files for a Xenograft."""

    def __init__(self) -> None:
        """Create a default instance of the class."""
        self.compressed = False
        self.uncompressed = False
        self.host = False
        self.graft = False

    def __repr__(self) -> str:
        """Return a human readable format of the data."""
        out = []
        if self.compressed:
            out.append("compressed")
        if self.uncompressed:
            out.append("uncompressed")
        if self.host:
            out.append("host")
        if self.graft:
            out.append("graft")

        return "[{}]".format(" ".join(out))


class XenomePreChecker:
    """A helper class to avoid unnecessary data processing.

    Performing Xenograft data analysis is time consuming. This class
    helps avoiding unnecessary steps, analysing the current availability
    and populating the `availability_for_sample_read` dict.
    """

    def __init__(
        self, analysis: Analysis, human_annotation: str, mouse_annotation: str
    ) -> None:
        """Create an instance of the class.

        The `availability_for_sample_read` is left empty. In order to
        populate it, the instance must be called using an executor.
        """
        self.analysis = analysis
        self.human_annotation = human_annotation
        self.mouse_annotation = mouse_annotation
        self.availability_for_sample_read: Dict[int, SampleFileAvailability] = {}

    def __call__(self, **kwargs: Any) -> None:
        """Populate the data analysing the available files.

        This function should be called using an executor.
        """
        input_filenames: SingleAnalysis = kwargs["input_filenames"]

        for filename in input_filenames:
            barcoded = filename.barcode
            assert barcoded.read_index is not None

            availability = self.availability_for_sample_read.setdefault(
                barcoded.read_index, SampleFileAvailability()
            )
            if barcoded.gzipped:
                assert not availability.uncompressed
                availability.compressed = True
            elif barcoded.organism:
                if barcoded.organism == self.human_annotation:
                    availability.graft = True
                elif barcoded.organism == self.mouse_annotation:
                    availability.host = True
                else:
                    raise AssertionError("invalid organism")
            else:
                assert not availability.compressed
                availability.uncompressed = True


class Xenome:
    """A class to handle the Xenome software.

    Xenome is able to analyse the original FASTQ files and to produce
    new ones with data split for host and graft.
    """

    RE_XENOME_OUTPUT_FORMAT = re.compile(
        r"^([^-]+-[^-]+-\d[0-9A-Za-z]-\d{3}-\d{3}(?:\.[^.]+)*?)_(\w{2}\d{1,2})_(\d)"
        r"\.fastq"
    )

    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        """Create the instance of the class.

        The BAM REPORTS directory is created if needed.
        """
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.xenome_command = "{} classify -T {} -P {}".format(
            analysis.config.xenome,
            analysis.config.xenome_threads,
            analysis.config.xenome_index,
        )
        reports_dir = os.path.join(self.analysis.get_bam_dir(), "REPORTS")
        os.makedirs(reports_dir, exist_ok=True)

        self.sample_base_out = os.path.join(reports_dir, self.analysis.sample)
        self.already_done: Dict[str, List[str]] = {}
        self.human_annotation = utils.get_human_annotation(analysis.config)
        self.mouse_annotation = utils.get_mouse_annotation(analysis.config)
        self.availability_for_sample_read: Dict[int, SampleFileAvailability] = {}
        self.analyzed_files: List[str] = []

    @staticmethod
    def check_correct_index_files(config: Config) -> None:
        """Check for needed files.

        Xenome needs some auxiliary files that must be created before
        running HaTSPiL. If some of these are missing, an error is
        printed and the process exists.
        """
        index_base = config.xenome_index

        non_existant_files: List[str] = []
        invalid_files: List[str] = []

        for suffix in ("-both", "-graft", "-host"):
            complete_filename = "{}{}.kmers.header".format(index_base, suffix)

            if not os.path.exists(complete_filename):
                non_existant_files.append(complete_filename)
                continue

            try:
                with open(complete_filename, "rb") as fd:
                    date_header = fd.read(4)
                    if date_header != b"\x25\x26\xed\x77":
                        invalid_files.append(complete_filename)
            except Exception:
                invalid_files.append(complete_filename)

        if not non_existant_files and not invalid_files:
            return

        if non_existant_files:
            print(
                "Not all Xenome index files are available.\n"
                "The following are missing:\n{}".format("\n".join(non_existant_files)),
                file=sys.stderr,
            )

        if invalid_files:
            print(
                "Some Xenome index files cannot be read or Xenome and index files are "
                "obsolete:\n{}\nUpdate Xenome and rebuild indices.".format(
                    "\n".join(invalid_files)
                ),
                file=sys.stderr,
            )
        exit(-1)

    def chdir(self) -> None:
        """Change current directory to the FASTQ folder."""
        os.chdir(self.fastq_dir)

    def check(self) -> None:
        """Check if Xenome must be run for the input file.

        It uses the `XenomePreChecker` with the `Executor`, and it
        updates the `availability_for_sample_read` according to the
        pre-checker.
        """
        self.analysis.logger.info("Checking if xenome must be performed")
        self.chdir()

        pre_checker = XenomePreChecker(
            self.analysis, self.human_annotation, self.mouse_annotation
        )

        executor = Executor(self.analysis)
        executor(pre_checker, input_split_reads=False, override_last_files=False)
        self.availability_for_sample_read = pre_checker.availability_for_sample_read

        self.analysis.logger.info("Checking complete")

    def xenome_must_run(self) -> bool:
        """Evaluate if Xenome must be run.

        It uses the content of `availability_for_sample_read` to
        determine if Xenome must be run or if all the files for host and
        graft are already available.
        """
        return bool(self.availability_for_sample_read) and any(
            filter(
                lambda availability: not cast(SampleFileAvailability, availability).host
                or not cast(SampleFileAvailability, availability).graft,
                self.availability_for_sample_read.values(),
            )
        )

    def decompress(self) -> None:
        """Decompress original file, if needed.

        The original files are generally compressed after Xenome is run.
        Therefore, if the analysis must be run again, they must be
        decompressed before being utilised.

        If `self.xenome_must_run()` is False no decompression is
        performed.
        """
        if not self.xenome_must_run():
            return

        if not any(
            filter(
                lambda availability: not cast(
                    SampleFileAvailability, availability
                ).uncompressed,
                self.availability_for_sample_read.values(),
            )
        ):
            return

        self.analysis.logger.info("Decompressing original files")
        self.chdir()

        def decompress_input_function(
            input_filename: Union[str, List[str]]
        ) -> Optional[str]:
            if not isinstance(input_filename, str):
                return None

            barcoded = BarcodedFilename(input_filename)
            if not barcoded.organism:
                return input_filename
            else:
                return None

        def decompress_output_function(input_filename: str) -> Iterable[str]:
            barcoded = BarcodedFilename(input_filename)
            if barcoded.gzipped:
                return [input_filename[:-3]]
            else:
                return [input_filename]

        def decompress_if_needed(**kwargs: SingleAnalysis) -> None:
            input_analyses = kwargs["input_filenames"]
            for input_analysis in input_analyses:
                input_filename = input_analysis.filename
                barcoded = BarcodedFilename(input_filename)
                if barcoded.gzipped:
                    try:
                        utils.gunzip(input_filename)
                    except Exception:
                        print(
                            "ERROR: cannot decompress file '{}'. "
                            "Please remove this file and try again.".format(
                                input_filename
                            ),
                            file=sys.stderr,
                        )
                        exit(-1)

        executor = Executor(self.analysis)
        executor(
            decompress_if_needed,
            input_function=decompress_input_function,
            output_function=decompress_output_function,
            split_by_organism=True,
            input_split_reads=True,
        )

        self.analysis.logger.info("Finished decompressing original files")

    def remove_gzipped_from_files(self) -> None:
        """Remove gzip files from the Analysis input files.

        If it is not necessary to run Xenome, the original fastq.gz
        files must not be used in the next step of the analysis. This
        helper function removes the original files from the `Analysis`.
        """
        operations_filenames: Dict[str, List[str]] = {}
        assert isinstance(self.analysis.last_operation_filenames, dict)
        for organism, filenames in self.analysis.last_operation_filenames.items():
            assert filenames
            organism_filenames = operations_filenames.setdefault(organism, [])
            for filename in filenames:
                try:
                    barcoded = BarcodedFilename(filename)
                except Exception:
                    continue

                if barcoded and not barcoded.gzipped:
                    organism_filenames.append(filename)
        self.analysis.last_operation_filenames = operations_filenames

    def xenome(self) -> None:
        """Run Xenome, if needed.

        Run Xenome, removing all the ambiguous and uninformative
        results. The parameters are sets automatically according to the
        type of type of NGS analysis (if single- or paired-end).

        In case all the output file are already available, Xenome is not
        run.
        """
        if not self.xenome_must_run():
            self.remove_gzipped_from_files()
            return

        self.analysis.logger.info("Running xenome")
        self.chdir()

        assert len(self.availability_for_sample_read) <= 2
        if len(self.availability_for_sample_read.keys()) == 2:
            pairing_parameter = "--pairs "
        else:
            pairing_parameter = ""

        temp_dir = os.path.join(os.path.sep, "tmp", "%s.xenome" % self.analysis.sample)
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir, exist_ok=True)

        def xenome_filtered_files(filenames: List[str]) -> List[str]:
            barcoded_filenames = [BarcodedFilename(filename) for filename in filenames]
            assert all(barcoded_filenames)
            return [
                filename
                for (filename, barcoded) in zip(filenames, barcoded_filenames)
                if not barcoded.gzipped and not barcoded.organism
            ]

        def xenome_input_function(filenames: Union[str, List[str]]) -> Optional[str]:
            assert isinstance(filenames, list)
            return " ".join(
                ["-i %s" % filename for filename in xenome_filtered_files(filenames)]
            )

        assert isinstance(self.analysis.last_operation_filenames, dict)
        for filenames in self.analysis.last_operation_filenames.values():
            self.analyzed_files += xenome_filtered_files(filenames)

        executor = Executor(self.analysis)
        executor(
            f"{self.xenome_command} {pairing_parameter}"
            f"--graft-name {self.human_annotation} "
            f"--host-name {self.mouse_annotation} "
            f"--output-filename-prefix {self.analysis.sample} "
            f"{{input_filename}} "
            f"--tmp-dir {temp_dir} "
            f'> "{self.sample_base_out}.xenome_summary.txt"',
            input_function=xenome_input_function,
            output_format=f"{self.analysis.sample}_{{}}_{{}}.fastq",
            output_function=lambda filename: [
                filename.format(organism, index)
                for organism, index in itertools.product(
                    [
                        self.human_annotation,
                        self.mouse_annotation,
                        "ambiguous",
                        "both",
                        "neither",
                    ],
                    [1, 2],
                )
            ],
            input_split_reads=False,
        )

        executor(
            lambda **kwargs: shutil.rmtree(temp_dir),
            output_format="",
            input_split_reads=False,
            override_last_files=False,
            allow_raw_filenames=True,
        )

        def get_fixed_filename(filename: str) -> str:
            basename = os.path.basename(filename)
            dirname = os.path.dirname(filename)
            new_basename = Xenome.RE_XENOME_OUTPUT_FORMAT.sub(
                r"\1.\2.R\3.fastq", basename
            )
            return os.path.join(dirname, new_basename)

        def rename_or_delete(**kwargs: SingleAnalysis) -> None:
            input_analyses = kwargs["input_filenames"]
            for input_analysis in input_analyses:
                input_filename = input_analysis.filename
                if (
                    "ambiguous" in input_filename
                    or "both" in input_filename
                    or "neither" in input_filename
                ):
                    os.unlink(input_filename)
                else:
                    os.rename(input_filename, get_fixed_filename(input_filename))

        def get_only_unambiguous(input_filename: str) -> Iterable[str]:
            if (
                "ambiguous" in input_filename
                or "both" in input_filename
                or "neither" in input_filename
            ):
                return []
            else:
                return [get_fixed_filename(input_filename)]

        executor(
            rename_or_delete,
            output_function=get_only_unambiguous,
            input_split_reads=False,
            split_by_organism=True,
            allow_raw_filenames=True,
        )

        self.analysis.logger.info("Finished xenome")

    def compress(self) -> None:
        """Compress the initial file.

        If Xenome is performed, the initial files can be compressed in
        order to save some disk space.
        """
        if self.analysis.run_fake or not self.analyzed_files:
            return

        self.analysis.logger.info("Compressing fastq files")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            lambda **kwargs: utils.gzip(kwargs["input_filename"].filename),
            input_filenames=self.analyzed_files,
            input_split_reads=True,
            override_last_files=False,
        )

        self.analysis.logger.info("Finished compressing fastq files")

    def cannot_unlink_results(self) -> None:
        """Set that the last files cannot be unlink."""
        self.analysis.can_unlink = False

    def run(self) -> None:
        """Perform all the steps to correctly handle Xenome.
        
        A preliminary check is performed in order to evaluate if Xenome
        must be run. In this case the run is performed and the
        original files are compressed. Moreover, the resulting FASTQ
        files are set to be not deletable.
        """
        self.check()
        self.decompress()
        self.xenome()
        self.compress()
        self.cannot_unlink_results()


class Disambiguate:
    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.aligned_fastq_files: Union[
            str, List[str], Dict[str, List[str]], None
        ] = None
        self.tempdir = os.path.join(
            analysis.get_bam_dir(), "disambiguate_{}".format(analysis.sample)
        )

    def create_avatar_links(self) -> None:
        self.analysis.logger.info("Creating links for avatar organism")

        def get_output_format(
            *args: AnalysisFileData, **kwargs: AnalysisFileData
        ) -> str:
            input_filename = kwargs["input_filename"].filename
            assert input_filename
            dirname = os.path.dirname(input_filename)
            barcode = BarcodedFilename(filename=input_filename)
            barcode.organism = utils.get_mouse_annotation(self.analysis.config)

            barcoded_filename = barcode.get_barcoded_filename()
            if not barcoded_filename:
                raise PipelineError(
                    f"'{input_filename}' is an invalid "
                    "filename that cannot be re-barcoded"
                )

            return os.path.join(dirname, barcoded_filename)

        def do_symlinks(*args: Any, **kwargs: Any) -> None:
            input_filename: str = kwargs["input_filename"].filename
            output_filenames: List[str] = kwargs["output_filename"]
            assert len(output_filenames) == 2

            os.symlink(input_filename, output_filenames[1])

        executor = Executor(self.analysis)
        executor(
            do_symlinks,
            input_split_reads=True,
            only_human=True,
            split_by_organism=True,
            output_format=["{input_filename}", get_output_format],
        )
        self.aligned_fastq_files = deepcopy(self.analysis.last_operation_filenames)
        self.analysis.can_unlink = False

    def unlink_symlinks(self) -> None:
        assert self.aligned_fastq_files
        self.analysis.logger.info("Removing fastq symlinks")
        last_operation_filenames = self.analysis.last_operation_filenames
        self.analysis.last_operation_filenames = self.aligned_fastq_files

        def unlink_files(*args: Any, **kwargs: Any) -> None:
            input_filename: AnalysisFileData = kwargs["input_filename"]
            assert input_filename.barcode
            organism = input_filename.barcode.organism
            if organism and organism.startswith("mm"):
                os.unlink(input_filename.filename)

        executor = Executor(self.analysis)
        executor(
            unlink_files,
            split_by_organism=True,
            input_split_reads=True,
            only_human=False,
            override_last_files=False,
        )

        self.analysis.last_operation_filenames = last_operation_filenames

    def disambiguate(self) -> None:
        self.analysis.logger.info("Running disambiguate")

        analysis = self.analysis
        config = analysis.config

        barcoded = BarcodedFilename.from_sample(self.analysis.sample)
        if barcoded.analyte == Analyte.RNASEQ:
            aligner = self.analysis.parameters["rnaseq_aligner"]
            if aligner == RnaSeqAligner.STAR:
                aligner_str = "star"
            else:
                raise PipelineError(
                    "cannot use aligner {} in combination "
                    "with disambiguate".format(aligner.name.lower())
                )
        else:
            aligner = self.analysis.parameters["aligner"]
            if aligner == GenericAligner.BWA:
                aligner_str = "bwa"
            else:
                raise PipelineError(
                    "cannot use aligner {} in combination "
                    "with disambiguate".format(aligner.name.lower())
                )

        executor = Executor(self.analysis)

        human_bam_index: Optional[int] = None

        def get_human_bam_index(*args: Any, **kwargs: Any) -> None:
            nonlocal human_bam_index

            analysis: SingleAnalysis = kwargs["input_filenames"]
            for index, analysis_file_data in enumerate(analysis):
                filename = analysis_file_data.filename
                try:
                    barcoded = BarcodedFilename(filename)
                except Exception:
                    continue

                organism = barcoded.organism
                if not organism or organism.startswith("hg"):
                    human_bam_index = index
                    return

        run_fake = self.analysis.run_fake
        self.analysis.run_fake = False
        executor(
            get_human_bam_index,
            split_by_organism=False,
            split_input_files=False,
            override_last_files=False,
        )
        self.analysis.run_fake = run_fake

        assert human_bam_index is not None

        executor(
            f"{config.disambiguate} "
            f"-s {analysis.sample} "
            f"-o {self.tempdir} "
            f"-a {aligner_str} "
            "{input_filename}",
            split_by_organism=False,
            split_input_files=False,
            output_format=[
                f"{os.path.join(self.tempdir, self.analysis.sample)}"
                ".disambiguatedSpeciesA.bam",
                f"{os.path.join(self.tempdir, self.analysis.sample)}"
                ".disambiguatedSpeciesB.bam",
            ],
            input_function=lambda filenames: " ".join(filenames),
            unlink_inputs=True,
        )

        bam_dir = self.analysis.get_bam_dir()
        out_prefix = os.path.join(bam_dir, self.analysis.basename)
        executor(
            lambda *args, **kwargs: os.rename(
                kwargs["input_filename"].filename, kwargs["output_filename"]
            ),
            split_by_organism=False,
            split_input_files=False,
            input_function=lambda filenames: cast(List[str], filenames)[
                cast(int, human_bam_index)
            ],
            output_format=f"{out_prefix}.disambiguated{{organism_str}}.bam",
        )

        if not self.analysis.run_fake:
            shutil.rmtree(self.tempdir)

        self.analysis.logger.info("Finished running disambiguate")

    def run(self) -> None:
        self.create_avatar_links()

        aligner = Aligner(self.analysis)
        aligner.only_human = False
        aligner.run()

        self.unlink_symlinks()

        aligner.sort_bam()

        self.disambiguate()


class Xenograft:
    CLASSIFIERS = [XenograftClassifier.DISAMBIGUATE, XenograftClassifier.XENOME]

    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.classifier = analysis.parameters["xenograft_classifier"]
        assert self.classifier is not None

    @staticmethod
    def get_available_classifier(
        cmdline_args: argparse.Namespace, config: Config
    ) -> Optional[XenograftClassifier]:

        classifier_name = cmdline_args.xenograft_classifier
        if classifier_name == "auto" or classifier_name is None:
            for classifier in Xenograft.CLASSIFIERS:
                classifier_exe = getattr(config, classifier.name.lower())
                if (
                    classifier_exe is not None
                    and shutil.which(classifier_exe) is not None
                ):
                    return classifier

            return None
        else:
            classifier_names = [
                classifier.name.lower() for classifier in Xenograft.CLASSIFIERS
            ]
            classifier_name_lower = classifier_name.lower()
            if classifier_name_lower not in classifier_names:
                return None

            classifier_exe = getattr(config, classifier_name_lower)
            if classifier_exe is not None and shutil.which(classifier_exe) is not None:
                classifier_index = classifier_names.index(classifier_name_lower)
                return Xenograft.CLASSIFIERS[classifier_index]
            else:
                return None

    def run(self) -> None:
        classifier: Union[Xenome, Disambiguate]

        if self.classifier == XenograftClassifier.XENOME:
            classifier = Xenome(self.analysis, self.fastq_dir)
        elif self.classifier == XenograftClassifier.DISAMBIGUATE:
            classifier = Disambiguate(self.analysis, self.fastq_dir)
        else:
            raise PipelineError("unexpected xenograft classifier")

        classifier.run()
