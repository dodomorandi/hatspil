import argparse
import itertools
import os
import re
import shutil
import subprocess
from enum import Enum, auto
from typing import Dict, List, Optional, Tuple, Union, cast, Any
from copy import deepcopy

from . import utils
from .analysis import Analysis
from .barcoded_filename import BarcodedFilename, Analyte
from .config import Config
from .exceptions import PipelineError
from .executor import AnalysisFileData, Executor, SingleAnalysis
from .aligner import Aligner, RnaSeqAligner, GenericAligner


class XenograftClassifier(Enum):
    XENOME = auto()
    DISAMBIGUATE = auto()


class Xenome:
    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.xenome_command = \
            f"{analysis.config.xenome} classify "\
            f"-T {analysis.config.xenome_threads} "\
            f"-P {analysis.config.xenome_index} --pairs"
        reports_dir = os.path.join(self.analysis.get_bam_dir(), "REPORTS")
        os.makedirs(reports_dir, exist_ok=True)

        self.sample_base_out = os.path.join(reports_dir, self.analysis.sample)
        self.skip_filenames: Dict[str, List[str]] = {}
        self.already_done: Dict[str, List[str]] = {}
        self.input_filenames = SingleAnalysis()
        self.human_annotation = utils.get_human_annotation(analysis.config)
        self.mouse_annotation = utils.get_mouse_annotation(analysis.config)

    def chdir(self) -> None:
        os.chdir(self.fastq_dir)

    def _check_lambda(self, **kwargs: SingleAnalysis) -> None:
        self.input_filenames = SingleAnalysis([
            cast(AnalysisFileData, filename)
            for filename in utils.flatten(kwargs["input_filenames"])
        ])
        valid_filenames = SingleAnalysis()
        for file_data in self.input_filenames:
            if file_data.barcode.is_xenograft():
                valid_filenames.append(file_data)
            else:
                organism = file_data.barcode.organism
                if not organism:
                    organism = utils.get_human_annotation(self.analysis.config)
                if organism not in self.skip_filenames:
                    self.skip_filenames[organism] = []
                self.skip_filenames[organism].append(
                    os.path.join(os.getcwd(), file_data.filename))

        self.input_filenames = valid_filenames
        self.input_filenames.sort(key=lambda file_data: file_data.filename)

        compressed = [
            file_data for file_data in self.input_filenames
            if file_data.filename.endswith(".gz")
        ]

        self.input_filenames = SingleAnalysis([
            file_data for file_data in self.input_filenames
            if not file_data.filename.endswith(".gz")
        ])

        for fastqgz_data in compressed:
            can_skip = True
            removed: Dict[str, List[str]] = {}
            barcoded_filename = BarcodedFilename(fastqgz_data.filename)
            for organism in (self.human_annotation, self.mouse_annotation):
                obtained_name = "%s.%s" % (barcoded_filename.get_barcode(),
                                           organism)

                if barcoded_filename.read_index:
                    obtained_name += ".R%s" % barcoded_filename.read_index
                obtained_name += ".fastq"

                removed[organism] = [os.path.join(os.getcwd(), obtained_name)]
                obtained_file_data = next(
                    (file_data
                     for file_data in self.input_filenames
                     if file_data.filename == obtained_name),
                    None)
                if not obtained_file_data:
                    can_skip = False
                else:
                    self.input_filenames.remove(obtained_file_data)

            if can_skip:
                self.already_done.update(removed)
            else:
                if not utils.check_gz(fastqgz_data.filename):
                    raise PipelineError(
                        "%s is corrupted. Fix the problem and retry" %
                        fastqgz_data.filename)
                utils.gunzip(fastqgz_data.filename)
                self.input_filenames.append(
                    AnalysisFileData(fastqgz_data.filename[:-3]))

    def check(self) -> bool:
        self.analysis.logger.info("Checking if xenome must be performed")
        self.chdir()
        retval = True

        executor = Executor(self.analysis)
        executor(self._check_lambda, input_split_reads=False)

        if len(self.input_filenames) == 0:
            retval = False

        self.analysis.logger.info("Checking complete")
        return retval

    def xenome(self) -> None:
        self.analysis.logger.info("Running xenome")
        self.chdir()

        temp_dir = os.path.join(os.path.sep, "tmp",
                                "%s.xenome" % self.analysis.sample)
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)

        executor = Executor(self.analysis)
        executor(
            f"{self.xenome_command} "
            f"--graft-name {self.human_annotation} "
            f"--host-name {self.mouse_annotation} "
            f"--output-filename-prefix {self.analysis.sample} "
            f"{{input_filename}} "
            f"--tmp-dir {temp_dir} "
            f"> \"{self.sample_base_out}.xenome_summary.txt\"",
            input_filenames=[
                file_data.filename for file_data in self.input_filenames],
            input_function=lambda filenames: " ".join(
                ["-i %s" % filename for filename in filenames]),
            output_format=f"{self.analysis.sample}_%s_%d.fastq",
            output_function=lambda filename: [
                filename % (organism, index) for organism, index in itertools.
                product([self.human_annotation, self.mouse_annotation], [1, 2])
            ],
            input_split_reads=False)
        shutil.rmtree(temp_dir)

        if self.analysis.run_fake:
            if self.analysis.last_operation_filenames is None:
                pass
            elif isinstance(self.analysis.last_operation_filenames, str):
                with open(self.analysis.last_operation_filenames, "a"):
                    pass
            elif isinstance(self.analysis.last_operation_filenames, list):
                for filename in self.analysis.last_operation_filenames:
                    with open(filename, "a"):
                        pass
            else:
                for filenames in self.analysis\
                        .last_operation_filenames.values():
                    for filename in filenames:
                        with open(filename, "a"):
                            pass

        self.analysis.logger.info("Finished xenome")

    def fix_fastq(self) -> None:
        match = re.match(R"^[^-]+-[^-]+-(\d)", self.analysis.sample)
        if not match or match.group(1) != "6":
            self.analysis.last_operation_filenames = self.skip_filenames
            return

        self.analysis.logger.info("Fixing xenome fastq")
        self.chdir()
        self.analysis.last_operation_filenames = {}
        re_fastq_filename = re.compile(
            R"^%s_((?:hg|mm)\d+)_([12])\.fastq$" % self.analysis.sample, re.I)
        fastq_files = [
            filename for filename in os.listdir()
            if re_fastq_filename.match(filename)
        ]
        for filename in fastq_files:
            match = re_fastq_filename.match(filename)
            assert match is not None

            organism = match.group(1)
            read_index = int(match.group(2))
            out_filename = "%s.%s.R%d.fastq" % (self.analysis.sample, organism,
                                                read_index)

            if not self.analysis.run_fake:
                with open(filename, "r") as in_fd,\
                        open(out_filename, "w") as out_fd:
                    for line_index, line in enumerate(in_fd):
                        if line_index % 4 == 0:
                            line = line.strip()
                            space_pos = line.find(" ")
                            if space_pos != -1:
                                space_pos = line.find(" ", space_pos + 1)

                            out_fd.write("@%s\n" % line[:space_pos])
                        elif line_index % 4 == 2:
                            out_fd.write("+\n")
                        else:
                            out_fd.write("%s\n" % line.strip())

            if organism not in self.analysis.last_operation_filenames:
                self.analysis.last_operation_filenames[organism] = []
            self.analysis.last_operation_filenames[organism].append(
                os.path.join(os.getcwd(), out_filename))

        other_fastq = [
            "%s_%s_%d.fastq" % cast(Tuple[str, str, int],
                                    ((self.analysis.sample, ) + combo))
            for combo in itertools.product(["ambiguous", "both", "neither"],
                                           range(1, 3))
        ]
        for filename in fastq_files + other_fastq:
            if os.path.exists(filename):
                os.unlink(filename)

        self.analysis.logger.info("Finished fixing xenome fastq")

    def compress(self) -> None:
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = [
            "%s.R%d.fastq" % (self.analysis.sample, index + 1)
            for index in range(2)
        ]
        for filename in fastq_files:
            utils.gzip(filename)

        self.analysis.logger.info("Finished compressing fastq files")

    def update_last_filenames(self) -> None:
        if self.analysis.last_operation_filenames is None:
            self.analysis.last_operation_filenames = {}

        assert isinstance(self.analysis.last_operation_filenames, dict)
        self.analysis.last_operation_filenames.update(self.skip_filenames)
        self.analysis.last_operation_filenames.update(self.already_done)

    def cannot_unlink_results(self) -> None:
        self.analysis.can_unlink = False

    def run(self) -> None:
        if self.check():
            self.xenome()
            self.fix_fastq()
            self.compress()
        self.update_last_filenames()
        self.cannot_unlink_results()


class Disambiguate:
    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.aligned_fastq_files: Union[str,
                                        List[str],
                                        Dict[str, List[str]],
                                        None] = None
        self.tempdir = os.path.join(analysis.get_bam_dir(),
                                    "disambiguate_{}".format(analysis.sample))

    def create_avatar_links(self) -> None:
        self.analysis.logger.info("Creating links for avatar organism")

        def get_output_format(*args: AnalysisFileData,
                              **kwargs: AnalysisFileData) -> str:
            input_filename = kwargs["input_filename"].filename
            assert(input_filename)
            dirname = os.path.dirname(input_filename)
            barcode = BarcodedFilename(filename=input_filename)
            barcode.organism = utils.get_mouse_annotation(
                self.analysis.config)

            barcoded_filename = barcode.get_barcoded_filename()
            if not barcoded_filename:
                raise PipelineError(f"'{input_filename}' is an invalid "
                                    "filename that cannot be re-barcoded")

            return os.path.join(dirname, barcoded_filename)

        def do_symlinks(*args: Any, **kwargs: Any) -> None:
            input_filename: str = kwargs["input_filename"].filename
            output_filenames: List[str] = kwargs["output_filename"]
            assert(len(output_filenames) == 2)

            os.symlink(input_filename, output_filenames[1])

        executor = Executor(self.analysis)
        executor(do_symlinks,
                 input_split_reads=True,
                 only_human=True,
                 split_by_organism=True,
                 output_format=["{input_filename}", get_output_format])
        self.aligned_fastq_files = deepcopy(
            self.analysis.last_operation_filenames)
        self.analysis.can_unlink = False

    def unlink_symlinks(self) -> None:
        assert(self.aligned_fastq_files)
        self.analysis.logger.info("Removing fastq symlinks")
        last_operation_filenames = self.analysis.last_operation_filenames
        self.analysis.last_operation_filenames = self.aligned_fastq_files

        def unlink_files(*args: Any, **kwargs: Any) -> None:
            input_filename: AnalysisFileData = kwargs["input_filename"]
            assert(input_filename.barcode)
            organism = input_filename.barcode.organism
            if organism and organism.startswith("mm"):
                os.unlink(input_filename.filename)

        executor = Executor(self.analysis)
        executor(unlink_files,
                 split_by_organism=True,
                 input_split_reads=True,
                 only_human=False,
                 override_last_files=False)

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
                raise PipelineError("cannot use aligner {} in combination "
                                    "with disambiguate"
                                    .format(aligner.name.lower()))
        else:
            aligner = self.analysis.parameters["aligner"]
            if aligner == GenericAligner.BWA:
                aligner_str = "bwa"
            else:
                raise PipelineError("cannot use aligner {} in combination "
                                    "with disambiguate"
                                    .format(aligner.name.lower()))

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

        executor(get_human_bam_index,
                 split_by_organism=False,
                 split_input_files=False,
                 override_last_files=False)

        assert(human_bam_index is not None)

        executor(f"{config.disambiguate} "
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
                     ".disambiguatedSpeciesB.bam"
                 ],
                 input_function=lambda filenames: " ".join(filenames),
                 unlink_inputs=True)

        bam_dir = self.analysis.get_bam_dir()
        out_prefix = os.path.join(bam_dir, self.analysis.basename)
        executor(lambda *args, **kwargs:
                 os.rename(kwargs["input_filename"].filename,
                           kwargs["output_filename"]),
                 split_by_organism=False,
                 split_input_files=False,
                 input_function=lambda filenames:
                 cast(List[str], filenames)[
                     cast(int, human_bam_index)],
                 output_format=f"{out_prefix}.disambiguated"
                 "{organism_str}.bam")

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
    CLASSIFIERS = [XenograftClassifier.DISAMBIGUATE,
                   XenograftClassifier.XENOME]

    def __init__(self, analysis: Analysis, fastq_dir: str) -> None:
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.classifier = analysis.parameters["xenograft_classifier"]
        assert(self.classifier is not None)

    @staticmethod
    def get_available_classifier(cmdline_args: argparse.Namespace,
                                 config: Config)\
            -> Optional[XenograftClassifier]:

        classifier_name = cmdline_args.xenograft_classifier
        if classifier_name == "auto" or classifier_name is None:
            for classifier in Xenograft.CLASSIFIERS:
                classifier_exe = getattr(config, classifier.name.lower())
                if classifier_exe is not None \
                        and shutil.which(classifier_exe) is not None:
                    return classifier

            return None
        else:
            classifier_names = [classifier.name.lower()
                                for classifier in Xenograft.CLASSIFIERS]
            classifier_name_lower = classifier_name.lower()
            if classifier_name_lower not in classifier_names:
                return None

            classifier_exe = getattr(config, classifier_name_lower)
            if classifier_exe is not None \
                    and shutil.which(classifier_exe) is not None:
                classifier_index = classifier_names.index(
                    classifier_name_lower)
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
