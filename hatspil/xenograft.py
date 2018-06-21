import itertools
import os
import re
import shutil
from typing import Dict, List, Tuple, cast

from . import utils
from .analysis import Analysis
from .barcoded_filename import BarcodedFilename, Tissue
from .exceptions import PipelineError
from .executor import AnalysisFileData, Executor, SingleAnalysis


class Xenograft:
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
                    organism = "hg19"
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
            for organism in ("hg19", "mm10"):
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
            f"{self.xenome_command} --graft-name hg19 --host-name mm10 "
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
                product(["hg19", "mm10"], [1, 2])
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
