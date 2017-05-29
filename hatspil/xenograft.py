from . import utils
from .executor import Executor
from .barcoded_filename import BarcodedFilename, Tissue
from .exceptions import PipelineError

from formatizer import f
import os
import re
import itertools
import shutil


class Xenograft:

    def __init__(self, analysis, fastq_dir):
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.xenome_command = f(
            "{analysis.config.xenome} classify "
            "-T {analysis.config.xenome_threads} "
            "-P {analysis.config.xenome_index} --pairs")
        reports_dir = os.path.join(self.fastq_dir, "REPORTS")
        try:
            os.makedirs(reports_dir, exist_ok=True)
        except OSError:
            pass

        self.sample_base_out = os.path.join(reports_dir, self.analysis.sample)
        self.skip_filenames = {}
        self.already_done = {}
        self.input_filenames = None

    def chdir(self):
        os.chdir(self.fastq_dir)

    def _check_lambda(self, **kwargs):
        self.input_filenames = [filename for filename
                                in utils.flatten(kwargs["input_filenames"]
                                                 .values())]
        valid_filenames = []
        for filename in self.input_filenames:
            barcoded_filename = BarcodedFilename(filename)
            if barcoded_filename.tissue == Tissue.PRIMARY_XENOGRAFT_TISSUE or \
                    barcoded_filename.tissue == Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE:
                valid_filenames.append(filename)
            else:
                organism = barcoded_filename.organism
                if not organism:
                    organism = "hg19"
                if organism not in self.skip_filenames:
                    self.skip_filenames[organism] = []
                self.skip_filenames[organism].append(os.path.join(os.getcwd(),
                                                                  filename))

        self.input_filenames = sorted(valid_filenames)

        compressed = [filename for filename
                      in self.input_filenames
                      if filename.endswith(".gz")]

        self.input_filenames = [filename for filename
                                in self.input_filenames
                                if not filename.endswith(".gz")]

        for fastqgz in compressed:
            can_skip = True
            removed = {}
            barcoded_filename = BarcodedFilename(fastqgz)
            for organism in ("hg19", "mm10"):
                obtained_name = "%s-%s-%d-%d%d%d-%d%d%d.%s" % (
                    barcoded_filename.project,
                    barcoded_filename.patient,
                    barcoded_filename.tissue,
                    barcoded_filename.molecule,
                    barcoded_filename.analyte,
                    barcoded_filename.kit,
                    barcoded_filename.biopsy,
                    barcoded_filename.sample,
                    barcoded_filename.sequencing,
                    organism
                )
                if barcoded_filename.read_index:
                    obtained_name += ".R%s" % barcoded_filename.read_index
                obtained_name += ".fastq"

                removed[organism] = [os.path.join(os.getcwd(), obtained_name)]
                if obtained_name not in self.input_filenames:
                    can_skip = False
                else:
                    self.input_filenames.remove(obtained_name)

            if can_skip:
                self.already_done.update(removed)
            else:
                if not utils.check_gz(fastqgz):
                    raise PipelineError(
                        "%s is corrupted. Fix the problem and retry" % fastqgz)
                utils.gunzip(fastqgz)
                self.input_filenames.append(fastqgz[:-3])

    def check(self):
        self.analysis.logger.info("Checking if xenome must be performed")
        self.chdir()
        retval = True

        executor = Executor(self.analysis)
        executor(self._check_lambda,
                 input_split_reads=False)

        if len(self.input_filenames) == 0:
            retval = False

        self.analysis.logger.info("Checking complete")
        return retval

    def xenome(self):
        self.analysis.logger.info("Running xenome")
        self.chdir()

        temp_dir = os.path.join(os.path.sep,
                                "tmp",
                                "%s.xenome" % self.analysis.sample)
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)

        executor = Executor(self.analysis)
        executor(f(
            "{self.xenome_command} --graft-name hg19 --host-name mm10 "
            "--output-filename-prefix {self.analysis.sample} "
            "{{input_filename}} "
            "--tmp-dir {temp_dir} "
            "> \"{self.sample_base_out}.xenome_summary.txt\""),
            input_filenames=self.input_filenames,
            input_function=lambda filenames: " ".join(["-i %s" % filename for filename in filenames]),
            output_format=f("{self.analysis.sample}_%s_%d.fastq"),
            output_function=lambda filename: [filename % (organism, index)
                                              for organism, index
                                              in itertools.product(["hg19", "mm10"],
                                                                   [1, 2])],
            input_split_reads=False
        )
        shutil.rmtree(temp_dir)

        self.analysis.logger.info("Finished xenome")

    def fix_fastq(self):
        match = re.match(R"^[^-]+-[^-]+-(\d{2})", self.analysis.sample)
        if not match or match.group(1) != "60":
            self.analysis.last_operation_filenames = self.skip_filenames
            return

        self.analysis.logger.info("Fixing xenome fastq")
        self.chdir()
        self.analysis.last_operation_filenames = {}
        re_fastq_filename = re.compile(R"^%s_((?:hg|mm)\d+)_([12])\.fastq$" % self.analysis.sample, re.I)
        fastq_files = [filename for filename in os.listdir() if re_fastq_filename.match(filename)]
        for filename in fastq_files:
            match = re_fastq_filename.match(filename)
            organism = match.group(1)
            read_index = int(match.group(2))
            out_filename = "%s.%s.R%d.fastq" % (self.analysis.sample, organism, read_index)
            with open(filename, "r") as in_fd,\
                    open(out_filename, "w") as out_fd:
                for line_index, line in enumerate(in_fd):
                    if line_index % 4 == 0:
                        splitted = line.strip().split(" ")
                        out_fd.write("@%s %s\n" % (splitted[0], splitted[1]))
                    elif line_index % 4 == 2:
                        out_fd.write("+\n")
                    else:
                        out_fd.write("%s\n" % line.strip())

            if not organism in self.analysis.last_operation_filenames:
                self.analysis.last_operation_filenames[organism] = []
            self.analysis.last_operation_filenames[organism].append(
                os.path.join(
                    os.getcwd(),
                    out_filename))

        other_fastq = [
            "%s_%s_%d.fastq" % tuple([self.analysis.sample] +
                                     list(combo))
            for combo in itertools.product(["ambiguous", "both", "neither"],
                                           range(1, 3))]
        for filename in fastq_files + other_fastq:
            if os.path.exists(filename):
                os.unlink(filename)

        self.analysis.logger.info("Finished fixing xenome fastq")

    def compress(self):
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = [
            "%s.R%d.fastq" % (self.analysis.sample, index + 1)
            for index in range(2)]
        for filename in fastq_files:
            utils.gzip(filename)

        self.analysis.logger.info("Finished compressing fastq files")

    def update_last_filenames(self):
        if self.analysis.last_operation_filenames is None:
            self.analysis.last_operation_filenames = {}
        self.analysis.last_operation_filenames.update(self.skip_filenames)
        self.analysis.last_operation_filenames.update(self.already_done)

    def run(self):
        if self.check():
            self.xenome()
            self.fix_fastq()
            self.compress()
        self.update_last_filenames()
