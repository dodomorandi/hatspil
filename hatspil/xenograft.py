from . import utils
from .executor import Executor

from formatizer import f
import os
import re
import itertools
import shutil
import gzip


class Xenograft:

    def __init__(self, analysis, fastq_dir):
        self.analysis = analysis
        self.fastq_dir = fastq_dir
        self.xenome_command = f(
            "{analysis.config.xenome} classify "
            "-T {analysis.config.xenome_threads} "
            "-P {analysis.config.xenome_index} --pairs")
        self.sample_base_out = os.path.join("REPORTS", self.analysis.sample)
        self.skip_filenames = {}

    def chdir(self):
        os.chdir(self.fastq_dir)

    def xenome(self):
        self.chdir()
        input_filenames = [
            filename
            for pair in utils.find_fastqs_by_organism(
                self.analysis.sample,
                ".").values()
            for filename, _ in pair]

        valid_filenames = []
        for filename in input_filenames:
            params = utils.get_params_from_filename(filename)
            if params[2] == 60 and not params[8]:
                valid_filenames.append(filename)
            elif params[2] != 60:
                organism = params[8]
                if not organism:
                    organism = "hg19"
                if organism not in self.skip_filenames:
                    self.skip_filenames[organism] = []
                self.skip_filenames[organism].append(os.path.join(os.getcwd(), filename))

        input_filenames = sorted(valid_filenames)

        if len(input_filenames) == 0:
            return

        self.analysis.logger.info("Running xenome")
        temp_dir = os.path.join(os.path.sep, "tmp", "%s.xenome" % self.analysis.sample)
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)

        executor = Executor(self.analysis)
        executor(f(
            "{self.xenome_command} --graft-name hg19 --host-name mm10 "
            "--output-filename-prefix {self.analysis.sample} "
            "{{input_filename}} "
            "--tmp-dir {temp_dir} "
            "> \"{self.sample_base_out}.xenome_summary.txt\""),
            input_filenames=input_filenames,
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

        self.analysis.last_operation_filenames.update(self.skip_filenames)
        self.analysis.logger.info("Finished fixing xenome fastq")

    def compress(self):
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = [
            "%s\.R%d.fastq" % (self.analysis.sample, index + 1)
            for index in range(2)]
        for filename in fastq_files:
            compressed_filename = filename + ".gz"
            with open(filename, "rb") as in_fd, \
                    gzip.open(compressed_filename, "wb") as out_fd:
                shutil.copyfileobj(in_fd, out_fd)
            os.unlink(filename)

        self.analysis.logger.info("Finished compressing fastq files")

    def run(self):
        self.xenome()
        self.fix_fastq()
        self.compress()
