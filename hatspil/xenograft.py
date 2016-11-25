from . import utils
from .exceptions import PipelineError

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

    def chdir(self):
        os.chdir(self.fastq_dir)

    def xenome(self):
        self.analysis.logger.info("Running xenome")
        self.chdir()
        retval = utils.run_and_log(f(
            "{self.xenome_command} --graft-name hg19 --host-name mm10 "
            "--output-filename-prefix {self.analysis.sample} "
            "-i \"{self.analysis.sample}_R1.fastq\" "
            "-i \"{self.analysis.sample}_R2.fastq\" "
            "> \"{self.sample_base_out}.xenome_summary.txt\""),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Xenome exited with status %d" % retval)
            raise PipelineError("xenome error")

        self.analysis.logger.info("Finished xenome")

    def fix_fastq(self):
        self.analysis.logger.info("Fixing xenome fastq")
        self.chdir()
        self.analysis.last_operation_filenames = {}
        re_fastq_filename = re.compile(R"^%s_((?:hg|mm)\d+)_([12])\.fastq$" % self.analysis.sample, re.I)
        fastq_files = [filename for filename in os.listdir() if re_fastq_filename.match(filename)]
        for filename in fastq_files:
            match = re_fastq_filename.match(filename)
            organism = match.group(1)
            read_index = int(match.group(2))
            out_filename = "%s_%s_R%d.fastq" % (self.analysis.sample, organism, read_index)
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
            os.unlink(filename)

        self.analysis.logger.info("Finished fixing xenome fastq")

    def compress(self):
        self.analysis.logger.info("Compressing fastq files")
        self.chdir()
        fastq_files = [
            "%s_R%d.fastq" % (self.analysis.sample, index + 1)
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
