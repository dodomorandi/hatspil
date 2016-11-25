from .exceptions import PipelineError
from . import utils

import os
import subprocess
from formatizer import f


class VarScan:

    def __init__(self, analysis):
        self.analysis = analysis

        self.minCov = 10
        self.minVar = 2
        self.minQ = 0
        self.freq = 0.05
        self.first_fifo = "/data/scratch/matteo/%s.fifo" % self.analysis.basename
        self.second_fifo = "/data/scratch/matteo/%s.indel.fifo" % self.analysis.basename

    def chdir(self):
        os.chdir(self.analysis.out_dir)

    def run(self):
        self.analysis.logger.info("Running VarScan")
        self.chdir()
        if len(self.analysis.bamfiles) == 0:
            self.analysis.logger.warning(f(
                "Empty BAM list. VarScan has nothing to do."
            ))
            return

        if not os.path.exists(self.first_fifo):
            os.mkfifo(self.first_fifo)
        if not os.path.exists(self.second_fifo):
            os.mkfifo(self.second_fifo)

        config = self.analysis.config
        input_files = self.analysis.bamfiles.items()

        for organism, bamfile in input_files:
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(config, organism)[0]
            pileup_process = subprocess.Popen(f(
                "{config.samtools} mpileup -d 8000 -f {genome_ref} "
                "{bamfile} | tee {self.first_fifo} > {self.second_fifo}"),
                shell=True)

            snp_process = subprocess.Popen(f(
                "{config.varscan} mpileup2snp {self.first_fifo} "
                "--min-coverage {self.minCov} --min-reads2 {self.minVar} "
                "--min-avg-qual {self.minQ} --min-var-freq {self.freq} "
                "--p-value 1 --output-vcf 1 > {self.analysis.basename}{organism_str}.varscan2.vcf"
                ), shell=True
            )

            indel_process = subprocess.Popen(f(
                "{config.varscan} mpileup2indel {self.second_fifo} "
                "--min-coverage {self.minCov} --min-reads2 {self.minVar} "
                "--min-avg-qual {self.minQ} --min-var-freq {self.freq} "
                "--p-value 1 --output-vcf 1 > {self.analysis.basename}{organism_str}.varscan2.indel.vcf"
                ), shell=True
            )

            self.analysis.logger.info("Waiting for Variant and InDel calls for id %s%s..." % (self.analysis.basename, organism_str))

            pileup_finished = False
            snp_finished = False
            indel_finished = False

            while not (pileup_finished and snp_finished and indel_finished):
                if not pileup_finished:
                    try:
                        pileup_retval = pileup_process.wait(5)
                        pileup_finished = True
                        self.analysis.logger.info("samtool pileup finished")
                    except:
                        pass

                if not snp_finished:
                    try:
                        snp_retval = snp_process.wait(5)
                        snp_finished = True
                        self.analysis.logger.info("VarScan mpileup2snp finished")
                    except:
                        pass

                if not indel_finished:
                    try:
                        indel_retval = indel_process.wait(5)
                        indel_finished = True
                        self.analysis.logger.info("VarScan mpileup2indel finished")
                    except:
                        pass

            if snp_retval != 0 or indel_retval != 0 or pileup_retval != 0:
                if snp_retval != 0:
                    self.analysis.logger.error("VarScan mpileup2snp failed")
                if indel_retval != 0:
                    self.analysis.logger.error("VarScan mpileup2indel failed")
                if pileup_process != 0:
                    self.analysis.logger.error("samtools pileup failed")

                raise PipelineError("varscan error")

        os.unlink(self.first_fifo)
        os.unlink(self.second_fifo)

        self.analysis.logger.info("Finished VarScan")
