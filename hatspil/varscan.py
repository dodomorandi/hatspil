from .exceptions import PipelineError
from .executor import Executor

import os
import subprocess
from formatizer import f


class VarScan:

    def __init__(self, analysis):
        self.analysis = analysis

        self.min_coverage_normal = 2
        self.min_coverage_tumor = 10
        self.minVar = 2
        self.minQ = 0
        self.min_var_frequency = 0.05
        if self.analysis.parameters["use_normals"]:
            self.first_fifo = "/data/scratch/matteo/%s.tumor.fifo" % self.analysis.basename
            self.second_fifo = "/data/scratch/matteo/%s.normal.fifo" % self.analysis.basename
        else:
            self.first_fifo = "/data/scratch/matteo/%s.fifo" % self.analysis.basename
            self.second_fifo = "/data/scratch/matteo/%s.indel.fifo" % self.analysis.basename

    def chdir(self):
        os.chdir(self.analysis.out_dir)

    def _run_varscan_normals(self, **kwargs):
        genome_ref = kwargs["genome_ref"]
        input_filenames = kwargs["input_filename"]
        organism_str = kwargs["organism_str"]

        config = self.analysis.config

        tumor_pileup_process = subprocess.Popen(f(
            "{config.samtools} mpileup -d 8000 -f {genome_ref} "
            "{input_filenames[\"tumor\"]} >{self.first_fifo}"),
            shell=True)

        normal_pileup_process = subprocess.Popen(f(
            "{config.samtools} mpileup -d 8000 -f {genome_ref} "
            "{input_filenames[\"normal\"]} >{self.second_fifo}"),
            shell=True)

        somatic_process = subprocess.Popen(
            f("{config.varscan} somatic {self.second_fifo} {self.first_fifo} "
              "--min-coverage-normal {self.min_coverage_normal} "
              "--min-coverage-tumor {self.min_coverage_tumor} "
              "--min-var-freq {self.min_var_frequency} "
              "--strand-filter 1 "
              "--output-snp {self.analysis.basename}{organism_str}.varscan2.snp.vcf "
              "--output-indel {self.analysis.basename}{organism_str}.varscan2.indel.vcf"),
            shell=True)

        self.analysis.logger.info("Waiting for SNP and InDel calls for id "
                                  "%s%s..." % (self.analysis.basename,
                                               organism_str))

        normal_pileup_finished = False
        tumor_pileup_finished = False
        somatic_finished = False

        while not (normal_pileup_finished and tumor_pileup_finished and somatic_finished):
            if not normal_pileup_finished:
                try:
                    normal_pileup_retval = normal_pileup_process.wait(5)
                    normal_pileup_finished = True
                    self.analysis.logger.info(
                        "samtool mpileup for normal sample finished")
                except:
                    pass

            if not tumor_pileup_finished:
                try:
                    tumor_pileup_retval = tumor_pileup_process.wait(5)
                    tumor_pileup_finished = True
                    self.analysis.logger.info(
                        "samtool mpileup for tumor sample finished")
                except:
                    pass

            if not somatic_finished:
                try:
                    somatic_retval = somatic_process.wait(5)
                    somatic_finished = True
                    self.analysis.logger.info("VarScan somatic finished")
                except:
                    pass

        if normal_pileup_retval != 0 or tumor_pileup_retval != 0 or somatic_retval != 0:
            if normal_pileup_retval != 0:
                self.analysis.logger.error(
                    "samtool mpileup for normal sample failed")
            if tumor_pileup_retval != 0:
                self.analysis.logger.error(
                    "samtool mpileup for tumor sample failed")
            if somatic_process != 0:
                self.analysis.logger.error("VarScan somatic failed")

            raise PipelineError("varscan error")

    def _run_varscan_no_normals(self, **kwargs):
        genome_ref = kwargs["genome_ref"]
        input_filename = kwargs["input_filename"]
        organism_str = kwargs["organism_str"]

        config = self.analysis.config

        pileup_process = subprocess.Popen(f(
            "{config.samtools} mpileup -d 8000 -f {genome_ref} "
            "{input_filename} | tee {self.first_fifo} > {self.second_fifo}"),
            shell=True)

        snp_process = subprocess.Popen(f(
            "{config.varscan} mpileup2snp {self.first_fifo} "
            "--min-coverage {self.min_coverage_tumor} "
            "--min-reads2 {self.minVar} "
            "--min-avg-qual {self.minQ} "
            "--min-var-freq {self.min_var_frequency} "
            "--p-value 1 --output-vcf 1 "
            ">{self.analysis.basename}{organism_str}.varscan2.vcf"
            ), shell=True
        )

        indel_process = subprocess.Popen(f(
            "{config.varscan} mpileup2indel {self.second_fifo} "
            "--min-coverage {self.min_coverage_tumor} "
            "--min-reads2 {self.minVar} "
            "--min-avg-qual {self.minQ} "
            "--min-var-freq {self.min_var_frequency} "
            "--p-value 1 --output-vcf 1 "
            ">{self.analysis.basename}{organism_str}.varscan2.indel.vcf"
            ), shell=True
        )

        self.analysis.logger.info("Waiting for Variant and InDel calls for id "
                                  "%s%s..." % (self.analysis.basename,
                                               organism_str))

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

    def run(self):
        self.analysis.logger.info("Running VarScan")
        self.chdir()

        if not os.path.exists(self.first_fifo):
            os.mkfifo(self.first_fifo)
        if not os.path.exists(self.second_fifo):
            os.mkfifo(self.second_fifo)

        executor = Executor(self.analysis)
        if self.analysis.parameters["use_normals"]:
            executor(
                self._run_varscan_normals,
                override_last_files=False,
                use_normals=True)
        else:
            executor(self._run_varscan_no_normals, override_last_files=False)

        os.unlink(self.first_fifo)
        os.unlink(self.second_fifo)

        self.analysis.logger.info("Finished VarScan")
