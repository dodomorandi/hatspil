from .exceptions import PipelineError
from .executor import Executor

import os
import subprocess


class VarScan:

    def __init__(self, analysis):
        self.analysis = analysis

        self.min_coverage_normal = 2
        self.min_coverage_tumor = 10
        self.minVar = 2
        self.minQ = 0
        self.min_var_frequency = 0.05
        self.normal_purity = 1
        self.tumor_purity = 1
        self.min_tumor_frequency = 0.01
        self.p_value_somatic = 0.05

        try:
            os.makedirs(analysis.get_out_dir(), exist_ok=True)
        except OSError:
            pass

        if self.analysis.using_normals:
            self.first_fifo = "/data/scratch/matteo/%s.fifo" % self.analysis.basename
        else:
            self.first_fifo = "/data/scratch/matteo/%s.fifo" % self.analysis.basename
            self.second_fifo = "/data/scratch/matteo/%s.indel.fifo" % self.analysis.basename

    def chdir(self):
        os.chdir(self.analysis.get_out_dir())

    def _run_varscan_normals(self, **kwargs):
        self.analysis.logger.info("Running VarScan Somatic")

        genome_ref = kwargs["genome_ref"]
        input_filenames = kwargs["input_filename"]
        organism_str = kwargs["organism_str"]

        config = self.analysis.config

        pileup_process = subprocess.Popen(
            f"{config.samtools} mpileup "
            f"-q 1 -d 6000 -f {genome_ref} "
            f"-B {input_filenames['normal']} "
            f"{input_filenames['tumor']} "
            "| awk '{if($4 >= 6) print $0}' "
            "| awk '{if($7 != 0) print $0}'" 
            f">{self.first_fifo}",
            shell=True)

        somatic_process = subprocess.Popen(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"somatic {self.first_fifo} {self.analysis.basename}.varscan2 "
            f"--mpileup 1 "
            f"--min-coverage-normal {self.min_coverage_normal} "
            f"--min-coverage-tumor {self.min_coverage_tumor} "
            f"--min-var-freq {self.min_var_frequency} "
            f"--strand-filter 1 "
            f"--normal-purity {self.normal_purity} "
            f"--tumor-purity {self.tumor_purity} "
            f"--output-vcf 1",
            shell=True)

        self.analysis.logger.info("Waiting for VarScan somatic for id "
                                  "%s%s..." % (self.analysis.basename,
                                               organism_str))

        pileup_finished = False
        somatic_finished = False

        while not (pileup_finished and somatic_finished):
            if not pileup_finished:
                try:
                    pileup_retval = pileup_process.wait(5)
                    pileup_finished = True
                    self.analysis.logger.info(
                        "samtool mpileup finished")
                except:
                    pass

            if not somatic_finished:
                try:
                    somatic_retval = somatic_process.wait(5)
                    somatic_finished = True
                    self.analysis.logger.info("VarScan somatic finished")
                except:
                    pass

        if pileup_retval != 0 or somatic_retval != 0:
            if pileup_retval != 0:
                self.analysis.logger.error(
                    "samtool mpileup failed")
            if somatic_process != 0:
                self.analysis.logger.error("VarScan somatic failed")

            raise PipelineError("varscan error")

        self.analysis.logger.info("Finished VarScan Somatic")

    def _run_varscan_no_normals(self, **kwargs):
        self.analysis.logger.info("Running VarScan mpileup2snp/indel")

        genome_ref = kwargs["genome_ref"]
        input_filename = kwargs["input_filename"]
        organism_str = kwargs["organism_str"]

        config = self.analysis.config

        pileup_process = subprocess.Popen(
            f"{config.samtools} mpileup -d 8000 -f {genome_ref} "
            f"{input_filename} | tee {self.first_fifo} > {self.second_fifo}",
            shell=True)

        snp_process = subprocess.Popen(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"mpileup2snp {self.first_fifo} "
            f"--min-coverage {self.min_coverage_tumor} "
            f"--min-reads2 {self.minVar} "
            f"--min-avg-qual {self.minQ} "
            f"--min-var-freq {self.min_var_frequency} "
            f"--p-value 1 --output-vcf 1 "
            f">{self.analysis.basename}{organism_str}.varscan2.vcf",
            shell=True
        )

        indel_process = subprocess.Popen(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"mpileup2indel {self.second_fifo} "
            f"--min-coverage {self.min_coverage_tumor} "
            f"--min-reads2 {self.minVar} "
            f"--min-avg-qual {self.minQ} "
            f"--min-var-freq {self.min_var_frequency} "
            f"--p-value 1 --output-vcf 1 "
            f">{self.analysis.basename}{organism_str}.varscan2.indel.vcf",
            shell=True
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

        self.analysis.logger.info("Finished VarScan mpileup2snp/indel")

    def process_somatic(self):
        self.analysis.logger.info("Running VarScan processSomatic")
        self.chdir()

        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"processSomatic "
            f"{{input_filename}} "
            f"--min-tumor-freq {self.min_tumor_freq} "
            f"--p-value {self.p_value_somatic}",
            shell=True,
            override_last_files=False,
            input_filename=[f"{self.analysis.basename}.varscan2.snp.vcf"],
            error_string="VarScan processSomatic exited with status {status}",
            exception_string="varscan processSomatic error")

        executor(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"processSomatic "
            f"{{input_filename}} "
            f"--min-tumor-freq {self.min_tumor_freq} "
            f"--p-value {self.p_value_somatic}",
            shell=True,
            override_last_files=False,
            input_filename=[f"{self.analysis.basename}.varscan2.snp.vcf"],
            error_string="VarScan processSomatic exited with status {status}",
            exception_string="varscan processSomatic error")

        self.analysis.logger.info("Finished VarScan processSomatic")

    def cnv(self):
        self.analysis.logger.info("Running VarScan copynumber")
        self.chdir()

        try:
            os.makedirs("CNV", exist_ok=True)
        except OSError:
            pass
        os.chdir("CNV")

        config = self.analysis.config

        executor = Executor(self.analysis)
        executor(
            f"{config.java} {config.varscan_jvm_args} -jar {config.varscan} "
            f"copynumber "
            f"{{input_filename}} "
            f"--min-tumor-freq {self.min_tumor_freq} "
            f"--p-value {self.p_value_somatic}",
            shell=True,
            override_last_files=False,
            input_filename=[f"{self.analysis.basename}.varscan2.cnv.vcf"],
            error_string="VarScan copynumber exited with status {status}",
            exception_string="varscan copynumber error")

    def run(self):
        self.analysis.logger.info("Running VarScan")
        self.chdir()

        if os.path.exists(self.first_fifo):
            os.unlink(self.first_fifo)
        os.mkfifo(self.first_fifo)

        if not self.analysis.using_normals:
            if os.path.exists(self.second_fifo):
                os.unlink(self.first_fifo)
            os.mkfifo(self.second_fifo)

        executor = Executor(self.analysis)
        if self.analysis.using_normals:
            executor(
                self._run_varscan_normals,
                override_last_files=False,
                use_normals=True)
        else:
            executor(self._run_varscan_no_normals, override_last_files=False)

        os.unlink(self.first_fifo)
        if not self.analysis.using_normals:
            os.unlink(self.second_fifo)

        self.analysis.logger.info("Finished VarScan")
