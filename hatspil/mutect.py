from .exceptions import PipelineError
from . import utils

import os
from formatizer import f


class Mutect:

    def __init__(self, analysis):
        self.analysis = analysis

    def chdir(self):
        os.chdir(self.analysis.out_dir)

    def run(self):
        self.analysis.logger.info("Running mutect")
        self.chdir()
        config = self.analysis.config

        retval = utils.run_and_log(f(
            '{config.java7} {config.mutect} --reference_sequence '
            '{config.genome_ref} --dbsnp {config.dbsnp138} --cosmic '
            '{config.cosmic} --input_file:tumor {self.analysis.bamfile} --out '
            '{self.analysis.basename}.mutect.1.17.vcf '
            # --input_file:normal {config.fnorm}
            ),
            self.analysis.logger
        )

        if retval != 0:
            self.analysis.logger.error("Mutect exited with status %d" % retval)
            raise PipelineError("mutect error")
