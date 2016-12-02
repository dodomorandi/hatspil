from .executor import Executor

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

        executor = Executor(self.analysis)
        executor(f(
            '{config.java7} {config.mutect} --reference_sequence '
            '{{genome_ref}} --dbsnp {{dbsnp}} --cosmic '
            '{{cosmic}} --input_file:tumor {{input_filename}} --out '
            '{self.analysis.basename}{{organism_str}}.mutect.1.17.vcf '),
            override_last_files=False,
            error_string="Mutect exited with status {status}",
            exception_string="mutect error"
        )
