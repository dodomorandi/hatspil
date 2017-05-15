from .executor import Executor

import os
from formatizer import f


class Mutect:

    def __init__(self, analysis):
        self.analysis = analysis

    def chdir(self):
        os.chdir(self.analysis.get_out_dir())

    def run(self):
        self.analysis.logger.info("Running mutect")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        if self.analysis.parameters["use_normals"]:
            executor(f(
                "{config.java7} {config.mutect_jvm_args} -jar {config.mutect} "
                '{config.mutect_args} --reference_sequence '
                '{{genome_ref}} --dbsnp {{dbsnp}} --cosmic '
                '{{cosmic}} --input_file:tumor {{input_filename["tumor"]}} '
                '--input_file:normal {{input_filename["normal"]}} --out '
                '{self.analysis.basename}{{organism_str}}.mutect.1.17.vcf '),
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error",
                use_normals=True
            )
        else:
            executor(f(
                "{config.java7} {config.mutect_jvm_args} -jar {config.mutect} "
                '{config.mutect_args} --reference_sequence '
                '{{genome_ref}} --dbsnp {{dbsnp}} --cosmic '
                '{{cosmic}} --input_file:tumor {{input_filename}} --out '
                '{self.analysis.basename}{{organism_str}}.mutect.1.17.vcf '),
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error"
            )

        self.analysis.logger.info("Finished mutect")
