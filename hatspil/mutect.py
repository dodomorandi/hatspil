import os

from .executor import Executor


class Mutect:
    def __init__(self, analysis):
        self.analysis = analysis

        os.makedirs(analysis.get_out_dir(), exist_ok=True)

    def chdir(self):
        os.chdir(self.analysis.get_out_dir())

    def run(self):
        self.analysis.logger.info("Running mutect")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        if self.analysis.using_normals:
            executor(
                f"{config.java7} {config.mutect_jvm_args} -jar {config.mutect} "
                f'{config.mutect_args} --reference_sequence '
                f'{{genome_ref}} --dbsnp {{dbsnp}} --cosmic '
                f'{{cosmic}} --input_file:tumor {{input_filename["tumor"]}} '
                f'--input_file:normal {{input_filename["normal"]}} --out '
                f'{self.analysis.basename}{{organism_str}}.mutect.1.17.vcf ',
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error",
                use_normals=True)
        else:
            executor(
                f"{config.java7} {config.mutect_jvm_args} -jar {config.mutect} "
                f'{config.mutect_args} --reference_sequence '
                f'{{genome_ref}} --dbsnp {{dbsnp}} --cosmic '
                f'{{cosmic}} --input_file:tumor {{input_filename}} --out '
                f'{self.analysis.basename}{{organism_str}}.mutect.1.17.vcf ',
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error")

        self.analysis.logger.info("Finished mutect")
