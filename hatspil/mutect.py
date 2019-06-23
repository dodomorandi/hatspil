"""The module to handle MuTect."""
import os

from .analysis import Analysis
from .executor import Executor


class Mutect:
    """Helper class to run MuTect on BAM files."""

    def __init__(self, analysis: Analysis) -> None:
        """Create an instance of the class."""
        self.analysis = analysis

        os.makedirs(analysis.get_out_dir(), exist_ok=True)

    def chdir(self) -> None:
        """Change current directory to the variant calling folder."""
        os.chdir(self.analysis.get_out_dir())

    def run(self) -> None:
        """Run MuTect according to the input files.

        The parameters are different depending on the presence of
        control samples.
        """
        self.analysis.logger.info("Running mutect")
        self.chdir()
        config = self.analysis.config

        executor = Executor(self.analysis)
        if self.analysis.using_normals:
            executor(
                f"{config.java7} {config.mutect_jvm_args} "
                f"-jar {config.mutect} "
                f"{config.mutect_args} --reference_sequence "
                f"{{genome_ref}} --dbsnp {{dbsnp}} --cosmic "
                f"{{cosmic}} --input_file:tumor {{input_filenames.sample}} "
                f"--input_file:normal {{input_filenames.control}} --out "
                f"{self.analysis.basename}.mutect.1.17{{organism_str}}.vcf ",
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error",
                use_normals=True,
                split_input_files=False,
            )
        else:
            executor(
                f"{config.java7} {config.mutect_jvm_args} "
                f"-jar {config.mutect} "
                f"{config.mutect_args} --reference_sequence "
                f"{{genome_ref}} --dbsnp {{dbsnp}} --cosmic "
                f"{{cosmic}} --input_file:tumor {{input_filename}} --out "
                f"{self.analysis.basename}.mutect.1.17{{organism_str}}.vcf ",
                override_last_files=False,
                error_string="Mutect exited with status {status}",
                exception_string="mutect error",
            )

        self.analysis.logger.info("Finished mutect")
