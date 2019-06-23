"""Module to handle Strelka software."""
import os
import shutil

from .analysis import Analysis
from .executor import Executor


class Strelka:
    """Class to help the multi-step process of running Strelka."""

    def __init__(self, analysis: Analysis) -> None:
        """Create an instance of the class.

        It creates the directory `Strelka` in the root directory in case
        it does not exist. If a directory with the basename of the
        sample already exists inside `Strelka`, it is automatically
        deleted.
        """
        self.analysis = analysis
        self.strelka_perl = os.path.join(
            self.analysis.config.strelka_basedir, "bin", "configureStrelkaWorkflow.pl"
        )
        self.strelka_dir = os.path.join("Strelka", self.analysis.basename)

        self.chdir()
        os.makedirs("Strelka", exist_ok=True)

        if os.path.exists(self.strelka_dir):
            shutil.rmtree(self.strelka_dir)

    def chdir(self) -> None:
        """Change current directory to the root directory."""
        os.chdir(self.analysis.root)

    def configure_strelka(self) -> None:
        """Run the configuration step of Strelka.

        The process uses the `configureStrelkaWorkflow.pl` tool.
        """
        self.analysis.logger.info("Configuring strelka")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            f"{self.strelka_perl} "
            f"--tumor={{input_filenames.sample}} "
            f"--normal={{input_filenames.control}} "
            f"--ref={{genome_ref}} "
            f"--config={self.analysis.config.strelka_config} "
            f"--output-dir={self.strelka_dir}",
            override_last_files=False,
            use_normals=True,
            split_input_files=False,
        )

        self.analysis.logger.info("Finished configuring strelka")

    def make(self) -> None:
        """Run `make` for `Strelka`.

        It changes the directory inside the Strelka directory for the
        current sample, and then it runs `make`.

        This function must be run only after `configure_strelka`.
        """
        self.analysis.logger.info("Running make for strelka")
        if not self.analysis.run_fake:
            self.chdir()
            os.chdir(self.strelka_dir)

        executor = Executor(self.analysis)
        executor(
            f"make -j {self.analysis.config.strelka_threads}",
            override_last_files=False,
            use_normals=True,
            split_input_files=False,
        )

        self.analysis.logger.info("Finished make for strelka")

    def run(self) -> None:
        """Run the Strelka software.

        It performs the configuration for Strelka, then it uses `make`
        in order to run the analysis for the current samples.

        It is important to note that Strelka needs a control sample to
        perform the analysis.
        """
        self.configure_strelka()
        self.make()
