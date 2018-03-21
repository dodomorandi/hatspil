import os
import shutil

from .executor import Executor


class Strelka:
    def __init__(self, analysis):
        self.analysis = analysis
        self.strelka_perl = os.path.join(self.analysis.config.strelka_basedir,
                                         "bin", "configureStrelkaWorkflow.pl")
        self.strelka_dir = os.path.join("Strelka", self.analysis.basename)

        self.chdir()
        os.makedirs("Strelka", exist_ok=True)

        if os.path.exists(self.strelka_dir):
            shutil.rmtree(self.strelka_dir)

    def chdir(self):
        os.chdir(self.analysis.root)

    def configure_strelka(self):
        self.analysis.logger.info("Configuring strelka")
        self.chdir()

        executor = Executor(self.analysis)
        executor(
            f"{self.strelka_perl} "
            f"--tumor={{input_filename[\"tumor\"]}} "
            f"--normal={{input_filename[\"normal\"]}} "
            f"--ref={{genome_ref}} "
            f"--config={self.analysis.config.strelka_config} "
            f"--output-dir={self.strelka_dir}",
            override_last_files=False,
            use_normals=True)

        self.analysis.logger.info("Finished configuring strelka")

    def make(self):
        self.analysis.logger.info("Running make for strelka")
        self.chdir()
        os.chdir(self.strelka_dir)

        executor = Executor(self.analysis)
        executor(
            f"make -j {self.analysis.config.strelka_threads}",
            override_last_files=False,
            use_normals=True)

        self.analysis.logger.info("Finished make for strelka")

    def run(self):
        self.configure_strelka()
        self.make()
