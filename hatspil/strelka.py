from .executor import Executor

from formatizer import f
import os
import shutil


class Strelka:

    def __init__(self, analysis):
        self.analysis = analysis
        self.strelka_perl = os.path.join(
            self.analysis.config.strelka_basedir,
            "bin", "configureStrelkaWorkflow.pl")
        self.strelka_dir = os.path.join("Strelka", self.analysis.basename)

        self.chdir()

        try:
            os.makedirs("Strelka", exist_ok=True)
        except OSError:
            pass

        if os.path.exists(self.strelka_dir):
            shutil.rmtree(self.strelka_dir)

    def chdir(self):
        os.chdir(self.analysis.root)

    def configure_strelka(self):
        self.analysis.logger.info("Configuring strelka")
        self.chdir()

        executor = Executor(self.analysis)
        executor(f(
            "{self.strelka_perl} "
            "--tumor={{input_filename[\"tumor\"]}} "
            "--normal={{input_filename[\"normal\"]}} "
            "--ref={{genome_ref}} "
            "--config={self.analysis.config.strelka_config} "
            "--output-dir={self.strelka_dir}"),
            override_last_files=False,
            use_normals=True
        )

        self.analysis.logger.info("Finished configuring strelka")

    def make(self):
        self.analysis.logger.info("Running make for strelka")
        self.chdir()
        os.chdir(self.strelka_dir)

        executor = Executor(self.analysis)
        executor(f(
            "make -j {self.analysis.config.strelka_threads}"),
            override_last_files=False,
            use_normals=True
        )

        self.analysis.logger.info("Finished make for strelka")

    def run(self):
        self.configure_strelka()
        self.make()
