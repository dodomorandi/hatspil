import os

from . import utils
from .executor import Executor


class Starter:
    def run(analysis, directory="."):
        input_filenames = [
            os.path.join(directory,
                         filename) for pair in utils.find_fastqs_by_organism(
                             analysis.sample, directory).values()
            for filename, _ in pair
        ]

        executor = Executor(analysis)
        executor(
            lambda **kwargs: None,
            input_filenames=input_filenames,
            output_format="{input_filename}")
        analysis.can_unlink = False
