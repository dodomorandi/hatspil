import os

from . import utils
from .analysis import Analysis
from .executor import Executor


class Starter:
    @staticmethod
    def run(analysis: Analysis, directory: str = ".") -> None:
        human_annotation = utils.get_human_annotation(analysis.config)
        input_filenames = [
            os.path.join(directory, filename)
            for pair
            in utils.find_fastqs_by_organism(
                analysis.sample,
                directory,
                human_annotation).values()
            for filename, _ in pair
        ]

        executor = Executor(analysis)
        executor(
            lambda **kwargs: None,
            input_filenames=input_filenames,
            output_format="{input_filename}")
        analysis.can_unlink = False
