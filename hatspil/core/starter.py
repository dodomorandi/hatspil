"""
Analysis bootstrapping module.

Each _execution module_ depends from the previous step, and it is
obvious that it must exist a special module that does not have
dependencies. This module contains this exact feature.
"""
import os

from . import utils
from .analysis import Analysis
from .executor import Executor


class Starter:
    """The class that can bootstrap an analysis.

    This class only contains the `run` static method, that can be used
    to start a new analysis.
    """

    @staticmethod
    def run(analysis: Analysis, directory: str = ".") -> None:
        """Bootstrap the analysis.

        This function does not perform any real execution task, but it
        uses the behaviour of the `Executor` class in order to update
        the `Analysis` instance to contain all the information for the
        real first execution step.

        Args:
            analysis: the instance of an `Analysis` class.
            directory: the directory where to search for FASTQ files.
        """
        human_annotation = utils.get_human_annotation(analysis.config)
        input_filenames = [
            os.path.join(directory, filename)
            for pair in utils.find_fastqs_by_organism(
                analysis.sample, directory, human_annotation
            ).values()
            for filename, _ in pair
        ]

        executor = Executor(analysis)
        executor(
            lambda **kwargs: None,
            input_filenames=input_filenames,
            output_format="{input_filename}",
        )
        analysis.can_unlink = False
