"""Module to ease the process of running an analysis.

Each analysis is run with a high degree of freedom, therefore it is
useful to have an abstraction layer to easily run a single analysis.
This module contains this layer of abstraction. 
"""
from multiprocessing.managers import SyncManager
from threading import Thread
from typing import Any, Dict, Iterable, List, Optional

from .config import Config
from .core.analysis import Analysis
from .core.barcoded_filename import Analyte, BarcodedFilename
from .core.starter import Starter
from .db import Db
from .mapping import Mapping
from .mutect import Mutect
from .reports.reports_generator import ReportsGenerator
from .strelka import Strelka
from .variant_calling import VariantCalling
from .varscan import VarScan


class Runner:
    """The class to ease the process of running an analysis.

    An instance of the class can be called to perform the main step of
    the analysis. Samples coupled with a control need a two-step
    analysis (both sample and control analyses have to be finished
    before the part that uses both) can be furtherly analysed using the
    `with_normals` method. The `generate_reports` can be run at the end
    of the analysis process in order to create a report file for the
    current analysis and a global report with all the information stored
    in the database.
    """

    def __init__(
        self,
        manager: SyncManager,
        root: str,
        config: Config,
        parameters: Dict[str, Any],
        fastq_dir: str,
    ) -> None:
        """Create an instance of the class.

        Args:
            manager: an instance of a multiprocessing Manager.
            root: the root directory for the analysis.
            config: the configuration of the analysis.
            parameters: the parameters obtained from command line.
            fastq_dir: the directory containing the FASTQ files that are
                       going to be analysed.
        """
        self.last_operations: Dict[str, Any] = manager.dict()
        self.root = root
        self.config = config
        self.parameters = parameters
        self.fastq_dir = fastq_dir

    def __call__(self, sample: str) -> None:
        """Run the analysis for the specified sample.

        If the sample has both tumor and control data, the mutations
        detection and the variant calling are not performed
        automatically. In this case it is necessary to call
        `Runner.with_normals` once both tumor and control files have
        been processed through this function.
        """
        analysis = Analysis(sample, self.root, self.config, self.parameters)

        barcoded_sample = BarcodedFilename.from_sample(sample)
        if self.config.use_mongodb:
            db = Db(analysis.config)
            db.store_barcoded(barcoded_sample)

        Starter.run(analysis, self.fastq_dir)

        filenames: List[str] = []
        if analysis.last_operation_filenames is not None:
            if isinstance(analysis.last_operation_filenames, str):
                filenames = [analysis.last_operation_filenames]
            elif isinstance(analysis.last_operation_filenames, dict):
                filenames = [
                    filename
                    for filenames_list in analysis.last_operation_filenames.values()
                    for filename in filenames_list
                ]
            else:
                filenames = analysis.last_operation_filenames

        mapping = Mapping(analysis, self.fastq_dir)
        mapping.run()

        if (
            not self.parameters["use_normals"]
            and not self.parameters["only_mapping"]
            and barcoded_sample.analyte != Analyte.RNASEQ
        ):
            self._run_mutation_analysis(analysis, False)

        self.last_operations[sample] = analysis.last_operation_filenames

    def with_normals(self, sample: str, tumor: str, normal: Optional[str]) -> None:
        """Perform the second part of the analysis using control data.

        Runs the mutation analysis and the variant calling pairing the
        tumor data with the controls.

        Args:
            sample: the name of the sample.
            tumor: the filename of the BAM obtained from the tumor.
            normal: the filename of the BAM obtained from the normal,
                    if any. In case `normal` is null, the analysis will
                    continue without a control for the tumor.
        """
        if not self.parameters["use_normals"] or self.parameters["only_mapping"]:
            return

        barcoded_filename = BarcodedFilename.from_sample(sample)
        if barcoded_filename.analyte == Analyte.RNASEQ:
            return

        analysis = Analysis(sample, self.root, self.config, self.parameters)

        db = Db(analysis.config)
        for filename in (tumor, normal):
            if filename is None:
                continue

            barcoded_filename = BarcodedFilename(filename)
            db.store_barcoded(barcoded_filename)

        if normal is None:
            analysis.last_operation_filenames = {"sample": [tumor], "control": []}

            self._run_mutation_analysis(analysis, False)
        else:
            analysis.last_operation_filenames = {"sample": [tumor], "control": [normal]}
            self._run_mutation_analysis(analysis, True)

    def generate_reports(self, samples: Iterable[str]) -> None:
        """Generate the analysis and the global reports.

        The behaviour of the function can be changed by the command line
        parameters `generate_report` and `generate_global_report`.

        Args:
            samples: the samples that are included in the analysis
                     report.
        """
        barcoded_samples = [BarcodedFilename.from_sample(sample) for sample in samples]
        reports_generator = ReportsGenerator(
            self.root, self.config, self.parameters, barcoded_samples
        )

        if self.parameters["generate_report"]:
            reports_generator.generate_analysis_reports()

        if self.parameters["generate_global_report"]:
            reports_generator.generate_global_reports()

    def _run_mutation_analysis(self, analysis: Analysis, use_strelka: bool) -> None:
        mutect = Mutect(analysis)
        varscan = VarScan(analysis)

        mutect.run()
        varscan.run()

        if use_strelka and analysis.using_normals:
            strelka = Strelka(analysis)
            strelka.run()

        variant_calling = VariantCalling(analysis)
        variant_calling.run()
