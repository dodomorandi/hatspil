from multiprocessing.managers import SyncManager
from threading import Thread
from typing import Any, Dict, Iterable, List, Optional

from .analysis import Analysis
from .barcoded_filename import Analyte, BarcodedFilename
from .config import Config
from .db import Db
from .mapping import Mapping
from .mutect import Mutect
from .reports_generator import ReportsGenerator
from .starter import Starter
from .strelka import Strelka
from .variant_calling import VariantCalling
from .varscan import VarScan


class Runner:
    def __init__(
        self,
        manager: SyncManager,
        root: str,
        config: Config,
        parameters: Dict[str, Any],
        fastq_dir: str,
    ) -> None:
        self.last_operations: Dict[str, Any] = manager.dict()
        self.root = root
        self.config = config
        self.parameters = parameters
        self.fastq_dir = fastq_dir

    def __call__(self, sample: str) -> None:
        analysis = Analysis(sample, self.root, self.config, self.parameters)
        Starter.run(analysis, self.fastq_dir)

        db = Db(analysis.config)
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

        for filename in filenames:
            barcoded_filename = BarcodedFilename(filename)
            db.store_barcoded(barcoded_filename)

        mapping = Mapping(analysis, self.fastq_dir)
        mapping.run()

        if (
            not self.parameters["use_normals"]
            and not self.parameters["only_mapping"]
            and barcoded_filename.analyte != Analyte.RNASEQ
        ):
            self._run_mutation_analysis(analysis, False)

        self.last_operations[sample] = analysis.last_operation_filenames

    def with_normals(self, sample: str, tumor: str, normal: Optional[str]) -> None:
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
