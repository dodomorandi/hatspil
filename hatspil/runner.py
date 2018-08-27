from multiprocessing.managers import SyncManager
from typing import Any, Dict, Optional

from .analysis import Analysis
from .barcoded_filename import BarcodedFilename, Tissue, Analyte
from .config import Config
from .mapping import Mapping
from .mutect import Mutect
from .starter import Starter
from .strelka import Strelka
from .variant_calling import VariantCalling
from .varscan import VarScan
from .xenograft import Xenograft


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
        barcoded_filename = BarcodedFilename.from_sample(sample)

        Starter.run(analysis, self.fastq_dir)

        if self.parameters["skip_mapping"]:
            analysis.run_fake = True
        mapping = Mapping(analysis, self.fastq_dir)
        mapping.run()
        analysis.run_fake = False

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
        if normal is None:
            analysis.last_operation_filenames = {"sample": [tumor], "control": []}

            self._run_mutation_analysis(analysis, False)
        else:
            analysis.last_operation_filenames = tumor
            self._run_mutation_analysis(analysis, True)

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
