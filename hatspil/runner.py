from .analysis import Analysis
from .xenograft import Xenograft
from .mapping import Mapping
from .mutect import Mutect
from .varscan import VarScan
from .strelka import Strelka
from .barcoded_filename import BarcodedFilename
from .starter import Starter


class Runner:

    def __init__(self, manager, root, config, parameters, fastq_dir):
        self.last_operations = manager.dict()
        self.root = root
        self.config = config
        self.parameters = parameters
        self.fastq_dir = fastq_dir

    def __call__(self, sample):
        analysis = Analysis(sample, self.root, self.config, self.parameters)
        barcoded_filename = BarcodedFilename.from_sample(sample)

        Starter.run(analysis, self.fastq_dir)

        if self.parameters["use_xenome"] and\
                barcoded_filename.tissue == 60:
            xenograft = Xenograft(analysis, self.fastq_dir)
            xenograft.run()

        mapping = Mapping(analysis, self.fastq_dir)
        mapping.run()

        if not self.parameters["use_normals"]:
            mutect = Mutect(analysis)
            varscan = VarScan(analysis)

            mutect.run()
            varscan.run()

        self.last_operations[sample] = analysis.last_operation_filenames

    def with_normals(self, sample, tumor, normal):
        if not self.parameters["use_normals"]:
            return

        analysis = Analysis(sample, self.root, self.config, self.parameters)
        analysis.last_operation_filenames = [tumor, normal]

        mutect = Mutect(analysis)
        varscan = VarScan(analysis)
        strelka = Strelka(analysis)

        mutect.run()
        varscan.run()
        strelka.run()
