from .analysis import Analysis
from .xenograft import Xenograft
from .mapping import Mapping
from .mutect import Mutect
from .varscan import VarScan
from .strelka import Strelka
from .variant_calling import VariantCalling
from .barcoded_filename import BarcodedFilename, Tissue
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

        if self.parameters["use_xenome"] and (
                barcoded_filename.tissue == Tissue.PRIMARY_XENOGRAFT_TISSUE or
                barcoded_filename.tissue == Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE):
            xenograft = Xenograft(analysis, self.fastq_dir)
            xenograft.run()

        if self.parameters["skip_mapping"]:
            analysis.run_fake = True
        mapping = Mapping(analysis, self.fastq_dir)
        mapping.run()
        analysis.run_fake = False

        if not self.parameters["use_normals"] and \
                not self.parameters["only_mapping"]:
            mutect = Mutect(analysis)
            varscan = VarScan(analysis)

            mutect.run()
            varscan.run()

            variant_calling = VariantCalling(analysis)
            variant_calling.run()

        self.last_operations[sample] = analysis.last_operation_filenames

    def with_normals(self, sample, tumor, normal):
        if not self.parameters["use_normals"] \
                or self.parameters["only_mapping"]:
            return

        analysis = Analysis(sample, self.root, self.config, self.parameters)
        analysis.last_operation_filenames = {
            "tumor": [tumor],
            "normal": [normal]
        }

        mutect = Mutect(analysis)
        varscan = VarScan(analysis)
        strelka = Strelka(analysis)

        mutect.run()
        varscan.run()
        strelka.run()

        variant_calling = VariantCalling(analysis)
        variant_calling.run()
