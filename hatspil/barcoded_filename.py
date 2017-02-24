import re
import os
from enum import IntEnum


class Molecule(IntEnum):
    DNA = 0,
    RNA = 1


class Analyte(IntEnum):
    WHOLE_EXOME = 0
    GENE_PANEL = 1


class Tissue(IntEnum):
    PRIMARY_SOLID_TUMOR = 1,
    RECURRENT_SOLID_TUMOR = 2,
    PRIMARY_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD = 3,
    RECURRENT_BLOOD_DERIVED_CANCER_BONE_MARROW = 4
    ADDITIONAL_NEW_PRIMARY = 5,
    METASTATIC = 6,
    ADDITIONAL_METASTATIC = 7
    HUMAN_TUMOR_ORGANIC_CELLS = 8,
    PRIMARY_BLOOD_DERIVED_CANCER_BONE_MARROW = 9,
    BLOOD_DERIVED_NORMAL = 10,
    SOLID_TISSUE_NORMAL = 11,
    BUCCAL_CELL_NORMAL = 12,
    EBV_IMMORTALIZED_NORMAL = 13,
    BONE_MARROW_NORMAL = 14,
    SAMPLE_TYPE_15 = 15,
    SAMPLE_TYPE_16 = 16,
    CONTROL_ANALYTE = 20,
    RECURRENT_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD = 40,
    CELL_LINES = 50,
    PRIMARY_XENOGRAFT_TISSUE = 60,
    CELL_LINE_DERIVED_XENOGRAFT_TISSUE = 61,
    SAMPLE_TYPE_99 = 99


class BarcodedFilename:
    re_filename = re.compile(R"^([^-]+)-([^-]+)-(\d{2})-(\d)(\d)(\d)-(\d)(\d)"
                             R"(?:\.[^.]+)*?(?:\.((?:hg|mm)\d+))?"
                             R"(?:\.R([12]))?(?:\.[^.]+)*?\.(\w+?)(\.gz)?$")
    re_sample = re.compile(R"^([^-]+)-([^-]+)-(\d{2})-(\d)(\d)(\d)-(\d)?(\d)?"
                           R"(?:\.[^.]+)*?(?:\.((?:hg|mm)\d+))?"
                           R"(?:\.R([12]))?(?:\.[^.]+)*?$")

    def __init__(self, filename=None):
        if filename is None:
            return

        filename = os.path.basename(filename)
        match = BarcodedFilename.re_filename.match(filename)
        if not match:
            raise RuntimeError(
                "Error parsing barcoded filename '%s'"
                % filename)

        self.project = match.group(1)
        self.patient = match.group(2)
        self.tissue = Tissue(int(match.group(3)))
        self.molecule = Molecule(int(match.group(4)))
        self.analyte = Analyte(int(match.group(5)))
        self.kit = int(match.group(6))
        self.biopsy = int(match.group(7))
        self.sample = int(match.group(8))
        self.organism = match.group(9)

        self.read_index = match.group(10)
        if self.read_index:
            self.read_index = int(self.read_index)

        self.extension = match.group(11)
        if match.group(12):
            self.gzipped = True
        else:
            self.gzipped = False

    def from_sample(sample):
        match = BarcodedFilename.re_sample.match(sample)
        if not match:
            raise RuntimeError(
                "Error parsing barcoded sampling '%s'"
                % sample)

        barcoded = BarcodedFilename()

        barcoded.project = match.group(1)
        barcoded.patient = match.group(2)
        barcoded.tissue = int(match.group(3))
        barcoded.molecule = int(match.group(4))
        barcoded.analyte = int(match.group(5))
        barcoded.kit = int(match.group(6))

        barcoded.biopsy = match.group(7)
        if barcoded.biopsy:
            barcoded.biopsy = int(barcoded.biopsy)

        barcoded.sample = match.group(8)
        if barcoded.sample:
            barcoded.sample = int(barcoded.sample)
        barcoded.organism = match.group(9)

        barcoded.read_index = match.group(10)
        if barcoded.read_index:
            barcoded.read_index = int(barcoded.read_index)

        barcoded.extension = None
        barcoded.gzipped = None

        return barcoded
