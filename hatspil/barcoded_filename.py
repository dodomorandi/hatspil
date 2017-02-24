import re
import os


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
        self.tissue = int(match.group(3))
        self.molecule = int(match.group(4))
        self.analyte = int(match.group(5))
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
