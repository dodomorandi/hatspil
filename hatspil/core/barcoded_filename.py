"""A module to handle HaTSPiL filename barcoding.

HaTSPiL highly relies on a filename barcoding in order to perform the
correct type of analysis. The barcode is heavily inspired from the one
used by TCGA, and it has been developed with the idea of being scalable
to a high variety of cases.

Here a brief example of a barcode:
lung-lc001-02-011-201.R1.fastq

lung - The name of the project. In this example, the project involves
       the analysis of lung cancer tissues.
lc001 - The name of the sample.
02 - The type of the tissue, in this case a recurrent solid tumor. See
     'Tissue' for more information.
0 - The molecule analysed, in this case DNA. See 'Molecule' for more
    information.
1 - The type of NGS experiment, in this case a gene panel assay. See
    'Analyte' for more information.
2 - The index of the biopsy, in this case the third (indices are
    0-based).
0 - The index of the sample, in this case the first (indices are
    0-based).
1 - The index of the sequencing, in this case the second (indices are
    0-based).

The R1 parameter is the standard way to annotate the index of the paired
index for paired-end analyses.
"""

import os
import re
from enum import IntEnum
from typing import Dict, Optional, Union

from .exceptions import BarcodeError


class Molecule(IntEnum):
    """The molecule used in the NGS assay.

    It specifies if the NGS library preparation started from DNA or RNA.
    """

    DNA = 0
    RNA = 1


class Analyte(IntEnum):
    """The type of NGS assay."""

    WHOLE_EXOME = 0
    GENE_PANEL = 1
    FUSION_PANEL = 2
    RNASEQ = 3


class Tissue(IntEnum):
    """The type of tissue.

    The values are taken from the specifications of TCGA, in order to
    provide some sort of compatibility.
    """

    PRIMARY_SOLID_TUMOR = 1
    RECURRENT_SOLID_TUMOR = 2
    PRIMARY_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD = 3
    RECURRENT_BLOOD_DERIVED_CANCER_BONE_MARROW = 4
    ADDITIONAL_NEW_PRIMARY = 5
    METASTATIC = 6
    ADDITIONAL_METASTATIC = 7
    HUMAN_TUMOR_ORGANIC_CELLS = 8
    PRIMARY_BLOOD_DERIVED_CANCER_BONE_MARROW = 9
    BLOOD_DERIVED_NORMAL = 10
    SOLID_TISSUE_NORMAL = 11
    BUCCAL_CELL_NORMAL = 12
    EBV_IMMORTALIZED_NORMAL = 13
    BONE_MARROW_NORMAL = 14
    SAMPLE_TYPE_15 = 15
    SAMPLE_TYPE_16 = 16
    CONTROL_ANALYTE = 20
    RECURRENT_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD = 40
    CELL_LINES = 50
    PRIMARY_XENOGRAFT_TISSUE = 60
    CELL_LINE_DERIVED_XENOGRAFT_TISSUE = 61
    SAMPLE_TYPE_99 = 99

    @staticmethod
    def create(value: Union[int, str]) -> "Tissue":
        """Create an instance from a raw value.

        When the tissue derives from a xenograft, some information about
        the animal is encoded in the tissue parameter. This function can
        be used to ease the conversion from any kind of partially-
        encoded tissue raw values to a valid 'Tissue' instance.
        """
        if isinstance(value, int):
            return Tissue(value)
        elif isinstance(value, str):
            if value[0] == "6":
                if value[1].isupper():
                    return Tissue(60)
                else:
                    return Tissue(61)
            else:
                return Tissue(int(value))
        else:
            raise BarcodeError("Unexpected value type for tissue")

    def is_normal(self) -> bool:
        """Check if the tissue is normal."""
        return self in (
            Tissue.BLOOD_DERIVED_NORMAL,
            Tissue.SOLID_TISSUE_NORMAL,
            Tissue.BUCCAL_CELL_NORMAL,
            Tissue.EBV_IMMORTALIZED_NORMAL,
            Tissue.BONE_MARROW_NORMAL,
        )

    def is_tumor(self) -> bool:
        """Check if the tissue is a tumor or tumor-related."""
        return self in (
            Tissue.PRIMARY_SOLID_TUMOR,
            Tissue.RECURRENT_SOLID_TUMOR,
            Tissue.PRIMARY_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD,
            Tissue.RECURRENT_BLOOD_DERIVED_CANCER_BONE_MARROW,
            Tissue.ADDITIONAL_NEW_PRIMARY,
            Tissue.METASTATIC,
            Tissue.ADDITIONAL_METASTATIC,
            Tissue.HUMAN_TUMOR_ORGANIC_CELLS,
            Tissue.PRIMARY_BLOOD_DERIVED_CANCER_BONE_MARROW,
            Tissue.RECURRENT_BLOOD_DERIVED_CANCER_PERIPHERAL_BLOOD,
            Tissue.CELL_LINES,
            Tissue.PRIMARY_XENOGRAFT_TISSUE,
            Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE,
        )

    def is_xenograft(self) -> bool:
        """Check if the tissue derives from a xenograft."""
        return self in (
            Tissue.PRIMARY_XENOGRAFT_TISSUE,
            Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE,
        )


class Xenograft:
    """A class to handle the encoded data of a xenograft sample.

    In case a tissue derives from a xenograft, additional information is
    encoded in the barcoded filename. This class allows the storage of
    this information and it provides helper functions to ease the
    encoding and the decoding.
    """

    def __init__(self, generation: int, parent: int, child: int) -> None:
        """Create a Xenograft from the generation, parent and child parameters.

        Args:
            generation: the index of the host generation, 0-based. If
                        the value is zero, then also ``parent`` must be
                        zero.
            parent: the index of the parent host, 0-based. It must be
                    less then 3.
            child: the index of host animal, 0-based. It must be less
                   then 3.
        """
        self.generation = generation
        self.parent = parent
        self.child = child

        if self.parent > 2 or self.child > 2:
            raise BarcodeError("max supported parents and children are 3 each")

        if self.generation == 0:
            if self.parent != 0:
                raise BarcodeError(
                    "the first generation of xenograft must "
                    "have the parent attribute set to 0 (ie "
                    "the sample field can only be 0, 1 or 2)"
                )

    @staticmethod
    def create(raw_tissue: str, raw_sample: Union[str, int]) -> Optional["Xenograft"]:
        """Create a Xenograft from the raw tissue and sample values.

        When a sample is a xenograft, additional information is encoded
        in the tissue and sample fields of the barcoded filename. This
        function decodes these information and eventually returns a
        Xenograft instance.

        Args:
            raw_tissue: the raw tissue field of a barcoded filename.
            raw_sample: the raw sample field of a barcoded filename.

        Returns:
            In case the tissue is a xenograft and the encoding is valid,
            a Xenograft instance is returned. None otherwise.

        """
        if raw_tissue[0] != "6":
            return None

        if not re.match("^6[A-Za-z]$", raw_tissue):
            raise BarcodeError(
                "old xenograft barcoding detected, unable to "
                "proceed.\n'6' must be followed by A-Z in case "
                "of primary xenograft tissue or a-z in case of "
                "cell line derived from xenograft tissue. The "
                "index of the letter corresponds to the "
                "generation. The tissue value must be adjusted "
                "as well."
            )

        generation = ord(raw_tissue[1].lower()) - ord("a")
        if isinstance(raw_sample, str):
            raw_sample = int(raw_sample)
        parent = int(raw_sample / 3)
        child = int(raw_sample % 3)

        return Xenograft(generation, parent, child)

    def get_sample_index(self) -> Optional[int]:
        """Get the encoded sample index."""
        if self.parent is not None and self.child is not None:
            return self.parent * 3 + self.child
        else:
            return None

    def to_dict(self) -> Dict[str, int]:
        """Transform the Xenograft instance into a dictionary.

        Transform the Xenograft instance into a dictionary with
        generation, parent and child properties.
        """
        return {
            "generation": self.generation,
            "parent": self.parent,
            "child": self.child,
        }

    @staticmethod
    def from_dict(dictionary: Dict[str, Union[int, str]]) -> "Xenograft":
        """Create a Xenograft instance from a dictionary.

        Create a Xenograft instance from a dictionary containing the
        generation, parent and child properties.
        """
        assert "generation" in dictionary
        assert "parent" in dictionary
        assert "child" in dictionary

        return Xenograft(
            int(dictionary["generation"]),
            int(dictionary["parent"]),
            int(dictionary["child"]),
        )

    def to_human(self) -> str:
        """Create a human-readable identification string.

        The human-readable format represents the host information with
        one number and two letter. In details, the first character is a
        number representing the generation, and it is 1-based. The
        second character is a lowercase letter, representing the index
        of the parent host. The third character is an uppercase letter
        and it represents the index of the host.

        Example:
            generation = 5
            parent = 2
            child = 0

            Human readable format = 6cA
        """
        if self.generation == 0:
            return "%d%s" % (self.generation + 1, chr(ord("A") + self.child))
        else:
            return "%d%s%s" % (
                self.generation + 1,
                chr(ord("a") + self.parent),
                chr(ord("A") + self.child),
            )


class BarcodedFilename:
    """The way to handle barcoded filenames.

    This class is responsible for handling the barcoded samples and
    filenames.
    """

    re_filename = re.compile(
        r"^([^-]+)-([^-]+)-(\d[0-9A-Za-z])-(\d)(\d)(\d)-"
        r"(\d)(\d)(\d)(?:\.[^.]+)*?(?:\.((?:hg|mm)\d+))?"
        r"(?:\.R([12]))?\.(\w+?)(\.gz)?$"
    )
    re_sample = re.compile(
        r"^([^-]+)-([^-]+)-(\d[0-9A-Za-z])-(\d)(\d)(\d)-"
        r"(\d)?(\d)?(\d)?(?:\.[^.]+)*?(?:\.((?:hg|mm)\d+))?"
        r"(?:\.R([12]))?$"
    )
    project: str
    patient: str
    tissue: Tissue
    molecule: Optional[Molecule]
    analyte: Optional[Analyte]
    kit: Optional[int]
    biopsy: Optional[int]
    xenograft: Optional[Xenograft]
    sample: Optional[int]
    sequencing: Optional[int]
    organism: Optional[str]
    read_index: Optional[int]
    extension: Optional[str]
    gzipped: Optional[bool]

    def __init__(self, filename: Optional[str] = None) -> None:
        """Create an instance from an optional filename

        An instance is created and, if a filename is provided, the
        properties are set according to the name of the file.

        Args:
            filename: a barcoded filename

        """
        if filename is None:
            return

        filename = os.path.basename(filename)
        match = BarcodedFilename.re_filename.match(filename)
        if not match:
            raise RuntimeError("Error parsing barcoded filename '%s'" % filename)

        self.project = match.group(1)
        self.patient = match.group(2)
        self.tissue = Tissue.create(match.group(3))
        self.molecule = Molecule(int(match.group(4)))
        self.analyte = Analyte(int(match.group(5)))
        self.kit = int(match.group(6))
        self.biopsy = int(match.group(7))

        self.xenograft = Xenograft.create(match.group(3), match.group(8))
        if self.xenograft is None:
            self.sample = int(match.group(8))
        else:
            self.sample = None

        self.sequencing = int(match.group(9))
        self.organism = match.group(10)

        read_index = match.group(11)
        if read_index:
            self.read_index = int(read_index)
        else:
            self.read_index = None

        self.extension = match.group(12)
        if match.group(13):
            self.gzipped = True
        else:
            self.gzipped = False

    @staticmethod
    def from_sample(sample: str) -> "BarcodedFilename":
        """Create an instance from a barcoded sample.

        Create a BarcodedFilename using the string of a sample instead
        of a filename. Some fields are allowed to be invalid.

        Args:
            sample: the barcoded sample

        Returns:
            An instance of BarcodedFilename.

        """
        match = BarcodedFilename.re_sample.match(sample)
        if not match:
            raise RuntimeError("Error parsing barcoded sampling '%s'" % sample)

        barcoded = BarcodedFilename()

        barcoded.project = match.group(1)
        barcoded.patient = match.group(2)
        barcoded.tissue = Tissue.create(match.group(3))
        barcoded.molecule = Molecule(int(match.group(4)))
        barcoded.analyte = Analyte(int(match.group(5)))
        barcoded.kit = int(match.group(6))

        biopsy = match.group(7)
        if biopsy:
            barcoded.biopsy = int(biopsy)
        else:
            barcoded.biopsy = None

        if match.group(8) is None:
            barcoded.xenograft = None
            barcoded.sample = None
        else:
            barcoded.xenograft = Xenograft.create(match.group(3), match.group(8))
            if barcoded.xenograft is None:
                barcoded.sample = int(match.group(8))
            else:
                barcoded.sample = None

        sequencing = match.group(9)
        if sequencing:
            barcoded.sequencing = int(sequencing)
        else:
            barcoded.sequencing = None
        barcoded.organism = match.group(10)

        read_index = match.group(11)
        if read_index:
            barcoded.read_index = int(read_index)
        else:
            barcoded.read_index = None

        barcoded.extension = None
        barcoded.gzipped = None

        return barcoded

    @staticmethod
    def from_parameters(
        project: str,
        patient: str,
        tissue: str,
        biopsy: Union[int, str],
        sample: Union[int, str, None],
        xenograft_generation: Union[int, str, None],
        xenograft_parent: Union[int, str, None],
        xenograft_child: Union[int, str, None],
        sequencing: Union[int, str],
        molecule: Union[int, str],
        analyte: Union[int, str],
        kit: Union[int, str],
    ) -> "BarcodedFilename":
        """Create a BarcodedFilename specifying each parameter.

        Create a instance, passing the value for each field.

        Args:
            project: the project field.
            patient: the patient field.
            tissue: the raw tissue string.
            biopsy: the index of the biopsy.
            sample: the raw sample string or None.
            xenograft_generation: the index of the generation of the
                                  host, or None.
            xenograft_parent: the index of the parent of the host, or
                              None.
            xenograft_child: the index of the host, or None.
            sequencing: the index of the sequencing.
            molecule: the raw value of the molecule.
            analyte: the raw value of the analyte.
            kit: the value of the kit.

        Returns:
            An instance of BarcodedFilename

        """
        barcoded = BarcodedFilename()
        barcoded.project = project
        barcoded.patient = patient
        barcoded.tissue = Tissue.create(tissue)
        barcoded.biopsy = int(biopsy)

        if xenograft_generation is not None:
            if xenograft_parent is None or xenograft_child is None:
                raise Exception(
                    "all three xenograft parameters (generation, parent and "
                    "child) must be specified"
                )

            if not barcoded.tissue.is_xenograft():
                raise Exception(
                    "xenograft parameters have been passed, but the tissue "
                    "is not a xenograft"
                )

            barcoded.xenograft = Xenograft(
                int(xenograft_generation), int(xenograft_parent), int(xenograft_child)
            )
            barcoded.sample = None
        else:
            if sample is None:
                raise Exception(
                    "'sample' or all xenograft parameters must be specified"
                )

            if barcoded.tissue.is_xenograft():
                raise Exception(
                    "'sample' parameter has been passed, but tissue is xenograft"
                )

            barcoded.sample = int(sample)
            barcoded.xenograft = None

        barcoded.sequencing = int(sequencing)
        barcoded.molecule = Molecule(int(molecule))
        barcoded.analyte = Analyte(int(analyte))
        barcoded.kit = int(kit)

        barcoded.organism = None
        barcoded.read_index = None
        barcoded.extension = None
        barcoded.gzipped = None

        return barcoded

    def get_barcode(self) -> str:
        """Get the barcode.

        Transform the instance into a barcode string. Invalid fields are
        represented with square chars.

        Returns:
            A string with the encoded barcode.

        """
        if self.biopsy is None:
            biopsy = "□"
        else:
            biopsy = str(self.biopsy)

        if self.sequencing is None:
            sequencing = "□"
        else:
            sequencing = str(self.sequencing)

        sample_index = self.get_sample_index()
        if sample_index is None:
            sample = "□"
        else:
            sample = str(sample_index)

        tissue = self.get_tissue_str_optional()
        if tissue is None:
            tissue_string = "□"
        else:
            tissue_string = tissue

        if self.molecule is None:
            molecule = "□"
        else:
            molecule = str(int(self.molecule))

        if self.analyte is None:
            analyte = "□"
        else:
            analyte = str(int(self.analyte))

        if self.kit is None:
            kit = "□"
        else:
            kit = str(int(self.kit))

        return "%s-%s-%02s-%s%s%s-%s%s%s" % (
            self.project,
            self.patient,
            tissue_string,
            molecule,
            analyte,
            kit,
            biopsy,
            sample,
            sequencing,
        )

    def get_barcoded_filename(self) -> Optional[str]:
        """Create a barcoded filename.

        Create a valid barcoded filename if all the fields are valid.

        Returns:
            A string representing a barcoded filename, or None if some
            required fields are None.

        """
        if self.biopsy is None or self.sequencing is None or self.extension is None:
            return None

        if self.read_index:
            read_str = f".R{self.read_index}"
        else:
            read_str = ""

        if self.gzipped:
            extension = f"{self.extension}.gz"
        else:
            extension = self.extension

        if self.organism:
            organism = f".{self.organism}"
        else:
            organism = ""

        return "{}-{}-{}-{}{}{}-{}{}{}{}{}.{}".format(
            self.project,
            self.patient,
            self.get_tissue_str(),
            self.molecule,
            self.analyte,
            self.kit,
            self.biopsy,
            self.get_sample_index(),
            self.sequencing,
            organism,
            read_str,
            extension,
        )

    def get_directory(self, analysis_dir: Optional[str] = None) -> str:
        """Return the directory for the analysis.

        Evaluate the output directory depending on the type of the
        analysis.

        Args:
            analysis_dir: the root directory for the analysis. If it is
                          None, a relative path is returned.

        Returns:
            The output directory depending on the type of analysis.

        """
        if analysis_dir is None:
            analysis_dir = ""

        if self.analyte == Analyte.WHOLE_EXOME:
            barcode_dir = "WXS"
        elif self.analyte == Analyte.GENE_PANEL:
            barcode_dir = "Panel"
        elif self.analyte == Analyte.FUSION_PANEL:
            barcode_dir = "Fusion"
        else:
            barcode_dir = ""

        return os.path.join(analysis_dir, barcode_dir)

    def __repr__(self) -> str:
        """Get a human readable representation of the instance."""
        repr = (
            "{project="
            + str(self.project)
            + " patient="
            + str(self.patient)
            + " tissue="
            + self.get_tissue_str()
            + " molecule="
            + str(self.molecule)
            + " analyte="
            + str(self.analyte)
            + " kit="
            + str(self.kit)
            + " biopsy="
            + str(self.biopsy)
        )

        if self.xenograft:
            repr += " xenograft=" + self.xenograft.to_human()

        repr += (
            " sample="
            + str(self.get_sample_index())
            + " sequencing="
            + str(self.sequencing)
        )

        if self.organism:
            repr += " organism=" + str(self.organism)
        if self.read_index:
            repr += " read_index=" + str(self.read_index)
        if self.extension:
            repr += " extension=" + self.extension
        if self.gzipped:
            " gzipped=" + str(self.gzipped)

        repr += "}"

        return repr

    def get_tissue_str(self) -> str:
        """Get the raw tissue string.

        The raw tissue string is created, including the xenograft data
        when appropriate. In case of invalid xenograft data combination,
        an exception is raised.
        """
        if self.tissue == Tissue.PRIMARY_XENOGRAFT_TISSUE:
            assert self.xenograft is not None
            return "6" + chr(ord("A") + self.xenograft.generation)
        elif self.tissue == Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE:
            assert self.xenograft is not None
            return "6" + chr(ord("a") + self.xenograft.generation)
        else:
            return "%02d" % self.tissue

    def get_tissue_str_optional(self) -> Optional[str]:
        """Get the raw tissue string.

        Similar to 'get_tissue_str', but in case of invalid data None is
        returned instead of raising an exception.
        """
        if self.tissue == Tissue.PRIMARY_XENOGRAFT_TISSUE:
            if self.xenograft is None:
                return None
            else:
                return "6" + chr(ord("A") + self.xenograft.generation)
        elif self.tissue == Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE:
            if self.xenograft is None:
                return None
            else:
                return "6" + chr(ord("a") + self.xenograft.generation)
        else:
            if self.tissue is None:
                return None
            else:
                return "%02d" % self.tissue

    def get_sample_index(self) -> Optional[int]:
        """Get the raw sample index.

        Return the index of the sample, evaluating whether this is a
        xenograft or not.
        """
        if self.xenograft is None:
            return self.sample
        else:
            return self.xenograft.get_sample_index()

    def is_xenograft(self) -> bool:
        """Check whether the sample is a xenograft."""
        return self.xenograft is not None

    def equals_without_tissue(self, other: "BarcodedFilename") -> bool:
        """Compare two barcodes ignoring the tissue.

        Return true if two barcodes has the same project, patient,
        molecule, analyte, biopsy and kit.
        """
        return (
            self.project == other.project
            and self.patient == other.patient
            and self.molecule == other.molecule
            and self.analyte == other.analyte
            and self.kit == other.kit
            and self.biopsy == other.biopsy
        )
