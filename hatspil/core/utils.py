"""A collection of utility function, shared across modules."""
import collections
import datetime
import gzip as gz
import logging
import os
import re
import shutil
import subprocess
from argparse import ArgumentTypeError
from copy import deepcopy
from logging import Logger
from typing import (Any, Callable, Dict, Generator, Iterable, List, Mapping,
                    Optional, Sequence, Tuple, TypeVar, Union, ValuesView,
                    cast)

from ..config import Config, KitData
from .barcoded_filename import BarcodedFilename
from .exceptions import AnnotationError, DataError


def get_current() -> str:
    """Get the current date in standard HaTSPiL format."""
    today = datetime.date.today()
    return "%04d_%02d_%02d" % (today.year, today.month, today.day)


def get_overridable_current_date(parameters: Dict[str, Any]) -> str:
    """Get an eventual overridden date.

    If the `parameters` dict contains a `use_date` value, return it.
    Otherwise return the result of `get_current`.
    """
    if parameters["use_date"] is None:
        return get_current()
    else:
        current_date = parameters["use_date"]
        assert isinstance(current_date, str)
        return current_date


def run_and_log(command: str, logger: Logger) -> int:
    """Run a command and log everything.

    Use `subprocess.Popen` to run a command. The standard output and the
    standard error are piped into the logger.

    Args:
        command: the command to run.
        logger: the logger.

    Returns:
        int: the exit status of the process.

    """
    logger.info("Running command: %s", command)
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        universal_newlines=True,
        bufsize=1,
    ) as process:
        (out, err) = process.communicate()

        for line in out.split("\n"):
            if line != "":
                logger.info(line)

        for line in err.split("\n"):
            if line != "":
                logger.warning(line)

        return process.wait()


def get_sample_filenames(
    obj: Union[Sequence[str], Mapping[str, List[str]], str],
    split_by_organism: bool = False,
) -> Union[List[str], Mapping[str, List[str]]]:
    """Return the filenames organised in a different way.

    Take a set of filenames in different possible shapes and reorganize
    them depending on the content and the value of `split_by_organism`.

    Args:
        obj: the filenames. It can be a string for one single filename,
             a list of filenames or a dict where each key is an organism
             code (i.e.: hg19) and the relative value is a list of
             filenames.
        split_by_organism: whether the filenames must be split by
                           organism or they must be returned all
                           together.
    Returns:
        The input filenames with the desired shape. There are different
        cases:
         * If `obj` is a list and its length is greater than 1 and
           `split_by_organism` is `True`, the organism for each file
           is obtained using `get_organism_from_filename`. A dict is
           created, where each organism maps to a list of filenames.
           If the dict contains more than one organism, it is returned,
           otherwise a list of the filenames is returned.
        * If `obj` is a list but its length is not greater than 1 or
          `split_by_organism` is `False`, a **copy** of `obj` is
          returned.
        * If `obj` is a dict and it contains more than one entry and
          `split_by_organism` is `True`, a **deep copy** of `obj` is
          returned.
        * If `obj` is a dict but it contains less than two entries or
          `split_by_organism` is `False`, a list of all the filenames
          in `obj` is returned.
        * If `obj` is a string and `split_by_organism` is `True`, the
          organism is obtained using `get_organism_from_filename`. If
          the organism is valid, a dict with the organism mapped to
          a list of one element, `obj`, is returned. Otherwise, if the
          organism is invalid (`None` or empty), a list of one element,
          `obj`, is returned.
        * If `obj` is a string but `split_by_organism` is `False`, a
          list of one element, `obj`, is returned.

    """
    if isinstance(obj, list):
        if split_by_organism and len(obj) > 1:
            filenames: Dict[str, List[str]] = {}
            for filename in obj:
                organism = get_organism_from_filename(filename)
                if organism is None:
                    organism = ""
                filenames.setdefault(organism, []).append(filename)
            if len(filenames) > 1:
                return filenames
            else:
                return list(next(iter(filenames.values())))
        else:
            return list(obj)
    elif isinstance(obj, dict):
        if split_by_organism and len(obj) > 1:
            return deepcopy(obj)
        else:
            values = obj.values()
            if not values:
                return []
            elif isinstance(next(iter(values)), list):
                return [filename for filenames in values for filename in filenames]
            elif isinstance(next(iter(values)), str):
                return list(cast(ValuesView[str], values))
            else:
                raise DataError("unexpected filenames type")
    else:
        assert isinstance(obj, str)
        if split_by_organism:
            organism = get_organism_from_filename(obj)
            if organism:
                return {organism: [obj]}
            else:
                return [obj]
        else:
            return [obj]


def get_organism_from_filename(filename: str) -> Optional[str]:
    """Get the organism from a filename.

    Try to analyse the barcode of a filename, and return the organism
    if available. Otherwise return `None`.
    """
    try:
        barcoded = BarcodedFilename(os.path.basename(filename))
        return barcoded.organism
    except Exception:
        return None


def get_samples_by_organism(
    obj: Union[List[str], Dict[str, List[str]], str], default_organism: str
) -> Dict[str, List[str]]:
    """Return the samples in a dict.

    Create a organism-samples dict.

    Args:
        obj: the samples that are collected.
        default_organism: when `obj` is not a dict, `default_organism`
                          is used as key for the output dict.
    Returns:
        A dictionary that maps organisms to lists of samples. If `obj`
        is a dict, a copy of `obj` is returned. If `obj` is a list,
        a dict with `default_organism` that maps to `obj` is returned.
        If `obj` is a string, a dict with `default_organism` that maps
        to a list of one element, `obj`, is returned.

    """
    if isinstance(obj, list):
        return {default_organism: obj}
    elif isinstance(obj, dict):
        return dict(obj)
    else:
        return {default_organism: [obj]}


def get_genome_ref_index_by_organism(config: Config, organism: str) -> Tuple[str, str]:
    """Return the reference file and the index file.

    Select the `config.*_ref` and `config.*_index` depending on
    `organism`.
    """
    if organism == "hg19":
        return (config.hg19_ref, config.hg19_index)
    elif organism == "hg38":
        return (config.hg38_ref, config.hg38_index)
    elif organism == "mm9":
        return (config.mm9_ref, config.mm9_index)
    elif organism == "mm10":
        return (config.mm10_ref, config.mm10_index)
    else:
        raise DataError("Invalid organism")


def get_dbsnp_by_organism(config: Config, organism: str) -> str:
    """Return the dbSNP filename.

    Select the `config.dbsnp_*` depending on `organism`.
    """
    if organism == "hg19":
        return config.dbsnp_hg19
    elif organism == "hg38":
        return config.dbsnp_hg38
    else:
        raise DataError("Invalid organism")


def get_cosmic_by_organism(config: Config, organism: str) -> str:
    """Return the cosmic DB filename.

    Select the `config.cosmic_*` depending on `organism`.
    """
    if organism == "hg19":
        return config.cosmic_hg19
    elif organism == "hg38":
        return config.cosmic_hg38
    else:
        raise DataError("Invalid organism")


def get_picard_max_records_string(max_records: str) -> str:
    """Get the max records string for Picard.

    Create the 'MAX_RECORDS_IN_RAM' parameter using `max_records`. If
    `max_records` is empty, an empty string is returned.
    """
    if max_records is None or max_records == "":
        return ""
    else:
        return " MAX_RECORDS_IN_RAM=%d" % int(max_records)


def find_fastqs_by_organism(
    sample: str, fastq_dir: str, default_organism: str
) -> Dict[str, List[Tuple[str, int]]]:
    """Search for FASTQ files and group them by organism.

    Find all the .fastq files inside `fastq_dir` that start with
    `sample` and have a valid suffix. Group all the files by organism.

    Args:
        sample: the barcoded sample as string.
        fastq_dir: the directory where the fastq files must be searched.
        default_organism: the organism to use in case the organism field
                          in a filename is absent.
    Returns:
        A dict that maps an organism to a list of fastq files.

    """
    re_fastq_filename = re.compile(
        r"^%s(?:\.((?:hg|mm)\d+))?\.R([12])\.fastq(?:\.gz)?$" % sample, re.I
    )
    fastq_files = [
        filename
        for filename in os.listdir(fastq_dir)
        if re_fastq_filename.match(filename)
    ]

    fastqs: Dict[str, List[Tuple[str, int]]] = {}
    for filename in fastq_files:
        match = re_fastq_filename.match(filename)
        assert match is not None

        organism = match.group(1)
        read_index = int(match.group(2))
        if organism is None or organism == "":
            organism = default_organism
        if organism in fastqs:
            fastqs[organism].append((filename, read_index))
        else:
            fastqs[organism] = [(filename, read_index)]

    return fastqs


def gzip(filename: str) -> None:
    """Compress a file with GZ compression."""
    compressed_filename = filename + ".gz"
    with open(filename, "rb") as in_fd, gz.open(
        compressed_filename, "wb", compresslevel=6
    ) as out_fd:
        shutil.copyfileobj(in_fd, out_fd)
    os.unlink(filename)


def gunzip(filename: str) -> None:
    """Decompress a GZ file."""
    decompressed_filename = filename[:-3]
    with open(decompressed_filename, "wb") as out_fd, gz.open(filename, "rb") as in_fd:
        shutil.copyfileobj(in_fd, out_fd)
    os.unlink(filename)


def check_gz(filename: str) -> bool:
    """Check if a GZ file is valid."""
    chunk_size = 2 ** 20
    with gz.open(filename, "rb") as fd:
        try:
            while fd.read(1):
                fd.seek(chunk_size, os.SEEK_CUR)

            return True
        except Exception:
            return False


def parsed_date(raw_date: str) -> str:
    """Parse a date in 'Y_M_D' format and return a std HaTSPiL date."""
    try:
        date = datetime.datetime.strptime(raw_date, "%Y_%m_%d")
    except ValueError:
        raise ArgumentTypeError("expected string in format YYYY_MM_DD")
    return "%04d_%02d_%02d" % (date.year, date.month, date.day)


def get_human_annotation(config: Config) -> str:
    """Get the best human genome annotation available in config."""
    if config.use_hg38:
        return "hg38"
    elif config.use_hg19:
        return "hg19"
    else:
        raise AnnotationError("no available human annotation in config")


def get_mouse_annotation(config: Config) -> str:
    """Get the best murine genome annotation available in config."""
    if config.use_mm10:
        return "mm10"
    elif config.use_mm9:
        return "mm9"
    else:
        raise AnnotationError("no available mouse annotation in config")


reFloat = re.compile(r"^(\d+\.\d*|\.\d+)$")
reInt = re.compile(r"^(\d+)$")


def parse_as_number(s: str) -> Union[int, float, str]:
    """Try to parse a string as number.

    If `s` matches a float format, a parsed float is returned. If `s`
    matches an int, a parset int is returned. Otherwise `s` is returned.
    """
    if reFloat.match(s):
        return float(s)
    elif reInt.match(s):
        return int(s)
    else:
        return s


T = TypeVar("T")
U = TypeVar("U")


def flatten(iterable: Iterable[Union[Iterable[T], Any]]) -> Generator[Any, None, None]:
    """Return a generator, flattening recusively an iterable object."""
    for element in iterable:
        if isinstance(element, collections.Iterable) and not isinstance(element, str):
            yield from flatten(element)
        else:
            yield element


def rfind_if(iterable: Sequence[T], fun: Callable[[T], bool]) -> Optional[int]:
    """Reverse find an object in an iterable that satisfies `fun`.

    Args:
        iterable: an iterable object.
        fun: a function that returns `True` when the item is found.
    Returns:
        The index of the first element for which `fun` returns `True`,
        performing the operation on the reversed iterable.

    """
    for index, element in enumerate(reversed(iterable)):
        if fun(element):
            return len(iterable) - index
    return None


def argmin(
    iterable: Iterable[T], key: Optional[Callable[[T], U]] = None
) -> Optional[int]:
    """Like `min`, but return the index of the element found."""
    best = min(
        ((index, element) for (index, element) in enumerate(iterable)),
        key=lambda x: key(x[1]) if key else x[1],
    )
    if best is not None:
        return best[0]
    else:
        return None


def create_logger(
    logger_name: str, handler: Optional[logging.FileHandler] = None
) -> Logger:
    """Create a named logger and add a handler to this."""
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    if handler:
        logger.addHandler(handler)

    return logger


def get_kit_from_barcoded(
    config: Config, barcoded: BarcodedFilename
) -> Optional[KitData]:
    """Get a kit from the config given a barcoded filename."""
    assert barcoded.kit is not None
    assert barcoded.analyte is not None

    return config.kits.get((barcoded.kit, barcoded.analyte))
