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

from .barcoded_filename import BarcodedFilename
from .config import Config, KitData
from .exceptions import AnnotationError, DataError


def get_current() -> str:
    today = datetime.date.today()
    return "%04d_%02d_%02d" % (today.year, today.month, today.day)


def get_overridable_current_date(parameters: Dict[str, Any]) -> str:
    if parameters["use_date"] is None:
        return get_current()
    else:
        current_date = parameters["use_date"]
        assert isinstance(current_date, str)
        return current_date


def run_and_log(command: str, logger: Logger) -> int:
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
    try:
        barcoded = BarcodedFilename(os.path.basename(filename))
        return barcoded.organism
    except Exception:
        return None


def get_samples_by_organism(
    obj: Union[List[str], Dict[str, List[str]], str], default_organism: str
) -> Dict[str, List[str]]:

    if isinstance(obj, list):
        return {default_organism: obj}
    elif isinstance(obj, dict):
        return dict(obj)
    else:
        return {default_organism: [obj]}


def get_genome_ref_index_by_organism(config: Config, organism: str) -> Tuple[str, str]:
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
    if organism == "hg19":
        return config.dbsnp138_hg19
    elif organism == "hg38":
        return config.dbsnp138_hg38
    else:
        raise DataError("Invalid organism")


def get_cosmic_by_organism(config: Config, organism: str) -> str:
    if organism == "hg19":
        return config.cosmic_hg19
    elif organism == "hg38":
        return config.cosmic_hg38
    else:
        raise DataError("Invalid organism")


def get_picard_max_records_string(max_records: str) -> str:
    if max_records is None or max_records == "":
        return ""
    else:
        return " MAX_RECORDS_IN_RAM=%d" % int(max_records)


def find_fastqs_by_organism(
    sample: str, fastq_dir: str, default_organism: str
) -> Dict[str, List[Tuple[str, int]]]:

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
    compressed_filename = filename + ".gz"
    with open(filename, "rb") as in_fd, gz.open(
        compressed_filename, "wb", compresslevel=6
    ) as out_fd:
        shutil.copyfileobj(in_fd, out_fd)
    os.unlink(filename)


def gunzip(filename: str) -> None:
    decompressed_filename = filename[:-3]
    with open(decompressed_filename, "wb") as out_fd, gz.open(filename, "rb") as in_fd:
        shutil.copyfileobj(in_fd, out_fd)
    os.unlink(filename)


def check_gz(filename: str) -> bool:
    chunk_size = 2 ** 20
    with gz.open(filename, "rb") as fd:
        try:
            while fd.read(1):
                fd.seek(chunk_size, os.SEEK_CUR)

            return True
        except Exception:
            return False


def parsed_date(raw_date: str) -> str:
    try:
        date = datetime.datetime.strptime(raw_date, "%Y_%m_%d")
    except ValueError:
        raise ArgumentTypeError("expected string in format YYYY_MM_DD")
    return "%04d_%02d_%02d" % (date.year, date.month, date.day)


def get_human_annotation(config: Config) -> str:
    if config.use_hg38:
        return "hg38"
    elif config.use_hg19:
        return "hg19"
    else:
        raise AnnotationError("no available human annotation in config")


def get_mouse_annotation(config: Config) -> str:
    if config.use_mm10:
        return "mm10"
    elif config.use_mm9:
        return "mm9"
    else:
        raise AnnotationError("no available mouse annotation in config")


reFloat = re.compile(r"^(\d+\.\d*|\.\d+)$")
reInt = re.compile(r"^(\d+)$")


def parse_as_number(s: str) -> Union[int, float, str]:
    if reFloat.match(s):
        return float(s)
    elif reInt.match(s):
        return int(s)
    else:
        return s


T = TypeVar("T")
U = TypeVar("U")


def flatten(iterable: Iterable[Union[Iterable[T], Any]]) -> Generator[Any, None, None]:
    for element in iterable:
        if isinstance(element, collections.Iterable) and not isinstance(element, str):
            yield from flatten(element)
        else:
            yield element


def rfind_if(iterable: Sequence[T], fun: Callable[[T], bool]) -> Optional[int]:
    for index, element in enumerate(reversed(iterable)):
        if fun(element):
            return len(iterable) - index
    return None


def argmin(
    iterable: Iterable[T], key: Optional[Callable[[T], U]] = None
) -> Optional[int]:
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
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    if handler:
        logger.addHandler(handler)

    return logger


def get_kit_from_barcoded(
    config: Config, barcoded: BarcodedFilename
) -> Optional[KitData]:
    assert barcoded.kit is not None
    assert barcoded.analyte is not None

    return config.kits.get((barcoded.kit, barcoded.analyte))