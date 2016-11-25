import datetime
import subprocess
import re
import os
from .exceptions import PipelineError


def get_current():
    today = datetime.date.today()
    return "%04d_%02d_%02d" % (today.year, today.month, today.day)


def run_and_log(command, logger):
    with subprocess.Popen(command,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True,
                          universal_newlines=True,
                          bufsize=1) as process:
        (out, err) = process.communicate()

        for line in out.split("\n"):
            if line != "":
                logger.info(line)

        for line in err.split("\n"):
            if line != "":
                logger.warning(line)

        return process.wait()


def get_read_index(filename):
    match = re.search(R"[._]R([12])[._]", filename)
    if not match:
        raise PipelineError("cannot find R1/R2 pattern in filename")
    return int(match.group(1))


def get_params_from_filename(filename, sample):
    filename = os.path.basename(filename)
    match = re.match(
        R"^%s(?:_((?:hg|mm)\d+))?.*[._]R([12])\.(fastq|bam|bai)$" % sample,
        filename,
        re.I)
    if not match:
        print(filename)
        raise RuntimeError(
            "Cannot parse filename. Filename in a strange format.")

    return (match.group(1), int(match.group(2)), match.group(3))


def get_sample_filenames(obj):
    if isinstance(obj, list):
        return obj
    elif isinstance(obj, dict):
        return [
            filename for filenames in obj.values() for filename in filenames]
    else:
        return [obj]


def get_samples_by_organism(obj, default_organism="hg19"):
    if isinstance(obj, list):
        return {default_organism: obj}
    elif isinstance(obj, dict):
        return obj
    else:
        return {default_organism: [obj]}


def get_genome_ref_index_by_organism(config, organism):
    if organism == "hg19":
        return (config.hg19_ref, config.hg19_index)
    elif organism == "hg38":
        return (config.hg38_ref, config.h38_index)
    elif organism == "mm9":
        return (config.mm9_ref, config.mm9_index)
    elif organism == "mm10":
        return (config.mm10_ref, config.mm10_index)
    else:
        raise RuntimeError("Invalid organism")


def get_dbsnp_by_organism(config, organism):
    if organism == "hg19":
        return config.dbsnp138_hg19
    elif organism == "hg38":
        return config.dbsnp138_hg38
    else:
        raise RuntimeError("Invalid organism")


def get_cosmic_by_organism(config, organism):
    if organism == "hg19":
        return config.cosmic_hg19
    elif organism == "hg38":
        return config.cosmic_hg38
    else:
        raise RuntimeError("Invalid organism")


def get_picard_max_records_string(max_records):
    if max_records is None or max_records == "":
        return ""
    else:
        return " MAX_RECORDS_IN_RAM=%d" % int(max_records)
