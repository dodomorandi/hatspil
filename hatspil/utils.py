import datetime
import subprocess
import re
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
