from .haloplex import Haloplex
from .mutect import Mutect
from .varscan import VarScan
from .analysis import Analysis
from .config import Config

import logging
import sys
import os
from multiprocessing import Pool
import functools
import smtplib
from email.mime.text import MIMEText
import traceback


def run(filename, config, root, fastq_dir):
    analysis = Analysis(filename, root, config)

    haloplex = Haloplex(analysis, fastq_dir)
    mutect = Mutect(analysis)
    varscan = VarScan(analysis)

    haloplex.run()
    mutect.run()
    varscan.run()


def main():
    usage = "Usage: %s list_file root_dir fastq_dir\n"\
            "       %s --configout config.ini"\
            % (sys.argv[0], sys.argv[0])

    if len(sys.argv) == 3:
        if sys.argv[1] == "--configout":
            config = Config()
            config.save(sys.argv[2])
            exit(0)
        else:
            print(usage)
            exit(-1)
    elif len(sys.argv) != 4:
        print(usage)
        exit(-1)

    if os.path.exists("config.ini"):
        config = Config("config.ini")
    else:
        config = Config()

    list_file = sys.argv[1]
    root = sys.argv[2]
    fastq_dir = sys.argv[3]

    if not os.path.exists(list_file):
        print("ERROR: list_file '%s' does not exist." % list_file)
        exit(-2)

    if not os.path.isdir(root):
        print("ERROR: root_dir '%s' is not a valid directory." % root)
        exit(-2)

    if not os.path.isdir(fastq_dir):
        print("ERROR: fastq_dir '%s' is not a valid directory." % fastq_dir)
        exit(-3)

    filenames = []
    with open(list_file) as fd:
        for line in fd:
            filenames.append(line.strip())

    logging.basicConfig(format="%(asctime)-15s %(message)s")
    runner = functools.partial(run, root=root, fastq_dir=fastq_dir, config=config)

    error_raised = False
    try:
        with Pool(5) as pool:
            pool.map(runner, filenames)
        msg = MIMEText("Pipeline for file list %s successfully completed." % list_file)
        msg["Subject"] = "Pipeline completed"
    except Exception:
        error_raised = True
        traceback.print_exc(file=sys.stdout)

        msg = MIMEText("Error while executing pipeline for file list %s.\nRaised exception:\n%s" % (list_file, traceback.format_exc()))
        msg["Subject"] = "Pipeline error"

    msg["From"] = "UV2000 Pipeline <pipeline@uv2000.hugef>"
    msg["To"] = config.mails

    smtp = smtplib.SMTP("localhost")
    smtp.send_message(msg)
    smtp.quit()

    if error_raised:
        exit(-1)
