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
import argparse
import re


def run(filename, config, root, fastq_dir):
    analysis = Analysis(filename, root, config)

    haloplex = Haloplex(analysis, fastq_dir)
    mutect = Mutect(analysis)
    varscan = VarScan(analysis)

    haloplex.run()
    mutect.run()
    varscan.run()


def get_parser():
    parser = argparse.ArgumentParser(
        description="Makes your life easier when performing some HTS "
                    "analysis.")
    parser.add_argument("--configout", action="store",
                        metavar="filename",
                        help="Dumps a default (almost empty configuration) in "
                        "a file.\nWhen this option is passed, any other "
                        "option will be ignored and the program will exit "
                        "after the file is being written.")
    parser.add_argument("--config", "-c", action="store",
                        metavar="config.ini",
                        help="Select the configuration file. If it is not "
                        "specified, the program will try to search for a file "
                        "called 'config.ini' in the current working "
                        "directory. If it is not available, an error will be "
                        "raised.")
    parser.add_argument("--no-mail", dest="mail", action="store_false",
                        help="Do not send emails.")

    list_file_group = parser.add_mutually_exclusive_group(required=False)
    list_file_group.add_argument("--list_file", action="store",
                                 help="The name of the file containing "
                                 "the name of the samples, one by line.")
    list_file_group.add_argument("--scan-samples", action="store_true",
                                 help="Scan for sample files instead of "
                                 "reading them from a file")
    parser.add_argument("--root_dir", action="store",
                        help="The root directory for the analysis",
                        )
    parser.add_argument("--fastq_dir", action="store",
                        help="The directory where the fastq files of the "
                        "samples are located.")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.configout is not None:
        config = Config()
        config.save(args.configout)
        print("Sample config written to '%s'" % args.configout)
        exit(0)
    elif (not args.list_file and not args.scan_samples) or \
            not args.root_dir or not args.fastq_dir:
        print("list_file (or --scan-samples), root_dir and fastq_dir are "
              "mandatory unless --configout is specified")
        parser.print_usage()
        exit(-1)

    if args.config is not None:
        if os.path.exists(args.config):
            config = Config(args.config)
        else:
            print("ERROR: config file '%s' does not exist." % args.config)
            exit(-1)
    else:
        if os.path.exists("config.ini"):
            config = Config("config.ini")
        else:
            print("ERROR: config file 'config.ini' does not exist in the "
                  "current directory.\nPlease use '--configout' option to "
                  "create a config file, then specify it with the '--config' "
                  "option or just name it 'config.ini' and put in the current "
                  "directory.")
            exit(-1)

    if args.list_file is not None and not os.path.exists(args.list_file):
        print("ERROR: list_file '%s' does not exist." % args.list_file)
        exit(-2)

    if not os.path.isdir(args.root_dir):
        print("ERROR: root_dir '%s' is not a valid directory." % args.root_dir)
        exit(-2)

    if not os.path.isdir(args.fastq_dir):
        print(
            "ERROR: fastq_dir '%s' is not a valid directory." %
            args.fastq_dir)
        exit(-3)

    if args.list_file:
        filenames = []
        with open(args.list_file) as fd:
            for line in fd:
                filenames.append(line.strip())
    elif args.scan_samples:
        fastq_files = [
            filename
            for filename in os.listdir
            (args.fastq_dir)
            if re.search(R"\.fastq$", filename, re.I)]
        filenames = list(
            set(
                [re.sub
                 (R"([._]R[12])?\.fastq$", "", filename, flags=re.I)
                    for filename in
                    fastq_files
                    if not re.search(R"trimmed|clipped", filename, re.I)]))
    else:
        raise RuntimeError("Unhandled condition")

    logging.basicConfig(format="%(asctime)-15s %(message)s")
    runner = functools.partial(
        run,
        root=args.root_dir,
        fastq_dir=args.fastq_dir,
        config=config)

    error_raised = False
    try:
        with Pool(5) as pool:
            pool.map(runner, filenames)

        if args.mail:
            msg = MIMEText(
                "Pipeline for file list %s successfully completed." %
                args.list_file)
            msg["Subject"] = "Pipeline completed"
    except Exception:
        error_raised = True
        traceback.print_exc(file=sys.stdout)

        if args.mail:
            msg = MIMEText(
                "Error while executing pipeline for file list %s.\n"
                "Raised exception:\n%s" %
                (args.list_file, traceback.format_exc()))
            msg["Subject"] = "Pipeline error"

    if args.mail:
        msg["From"] = "UV2000 Pipeline <pipeline@uv2000.hugef>"
        msg["To"] = config.mails

        smtp = smtplib.SMTP("localhost")
        smtp.send_message(msg)
        smtp.quit()

    if error_raised:
        exit(-1)
