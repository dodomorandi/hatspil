from .config import Config
from . import utils
from .barcoded_filename import BarcodedFilename, Tissue
from .runner import Runner

import logging
import sys
import os
from multiprocessing import Pool, Manager
import smtplib
from email.mime.text import MIMEText
import traceback
import argparse
import re


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
    parser.add_argument("--no-cutadapt", dest="use_cutadapt",
                        action="store_false", help="Skips cutadapt.")
    parser.add_argument("--no-tdf", dest="use_tdf",
                        action="store_false", help="Skips tdf generation.")
    parser.add_argument("--no-R-checks", dest="r_checks",
                        action="store_false",
                        help="Skips some R dependency checks. If omitted, "
                        "the program will check some R depencencies and, "
                        "if some packages are found missing, it will try to "
                        "install them.")
    parser.add_argument("--dont-use-normals", action="store_false",
                        dest="use_normals",
                        help="Normally, whenever a normal sample is found, it is used. "
                        "In this case many phases of the analysis are "
                        "performed using different parameters. If this option "
                        "is passed, these phases are skipped.")
    parser.add_argument("--dont-mark-duplicates", action="store_false",
                        dest="mark_duplicates",
                        help="Do not mark PCR duplicates during mapping phase "
                        "for xenograft tissues.")
    parser.add_argument("--no-xenome", action="store_false", dest="use_xenome",
                        help="Do not use xenograft analysis using xenome even "
                        "if the tissue has the xenograft flag")
    parser.add_argument("--post-recalibration", action="store_true",
                        help="Perform a post-recalibration analysis after the "
                        "basic recalibration.")
    parser.add_argument("--compress-fastq", action="store_true",
                        help="If set, the fastqs files are compressed at the "
                        "end of the mapping phase.")
    parser.add_argument("--trim-5", action="store", type=int, metavar="n_bases",
                        default=5,
                        help="Bases that will be trimmed at 5' (default=5)")
    parser.add_argument("--trim-3", action="store", type=int, metavar="n_bases",
                        default=None, help="Bases that will be trimmed at 3' "
                        "(default=10 when --xenograft is passed, 0 otherwise)")
    parser.add_argument("--gatk-threads", metavar="n",
                        action="store", type=int, default=20,
                        help="Number of threads to be used for GATK. "
                        "Default=20.")
    parser.add_argument("--picard-max-records", action="store",
                        help="Sets the MAX_RECORDS_IN_RAM for Picard. "
                        "If unspecified, the parameter is not passed.")
    parser.add_argument("--use-date", action="store", default=None,
                        type=utils.parsed_date, metavar="YYYY_MM_DD",
                        dest="use_date",
                        help="Use the specified date instead of the current "
                        "one")
    parser.add_argument("--skip-mapping", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--only-mapping", action="store_true",
                        help=argparse.SUPPRESS)

    list_file_group = parser.add_mutually_exclusive_group(required=False)
    list_file_group.add_argument("--list-file", action="store",
                                 help="The name of the file containing "
                                 "the name of the samples, one by line.")
    list_file_group.add_argument("--scan-samples", action="store_true",
                                 help="Scan for sample files instead of "
                                 "reading them from a file")
    parser.add_argument("--root-dir", action="store",
                        help="The root directory for the analysis",
                        )
    parser.add_argument("--fastq-dir", action="store",
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
        print("--list-file (or --scan-samples), --root-dir and --fastq-dir "
              "are mandatory unless --configout is specified")
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

    if not config.check_programs():
        print("ERROR: not every program inside configuration is correctly "
              "set. Please fix the problem and try again.")
        exit(-1)

    if not config.check_files():
        print("ERROR: not every file inside configuration is correctly "
              "set. Please fix the problem and try again.")
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

    if config.use_mongodb:
        try:
            from pymongo import MongoClient
        except:
            print("ERROR: use_mongodb is set to true inside config file but "
                  "pymongo is not installed. Please install it using pip3.")
            exit(-4)

        try:
            client = MongoClient(config.mongodb_host, config.mongodb_port)
        except:
            print("ERROR: cannot connect to MongoDB. Please check the config "
                  "file under the section MONGODB")
            exit(-4)

        db = client[config.mongodb_database]
        try:
            db.authenticate(config.mongodb_username,
                            config.mongodb_password)
        except:
            print("ERROR: MongoDB authentication failed. Please check the "
                  "config file under the section MONGODB")
            exit(-4)

    if args.list_file:
        re_pattern = re.compile(R"^([^-]+)(?:-([^-]+)(?:-(\d[0-9A-Za-z]|\*)"
                                R"(?:-(\d|\*)(?:(\d|\*)(\d|\*)?)?"
                                R"(?:-(\d|\*)(\d|\*)(\d|\*)?)?)?)?)?$")

        all_filenames = os.listdir(args.fastq_dir)
        filenames = set()
        with open(args.list_file) as fd:
            for line_index, line in enumerate(fd):
                match = re_pattern.match(line.strip())
                if not match:
                    print("Invalid file at line %d of file list"
                          % (line_index + 1))
                    exit(-1)

                current_pattern = R"^("
                current_pattern += match.group(1).replace("*", R"[^-]")
                group = match.group(2)
                if group:
                    current_pattern += "-"
                    current_pattern += group.replace("*", R"[^-]")

                    group = match.group(3)
                    if group:
                        current_pattern += "-"
                        current_pattern += group.replace("*", R"\d{2}")

                        group = match.group(4)
                        if group:
                            current_pattern += "-"
                            current_pattern += match.group(4).replace("*",
                                                                      R"\d")

                            group = match.group(5)
                            if group:
                                current_pattern += group.replace("*", R"\d")

                                group = match.group(6)
                                if group:
                                    current_pattern += group.replace("*",
                                                                     R"\d")

                                    group = match.group(7)
                                    if group:
                                        current_pattern += "-"
                                        current_pattern += group.replace("*", R"\d")

                                        group = match.group(8)
                                        if group:
                                            current_pattern += group.replace("*", R"\d")

                                            group = match.group(9)
                                            if group:
                                                current_pattern += group.replace("*", R"\d")
                                            else:
                                                current_pattern += "R\d"
                                        else:
                                            current_pattern += R"\d{2}"
                                    else:
                                        current_pattern += R"\d{3}"
                                else:
                                    current_pattern += R"\d-\d{3}"
                            else:
                                current_pattern += R"\d{2}-\d{3}"
                        else:
                            current_pattern += R"-\d{3}-\d{3}"
                    else:
                        current_pattern += R"-\d[0-9A-Za-z]-\d{3}-\d{3}"
                else:
                    current_pattern += R"-[^-]+-\d[0-9A-Za-z]-\d{3}-\d{3}"
                current_pattern += R")(?:\.(?:hg|mm)\d+)?(?:\.R[12])?\.fastq(\.gz)?$"

                re_current_pattern = re.compile(current_pattern, re.I)
                added_files = 0
                for filename in all_filenames:
                    match = re_current_pattern.match(
                        os.path.basename(filename))
                    if match:
                        barcoded_filename = BarcodedFilename(filename)
                        if (barcoded_filename.tissue != Tissue.PRIMARY_XENOGRAFT_TISSUE and
                            barcoded_filename.tissue != Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE) or \
                                barcoded_filename.organism is None:
                            filenames.add(match.group(1))
                            added_files += 1

                if added_files == 0:
                    print("ERROR: cannot find any file for sample %s"
                          % line.strip())
                    exit(-1)

    elif args.scan_samples:
        fastq_files = [
            filename
            for filename in os.listdir
            (args.fastq_dir)
            if re.search(R"\.fastq$", filename, re.I)]
        filenames = list(
            set(
                [re.sub
                 (R"(\.R[12])?\.fastq$", "", filename, flags=re.I)
                    for filename in
                    fastq_files
                    if not re.search(R"trimmed|clipped", filename, re.I)]))
    else:
        raise RuntimeError("Unhandled condition")

    parameters = {
        "use_xenome": args.use_xenome,
        "use_cutadapt": args.use_cutadapt,
        "mark_duplicates": args.mark_duplicates,
        "run_post_recalibration": args.post_recalibration,
        "compress_fastq": args.compress_fastq,
        "gatk_threads": args.gatk_threads,
        "picard_max_records": args.picard_max_records,
        "use_normals": args.use_normals,
        "trim_5": args.trim_5,
        "trim_3": args.trim_3,
        "use_date": args.use_date,
        "skip_mapping": args.skip_mapping,
        "only_mapping": args.only_mapping,
        "use_tdf": args.use_tdf
    }
    logging.basicConfig(format="%(asctime)-15s %(message)s")

    if args.r_checks and args.post_recalibration:
        try:
            import rpy2.robjects.packages as rpackages
            from rpy2.robjects.vectors import StrVector
        except:
            rpackages = None
            print("cannot correctly import rpy2. "
                  "R checks are skipped.")

        if rpackages:
            print("Checking R packages and, eventually, performing "
                  "installations")
            dependencies = ("ggplot2", "gplots", "reshape", "grid", "tools",
                            "gsalib")
            rutils = rpackages.importr("utils")
            base = rpackages.importr("base")
            rutils.chooseCRANmirror(ind=1)
            installed_packages = rutils.installed_packages().rx(True, 1)
            for package in dependencies:
                if not base.is_element(package, installed_packages)[0]:
                    sys.stdout.write("Installing R package %s..." % package)
                    rutils.install_packages(
                        StrVector(dependencies),
                        quiet=True)
                    print(" done.")

            installed_packages = rutils.installed_packages().rx(True, 1)
            for package in dependencies:
                if not base.is_element(package, installed_packages)[0]:
                    print("Package %s has not been correctly installed. "
                          "Try again or check this dependency manually."
                          % package)
                    exit(-1)
            print("Done with R packages")

    manager = Manager()
    runner = Runner(
        manager,
        root=args.root_dir,
        fastq_dir=args.fastq_dir,
        parameters=parameters,
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

    if args.use_normals and not args.only_mapping:
        samples = {}
        normals = {}
        last_operations = {}
        for sample, last_operation in runner.last_operations.items():
            last_operations[sample] = last_operation
            for filename in utils.get_sample_filenames(last_operation):
                barcoded_filename = BarcodedFilename(filename)
                fake_sample = "%s-%s-%d%d%d-%d%d%d" % (
                    barcoded_filename.project,
                    barcoded_filename.patient,
                    barcoded_filename.molecule,
                    barcoded_filename.analyte,
                    barcoded_filename.kit,
                    barcoded_filename.biopsy,
                    barcoded_filename.get_sample_index(),
                    barcoded_filename.sequencing)

                if barcoded_filename.tissue == Tissue.BONE_MARROW_NORMAL or \
                        barcoded_filename.tissue == Tissue.BUCCAL_CELL_NORMAL or \
                        barcoded_filename.tissue == Tissue.SOLID_TISSUE_NORMAL or \
                        barcoded_filename.tissue == Tissue.BLOOD_DERIVED_NORMAL or \
                        barcoded_filename.tissue == Tissue.EBV_IMMORTALIZED_NORMAL:
                    sample_type = "normal"
                else:
                    sample_type = "tumor"

                if fake_sample not in samples:
                    samples[fake_sample] = {}

                if sample_type == "normal":
                    samples[fake_sample]["normal"] = filename
                    normals[filename] = barcoded_filename
                else:
                    if "tumor" not in samples[fake_sample]:
                        samples[fake_sample]["tumor"] = []
                    samples[fake_sample]["tumor"].append((filename, sample))

        samples_with_no_normal = {sample: filenames["tumor"][0][0]
                                  for sample, filenames
                                  in samples.items()
                                  if "normal" not in filenames and
                                  len(filenames["tumor"]) > 0}

        for sample, filename in samples_with_no_normal.items():
            sample_barcode = BarcodedFilename(filename)
            candidates = {filename: barcode
                          for filename, barcode
                          in normals.items()
                          if barcode.project == sample_barcode.project and
                          barcode.patient == sample_barcode.patient and
                          barcode.molecule == sample_barcode.molecule and
                          barcode.analyte == sample_barcode.analyte and
                          barcode.kit == sample_barcode.kit and
                          barcode.biopsy == sample_barcode.biopsy and
                          barcode.sample == sample_barcode.sample}

            if len(candidates) == 0:
                candidates = {filename: barcode
                              for filename, barcode
                              in normals.items()
                              if barcode.project == sample_barcode.project and
                              barcode.patient == sample_barcode.patient and
                              barcode.molecule == sample_barcode.molecule and
                              barcode.analyte == sample_barcode.analyte and
                              barcode.kit == sample_barcode.kit and
                              barcode.biopsy == sample_barcode.biopsy}

            if len(candidates) == 0:
                candidates = {filename: barcode
                              for filename, barcode
                              in normals.items()
                              if barcode.project == sample_barcode.project and
                              barcode.patient == sample_barcode.patient and
                              barcode.molecule == sample_barcode.molecule and
                              barcode.analyte == sample_barcode.analyte and
                              barcode.kit == sample_barcode.kit}

            if len(candidates) == 0:
                candidates = {filename: barcode
                              for filename, barcode
                              in normals.items()
                              if barcode.project == sample_barcode.project and
                              barcode.patient == sample_barcode.patient and
                              barcode.molecule == sample_barcode.molecule and
                              barcode.analyte == sample_barcode.analyte}

            if len(candidates) == 0:
                candidates = {filename: barcode
                              for filename, barcode
                              in normals.items()
                              if barcode.project == sample_barcode.project and
                              barcode.patient == sample_barcode.patient and
                              barcode.molecule == sample_barcode.molecule}

            if len(candidates) == 1:
                samples[sample]["normal"] = list(candidates.items())[0][0]
            elif len(candidates) > 1:
                candidates = list(candidates.items())
                candidates.sort(key=lambda x: os.stat(x[0]).st_size)
                samples[sample]["normal"] = candidates[-1][0]

        triplets = []
        for sample, values in samples.items():
            if len(values) < 2:
                continue
            for tumor in values["tumor"]:
                triplets.append((tumor[1], tumor[0], values["normal"]))

        error_raised = False
        try:
            with Pool(5) as pool:
                pool.starmap(runner.with_normals, triplets)

            if args.mail:
                msg = MIMEText(
                    "Pipeline for file list %s successfully completed." %
                    args.list_file)
                msg["Subject"] = "Pipeline completed"

        except RuntimeError:
            error_raised = True
            traceback.print_exc(file=sys.stdout)

            if args.mail:
                msg = MIMEText(
                    "Error while executing pipeline for file list %s.\n"
                    "Raised exception:\n%s" %
                    (args.list_file, traceback.format_exc()))
                msg["Subject"] = "Pipeline error"

    if args.mail and len(config.mails) > 0:
        msg["From"] = "UV2000 Pipeline <pipeline@uv2000.hugef>"
        msg["To"] = config.mails

        smtp = smtplib.SMTP("localhost")
        smtp.send_message(msg)
        smtp.quit()

    if error_raised:
        exit(-1)
