import argparse
import itertools
import logging
import os
import re
import shutil
import smtplib
import sys
import traceback
from email.mime.text import MIMEText
from enum import Enum
from multiprocessing import Manager, Pool
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    MutableMapping,
    Optional,
    Tuple,
    Union,
    cast,
)

from . import utils
from .barcoded_filename import Analyte, BarcodedFilename, Tissue
from .config import Config
from .aligner import GenericAligner, RnaSeqAligner
from .runner import Runner
from .xenograft import XenograftClassifier, Xenograft


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Makes your life easier when performing some HTS analysis."
    )
    parser.add_argument(
        "--aligner",
        action="store",
        dest="aligner",
        choices=[aligner.name.lower() for aligner in GenericAligner] + ["auto"],
        default="auto",
        help="Select the aligner. When this option is set to "
        "'auto' will be used the first aligner available",
    )
    parser.add_argument(
        "--rnaseq-aligner",
        action="store",
        dest="rnaseq_aligner",
        choices=[aligner.name.lower() for aligner in RnaSeqAligner] + ["auto"],
        default="auto",
        help="Select the aligner for RNA-Seq data. When this option is set to "
        "'auto' will be used the first aligner available",
    )
    parser.add_argument(
        "--xenograft-classifier",
        action="store",
        dest="xenograft_classifier",
        choices=[classifier.name.lower() for classifier in XenograftClassifier]
        + ["auto"],
        default="auto",
        help="Select the xenograft classifier, the software used to "
        "distinguish the reads belonging to different organisms. When this "
        "option is set to 'auto' will be used the first classifier "
        "available",
    )
    parser.add_argument(
        "--configout",
        action="store",
        metavar="filename",
        help="Dumps a default (almost empty configuration) in "
        "a file.\nWhen this option is passed, any other "
        "option will be ignored and the program will exit "
        "after the file is being written.",
    )
    parser.add_argument(
        "--config",
        "-c",
        action="store",
        metavar="config.ini",
        help="Select the configuration file. If it is not "
        "specified, the program will try to search for a file "
        "called 'config.ini' in the current working "
        "directory. If it is not available, an error will be "
        "raised.",
    )
    parser.add_argument(
        "--no-mail", dest="mail", action="store_false", help="Do not send emails."
    )
    parser.add_argument(
        "--no-cutadapt",
        dest="use_cutadapt",
        action="store_false",
        help="Skips cutadapt.",
    )
    parser.add_argument(
        "--no-tdf", dest="use_tdf", action="store_false", help="Skips tdf generation."
    )
    parser.add_argument(
        "--no-R-checks",
        dest="r_checks",
        action="store_false",
        help="Skips some R dependency checks. If omitted, "
        "the program will check some R depencencies and, "
        "if some packages are found missing, it will try to "
        "install them.",
    )
    parser.add_argument(
        "--dont-use-normals",
        action="store_false",
        dest="use_normals",
        help="Normally, whenever a normal sample is found, it is used. "
        "In this case many phases of the analysis are "
        "performed using different parameters. If this option "
        "is passed, these phases are skipped.",
    )
    parser.add_argument(
        "--dont-mark-duplicates",
        action="store_false",
        dest="mark_duplicates",
        help="Do not mark PCR duplicates during mapping phase "
        "for xenograft tissues.",
    )
    parser.add_argument(
        "--skip-xenograft-classifier",
        action="store_false",
        dest="use_xenograft_classifier",
        help="Do not use xenograft classifiers even with xenograft tissues",
    )
    parser.add_argument(
        "--post-recalibration",
        action="store_true",
        help="Perform a post-recalibration analysis after the basic recalibration.",
    )
    parser.add_argument(
        "--compress-fastq",
        action="store_true",
        help="If set, the fastqs files are compressed at the "
        "end of the mapping phase.",
    )
    parser.add_argument(
        "--trim-5",
        action="store",
        type=int,
        metavar="n_bases",
        default=5,
        help="Bases that will be trimmed at 5' (default=5)",
    )
    parser.add_argument(
        "--trim-3",
        action="store",
        type=int,
        metavar="n_bases",
        default=None,
        help="Bases that will be trimmed at 3' "
        "(default=10 when --xenograft is passed, 0 otherwise)",
    )
    parser.add_argument(
        "--gatk-threads",
        metavar="n",
        action="store",
        type=int,
        default=20,
        help="Number of threads to be used for GATK. Default=20.",
    )
    parser.add_argument(
        "--picard-max-records",
        action="store",
        help="Sets the MAX_RECORDS_IN_RAM for Picard. "
        "If unspecified, the parameter is not passed.",
    )
    parser.add_argument(
        "--use-date",
        action="store",
        default=None,
        type=utils.parsed_date,
        metavar="YYYY_MM_DD",
        dest="use_date",
        help="Use the specified date instead of the current one",
    )
    parser.add_argument("--skip-mapping", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--only-mapping", action="store_true", help=argparse.SUPPRESS)

    list_file_group = parser.add_mutually_exclusive_group(required=False)
    list_file_group.add_argument(
        "--list-file",
        action="store",
        help="The name of the file containing the name of the samples, one by line.",
    )
    list_file_group.add_argument(
        "--scan-samples",
        action="store_true",
        help="Scan for sample files instead of reading them from a file",
    )
    parser.add_argument(
        "--root-dir", action="store", help="The root directory for the analysis"
    )
    parser.add_argument(
        "--fastq-dir",
        action="store",
        help="The directory where the fastq files of the samples are located.",
    )

    return parser


def set_aligner_param(
    args: argparse.Namespace,
    parameters: MutableMapping[str, Any],
    param_name: str,
    available_aligners: Iterable[Enum],
    files_checkers: Iterable[Optional[Callable[[Config], bool]]],
    config: Config,
) -> None:

    type_aligner = getattr(args, param_name)
    if type_aligner == "auto" or type_aligner is None:
        for aligner, files_checker in zip(available_aligners, files_checkers):
            aligner_exec = getattr(config, aligner.name.lower())
            if aligner_exec is not None and shutil.which(aligner_exec) is not None:
                if files_checker is not None and not files_checker(config):
                    continue

                parameters[param_name] = aligner
                break

        if param_name not in parameters:
            print(
                "No valid aligner is available. "
                "Please check your configuration file.\n"
                "One of the following parameters must be valid: "
                + " ".join(aligner.name.lower() for aligner in available_aligners),
                file=sys.stderr,
            )
            exit(-5)
    else:
        aligner_names = [aligner.name.lower() for aligner in available_aligners]
        aligner_lower = type_aligner.lower()
        if aligner_lower not in aligner_names:
            print(
                "The aligner %s is not present in the available ones" % type_aligner,
                file=sys.stderr,
            )
            exit(-5)
        aligner_index = aligner_names.index(aligner_lower)
        files_checker = next(itertools.islice(files_checkers, aligner_index, None))

        if files_checker is not None and not files_checker(config):
            print(
                "The aligner %s does not have all the needed config "
                "parameters correctly set" % type_aligner
            )
            exit(-5)

        aligner_exec = getattr(config, aligner_lower)
        if aligner_exec is not None and shutil.which(aligner_exec) is not None:
            if param_name == "aligner":
                parameters[param_name] = GenericAligner[type_aligner.upper()]
            else:
                parameters[param_name] = RnaSeqAligner[type_aligner.upper()]
        else:
            print("The chosen aligner is not executable", file=sys.stderr)
            exit(-5)


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    if args.configout is not None:
        config = Config()
        config.save(args.configout)
        print("Sample config written to '%s'" % args.configout)
        exit(0)
    elif (
        (not args.list_file and not args.scan_samples)
        or not args.root_dir
        or not args.fastq_dir
    ):
        print(
            "--list-file (or --scan-samples), --root-dir and --fastq-dir "
            "are mandatory unless --configout is specified"
        )
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
            print(
                "ERROR: config file 'config.ini' does not exist in the "
                "current directory.\nPlease use '--configout' option to "
                "create a config file, then specify it with the '--config' "
                "option or just name it 'config.ini' and put in the current "
                "directory."
            )
            exit(-1)

    if not config.check_programs():
        print(
            "ERROR: not every program inside configuration is correctly "
            "set. Please fix the problem and try again."
        )
        exit(-1)

    if not config.check_files():
        print(
            "ERROR: not every file inside configuration is correctly "
            "set. Please fix the problem and try again."
        )
        exit(-1)

    if args.list_file is not None and not os.path.exists(args.list_file):
        print("ERROR: list_file '%s' does not exist." % args.list_file)
        exit(-2)

    if not os.path.isdir(args.root_dir):
        print("ERROR: root_dir '%s' is not a valid directory." % args.root_dir)
        exit(-2)

    if not os.path.isdir(args.fastq_dir):
        print("ERROR: fastq_dir '%s' is not a valid directory." % args.fastq_dir)
        exit(-3)

    if config.use_mongodb:
        try:
            from pymongo import MongoClient
        except ImportError:
            print(
                "ERROR: use_mongodb is set to true inside config file but "
                "pymongo is not installed. Please install it using pip3."
            )
            exit(-4)

        try:
            client = MongoClient(config.mongodb_host, config.mongodb_port)
        except Exception:
            print(
                "ERROR: cannot connect to MongoDB. Please check the config "
                "file under the section MONGODB"
            )
            exit(-4)

        db = client[config.mongodb_database]
        try:
            db.authenticate(config.mongodb_username, config.mongodb_password)
        except Exception:
            print(
                "ERROR: MongoDB authentication failed. Please check the "
                "config file under the section MONGODB"
            )
            exit(-4)

    args.fastq_dir = os.path.abspath(args.fastq_dir)
    args.root_dir = os.path.abspath(args.root_dir)

    if args.list_file:
        re_pattern = re.compile(
            r"^([^-]+)(?:-([^-]+)(?:-(\d[0-9A-Za-z]|\*)"
            r"(?:-(\d|\*)(?:(\d|\*)(\d|\*)?)?"
            r"(?:-(\d|\*)(\d|\*)(\d|\*)?)?)?)?)?$"
        )

        all_filenames = os.listdir(args.fastq_dir)
        filenames_set = set()
        with open(args.list_file) as fd:
            for line_index, line in enumerate(fd):
                match = re_pattern.match(line.strip())
                if not match:
                    print("Invalid file at line %d of file list" % (line_index + 1))
                    exit(-1)

                current_pattern = r"^("
                current_pattern += match.group(1).replace("*", r"[^-]+")
                group = match.group(2)
                if group:
                    current_pattern += "-"
                    current_pattern += group.replace("*", r"[^-]+")

                    group = match.group(3)
                    if group:
                        current_pattern += "-"
                        current_pattern += group.replace("*", r"\d[0-9A-Za-z]")

                        group = match.group(4)
                        if group:
                            current_pattern += "-"
                            current_pattern += match.group(4).replace("*", r"\d")

                            group = match.group(5)
                            if group:
                                current_pattern += group.replace("*", r"\d")

                                group = match.group(6)
                                if group:
                                    current_pattern += group.replace("*", r"\d")

                                    group = match.group(7)
                                    if group:
                                        current_pattern += "-"
                                        current_pattern += group.replace("*", r"\d")

                                        group = match.group(8)
                                        if group:
                                            current_pattern += group.replace("*", r"\d")

                                            group = match.group(9)
                                            if group:
                                                current_pattern += group.replace(
                                                    "*", r"\d"
                                                )
                                            else:
                                                current_pattern += r"\d"
                                        else:
                                            current_pattern += r"\d{2}"
                                    else:
                                        current_pattern += r"-\d{3}"
                                else:
                                    current_pattern += r"\d-\d{3}"
                            else:
                                current_pattern += r"\d{2}-\d{3}"
                        else:
                            current_pattern += r"-\d{3}-\d{3}"
                    else:
                        current_pattern += r"-\d[0-9A-Za-z]-\d{3}-\d{3}"
                else:
                    current_pattern += r"-[^-]+-\d[0-9A-Za-z]-\d{3}-\d{3}"
                current_pattern += r")(?:\.(?:hg|mm)\d+)?(?:\.R[12])?\.fastq(\.gz)?$"

                re_current_pattern = re.compile(current_pattern, re.I)
                added_files = 0
                for filename in all_filenames:
                    match = re_current_pattern.match(os.path.basename(filename))
                    if match:
                        barcoded_filename = BarcodedFilename(filename)
                        if (
                            barcoded_filename.tissue != Tissue.PRIMARY_XENOGRAFT_TISSUE
                            and barcoded_filename.tissue
                            != Tissue.CELL_LINE_DERIVED_XENOGRAFT_TISSUE
                        ) or barcoded_filename.organism is None:
                            filenames_set.add(match.group(1))
                            added_files += 1

                if added_files == 0:
                    print("ERROR: cannot find any file for sample %s" % line.strip())
                    exit(-1)
        filenames = list(filenames_set)

    elif args.scan_samples:
        fastq_files = [
            filename
            for filename in os.listdir(args.fastq_dir)
            if re.search(r"\.fastq$", filename, re.I)
        ]
        filenames = list(
            set(
                [
                    re.sub(r"(\.R[12])?\.fastq$", "", filename, flags=re.I)
                    for filename in fastq_files
                    if not re.search(r"trimmed|clipped", filename, re.I)
                ]
            )
        )
    else:
        raise RuntimeError("Unhandled condition")

    parameters = {
        "use_xenograft_classifier": args.use_xenograft_classifier,
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
        "use_tdf": args.use_tdf,
    }

    logging.basicConfig(format="%(asctime)-15s %(message)s")

    aligner_is_needed = False
    rnaseq_aligner_is_needed = False
    xenograft_classifier_is_needed = False
    for sample in filenames:
        barcoded_filename = BarcodedFilename.from_sample(sample)
        if barcoded_filename.analyte == Analyte.RNASEQ:
            rnaseq_aligner_is_needed = True
        else:
            aligner_is_needed = True

        if barcoded_filename.is_xenograft():
            xenograft_classifier_is_needed = True

    if rnaseq_aligner_is_needed:
        set_aligner_param(
            args,
            parameters,
            "rnaseq_aligner",
            [RnaSeqAligner.STAR],
            [Config.check_star_files],
            config,
        )
    if aligner_is_needed:
        set_aligner_param(
            args,
            parameters,
            "aligner",
            [GenericAligner.NOVOALIGN, GenericAligner.BWA],
            [None, None],
            config,
        )

    if xenograft_classifier_is_needed:
        xenograft_classifier = Xenograft.get_available_classifier(args, config)
        if xenograft_classifier is None:
            if args.xenograft_classifier == "auto":
                print(
                    "Xenograft classifier is needed but none is available. "
                    "Check config.",
                    file=sys.stderr,
                )
            else:
                print(
                    f"Xenograft classifier '{args.xenograft_classifier}' is "
                    "selected, but it cannot be used. Check config.",
                    file=sys.stderr,
                )
            exit(-1)
        parameters["xenograft_classifier"] = xenograft_classifier

    if args.r_checks and args.post_recalibration:
        try:
            import rpy2.robjects.packages as rpackages
            from rpy2.robjects.vectors import StrVector
        except ImportError:
            rpackages = None
            print("cannot correctly import rpy2. R checks are skipped.")

        if rpackages:
            print("Checking R packages and, eventually, performing installations")
            dependencies = ("ggplot2", "gplots", "reshape", "grid", "tools", "gsalib")
            rutils = rpackages.importr("utils")
            base = rpackages.importr("base")
            rutils.chooseCRANmirror(ind=1)
            installed_packages = rutils.installed_packages().rx(True, 1)
            for package in dependencies:
                if not base.is_element(package, installed_packages)[0]:
                    sys.stdout.write("Installing R package %s..." % package)
                    rutils.install_packages(StrVector(dependencies), quiet=True)
                    print(" done.")

            installed_packages = rutils.installed_packages().rx(True, 1)
            for package in dependencies:
                if not base.is_element(package, installed_packages)[0]:
                    print(
                        "Package %s has not been correctly installed. "
                        "Try again or check this dependency manually." % package
                    )
                    exit(-1)
            print("Done with R packages")

    manager = Manager()
    runner = Runner(
        manager,
        root=args.root_dir,
        fastq_dir=args.fastq_dir,
        parameters=parameters,
        config=config,
    )

    error_raised = False
    try:
        with Pool(5) as pool:
            pool.map(runner, filenames)

        if args.mail:
            msg = MIMEText(
                "Pipeline for file list %s successfully completed." % args.list_file
            )
            msg["Subject"] = "Pipeline completed"
    except Exception:
        error_raised = True
        traceback.print_exc(file=sys.stdout)

        if args.mail:
            msg = MIMEText(
                "Error while executing pipeline for file list %s.\n"
                "Raised exception:\n%s" % (args.list_file, traceback.format_exc())
            )
            msg["Subject"] = "Pipeline error"

    if args.use_normals and not args.only_mapping:
        samples: Dict[str, Dict[str, Union[str, List[Tuple[str, str]]]]] = {}
        normals: Dict[str, BarcodedFilename] = {}
        last_operations = {}
        for sample, last_operation in runner.last_operations.items():
            last_operations[sample] = last_operation
            for filename in utils.get_sample_filenames(last_operation):
                barcoded_filename = BarcodedFilename(filename)
                assert barcoded_filename.biopsy is not None
                assert barcoded_filename.sequencing is not None

                fake_sample = "%s-%s-%d%d%d-%d%d%d" % (
                    barcoded_filename.project,
                    barcoded_filename.patient,
                    barcoded_filename.molecule,
                    barcoded_filename.analyte,
                    barcoded_filename.kit,
                    barcoded_filename.biopsy,
                    barcoded_filename.get_sample_index(),
                    barcoded_filename.sequencing,
                )

                if barcoded_filename.tissue in (
                    Tissue.BONE_MARROW_NORMAL,
                    Tissue.BUCCAL_CELL_NORMAL,
                    Tissue.SOLID_TISSUE_NORMAL,
                    Tissue.BLOOD_DERIVED_NORMAL,
                    Tissue.EBV_IMMORTALIZED_NORMAL,
                ):
                    sample_type = "control"
                else:
                    sample_type = "sample"

                if fake_sample not in samples:
                    samples[fake_sample] = {}

                if sample_type == "control":
                    samples[fake_sample]["control"] = filename
                    normals[filename] = barcoded_filename
                else:
                    if "sample" not in samples[fake_sample]:
                        samples[fake_sample]["sample"] = []

                    cast(List[Tuple[str, str]], samples[fake_sample]["sample"]).append(
                        (filename, sample)
                    )

        samples_with_no_normal = {
            sample: filenames["sample"][0][0]
            for sample, filenames in samples.items()
            if "control" not in filenames and len(filenames["sample"]) > 0
        }

        for sample, filename in samples_with_no_normal.items():
            sample_barcode = BarcodedFilename(filename)
            candidates = {
                filename: barcode
                for filename, barcode in normals.items()
                if barcode.project == sample_barcode.project
                and barcode.patient == sample_barcode.patient
                and barcode.molecule == sample_barcode.molecule
                and barcode.analyte == sample_barcode.analyte
                and barcode.kit == sample_barcode.kit
                and barcode.biopsy == sample_barcode.biopsy
                and barcode.sample == sample_barcode.sample
            }

            if len(candidates) == 0:
                candidates = {
                    filename: barcode
                    for filename, barcode in normals.items()
                    if barcode.project == sample_barcode.project
                    and barcode.patient == sample_barcode.patient
                    and barcode.molecule == sample_barcode.molecule
                    and barcode.analyte == sample_barcode.analyte
                    and barcode.kit == sample_barcode.kit
                    and barcode.biopsy == sample_barcode.biopsy
                }

            if len(candidates) == 0:
                candidates = {
                    filename: barcode
                    for filename, barcode in normals.items()
                    if barcode.project == sample_barcode.project
                    and barcode.patient == sample_barcode.patient
                    and barcode.molecule == sample_barcode.molecule
                    and barcode.analyte == sample_barcode.analyte
                    and barcode.kit == sample_barcode.kit
                }

            if len(candidates) == 0:
                candidates = {
                    filename: barcode
                    for filename, barcode in normals.items()
                    if barcode.project == sample_barcode.project
                    and barcode.patient == sample_barcode.patient
                    and barcode.molecule == sample_barcode.molecule
                    and barcode.analyte == sample_barcode.analyte
                }

            if len(candidates) == 0:
                candidates = {
                    filename: barcode
                    for filename, barcode in normals.items()
                    if barcode.project == sample_barcode.project
                    and barcode.patient == sample_barcode.patient
                    and barcode.molecule == sample_barcode.molecule
                }

            if len(candidates) == 1:
                samples[sample]["control"] = list(candidates.items())[0][0]
            elif len(candidates) > 1:
                candidates_list = list(candidates.items())
                del candidates
                candidates_list.sort(key=lambda x: os.stat(x[0]).st_size)
                samples[sample]["control"] = candidates_list[-1][0]

        triplets: List[Tuple[str, str, Optional[str]]] = []
        for sample, values in samples.items():
            if "sample" not in values:
                continue

            for tumor in values["sample"]:
                if "control" in values:
                    triplets.append(
                        cast(
                            Tuple[str, str, str],
                            (tumor[1], tumor[0], values["control"]),
                        )
                    )
                else:
                    triplets.append((tumor[1], tumor[0], None))

        error_raised = False
        try:
            with Pool(5) as pool:
                pool.starmap(runner.with_normals, triplets)

            if args.mail:
                msg = MIMEText(
                    "Pipeline for file list %s successfully completed." % args.list_file
                )
                msg["Subject"] = "Pipeline completed"

        except RuntimeError:
            error_raised = True
            traceback.print_exc(file=sys.stdout)

            if args.mail:
                msg = MIMEText(
                    "Error while executing pipeline for file list %s.\n"
                    "Raised exception:\n%s" % (args.list_file, traceback.format_exc())
                )
                msg["Subject"] = "Pipeline error"

    if args.mail and len(config.mails) > 0:
        msg["From"] = "UV2000 Pipeline <pipeline@uv2000.hugef>"
        msg["To"] = config.mails

        smtp = smtplib.SMTP("localhost")
        smtp.send_message(msg)
        smtp.quit()

    if error_raised:
        exit(-1)
