"""The configuration module of HaTSPiL."""

import os
import subprocess
import sys
from configparser import ConfigParser, SectionProxy
from itertools import chain
from typing import Any, Dict, Optional, Tuple

from .core.barcoded_filename import Analyte


class KitData:
    """The information available for each kit.

    A configuration file must include at least one kit, which defines a
    set of properties for a subset of possible barcode values. This
    class stores all the information for a kit.
    """

    def __init__(self) -> None:
        """Create a default kit.

        The values are just placeholder and they should be overridden
        by the content of a `config.ini` file.
        """
        self.name = "Kit"
        self.target_list = "04818-1466508813_target.interval_list"
        self.bait_list = "04818-1466508813_amplicons.interval_list"
        self.indels_hg19 = [
            "1000G_phase1.indels.hg19.sites.vcf",
            "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
        ]
        self.indels_hg38 = [
            "1000G_phase1.indels.hg38.sites.vcf",
            "Mills_and_1000G_gold_standard.indels.hg38.sites.vcf",
        ]
        self.amplicons = "amplicons.bed"
        self.cancer_site = "soft_tissue"
        self.adapter_r1 = ""
        self.adapter_r2 = ""
        self.mean_len_library = 200
        self.sd_len_library = 100


class Config:
    """The configuration of HaTSPiL.

    This class is responsible for reading the configuration ini file and
    check for missing parameters. A plausible default value is available
    for most of the fields, in order to create a template config file
    and ease the process of creating a brand new ini file.
    """

    executables = ("java", "java7", "perl", "seqtk", "fastqc", "samtools")
    optional_executables = ("novoalign", "bwa", "star", "xenome", "disambiguate")
    jars = ("picard", "varscan", "gatk", "mutect", "bam2tdf")
    files = (
        "strelka_basedir",
        "strelka_config",
        "hg19_ref",
        "hg19_index",
        "hg38_ref",
        "hg38_index",
        "mm9_ref",
        "mm9_index",
        "mm10_ref",
        "mm10_index",
        "cosmic_hg19",
        "cosmic_hg38",
        "dbsnp_hg19",
        "dbsnp_hg38",
        "annovar_basedir",
        "temporary_dir",
    )
    optional_files = (
        "star_index_hg19",
        "star_index_hg38",
        "star_index_mm9",
        "star_index_mm10",
        "features_hg19",
        "features_hg38",
        "features_mm9",
        "features_mm10",
    )
    parameters = (
        "xenome_index",
        "xenome_threads",
        "strelka_threads",
        "use_hg19",
        "use_hg38",
        "use_mm9",
        "use_mm10",
        "mails",
        "use_mongodb",
        "picard_jvm_args",
        "varscan_jvm_args",
        "gatk_jvm_args",
        "mutect_jvm_args",
        "mutect_args",
    )
    kits_files = ("target_list", "bait_list", "indels_hg19", "indels_hg38", "amplicons")
    kits_parameters = (
        "name",
        "cancer_site",
        "adapter_r1",
        "adapter_r2",
        "mean_len_library",
        "sd_len_library",
    )
    mongodb = ("database", "host", "port", "username", "password")

    def __init__(self, filename: Optional[str] = None) -> None:
        """Create a template config and fill it from a file content.

        A config file is created with a set of plausible default values.
        If a `filename` is specified, the parameters are filled with the
        content of the specified config file.
        """
        self.java = "java"
        self.java7 = "java"
        self.perl = "perl"
        self.star = "STAR"
        self.star_index_hg19 = "star_index_hg19"
        self.star_index_hg38 = "star_index_hg38"
        self.star_index_mm9 = "star_index_mm9"
        self.star_index_mm10 = "star_index_mm10"
        self.features_hg19 = "features_hg19.gtf"
        self.features_hg38 = "features_hg38.gtf"
        self.features_mm9 = "features_mm9.gtf"
        self.features_mm10 = "features_mm10.gtf"
        self.bwa = "bwa"
        self.novoalign = "novoalign"
        self.picard = "picard.jar"
        self.picard_jvm_args = "-Xmx128g"
        self.varscan = "VarScan.jar"
        self.varscan_jvm_args = "-Xmx36g"
        self.gatk = "GenomeAnalysisTK.jar"
        self.gatk_jvm_args = "-Xmx128g"
        self.seqtk = "seqtk"
        self.fastqc = "fastqc"
        self.mutect = "mutect.jar "
        self.mutect_jvm_args = "-Xmx128g"
        self.mutect_args = "--analysis_type MuTect --enable_extended_output"
        self.samtools = "samtools"
        self.xenome = "xenome"
        self.xenome_index = "xenome_idx"
        self.xenome_threads = 1
        self.disambiguate = "disambiguate"
        self.strelka_basedir = "/usr/share/strelka"
        self.strelka_config = "/usr/share/strelka/config.ini"
        self.strelka_threads = 1
        self.annovar_basedir = "/usr/share/annovar"
        self.temporary_dir = "/tmp"
        self.bam2tdf = "bam2tdf.jar"
        self.hg19_ref = "ucsc.hg19.fasta"
        self.hg19_index = "ucsc.hg19"
        self.hg38_ref = "ucsc.hg38.fasta"
        self.hg38_index = "ucsc.hg38"
        self.mm9_ref = "ucsc.mm9.fasta"
        self.mm9_index = "ucsc.mm9"
        self.mm10_ref = "ucsc.mm10.fasta"
        self.mm10_index = "ucsc.mm10"
        self.cosmic_hg19 = "Cosmic.hg19.vcf"
        self.cosmic_hg38 = "Cosmic.hg38.vcf"
        self.dbsnp_hg19 = "dbsnp.hg19.vcf"
        self.dbsnp_hg38 = "dbsnp.hg38.vcf"
        self.use_hg19 = True
        self.use_hg38 = True
        self.use_mm9 = True
        self.use_mm10 = True
        self.kits: Dict[Tuple[int, Analyte], KitData] = {}
        self.mails = ""
        self.use_mongodb = True
        self.mongodb_database = "hatspil"
        self.mongodb_username = "hatspil"
        self.mongodb_password = "hatspil"
        self.mongodb_host = "localhost"
        self.mongodb_port = 27017

        if filename:
            self._check_after_init(filename)

    def _check_after_init(self, filename: str) -> None:
        VALID_SECTIONS = ("EXECUTABLES", "JARS", "FILES", "PARAMETERS", "MONGODB")

        parser = ConfigParser()
        parser.read(filename)
        for section_name in VALID_SECTIONS:
            if section_name not in parser:
                sys.stderr.write(
                    "WARNING: %s section in config not "
                    "found. Values are set to default\n" % section_name
                )
                continue

            section = parser[section_name]
            for key in section:
                if section_name == "MONGODB":
                    setattr(self, "mongodb_" + key, section[key])
                else:
                    setattr(self, key, section[key])

        if "PARAMETERS" in parser:
            for param in (
                "xenome_threads",
                "strelka_threads",
                "mean_len_library",
                "sd_len_library",
            ):
                if param in parser["PARAMETERS"]:
                    setattr(self, param, parser["PARAMETERS"].getint(param))
            for param in ("use_hg19", "use_hg38", "use_mm9", "use_mm10", "use_mongodb"):
                if param in parser["PARAMETERS"]:
                    setattr(self, param, parser["PARAMETERS"].getboolean(param))

        for section_name, section_params in parser.items():
            if section_name.startswith("KIT"):
                self._check_kit_section(section_name, section_params)
            elif section_name not in chain(VALID_SECTIONS, ("DEFAULT",)):
                sys.stderr.write("WARNING: '%s' section is invalid\n" % section_name)

        if "MONGODB" in parser:
            if "port" in parser["MONGODB"]:
                self.mongodb_port = parser["MONGODB"].getint("port")

    def _check_kit_section(
        self, section_name: str, section_params: SectionProxy
    ) -> None:
        splitted_section_name = section_name.split(" ")
        valid_section = True
        if len(splitted_section_name) != 3:
            valid_section = False

        if valid_section:
            try:
                kit_index = int(splitted_section_name[1])
                analyte = Analyte.__members__[splitted_section_name[2]]
            except Exception:
                try:
                    kit_index = int(splitted_section_name[2])
                    analyte = Analyte.__members__[splitted_section_name[1]]
                except Exception:
                    valid_section = False

        if valid_section:
            current_kit = self.kits.setdefault((kit_index, analyte), KitData())
        else:
            analytes_str = "|".join(Analyte.__members__.keys())
            sys.stderr.write(
                "WARNING: '{}' kit section is invalid. "
                "The format is as following:\n"
                "KIT n [{}]\nor\nKIT [{}] n\n".format(
                    section_name, analytes_str, analytes_str
                )
            )
            return

        for param_name in (
            "name",
            "cancer_site",
            "adapter_r1",
            "adapter_r2",
            "target_list",
            "bait_list",
            "indels_hg19",
            "indels_hg38",
            "amplicons",
        ):
            if param_name in section_params:
                param_value = section_params[param_name]
                if param_name.startswith("indels_"):
                    param_value = [
                        filename.strip() for filename in param_value.split(",")
                    ]
                setattr(current_kit, param_name, param_value)

        for param_name in ("mean_len_library", "sd_len_library"):
            if param_name in section_params:
                setattr(current_kit, param_name, section_params.getint(param_name))

    def save(self, filename: str) -> None:
        """Save the object into a config file."""
        config = ConfigParser()

        executables = {}
        for executable in Config.executables:
            executables[executable] = getattr(self, executable)
        config["EXECUTABLES"] = executables

        jars = {}
        for jar in Config.jars:
            jars[jar] = getattr(self, jar)
        config["JARS"] = jars

        files = {}
        for filepath in Config.files:
            files[filepath] = getattr(self, filepath)
        config["FILES"] = files

        parameters = {}
        for parameter in Config.parameters:
            parameters[parameter] = str(getattr(self, parameter))
        config["PARAMETERS"] = parameters

        available_kits = dict(self.kits)
        if not available_kits:
            available_kits[(0, Analyte.WHOLE_EXOME)] = KitData()

        for kit_name, kit_data in available_kits.items():
            section_name = "KIT {} {}".format(kit_name[0], kit_name[1].name)
            kit = {}
            for param in chain(Config.kits_files, Config.kits_parameters):
                param_value = getattr(kit_data, param)
                if param.startswith("indel_"):
                    param_value = ",".join(param_value)
                kit[param] = param_value
            config[section_name] = kit

        mongodb_params = {}
        for mongodb_param in Config.mongodb:
            mongodb_params["mongodb_" + mongodb_param] = str(
                getattr(self, "mongodb_" + mongodb_param)
            )
        config["MONGODB"] = mongodb_params

        with open(filename, "w") as fd:
            config.write(fd)

    def check_programs(self) -> bool:
        """Check if the mandatory software are available."""
        ok = True

        for param in Config.jars:
            jar = getattr(self, param)
            if not os.access(jar, os.R_OK):
                sys.stderr.write(
                    "ERROR: %s jar file cannot be found. Check config.\n" % param
                )
                ok = False

        for param in Config.executables:
            executable = getattr(self, param)
            if (
                subprocess.call(
                    executable,
                    shell=True,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                == 127
            ):
                sys.stderr.write(
                    "ERROR: %s cannot be executed. Check config.\n" % param
                )
                ok = False

        for param in Config.optional_executables:
            executable = getattr(self, param)
            if (
                subprocess.call(
                    executable,
                    shell=True,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                == 127
            ):
                setattr(self, param, None)

        return ok

    def _check_valid_annotation(self, param: str) -> bool:
        if "hg19" in param:
            return self.use_hg19
        elif "hg38" in param:
            return self.use_hg38
        elif "mm9" in param:
            return self.use_mm9
        elif "mm10" in param:
            return self.use_mm10

        return True

    @staticmethod
    def _check_readable_filename(filename: str, param: str) -> bool:
        if not os.access(filename, os.R_OK):
            sys.stderr.write(
                "ERROR: {} for param {} cannot be read. Check config.\n".format(
                    filename, param
                )
            )
            return False
        else:
            return True

    def _check_file_param(self, obj: Any, param: str) -> Optional[bool]:
        if not self._check_valid_annotation(param):
            return None

        filepath = getattr(obj, param)
        return Config._check_readable_filename(filepath, param)

    def check_files(self) -> bool:
        """Check if all files are available."""
        ok = True

        for param in Config.files:
            param_is_valid = self._check_file_param(self, param)
            if param_is_valid is None:
                continue

            if not param_is_valid:
                ok = False

        for kit_data in self.kits.values():
            for param in Config.kits_files:
                if not getattr(kit_data, param):
                    continue

                if param.startswith("indels_"):
                    if not self._check_valid_annotation(param):
                        continue

                    for filename in getattr(kit_data, param):
                        if not Config._check_readable_filename(filename, param):
                            ok = False
                else:
                    param_is_valid = self._check_file_param(kit_data, param)
                    if param_is_valid is None:
                        continue

                    if not param_is_valid:
                        ok = False

        return ok

    def check_star_files(self) -> bool:
        """Check if all files needed by start are available.

        This check is separated from `check_files` because the
        availability of STAR is not mandatory.
        """
        ok = True
        for param in ("star_index", "features"):
            for organism in ("hg19", "hg38", "mm9", "mm10"):
                param_is_valid = self._check_file_param(self, f"{param}_{organism}")
                if param_is_valid is None:
                    continue

                if not param_is_valid:
                    ok = False

        return ok
