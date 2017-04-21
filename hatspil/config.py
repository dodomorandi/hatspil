from configparser import ConfigParser

import os
import sys
import subprocess


class Config:
    executables = ("java", "perl", "novoalign", "seqtk", "fastqc", "samtools", "xenome")
    jars = ("picard", "varscan", "gatk", "mutect", "bam2tdf")
    files = ("strelka_basedir", "strelka_config", "hg19_ref",
             "hg19_index", "hg38_ref", "hg38_index", "mm9_ref", "mm9_index",
             "mm10_ref", "mm10_index", "cosmic_hg19", "cosmic_hg38",
             "dbsnp138_hg19", "dbsnp138_hg38", "target_list", "bait_list",
             "indel_1", "indel_2", "annovar_basedir", "annotations")
    parameters = ("xenome_index", "xenome_threads", "strelka_threads",
                  "mean_len_library", "sd_len_library", "use_hg19", "use_hg38",
                  "use_mm9", "use_mm10", "kit", "mails")

    def __init__(self, filename=None):
        self.java = "java"
        self.perl = "perl"
        self.novoalign = "novoalign"
        self.picard = 'picard.jar'
        self.picard_jvm_args = '-Xmx128g'
        self.varscan = 'VarScan.jar'
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
        self.strelka_basedir = "/usr/share/strelka"
        self.strelka_config = "/usr/share/strelka/config.ini"
        self.strelka_threads = 1
        self.annovar_basedir = "/usr/share/annovar"
        self.bam2tdf = "bam2tdf.jar"
        self.hg19_ref = 'ucsc.hg19.fasta'
        self.hg19_index = "ucsc.hg19"
        self.hg38_ref = 'ucsc.hg38.fasta'
        self.hg38_index = "ucsc.hg38"
        self.mm9_ref = 'ucsc.mm9.fasta'
        self.mm9_index = "ucsc.mm9"
        self.mm10_ref = 'ucsc.mm10.fasta'
        self.mm10_index = "ucsc.mm10"
        self.cosmic_hg19 = "Cosmic.hg19.vcf"
        self.cosmic_hg38 = "Cosmic.hg38.vcf"
        self.dbsnp138_hg19 = 'dbsnp_138.hg19.vcf'
        self.dbsnp138_hg38 = 'dbsnp_138.hg38.vcf'
        self.use_hg19 = True
        self.use_hg38 = True
        self.use_mm9 = True
        self.use_mm10 = True
        self.mean_len_library = 200
        self.sd_len_library = 100
        self.kit = "Kit"
        self.target_list = "04818-1466508813_target.interval_list"
        self.bait_list = "04818-1466508813_amplicons.interval_list"
        self.indel_1 = "1000G_phase1.indels.hg19.sites.vcf"
        self.indel_2 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
        self.annotations = "48379-1473715058_Amplicons.bed"
        self.mails = ""

        if filename:
            parser = ConfigParser()
            parser.read(filename)
            for section_name in "EXECUTABLES", "JARS", "FILES", "PARAMETERS":
                if section_name not in parser:
                    sys.stderr.write("WARNING: %s section in config not "
                                     "found. Values are set to default\n" %
                                     section_name)
                    continue

                section = parser[section_name]
                for key in section:
                    setattr(self, key, section[key])

            self.xenome_threads = parser["PARAMETERS"].getint("xenome_threads")
            self.strelka_threads = parser["PARAMETERS"].getint("strelka_threads")
            self.mean_len_library = parser["PARAMETERS"].getint("mean_len_library")
            self.sd_len_library = parser["PARAMETERS"].getint("sd_len_library")
            self.use_hg19 = parser["PARAMETERS"].getboolean("use_hg19")
            self.use_hg38 = parser["PARAMETERS"].getboolean("use_hg38")
            self.use_mm9 = parser["PARAMETERS"].getboolean("use_mm9")
            self.use_mm10 = parser["PARAMETERS"].getboolean("use_mm10")

    def save(self, filename):
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

        with open(filename, "w") as fd:
            config.write(fd)

    def check_programs(self):
        ok = True

        for param in Config.jars:
            jar = getattr(self, param)
            if not os.access(jar, os.R_OK):
                sys.stderr.write("ERROR: %s jar file cannot be found. "
                                 "Check config.\n" % param)
                ok = False

        for param in Config.executables:
            executable = getattr(self, param)
            if subprocess.call(executable,
                               shell=True,
                               stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL) == 127:
                sys.stderr.write("ERROR: %s cannot be executed. "
                                 "Check config.\n" % param)
                ok = False

        return ok


    def check_files(self):
        ok = True

        for param in Config.files:
            if "hg19" in param:
                if not self.use_hg19:
                    continue
            elif "hg38" in param:
                if not self.use_hg38:
                    continue
            elif "mm9" in param:
                if not self.use_mm9:
                    continue
            elif "mm10" in param:
                if not self.use_mm10:
                    continue

            filepath = getattr(self, param)
            if not os.access(filepath, os.R_OK):
                sys.stderr.write("ERROR: %s cannot be read. "
                                 "Check config.\n" % param)
                ok = False

        return ok
