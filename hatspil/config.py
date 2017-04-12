from configparser import ConfigParser

import os
import sys
import subprocess


class Config:
    executables = ("novoalign", "seqtk", "fastqc", "samtools", "xenome")
    jars = ("picard", "varscan", "gatk", "mutect", "bam2tdf")

    def __init__(self, filename=None):
        self.java = "java"
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
        self.mean_len_library = 200
        self.sd_len_library = 100
        self.kit = "Kit"
        self.target_list = "04818-1466508813_target.interval_list"
        self.bait_list = "04818-1466508813_amplicons.interval_list"
        self.indel_1 = "1000G_phase1.indels.hg19.sites.vcf"
        self.indel_2 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
        self.mails = ""

        if filename:
            parser = ConfigParser()
            parser.read(filename)
            for section_name in "EXECUTABLES", "PARAMETERS", "JARS":
                if section_name not in parser:
                    sys.stderr.write("WARNING: %s section in config not "
                                     "found. Values are set to default\n" %
                                     section_name)
                    continue

                section = parser[section_name]
                for key in section:
                    setattr(self, key, section[key])

            self.xenome_threads = int(self.xenome_threads)
            self.strelka_threads = int(self.strelka_threads)
            self.mean_len_library = int(self.mean_len_library)
            self.sd_len_library = int(self.sd_len_library)

    def save(self, filename):
        config = ConfigParser()

        executables = config["EXECUTABLES"]
        for executable in Config.executables:
            executables[executable] = getattr(self, executable)

        jars = config["JARS"]
        for jar in Config.jars:
            jars[jar] = getattr(self, jar)

        parameters = config["PARAMETERS"]
        parameters["xenome_index"] = self.xenome_index
        parameters["xenome_threads"] = str(self.xenome_threads)
        parameters["strelka_basedir"] = self.strelka_basedir
        parameters["strelka_config"] = self.strelka_config
        parameters["strelka_threads"] = str(self.strelka_threads)
        parameters["hg19_ref"] = self.hg19_ref
        parameters["hg19_index"] = self.hg19_index
        parameters["hg38_ref"] = self.hg38_ref
        parameters["hg38_index"] = self.hg38_index
        parameters["mm9_ref"] = self.mm9_ref
        parameters["mm9_index"] = self.mm9_index
        parameters["mm10_ref"] = self.mm10_ref
        parameters["mm10_index"] = self.mm10_index
        parameters["cosmic_hg19"] = self.cosmic_hg19
        parameters["cosmic_hg38"] = self.cosmic_hg38
        parameters["dpsnp138_hg19"] = self.dbsnp138_hg19
        parameters["dpsnp138_hg38"] = self.dbsnp138_hg38
        parameters["mean_len_library"] = str(self.mean_len_library)
        parameters["sd_len_library"] = str(self.sd_len_library)
        parameters["kit"] = self.kit
        parameters["target_list"] = self.target_list
        parameters["bait_list"] = self.bait_list
        parameters["indel_1"] = self.indel_1
        parameters["indel_2"] = self.indel_2
        parameters["mails"] = self.mails

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
            if not os.access(executable.split(" ")[0], os.X_OK) or \
                    subprocess.call(executable,
                                    shell=True,
                                    stdin=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL) == 127:
                sys.stderr.write("ERROR: %s cannot be executed. "
                                 "Check config.\n" % param)
                ok = False

        return ok
