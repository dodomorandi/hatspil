from configparser import ConfigParser


class Config:

    def __init__(self, filename=None):
        self.novoalign = "novoalign"
        self.picard = 'java -Xmx128g -jar picard.jar'
        self.varscan = 'java -Xmx36g -jar VarScan.jar'
        self.gatk = "java -Xmx128g -jar GenomeAnalysisTK.jar"
        self.seqtk = "seqtk"
        self.fastqc = "fastqc"
        self.java7 = "java"
        self.mutect = "-Xmx128g -jar mutect.jar "\
            "--analysis_type MuTect --enable_extended_output"
        self.samtools = "samtools"
        self.xenome = "xenome"
        self.xenome_index = "xenome_idx"
        self.xenome_threads = 1
        self.strelka_basedir = "/usr/share/strelka"
        self.strelka_config = "/usr/share/strelka/config.ini"
        self.strelka_threads = 1
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
            default_section = parser["DEFAULT"]
            for key in default_section:
                setattr(self, key, default_section[key])

            self.xenome_threads = int(self.xenome_threads)
            self.strelka_threads = int(self.strelka_threads)
            self.mean_len_library = int(self.mean_len_library)
            self.sd_len_library = int(self.sd_len_library)

    def save(self, filename):
        config = ConfigParser()
        default = config["DEFAULT"]
        default["novoalign"] = self.novoalign
        default["picard"] = self.picard
        default["varscan"] = self.varscan
        default["gatk"] = self.gatk
        default["seqtk"] = self.seqtk
        default["fastqc"] = self.fastqc
        default["java7"] = self.java7
        default["mutect"] = self.mutect
        default["samtools"] = self.samtools
        default["xenome"] = self.xenome
        default["xenome_index"] = self.xenome_index
        default["xenome_threads"] = str(self.xenome_threads)
        default["strelka_basedir"] = self.strelka_basedir
        default["strelka_config"] = self.strelka_config
        default["strelka_threads"] = self.strelka_threads
        default["hg19_ref"] = self.hg19_ref
        default["hg19_index"] = self.hg19_index
        default["hg38_ref"] = self.hg38_ref
        default["hg38_index"] = self.hg38_index
        default["mm9_ref"] = self.mm9_ref
        default["mm9_index"] = self.mm9_index
        default["mm10_ref"] = self.mm10_ref
        default["mm10_index"] = self.mm10_index
        default["cosmic_hg19"] = self.cosmic_hg19
        default["cosmic_hg38"] = self.cosmic_hg38
        default["dpsnp138_hg19"] = self.dbsnp138_hg19
        default["dpsnp138_hg38"] = self.dbsnp138_hg38
        default["mean_len_library"] = str(self.mean_len_library)
        default["sd_len_library"] = str(self.sd_len_library)
        default["kit"] = self.kit
        default["target_list"] = self.target_list
        default["bait_list"] = self.bait_list
        default["indel_1"] = self.indel_1
        default["indel_2"] = self.indel_2
        default["mails"] = self.mails

        with open(filename, "w") as fd:
            config.write(fd)
