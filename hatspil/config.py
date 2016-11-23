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
        self.genome_ref = 'ucsc.hg19.fasta'
        self.hg19_index = "ucsc.hg19"
        self.cosmic = "Cosmic.hg19.vcf"
        self.dbsnp138 = 'dbsnp_138.hg19.vcf'
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
        default["genome_ref"] = self.genome_ref
        default["hg19_index"] = self.hg19_index
        default["cosmic"] = self.cosmic
        default["dpsnp138"] = self.dbsnp138
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
