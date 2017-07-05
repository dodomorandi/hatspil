from .executor import Executor
from .ranges import GenomicRange, GenomicRanges
from .barcoded_filename import BarcodedFilename

import os
import glob
import pandas as pd
import vcf
import math
from formatizer import f
import re


class VariantCalling:
    min_allele_frequency = 0.01
    min_cov_position = 10
    build_version = "hg19"
    cancer_site = "soft_tissue"
    dataset_filename = os.path.join(os.path.dirname(__file__), "data.hdf")
    medium_damage = 5
    high_damage = 20

    def __init__(self, analysis):
        self.analysis = analysis

        try:
            os.makedirs(analysis.get_out_dir(), exist_ok=True)
        except OSError:
            pass

        self.mutect_filenames = glob.glob(os.path.join(analysis.get_out_dir(),
                                                       self.analysis.basename)
                                          + "*.mutect*.vcf")
        self.varscan_filenames = {}
        for varscan_type in ("snp", "indel"):
            self.varscan_filenames[varscan_type] =\
                glob.glob(os.path.join(analysis.get_out_dir(), self.analysis.basename)
                          + "*.varscan2." + varscan_type + ".vcf")
        self.annovar_dirname = os.path.join(analysis.get_out_dir(),
                                            self.analysis.basename + "_annovar")
        try:
            os.makedirs(self.annovar_dirname, exist_ok=True)
        except:
            pass

        self.annovar_file = os.path.join(self.annovar_dirname, self.analysis.basename + "_annovar_input")
        self.variants_filename = os.path.join(self.annovar_dirname, self.analysis.basename + "_variants.csv")
        self.multianno_filename = os.path.join("%s.%s_multianno.txt" % (self.annovar_file, VariantCalling.build_version))

    def chdir(self):
        os.chdir(self.analysis.get_out_dir())

    def prepare_for_annovar(self):
        self.analysis.logger.info("Starting data preparation for ANNOVAR")

        self.variants = None

        mutect_data = None
        varscan_data = None

        if len(self.mutect_filenames) > 0:
            for mutect_filename in self.mutect_filenames:
                if mutect_data is None:
                    mutect_data = pd.read_table(mutect_filename, comment="#")
                else:
                    mutect_data = mutect_data.join(pd.read_table(mutect_filename, comment="#"))

            mutect_data.insert(0, "key", mutect_data.apply(lambda row: "%s:%d-%d_%s_%s" % (row.contig, row.position, row.position, row.ref_allele, row.alt_allele), axis=1))

            mutect_data.drop(mutect_data[mutect_data.tumor_f < VariantCalling.min_allele_frequency].index, inplace=True)
            mutect_data.drop(mutect_data[(mutect_data.judgement == "REJECT") & (mutect_data.failure_reasons != "possible_contamination")].index, inplace=True)
            mutect_data.drop(mutect_data[
                    mutect_data.strand_bias_counts.str.split(R"\W")
                        .map(lambda values: [int(value) for value in values if value != ""])
                        .map(lambda values: values[2] <= 0 or values[3] <= 0)].index,
                    inplace=True)
            mutect_data["tot_cov"] = mutect_data.t_ref_count + mutect_data.t_alt_count
            mutect_data.drop(mutect_data[mutect_data.tot_cov <= VariantCalling.min_cov_position].index, inplace=True)

        for varscan_type in ("snp", "indel"):
            filenames = self.varscan_filenames[varscan_type]
            if len(filenames) > 0:
                varscan_data_list = []
                for varscan_filename in filenames:
                    for record in vcf.Reader(filename=varscan_filename):
                        for sample in record.samples:
                            varscan_data_list.append({
                                "chr": record.CHROM,
                                "start": record.start + 1,
                                "end": record.end,
                                "width": record.end - record.start,
                                "ref": str(record.REF),
                                "alt": ",".join([str(alt) for alt in record.ALT]),
                                "totalDepth": sample.data.DP,
                                "refDepth": sample.data.RD,
                                "altDepth": sample.data.AD,
                                "sampleNames": sample.sample,
                                "indelError": "indelError" not in record.FILTER,
                                "QUAL": record.QUAL,
                                "DP": record.INFO["DP"],
                                "SOMATIC": "SOMATIC" in record.INFO and record.INFO["SOMATIC"],
                                "SS": record.INFO["SS"],
                                "SSC": record.INFO["SSC"],
                                "GPV": record.INFO["GPV"],
                                "SPV": record.INFO["SPV"],
                                "GT": sample.data.GT,
                                "GQ": sample.data.GQ,
                                "RD": sample.data.RD,
                                "FREQ": float(sample.data.FREQ[:-1]) / 100,
                                "DP4": int(sample.data.DP4.split(",")[0])
                            })

                current_data = pd.DataFrame(data=varscan_data_list, columns=[
                    "chr", "start", "end", "width", "ref", "alt", "totalDepth",
                    "refDepth", "altDepth", "sampleNames", "indelError", "QUAL", "DP",
                    "SOMATIC", "SS", "SSC", "GPV", "SPV", "GT", "GQ", "RD", "FREQ",
                    "DP4"])
                current_data.drop_duplicates(inplace=True)
                current_data["key"] = current_data.apply(lambda row: "%s:%d-%d_%s_%s" % (row.chr, row.start, row.end, row.ref, row.alt), axis=1)
                tumor = current_data[current_data.sampleNames == "TUMOR"].copy()
                normal = current_data[current_data.sampleNames == "NORMAL"]

                normal.reset_index(inplace=True, drop=True)
                tumor.reset_index(inplace=True, drop=True)
                normals = normal.loc[pd.match(tumor.key, normal.key)][["totalDepth", "refDepth", "altDepth"]]
                tumor["totalDepth.NORMAL"] = normals["totalDepth"].values
                tumor["refDepth.NORMAL"] = normals["refDepth"].values
                tumor["altDepth.NORMAL"] = normals["altDepth"].values

                tumor.drop(tumor[(~tumor.indelError) | (~tumor.SOMATIC) | (tumor.DP < VariantCalling.min_cov_position) | (tumor.FREQ < VariantCalling.min_allele_frequency)].index, inplace=True)

                if varscan_data is None:
                    varscan_data = tumor
                else:
                    varscan_data = pd.concat((varscan_data, tumor))

        self.variants = pd.DataFrame(columns=("key", "DP", "FREQ", "method"))

        if varscan_data is None:
            varscan_data = pd.DataFrame(columns=("key", "DP", "FREQ"))
        else:
            varscan_data = varscan_data[["key", "DP", "FREQ"]].copy()
            self.variants = pd.concat([self.variants, varscan_data])

        if mutect_data is None:
            mutect_data = pd.DataFrame(columns=("key", "DP", "FREQ"))
        else:
            mutect_data = mutect_data.rename(index=str, columns={"tot_cov": "DP", "tumor_f": "FREQ"})[["key", "DP", "FREQ"]].copy()
            self.variants = pd.concat([self.variants, mutect_data])

        self.variants.drop_duplicates(("key",), inplace=True)
        self.variants.loc[(self.variants.key.isin(mutect_data.key)) & (~self.variants.key.isin(varscan_data.key)), "method"] = "Mutect1.17"
        self.variants.loc[(~self.variants.key.isin(mutect_data.key)) & (self.variants.key.isin(varscan_data.key)), "method"] = "VarScan2"
        self.variants.loc[(self.variants.key.isin(mutect_data.key)) & (self.variants.key.isin(varscan_data.key)), "method"] = "Mutect1.17:VarScan2"

        with open(self.annovar_file, "w") as fd:
            for key in self.variants.key:
                splitted = re.split(R"[:_-]", key)
                if len(splitted) != 5:
                    continue
                fd.write("%s\n" % "\t".join(splitted))

        self.variants.to_csv(self.variants_filename, index=False)

        self.analysis.logger.info("Finished data preparation for ANNOVAR")

    def annovar(self):
        self.analysis.logger.info("Running ANNOVAR")
        config = self.analysis.config

        os.chdir(config.annovar_basedir)
        executor = Executor(self.analysis)

        executor(f(
            "{config.perl} table_annovar.pl "
            "{{input_filename}} humandb/ "
            "-buildver {VariantCalling.build_version} "
            "-protocol refGene,snp138,cosmic70,clinvar_20160302,popfreq_all_20150413,dbnsfp30a,cadd13 "
            "-operation g,f,f,f,f,f,f -nastring NA -remove -v"),
            input_filenames=self.annovar_file,
            override_last_files=False,
            error_string="ANNOVAR exited with status {status}",
            exception_string="annovar error",
            only_human=True
        )

        self.analysis.logger.info("Finished running ANNOVAR")

    def collect_annotated_variants(self):
        self.analysis.logger.info("Collecting annotated variants from ANNOVAR")

        annotation = pd.read_table(self.multianno_filename)
        annotation = annotation[annotation.Chr.str.match(R"chr(?:\d{1,2}|[xXyY])")]
        annotation.insert(0, "id", annotation.apply(lambda row: "%s:%d-%d_%s_%s" % (row.Chr, row.Start, row.End, row.Ref, row.Alt), axis=1))

        cancer_genes = pd.read_hdf(VariantCalling.dataset_filename, "cancer_genes")
        panel_drug = pd.read_hdf(VariantCalling.dataset_filename, "panel_drug")
        gene_info = pd.read_hdf(VariantCalling.dataset_filename, "gene_info")

        selected_cancer_genes = cancer_genes[cancer_genes.cancer_site == VariantCalling.cancer_site]

        annotation.reset_index(inplace=True, drop=True)
        gene_info.reset_index(inplace=True, drop=True)
        annotation["gene_type"] = gene_info.loc[pd.match(annotation["Gene.refGene"], gene_info.symbol)].cancer_type.values
        annotation.loc[annotation.gene_type == "rst", "gene_type"] = float("nan")
        annotation["cancer_gene_site"] = float("nan")
        annotation.loc[annotation["Gene.refGene"].isin(selected_cancer_genes.symbol), "cancer_gene_site"] = VariantCalling.cancer_site

        annotation.loc[annotation["Func.refGene"] == "splicing", "ExonicFunc.refGene"] = "splicing"
        annotation.loc[annotation["Func.refGene"] == "splicing", "AAChange.refGene"] = annotation[annotation["Func.refGene"] == "splicing"]["GeneDetail.refGene"]
        annotation.loc[annotation["Func.refGene"] == "UTR3", "AAChange.refGene"] = annotation[annotation["Func.refGene"] == "UTR3"]["GeneDetail.refGene"]
        annotation.loc[annotation["Func.refGene"] == "UTR5", "AAChange.refGene"] = annotation[annotation["Func.refGene"] == "UTR5"]["GeneDetail.refGene"]

        annotation["Gene.refGene"] = annotation["Gene.refGene"].map(lambda value: re.split(R"[^A-Za-z0-9]", value)[0])

        annotation.loc[annotation.CADD13_PHRED > VariantCalling.medium_damage, "damaging"] = "Medium"
        annotation.loc[(annotation.CADD13_PHRED > VariantCalling.high_damage) |
                      (annotation["ExonicFunc.refGene"].isin(("frameshift substitution","nonframeshift substitution","splicing"))), "damaging"] = "High"

        hgnc = pd.read_hdf(VariantCalling.dataset_filename, "hgnc")
        hgnc.dropna(subset=("chromosome",), inplace=True)

        annotation.sort_values(["Chr", "Start", "End"], inplace=True)
        hgnc.sort_values(["chromosome", "tstart", "tend"], inplace=True)
        annotation.reset_index(drop=True, inplace=True)
        hgnc.reset_index(drop=True, inplace=True)

        annotation_ranges = GenomicRanges([GenomicRange(row.Chr, row.Start - 1, row.End) for _, row in annotation.iterrows()])
        hgnc_ranges = GenomicRanges([GenomicRange(row.chromosome, row.tstart - 1, row.tend) for _, row in hgnc.iterrows()])

        if annotation_ranges.resorted:
            raise Exception("annotation ranges expected to be not resorted")
        if hgnc_ranges.resorted:
            raise Exception("hgnc ranges expected to be not resorted")
        overlaps = annotation_ranges.overlaps(hgnc_ranges)
        annotation.loc[[index[0] for index in overlaps],
                      "hgnc_refseq_accession"] =\
            hgnc.loc[[index[1] for index in overlaps]].refseq_accession.tolist()

        annotation["hgnc_canonical_refseq"] = annotation.apply(
            lambda row : "|".join([refgene
                                   for refgene
                                   in str(row["AAChange.refGene"]).split(",")
                                   if re.search(str(row.hgnc_refseq_accession), refgene)]), axis=1)
        empty_canonical_refseqs = annotation.hgnc_canonical_refseq == ""
        annotation.loc[empty_canonical_refseqs, "alternative_refseq"] = annotation.loc[empty_canonical_refseqs, "AAChange.refGene"].tolist()
        annotation.loc[empty_canonical_refseqs, "hgnc_canonical_refseq"] = float("nan")

        haloplex = pd.read_table(self.analysis.config.annotations,
                                 header=None,
                                 names=("chrom", "start", "end", "name", "length", "strand"),
                                 skiprows=2,
                                 delim_whitespace=True)

        haloplex.sort_values(["chrom", "start", "end"], inplace=True)
        haloplex_ranges = GenomicRanges([GenomicRange(row.chrom, row.start, row.end) for _, row in haloplex.iterrows()])
        overlaps = annotation_ranges.overlaps(haloplex_ranges)
        annotation["in_gene_panel"] = False
        annotation.loc[set([index[0] for index in overlaps]), "in_gene_panel"] = True
        annotation["druggable"] = annotation["Gene.refGene"].isin(panel_drug.gene_symbol)

        self.variants = pd.read_csv(self.variants_filename, index_col=False)
        annotation.reset_index(inplace=True, drop=True)
        self.variants.reset_index(inplace=True, drop=True)

        annotation_all = pd.merge(annotation,
                                  self.variants.loc[
                                      pd.match(annotation.id,
                                               self.variants.key)]
                                  .rename(columns={"key": "id"}),
                                  on="id",
                                  how="inner")
        annotation_all.to_csv(os.path.join(
            self.analysis.get_out_dir(),
            self.analysis.basename + "_variants.csv"), index=False)
        annotation_all[[
            "id", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
            "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene",
            "AAChange.refGene", "snp138", "cosmic70", "CLINSIG", "CLNDBN",
            "CLNACC", "CLNDSDB", "CLNDSDBID", "1000G_ALL", "CADD13_PHRED",
            "gene_type", "damaging", "hgnc_refseq_accession",
            "hgnc_canonical_refseq", "alternative_refseq", "druggable", "DP",
            "FREQ", "method"
        ]].to_csv(os.path.join(self.analysis.get_out_dir(),
                               self.analysis.basename +
                               "_variants_pretty.csv"),
                  index=False)

        del annotation_all

        config = self.analysis.config
        if config.use_mongodb:
            from pymongo import MongoClient, ReturnDocument

            mongo = MongoClient(config.mongodb_host, config.mongodb_port)
            db = mongo[config.mongodb_database]
            db.authenticate(config.mongodb_username, config.mongodb_password)

            barcoded_sample = BarcodedFilename.from_sample(self.analysis.sample)

            self.variants = pd.read_csv(self.variants_filename, index_col=False)
            variants = [{key.replace(".", " "): value
                         for key, value in variant.items()
                         if type(value) != float or not math.isnan(value)}
                        for variant in self.variants.to_dict("records")]
            annotations = [{key.replace(".", " "): value
                            for key, value in record.items()
                            if type(value) != float or not math.isnan(value)}
                           for record in annotation.to_dict("records")]

            def find_or_insert(collection, data, new_data=None):
                set_data = dict(data)
                if new_data is not None:
                    set_data.update(new_data)
                return collection.find_one_and_update(
                    data, {"$set": set_data}, upsert=True,
                    return_document=ReturnDocument.AFTER)

            project_obj = find_or_insert(db.projects, {
                "name": barcoded_sample.project
            })

            patient_obj = find_or_insert(db.patients, {
                "project": project_obj["_id"],
                "name": barcoded_sample.patient
            })

            biopsy_obj = find_or_insert(db.biopsies, {
                "patient": patient_obj["_id"],
                "index": barcoded_sample.biopsy,
                "tissue": barcoded_sample.tissue
            })

            sample_obj = find_or_insert(db.samples, {
                "biopsy": biopsy_obj["_id"],
                "index": barcoded_sample.sample
            })

            sequencing_obj = find_or_insert(db.sequencings, {
                "sample": sample_obj["_id"],
                "index": barcoded_sample.sequencing
            }, {
                "molecule": barcoded_sample.molecule,
                "analyte": barcoded_sample.analyte,
                "kit": barcoded_sample.kit
            })

            find_or_insert(db.analyses, {
                "sequencing": sequencing_obj["_id"],
                "date": self.analysis.current
            }, {
                "variants": variants,
                "annotations": annotations,
            })

        self.analysis.logger.info("Finished collecting annotated variants "
                                  "from ANNOVAR")

    def run(self):
        if len(self.mutect_filenames) + sum([len(values) for values in self.varscan_filenames.values()]) > 0:
            self.prepare_for_annovar()
            self.annovar()
            self.collect_annotated_variants()
