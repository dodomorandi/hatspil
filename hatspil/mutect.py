from .exceptions import PipelineError
from . import utils

import os
from formatizer import f


class Mutect:

    def __init__(self, analysis):
        self.analysis = analysis

    def chdir(self):
        os.chdir(self.analysis.out_dir)

    def run(self):
        self.analysis.logger.info("Running mutect")
        self.chdir()
        config = self.analysis.config

        input_files = self.analysis.bamfiles.items()
        for key in self.analysis.bamfiles.keys():
            if not key.startswith("hg"):
                del input_files[key]

        for organism, bamfile in input_files:
            if len(input_files) == 1:
                organism_str = ""
            else:
                organism_str = "_%s" % organism

            genome_ref = utils.get_genome_ref_index_by_organism(config, organism)[0]
            dbsnp = utils.get_dbsnp_by_organism(config, organism)
            cosmic = utils.get_cosmic_by_organism(config, organism)
            retval = utils.run_and_log(f(
                '{config.java7} {config.mutect} --reference_sequence '
                '{genome_ref} --dbsnp {dbsnp} --cosmic '
                '{cosmic} --input_file:tumor {bamfile} --out '
                '{self.analysis.basename}{organism_str}.mutect.1.17.vcf '
                ),
                self.analysis.logger
            )

            if retval != 0:
                self.analysis.logger.error("Mutect exited with status %d" % retval)
                raise PipelineError("mutect error")
