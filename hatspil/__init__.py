"""
The HaTSPiL package.

Welcome to HaTSPiL documentation! Here you can find all the information
you need to extend and improve the tool. You can also use the provided
information to use the APIs in order to extend the features of other
software.

Here a brief description of all the modules and packages available:
    * hatspil - the starting point of the executable. Argument parsing
                and initial checks are perfomed in this module
    * core - A package with the main core features of HaTSPiL.
    * config - The module responsible for handling the configuration
               file data.
    * runner - Where the analysis of a single sample starts. This module
               starts all the macro-steps required to analyse each
               sample.
    * mapping - This module handles the raw fastq files in order to
                produce aligned BAM files. It depends on the other two
                modules 'aligner' and 'xenograft'.
    * aligner - The module responsible for the alignment of the fastq
                files.
    * xenograft - When a sample derives from xenotransplanted tissues,
                  this module is used in order to split the data between
                  graft and host.
    * mutect - The module responsible for running the MuTect software.
    * varscan - The module responsible for running the VarScan software.
    * strelka - The module responsible for running the Strelka software.
    * variant_calling - the previous three modules produce results that
                        are integrated in this module, in order to
                        produce better and more meaningful data.
    * db - A package containing all the tools to store the results into
           a MongoDB.
    * reports - A packages with two modules related to reports
                generation.
"""

from .core import *
from .reports import *
