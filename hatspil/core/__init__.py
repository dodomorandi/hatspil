"""
The core of HaTSPiL.

The modules contained in this package are the hard core of HaTSPiL. Indeed, the functionalities provided are mostly meant to be used from other modules or other software, because extending HaTPSiL should not involve changing the modules of 'core'.

The modules contained in this package are the following:
    * barcoded_filename - Everything related to encoding and decoding the information of a filename barcode is in this module.
    * analysis - The module containing the information that needs to be shared between different steps of the analysis for a sample.
    * starter - A helper module needed to bootstrap an analysis.
    * executor - The real core of HaTSPiL. This module contains the class responsible to run each command of the analyses.
    * utils - A set of utility functions.
    * ranges - A module useful to perform genomic ranges operations.
    * exceptions - A small set of custom exceptions.
"""
