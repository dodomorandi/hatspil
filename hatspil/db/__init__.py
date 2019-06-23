"""A package to easily handle a MongoDB for the results.

Saving the information obtained during the steps of an analysis is very
crucial, and doing it in a way that the data can be retrieved easily
from a database generally requires an abstraction layer.

This package aims to be this layer. All the serialization and
deserialization from/to the MongoDB are handled by modules inside this
package.

There are 4 modules available:
* db: the main and common database helper features;
* collection: the abstration layer for the collections of the DB;
* cutadapt: the easy way to handle results from Cutadapt;
* picard_metrics: the easy way to handle the metrics files obtained from
                  Picard.
"""
__all__ = ["db", "cutadapt", "picard_metrics"]

from .db import Db
