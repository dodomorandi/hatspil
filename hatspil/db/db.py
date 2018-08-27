from typing import Any, Dict, Optional

from hatspil.barcoded_filename import BarcodedFilename
from hatspil.config import Config
from pymongo import MongoClient

from .collection import Collection


class Db:
    _COLLECTIONS = [
        "projects",
        "patients",
        "biopsies",
        "samples",
        "sequencings",
        "annotations",
        "variants",
        "analyses",
        "cutadapt",
    ]
    projects: Collection
    patients: Collection
    biopsies: Collection
    samples: Collection
    sequencings: Collection
    annotations: Collection
    variants: Collection
    analyses: Collection
    cutadapt: Collection

    def __init__(self, config: Config) -> None:
        self.config = config

        if config.use_mongodb:
            mongo = MongoClient(config.mongodb_host, config.mongodb_port)
            self.db = mongo[config.mongodb_database]
            self.db.authenticate(config.mongodb_username, config.mongodb_password)
        else:
            self.db = None

        for collection_name in Db._COLLECTIONS:
            setattr(self, collection_name, Collection(self, collection_name))

    def store_barcoded(self, barcoded: BarcodedFilename) -> Optional[Dict[str, Any]]:
        if not self.config.use_mongodb:
            return None

        project = self.projects.find_or_insert({"name": barcoded.project})
        if not project:
            return None

        patient = self.patients.find_or_insert(
            {"project": project["_id"], "name": barcoded.patient}
        )
        if not patient:
            return None

        biopsy = self.biopsies.find_or_insert(
            {
                "patient": patient["_id"],
                "index": barcoded.biopsy,
                "tissue": int(barcoded.tissue),
            }
        )
        if not biopsy:
            return None

        sample_data = {"biopsy": biopsy["_id"]}
        if barcoded.xenograft is None:
            assert not barcoded.tissue.is_xenograft()
            sample_data["index"] = barcoded.sample
        else:
            assert barcoded.tissue.is_xenograft()
            sample_data["xenograft"] = barcoded.xenograft.to_dict()
        sample = self.samples.find_or_insert(sample_data)
        if not sample:
            return None

        sequencing = self.sequencings.find_or_insert(
            {"sample": sample["_id"], "index": barcoded.sequencing},
            {
                "molecule": int(barcoded.molecule),
                "analyte": int(barcoded.analyte),
                "kit": barcoded.kit,
            },
        )
        if not sequencing:
            return None

        return {
            "project": project,
            "patient": patient,
            "biopsy": biopsy,
            "sample": sample,
            "sequencing": sequencing,
        }

    def from_barcoded(
        self, barcoded: BarcodedFilename
    ) -> Optional[Dict[str, Dict[str, Any]]]:
        if not self.config.use_mongodb:
            return None

        project = self.projects.find({"name": barcoded.project})
        if not project:
            return None

        patient = self.patients.find(
            {"project": project["_id"], "name": barcoded.patient}
        )
        if not patient:
            return None

        biopsy = self.biopsies.find(
            {
                "patient": patient["_id"],
                "index": barcoded.biopsy,
                "tissue": int(barcoded.tissue),
            }
        )
        if not biopsy:
            return None

        sample_data = {"biopsy": biopsy["_id"]}
        if barcoded.xenograft is None:
            assert not barcoded.tissue.is_xenograft()
            sample_data["index"] = barcoded.sample
        else:
            assert barcoded.tissue.is_xenograft()
            sample_data["xenograft"] = barcoded.xenograft.to_dict()
        sample = self.samples.find(sample_data)
        if not sample:
            return None

        sequencing = self.sequencings.find(
            {"sample": sample["_id"], "index": barcoded.sequencing}
        )
        if not sequencing:
            return None

        return {
            "project": project,
            "patient": patient,
            "biopsy": biopsy,
            "sample": sample,
            "sequencing": sequencing,
        }
