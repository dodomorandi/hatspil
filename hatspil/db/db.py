from typing import Any, Dict, Optional

from bson import ObjectId
from hatspil.barcoded_filename import BarcodedFilename, Xenograft
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
        "picard_metrics",
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
    picard_metrics: Collection

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

        if not barcoded.molecule or not barcoded.analyte:
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

        if not barcoded.molecule or not barcoded.analyte:
            return None

        sequencing = self.sequencings.find(
            {
                "sample": sample["_id"],
                "index": barcoded.sequencing,
                "analyte": int(barcoded.analyte),
                "molecule": int(barcoded.molecule),
            }
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

    def from_sequencing_id(
        self, sequencing_id: ObjectId
    ) -> Optional[Dict[str, Dict[str, Any]]]:
        if not self.config.use_mongodb:
            return None

        sequencing = self.sequencings.find({"_id": sequencing_id})
        if not sequencing:
            return None

        sample = self.samples.find({"_id": sequencing["sample"]})
        if not sample:
            return None

        biopsy = self.biopsies.find({"_id": sample["biopsy"]})
        if not biopsy:
            return None

        patient = self.patients.find({"_id": biopsy["patient"]})
        if not patient:
            return None

        project = self.projects.find({"_id": patient["project"]})
        if not project:
            return None

        return {
            "project": project,
            "patient": patient,
            "biopsy": biopsy,
            "sample": sample,
            "sequencing": sequencing,
        }

    @staticmethod
    def to_barcoded(data: Dict[str, Any]) -> Optional[BarcodedFilename]:
        project = data["project"]
        assert project
        project_name = project["name"]

        patient = data["patient"]
        assert patient
        patient_name = patient["name"]

        biopsy = data["biopsy"]
        assert biopsy
        biopsy_index = biopsy["index"]
        tissue = biopsy["tissue"]

        sample = data["sample"]
        assert sample
        xenograft_generation: Optional[int]
        xenograft_parent: Optional[int]
        xenograft_child: Optional[int]
        sample_index: Optional[str]
        if "xenograft" in sample:
            xenograft = Xenograft.from_dict(sample["xenograft"])
            assert xenograft
            xenograft_generation = xenograft.generation
            xenograft_parent = xenograft.parent
            xenograft_child = xenograft.child
            sample_index = None

            assert xenograft_generation is not None
            assert xenograft_parent is not None
            assert xenograft_child is not None
        else:
            sample_index = sample["index"]
            xenograft_generation = None
            xenograft_parent = None
            xenograft_child = None

            assert sample_index is not None and sample_index != ""

        sequencing = data["sequencing"]
        assert sequencing
        sequencing_index = sequencing["index"]
        molecule = sequencing["molecule"]
        analyte = sequencing["analyte"]
        # At to date, kit  is unused and can be inexistent
        if "kit" in sequencing:
            kit = sequencing["kit"]
        else:
            kit = 0

        assert project_name
        assert patient_name
        assert biopsy_index is not None
        assert tissue
        assert molecule is not None
        assert analyte is not None
        assert kit is not None

        return BarcodedFilename.from_parameters(
            project_name,
            patient_name,
            tissue,
            biopsy_index,
            sample_index,
            xenograft_generation,
            xenograft_parent,
            xenograft_child,
            sequencing_index,
            molecule,
            analyte,
            kit,
        )
