from pymongo import MongoClient, ReturnDocument

from .config import Config


class Collection:
    def __init__(self, db, collection_name):
        assert isinstance(db, Db)

        self.collection_name = collection_name
        if not db.config.use_db:
            self.collection = None
            return

        self.db = db
        self.collection = db.db[collection_name]

    def find_or_insert(self, data, new_data=None):
        if not self.db.config.use_mongodb:
            return None

        assert self.collection is not None

        set_data = dict(data)
        if new_data is not None:
            set_data.update(new_data)
        return self.collection.find_one_and_update(
            data, {"$set": set_data},
            upsert=True,
            return_document=ReturnDocument.AFTER)

    def find(self, data):
        if not self.db.config.use_mongodb:
            return None

        assert self.collection is not None

        set_data = dict(data)
        return self.collection.find_one(data, {"$set": set_data})


class Db:
    _collections = [
        "projects", "patients", "biopsies", "samples", "sequencings",
        "annotations", "variants", "analyses"
    ]

    def __init__(self, config):
        assert isinstance(config, Config)

        self.config = config
        for collection_name in Db._collections:
            setattr(self, f"_{collection_name}",
                    Collection(self, collection_name))

        if config.use_mongodb:
            mongo = MongoClient(config.mongodb_host, config.mongodb_port)
            self.db = mongo[config.mongodb_database]
            self.db.authenticate(config.mongodb_username,
                                 config.mongodb_password)
        else:
            self.db = None

    def store_barcoded(self, barcoded):
        if not self.config.use_mongodb:
            return None

        project = self.projects.find_or_insert({"name": barcoded.project})

        patient = self.patient.find_or_insert({
            "project": project["_id"],
            "name": barcoded.patient
        })

        biopsy = self.biopsies.find_or_insert({
            "patient": patient["_id"],
            "index": barcoded.biopsy,
            "tissue": barcoded.tissue
        })

        sample_data = {"biopsy": biopsy["_id"]}
        if barcoded.xenograft is None:
            sample_data["index"] = barcoded.sample
        else:
            sample_data["xenograft"] = barcoded.xenograft.to_dict()
        sample = self.samples.find_or_insert(sample_data)

        sequencing = self.sequencings.find_or_insert(
            {
                "sample": sample["_id"],
                "index": barcoded.sequencing
            }, {
                "molecule": barcoded.molecule,
                "analyte": barcoded.analyte,
                "kit": barcoded.kit
            })

        return {
            "project": project,
            "patient": patient,
            "biopsy": biopsy,
            "sample": sample,
            "sequencing": sequencing
        }

    def from_barcoded(self, barcoded):
        if not self.config.use_mongodb:
            return None

        project = self.projects.find({"name": barcoded.project})
        patient = self.patient.find({
            "project": project["_id"],
            "name": barcoded.patient
        })
        biopsy = self.biopsies.find({
            "patient": patient["_id"],
            "index": barcoded.biopsy
        })
        sample = self.samples.find({
            "biopsy": biopsy["_id"],
            "index": barcoded.sample
        })
        sequencing = self.sequencings.find({
            "sample": sample["_id"],
            "index": barcoded.sequencing
        })

        return {
            "project": project,
            "patient": patient,
            "biopsy": biopsy,
            "sample": sample,
            "sequencing": sequencing
        }
