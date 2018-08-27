import hatspil.db

from typing import Dict
from pymongo import ReturnDocument


class Collection:
    def __init__(self, db: "hatspil.db.Db", collection_name: str) -> None:
        self.collection_name = collection_name
        if db.config.use_mongodb:
            self.db = db
            self.collection = db.db[collection_name]
        else:
            self.db = None
            self.collection = None

    def find_or_insert(self, data: Dict, new_data: Dict = None):
        if not self.db.config.use_mongodb:
            return None

        assert self.collection is not None

        set_data = dict(data)
        if new_data is not None:
            set_data.update(new_data)
        return self.collection.find_one_and_update(
            data, {"$set": set_data}, upsert=True, return_document=ReturnDocument.AFTER
        )

    def find(self, data: Dict):
        if not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        return self.collection.find_one(data)
