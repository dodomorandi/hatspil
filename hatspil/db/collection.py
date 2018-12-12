from typing import Any, Dict, List, Optional, cast

import hatspil.db

try:
    from pymongo import ReturnDocument
    from pymongo.cursor import Cursor
except Exception:
    pass


class Collection:
    db: Optional["hatspil.db.Db"]

    def __init__(self, db: "hatspil.db.Db", collection_name: str) -> None:
        self.collection_name = collection_name
        if db.config.use_mongodb:
            self.db = db
            self.collection = db.db[collection_name]
        else:
            self.db = None
            self.collection = None

    def find_or_insert(
        self, data: Dict[Any, Any], new_data: Optional[Dict[Any, Any]] = None
    ) -> Optional[Dict[str, Any]]:
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None

        set_data = dict(data)
        if new_data is not None:
            set_data.update(new_data)
        retval = self.collection.find_one_and_update(
            data, {"$set": set_data}, upsert=True, return_document=ReturnDocument.AFTER
        )
        return cast(Optional[Dict[str, Any]], retval)

    def find(self, data: Dict[Any, Any]) -> Optional[Dict[str, Any]]:
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        retval = self.collection.find_one(data)
        return cast(Optional[Dict[str, Any]], retval)

    def find_all(
        self, data: Optional[Dict[Any, Any]] = None
    ) -> Optional[List[Dict[str, Any]]]:
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        return [cast(Dict[str, Any], element) for element in self.collection.find(data)]

    def iter(self, data: Optional[Dict[Any, Any]] = None) -> Optional["Cursor"]:
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        return self.collection.find(data)
