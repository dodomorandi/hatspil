"""The module for the abstraction of a MongoDB collection."""
from typing import Any, Dict, List, Optional, cast

import hatspil.db

try:
    from pymongo import ReturnDocument
    from pymongo.cursor import Cursor
except Exception:
    pass


class Collection:
    """A class to help with MondoDB collections.

    Working directly with MongoDB data and functions can be
    overwhelming. This class helps with the handling of MongoDB
    collections.
    """

    db: Optional["hatspil.db.Db"]

    def __init__(self, db: "hatspil.db.Db", collection_name: str) -> None:
        """Create a new instance of the class.

        Args:
            db: an instance of `Db`.
            collection_name: the name of the collection in the database.
        """
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
        """Find an element and insert it if not found.

        Search for some data inside the collection. If it is found, the
        content is returned, otherwise a new element is inserted into
        the collection and it is returned.

        Args:
            data: the things to match inside the collection. All the
                  MongoDB special syntax can be used.
            new_data: the data that needs to be inserted into the
                      database if `data` is not found. It does not need
                      to include the parameters already in `data`. If
                      `None`, only the properties of `data` are
                      inserted.

        Returns:
            The content in the database, after the eventual insertion.
            `None` if the MongoDB cannot be used.

        """
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
        """Find an element in the collection.

        Searches for an element in the collection and returns it in case
        it is found.

        Args:
            data: the things to match inside the collection. All the
                  MongoDB special syntax can be used.

        Returns:
            The element found in the database or `None`.
            `None` if the MongoDB cannot be used.

        """
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        retval = self.collection.find_one(data)
        return cast(Optional[Dict[str, Any]], retval)

    def find_all(
        self, data: Optional[Dict[Any, Any]] = None
    ) -> Optional[List[Dict[str, Any]]]:
        """Find all the matching elements in the collection.

        Searches for all the elements in the collection that match the
        input data. See `Collection.iter` for a version that returns a
        cursor (an iterable object).

        Args:
            data: the things to match inside the collection. All the
                  MongoDB special syntax can be used.

        Returns:
            A list of the elements found in the collection.
            `None` if the MongoDB cannot be used.

        """
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        return [cast(Dict[str, Any], element) for element in self.collection.find(data)]

    def iter(self, data: Optional[Dict[Any, Any]] = None) -> Optional["Cursor"]:
        """Iterate among all the matching elements in the collection.

        This is the version of `Collection.find_all` that returns a
        cursor instead of a list. Use this function if many results are
        expected.

        See `Collection.find_all` for more information.
        """
        if not self.db or not self.db.config.use_mongodb:
            return None

        assert self.collection is not None
        return self.collection.find(data)
