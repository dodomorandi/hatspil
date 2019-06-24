"""Module to help deserializing Picard results into a MongoDB."""
from collections import OrderedDict
from enum import Enum, auto
from typing import Any, Dict, Union, cast

import hatspil.db

from ..analysis import Analysis
from ..barcoded_filename import BarcodedFilename
from ..utils import parse_as_number
from .collection import Collection


class PicardMetricsType(Enum):
    """The type of a Picard Metrics file."""

    hs = auto()
    gcbias = auto()
    marked_duplicates = auto()
    umi = auto()
    no_duplicates = auto()


class PicardMetrics:
    """A class to help deserializing the Picard Metrics file into a DB.

    Using this class it is possible to store the different kind of
    Picard metrics files into the MongoDB using the same abstraction.
    """

    def __init__(self, db: "hatspil.db.Db") -> None:
        """Create an instance of the class."""
        self.db = db
        self.collection = Collection(db, "picard_metrics")

    def store_from_file(
        self, analysis: Analysis, filename: str, metrics_type: PicardMetricsType
    ) -> None:
        """Deserialize the content of a Picard Metrics file into the DB.

        The type of the metrics is used in the database to store
        multiple metrics data for the same sample.

        If the MongoDB cannot be used, no operation is performed.

        Args:
            analysis: the current analysis.
            filename: the name of the Picard Metrics file.
            metrics_type: the type of the Picard Metrics stored in the
                          file.
        """
        if not analysis.config.use_mongodb or analysis.run_fake:
            return

        data = cast(
            Dict[str, Union[str, Dict[str, Any]]],
            PicardMetrics.from_file_to_dict(filename),
        )
        barcoded = BarcodedFilename.from_sample(analysis.sample)
        db_from_barcoded = self.db.from_barcoded(barcoded)
        assert db_from_barcoded
        sequencing = db_from_barcoded["sequencing"]
        data["sequencing"] = sequencing["_id"]
        data["date"] = analysis.current
        data["type"] = metrics_type.name

        self.collection.find_or_insert(
            {
                "sequencing": sequencing["_id"],
                "date": analysis.current,
                "type": metrics_type.name,
            },
            data,
        )

    @staticmethod
    def from_file_to_dict(picard_generated_filename: str) -> Dict[str, Dict[str, Any]]:
        """Deserialize a Picard Metrics file into a dict."""
        data: Dict[str, Dict[str, Any]] = {}

        with open(picard_generated_filename) as fd:
            line = next(fd)
            # Search the first empty line
            while line.strip():
                line = next(fd)
                continue

            file_finished = False
            while not file_finished:
                try:
                    line = next(fd).strip()
                except StopIteration:
                    break

                if not line:
                    continue

                assert line[:3] == "## "
                splitted_line = line[3:].split("\t")
                section_name = splitted_line[0].lower()
                data[section_name] = OrderedDict()
                section = data[section_name]

                file_finished = True
                for line in fd:
                    line = line.strip()
                    if not line:
                        file_finished = False
                        break

                    if line.startswith("#"):
                        continue

                    splitted_line = line.split("\t")

                    if not section:
                        content_header = [param.lower() for param in splitted_line]
                        data[section_name] = OrderedDict(
                            ((param, []) for param in content_header)
                        )
                        section = data[section_name]
                    else:
                        values = [parse_as_number(value) for value in splitted_line]
                        for param_name, param_value in zip(content_header, values):
                            section[param_name].append(param_value)

                if all(map(lambda l: len(l) <= 1, section.values())):
                    data[section_name] = OrderedDict(
                        (key, values[0] if values else "")
                        for key, values in section.items()
                    )

        return data
