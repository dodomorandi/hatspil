import re
from typing import Any, Dict, List, TextIO, Tuple

import hatspil.db
from hatspil.analysis import Analysis
from hatspil.barcoded_filename import BarcodedFilename

from .collection import Collection


class Cutadapt:
    RE_READ_WITH_ADAPTER = re.compile(R"\s*Read (\d) with adapter:\s*(.*)")
    RE_SHORT_PAIRS = re.compile(R"\s*Pairs that were too short:\s*(.*)")
    RE_PAIRS_WRITTEN = re.compile(
        R"\s*Pairs written (passing filters):\s*(.*)")
    RE_TOTAL_BP_PROCESSED = re.compile(
        R"\s*Total basepairs processed:\s*(.*?) bp")
    RE_READ = re.compile(R"\s*Read (\d):\s*(.*?) bp.*")
    RE_TOTAL_WRITTEN = re.compile(R"\s*Total written (filtered):\s*(.*)")
    RE_ADAPTER_HEADER = re.compile(
        R"=== (First|Second) read: Adapter (\d) ===")
    RE_ADAPTER_INFO = re.compile(
        R"Sequence: ([^;]+); Type: ([^;]+); Length: (\d+); "
        R"Trimmed: (\d+) times.")
    RE_ADAPTER_ERROR = re.compile(R"(\d+)(?:-(\d+))? bp: (\d+)")

    def __init__(self, db: "hatspil.db.Db") -> None:
        self.db = db
        self.collection = Collection(db, "cutadapt")

    def store_from_file(self, analysis: Analysis,
                        cutadapt_report_filename: str) -> None:
        if not analysis.config.use_mongodb:
            return

        data = Cutadapt.from_file_to_dict(cutadapt_report_filename)
        barcoded = BarcodedFilename.from_sample(analysis.sample)
        sequencing = self.db.from_barcoded(barcoded)["sequencing"]
        data["sequencing"] = sequencing["_id"]
        data["date"] = analysis.current

        self.collection.find_or_insert({
            "sequencing": sequencing["_id"],
            "date": analysis.current},
            data)

    @staticmethod
    def _parse_summary(fd: TextIO, data: Dict[str, Any]) -> None:
        next(fd)  # Empty line

        line = next(fd)
        data["processed_read_pairs"] =\
            Cutadapt._parse_int_with_comas(line[:27].strip())

        while True:
            line = next(fd)
            match = Cutadapt.RE_READ_WITH_ADAPTER.match(line)
            if not match:
                break

            adapter = int(match.group(1))
            count = Cutadapt._parse_int_with_comas(
                match.group(2).split(" ")[0])

            data[f"processed_read_pairs_adapter_{adapter}"] = count

        match = Cutadapt.RE_SHORT_PAIRS.match(line)
        assert match
        data["pairs_too_short"] = Cutadapt._parse_int_with_comas(
            match.group(1).split(" ")[0])

        line = next(fd)
        match = Cutadapt.RE_PAIRS_WRITTEN.match(line)
        assert match
        data["pairs_written"] = Cutadapt._parse_int_with_comas(
            match.group(1).split(" ")[0])

        next(fd)  # Empty line

        line = next(fd)
        match = Cutadapt.RE_TOTAL_BP_PROCESSED.match(line)
        assert match
        data["total_bp_processed"] = \
            Cutadapt._parse_int_with_comas(match.group(1))

        while True:
            line = next(fd)
            match = Cutadapt.RE_READ.match(line)
            if not match:
                break

            read_index = int(match.group(1))
            count = Cutadapt._parse_int_with_comas(match.group(2))
            data[f"read_{read_index}_bp_processed"] = count

        match = Cutadapt.RE_TOTAL_WRITTEN.match(line)
        assert match
        data["total_bp_written"] = \
            Cutadapt._parse_int_with_comas(
                match.group(1).split(" ")[0])

        while True:
            line = next(fd)
            match = Cutadapt.RE_READ.match(line)
            if not match:
                break

            read_index = int(match.group(1))
            count = Cutadapt._parse_int_with_comas(match.group(2))
            data[f"read_{read_index}_bp_written"] = count

        # After this case there is an empty line, `line` content
        # can be ignored

    @staticmethod
    def _parse_adapter_entry(fd: TextIO, data: Dict[str, Any],
                             adapter_index: int) -> None:
        prefix = f"adapter_{adapter_index}"
        next(fd)  # Empty line

        line = next(fd)
        match = Cutadapt.RE_ADAPTER_INFO.match(line)

        data[f"{prefix}_sequence"] = match.group(1)
        data[f"{prefix}_type"] = match.group(2)
        data[f"{prefix}_length"] = int(match.group(3))
        data[f"{prefix}_trimmed"] = int(match.group(4))

        next(fd)  # Empty line
        next(fd)  # No. of allowed errors

        allowed_errors: Dict[str, int] = {}
        line = next(fd)
        for adapter_error in [x.strip() for x in line.split(";")]:
            match = Cutadapt.RE_ADAPTER_ERROR.match(adapter_error)
            assert match

            if match.group(2):
                interval = (int(match.group(1)), int(match.group(2)))
            else:
                base = int(match.group(1))
                interval = (base, base)
            errors = int(match.group(3))
            allowed_errors["%d-%d" % interval] = errors
        data[f"{prefix}_allowed_errors"] = allowed_errors

        next(fd)  # Empty line
        next(fd)  # Bases preceding removed adapters

        bases_preceding_removed_adapter: Dict[str, float] = {}
        for _ in range(5):
            line = next(fd)
            base, fraction = [x.strip() for x in line.split(":")]
            if len(base) != 1:
                base = "other"
            bases_preceding_removed_adapter[base] = float(fraction[:-1]) / 100
        data[f"{prefix}_preceding_bases"] = \
            bases_preceding_removed_adapter

        next(fd)  # Empty line
        next(fd)  # Overview of removed sequences

        removed_sequences: Dict[str, List[Any]] = {"header": [], "data": []}
        removed_sequences["header"] = next(fd).strip().split("\t")

        while True:
            line = next(fd).strip()
            if not line:
                break

            params = line.split("\t")
            length = int(params[0])
            count = int(params[1])
            expect = float(params[2])
            max_err = int(params[3])
            error_counts = [int(x) for x in params[4].split(" ")]

            removed_sequences["data"].append(
                [length, count, expect, max_err, error_counts])

        data[f"{prefix}_removed_sequences"] = removed_sequences

    @staticmethod
    def from_file_to_dict(cutadapt_report_filename: str) -> Dict[str, Any]:
        data: Dict[str, Any] = {}

        with open(cutadapt_report_filename) as fd:
            for line in fd:
                if line == "=== Summary ===":
                    Cutadapt._parse_summary(fd, data)
                    continue

                match = Cutadapt.RE_ADAPTER_HEADER.match(line)
                if match:
                    adapter_index = int(match.group(2))
                    Cutadapt._parse_adapter_entry(fd, data, adapter_index)

        return data

    @staticmethod
    def _parse_int_with_comas(raw: str) -> int:
        return int(raw.replace(",", ""))
