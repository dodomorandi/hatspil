from collections import OrderedDict
from enum import Enum
from typing import (Any, Callable, Iterable, List, Optional, Sequence,
                    Union)


class ReportTableSortingType(Enum):
    ASCENDING = "asc"
    DESCENDING = "des"


class ReportTableColumnOrdering:
    def __init__(self, column_index: int, sorting_type: ReportTableSortingType) -> None:
        self.column_index = column_index
        self.sorting_type = sorting_type

    def html(self) -> str:
        return "[{}, {}]".format(self.column_index, str(self.sorting_type))


class ReportTable:
    Entry = Union[str, int, float]
    columns: OrderedDict[str, List[Entry]]

    def __init__(self, html_id: str, *column_names: str) -> None:
        self.html_id = html_id
        self._order: Optional[List[ReportTableColumnOrdering]] = None
        self._row_modifier: Optional[Callable[[OrderedDict], str]] = None
        self._table_class: Optional[str] = None
        self.style: Optional[str] = None

        self.set_columns(*column_names)

    def set_columns(self, column_names: Iterable[str]) -> None:
        self.columns = OrderedDict()
        for column_name in column_names:
            self.columns[column_name] = []

    def add_row(self, row: Sequence[Entry]) -> None:
        assert all(map(lambda column: not column, self.columns.values()))
        assert len(row) == len(self.columns)

        for column, entry in zip(self.columns.values(), row):
            column.append(entry)

    @property
    def order(self) -> Optional[List[ReportTableColumnOrdering]]:
        return self.order

    @order.setter
    def order(
        self,
        order: Union[Iterable[ReportTableColumnOrdering], ReportTableColumnOrdering],
    ) -> None:
        if isinstance(order, ReportTableColumnOrdering):
            self._order = [order]
        else:
            self._order = list(order)

            if not self.order:
                self._order = None

    @property
    def row_modifier(self) -> Optional[Callable[[OrderedDict[str, Entry]], str]]:
        return self._row_modifier

    @row_modifier.setter
    def row_modifier(
        self, modifier: Optional[Callable[[OrderedDict[str, Entry]], str]]
    ) -> None:
        self._row_modifier = modifier

    @property
    def table_class(self) -> Optional[str]:
        return self._table_class

    @table_class.setter
    def table_class(self, table_class: Optional[str]) -> None:
        self._table_class = table_class

    def html(self) -> str:
        if not self.columns:
            return ""

        def default_row_modifier(_: Any) -> str:
            return ""

        if self._row_modifier:
            row_modifier = self._row_modifier
        else:
            row_modifier = default_row_modifier

        html_rows: List[str] = []
        for row_index in range(len(next(iter(self.columns)))):
            row = self._get_row(row_index)
            html_rows.append(
                "<tr{}>{}</tr>".format(
                    row_modifier(row),
                    "".join(f"<td>{value}</td>" for value in row.values()),
                )
            )

        return (
            '<table id="{}"{}>'
            "<thead><tr>{}</tr></thead>"
            "{}</table>"
            "<script>$(document).ready( function () {{"
            "$('#{}').DataTable({});"
            "}} );"
            "</script>"
        ).format(
            self.html_id,
            f' class="{self.table_class}"' if self.table_class else "",
            "".join([f"<th>{col_name}</th>" for col_name in self.columns.keys()]),
            "".join(html_rows),
            self.html_id,
            "{{'order': [{}]}}".format(
                ",".join([order.html() for order in self._order])
            ) if self._order else "",
        )

    def _get_row(self, row_index: int) -> OrderedDict[str, Entry]:
        assert self.columns
        assert row_index >= 0
        assert row_index < len(next(iter(self.columns)))

        row_content: OrderedDict[str, ReportTable.Entry] = OrderedDict()
        for col_name, entries in self.columns.items():
            row_content[col_name] = entries[row_index]
        return row_content
        assert row_index < len(self.columns)
