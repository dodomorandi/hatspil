"""An abstraction layer to handle tables in a report.

Reports need some HTML, therefore an abstraction layer is extremely
useful to create and customize a table before the HTML code is
generated.

This module contains some classes that ease the process.
"""
from collections import OrderedDict
from enum import Enum
from typing import (Any, Callable, Iterable, List, MutableMapping, Optional,
                    Sequence, Union)


class ReportTableSortingType(Enum):
    """The sorting of a table column."""

    ASCENDING = "asc"
    DESCENDING = "des"


class ReportTableColumnOrdering:
    """Helper class to store the sorting for a specific table column."""

    def __init__(self, column_index: int, sorting_type: ReportTableSortingType) -> None:
        """Create an instance of the class."""
        self.column_index = column_index
        self.sorting_type = sorting_type

    def html(self) -> str:
        """Create the configuration string to put inside the HTML."""
        return '[{}, "{}"]'.format(self.column_index, self.sorting_type.value)


class ReportTable:
    """A helper class to handle a table inside a report."""

    Entry = Union[str, int, float, bool]
    # This is an OrderedDict, but the annotation is not accepted in Python 3.6
    columns: MutableMapping[str, List[Entry]]

    def __init__(self, html_id: str, *column_names: str) -> None:
        """Create an instance of the class.

        Args:
            html_id: the 'id' of the html table that will be created.
            column_names: the name for each column that will be created
                          in the table.
        """
        self.html_id = html_id
        self._order: Optional[List[ReportTableColumnOrdering]] = None
        self._row_modifier: Optional[
            Callable[[MutableMapping[str, ReportTable.Entry]], str]
        ] = None
        self._table_class: Optional[str] = None
        self.style: Optional[str] = None

        self.set_columns(column_names)

    def set_columns(self, column_names: Iterable[str]) -> None:
        """Set the name of the columns.

        This function overrides any previous names, therefore it is a
        good idea to call this function once per instance.
        """
        self.columns: MutableMapping[str, ReportTable.Entry] = OrderedDict()
        for column_name in column_names:
            self.columns[column_name] = []

    def add_row(self, row: Sequence[Entry]) -> None:
        """Add a row to the table.

        The number of elements in the row must be the same as the length
        of the columns. It is necessary to call
        `ReportTable.set_columns` beforehand to perform a useful
        operation.

        Args:
            row: a sequence of entries. Each entry can be a str, an int,
                 a float or a bool.
        """
        assert all(map(lambda column: column is not None, self.columns.values()))
        assert len(row) == len(self.columns)

        for column, entry in zip(self.columns.values(), row):
            column.append(entry)

    @property
    def order(self) -> Optional[List[ReportTableColumnOrdering]]:
        """Get the order of the columns.

        The value is set using `ReportTable.set_order`.
        """
        return self._order

    def set_order(
        self,
        order: Union[Iterable[ReportTableColumnOrdering], ReportTableColumnOrdering],
    ) -> None:
        """Set the order of the columns.

        The arguments will be always converted into a list of ordering
        elements.

        Args:
            order: a single column ordering or an iterable of column
                   ordering objects. In both cases a list is saved
                   internally.
        """
        if isinstance(order, ReportTableColumnOrdering):
            self._order = [order]
        else:
            self._order = list(order)

            if not self._order:
                self._order = None

    @property
    def row_modifier(self) -> Optional[Callable[[MutableMapping[str, Entry]], str]]:
        """Get the row modifier function, if any.

        See the documentation of the setter for more information.
        """
        return self._row_modifier

    @row_modifier.setter
    def row_modifier(
        self, modifier: Optional[Callable[[MutableMapping[str, Entry]], str]]
    ) -> None:
        """Set the row modifier function.

        The `modifier` function is used to add some HTML properties to
        each 'tr' element depending on the content of the row. This is
        mostly useful to change the style of the row when one or more
        properties have interesting values.
        """
        self._row_modifier = modifier

    @property
    def table_class(self) -> Optional[str]:
        """Get the table HTML class."""
        return self._table_class

    @table_class.setter
    def table_class(self, table_class: Optional[str]) -> None:
        """Set the table HTML class."""
        self._table_class = table_class

    def html(self) -> str:
        """Create the HTML representation of the table.

        In case no columns have been set, an empty str is returned.
        """
        if not self.columns:
            return ""

        def default_row_modifier(_: Any) -> str:
            return ""

        if self._row_modifier:
            row_modifier = self._row_modifier
        else:
            row_modifier = default_row_modifier

        html_rows: List[str] = []
        if self.columns:
            first_column_len = len(next(iter(self.columns.values())))
            assert all(
                map(
                    lambda column: len(column) == first_column_len,
                    self.columns.values(),
                )
            )

        for row_index in range(len(next(iter(self.columns.values())))):
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
            )
            if self._order
            else "",
        )

    def _get_row(self, row_index: int) -> MutableMapping[str, Entry]:
        assert self.columns
        assert row_index >= 0
        assert row_index < len(next(iter(self.columns.values())))

        row_content: MutableMapping[str, ReportTable.Entry] = OrderedDict()
        for col_name, entries in self.columns.items():
            row_content[col_name] = entries[row_index]
        return row_content
        assert row_index < len(self.columns)
