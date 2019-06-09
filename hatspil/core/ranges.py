"""Easily handle genomic ranges."""
from typing import (Iterable, List, Optional, Sequence, Tuple, TypeVar, Union,
                    cast)

RangeType = TypeVar("RangeType", bound="Range")
RangesType = TypeVar("RangesType", bound="Ranges")


class Range:
    """A generic range, with a start and an end."""

    def __init__(self, start: Optional[int], end: Optional[int]) -> None:
        """Create a range given start and end, optionally."""
        self.start = start
        self.end = end

    def __len__(self) -> int:
        """Return the len of a valid Range."""
        assert self.start is not None
        assert self.end is not None

        return self.end - self.start

    def intersect(self: RangeType, other: RangeType) -> RangeType:
        """Return the intersection between two valid Range objects."""
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)

        if new_end - new_start > max(len(self), len(other)):
            new_end = 0
            new_start = 0

        return cast(RangeType, Range(new_start, new_end))

    def union(self: RangeType, other: RangeType) -> Union[RangeType, List[RangeType]]:
        """Return the union between two valid Range objects.

        The result type depends on the fact that the two initial ranges
        overlap or not. If the ranges overlap, a single Range is
        returned. Otherwise, a list of two separated ranges is returned.
        """
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        if (self.start < other.start and self.end > other.start) or (
            other.start < self.start and other.end > self.start
        ):
            return cast(
                RangeType, Range(min(self.start, other.start), max(self.end, other.end))
            )
        else:
            return Ranges([self, other])

    def valid(self) -> bool:
        """Check if the range is valid."""
        return self.start != self.end

    def __repr__(self) -> str:
        """Return a string representation for the range."""
        assert self.start is not None
        assert self.end is not None
        return "%d-%d" % (self.start, self.end)

    def __lt__(self: RangeType, other: RangeType) -> bool:
        """Check if a range is less than another.

        A range `a` is less than a range `b` if `a.start < b.start` or
        `a.start == b.start` and `a.end < b.end`.
        """
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        return self.start < other.start or (
            self.start == other.start and self.end < other.end
        )

    def __le__(self: RangeType, other: RangeType) -> bool:
        """Check if a range is less or equal than another.

        A range `a` is less or equal than a range `b` if
        `a.start < b.start` or `a.start == b.start` and
        `a.end <= b.end`.
        """
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        return self.start < other.start or (
            self.start == other.start and self.end <= other.end
        )


class Ranges(List[RangeType]):
    """A useful list of Range."""

    def __init__(self, ranges: Sequence[RangeType] = []) -> None:
        """Create a `Ranges` from a list of Range.

        If the ranges are not not sorted, they will be sorted
        automatically.
        """
        super().__init__(ranges)
        self.resorted = not all(
            self[index] <= self[index + 1] for index in range(len(self) - 1)
        )
        if self.resorted:
            super().sort()

    def overlaps(self: RangesType, other: RangesType) -> List[Tuple[int, int]]:
        """Find the overlaps between two `Ranges`.

        Return a list with the pair of indices. Each pair of indices
        identify the `Range` object that overlaps between the two
        `Ranges`.
        """
        overlaps_indices = []
        self_iter: Iterable[Tuple[int, RangeType]] = iter(enumerate(self))
        for self_index, self_range in self_iter:
            assert self_range.start is not None
            assert self_range.end is not None

            other_iter: Iterable[Tuple[int, RangeType]] = iter(enumerate(other))
            for other_index, other_range in other_iter:
                assert other_range.start is not None
                assert other_range.end is not None

                if other_range.end < self_range.start:
                    continue
                elif self_range.end < other_range.start:
                    break

                if (
                    self_range.start < other_range.start
                    and self_range.end >= other_range.start
                ) or (
                    other_range.start < self_range.start
                    and other_range.end >= self_range.start
                ):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices


class GenomicRange(Range):
    """A range object for genomic data.

    Like a `Range` object, with the information about the chromosome and
    the strand.
    """

    def __init__(
        self,
        chrom: Optional[str],
        start: Optional[int],
        end: Optional[int],
        strand: str = "*",
    ) -> None:
        """Create a genomic range."""
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

    def intersect(self, other: "GenomicRange") -> "GenomicRange":
        """Return the intersection between two genomic ranges."""
        if self.chrom != other.chrom:
            return GenomicRange(None, None, None)

        new_range = super().intersect(other)
        if self.strand == other.strand:
            strand = self.strand
        else:
            strand = "*"

        return GenomicRange(self.chrom, new_range.start, new_range.end, strand)

    def union(self, other: "GenomicRange") -> Union["GenomicRange", "GenomicRanges"]:
        """Return the union between two genomic ranges.

        The result type depends on the fact the two genomic ranges
        intersect or not. See `Range.union` for more information.
        """
        if self.chrom == other.chrom:
            union_result = super().union(other)
            if isinstance(union_result, Range):
                if self.strand == other.strand:
                    strand = self.strand
                else:
                    strand = "*"

                return GenomicRange(
                    self.chrom, union_result.start, union_result.end, strand
                )
            else:
                return GenomicRanges([self, other])
        else:
            return GenomicRanges([self, other])

    def __repr__(self) -> str:
        """Return the string representation for the genomic range."""
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None

        out = "%s:%d-%d" % (self.chrom, self.start, self.end)
        if self.strand != "*":
            out += "(%s)" % self.strand
        return out

    def __lt__(self, other: "GenomicRange") -> bool:
        """Check whether a genomic range is less than another.

        A genomic range `a` is less than the genomic range `b` if
        `a.chrom < b.chrom` or `a.chrom == b.chrom` and
        `Range(a) < Range(b)`. See `Range.__lt__` for more information.
        """
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None
        assert other.chrom is not None
        assert other.start is not None
        assert other.end is not None

        return self.chrom < other.chrom or (
            self.chrom == other.chrom and super().__lt__(other)
        )

    def __le__(self, other: "GenomicRange") -> bool:
        """Check whether a genomic range is less than or equal another.

        A genomic range `a` is less than or equal the genomic range `b`
        if `a.chrom < b.chrom` or `a.chrom == b.chrom` and
        `Range(a) <= Range(b)`. See `Range.__le__` for more information.
        """
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None
        assert other.chrom is not None
        assert other.start is not None
        assert other.end is not None

        return self.chrom < other.chrom or (
            self.chrom == other.chrom and super().__le__(other)
        )


class GenomicRanges(Ranges):
    """A useful list of GenomicRange."""

    def __init__(self, ranges: Sequence[GenomicRange] = []) -> None:
        """Create a GenomicRanges from a raw list."""
        super().__init__(ranges)

    def overlaps(self, other: "GenomicRanges") -> List[Tuple[int, int]]:
        """Find the overlaps between two GenomicRanges.

        See `Ranges.overlaps` for more information.
        """
        overlaps_indices = []
        for self_index, self_range in enumerate(self):
            assert self_range.chrom is not None
            assert self_range.start is not None
            assert self_range.end is not None

            for other_index, other_range in enumerate(other):
                assert other_range.chrom is not None
                assert other_range.start is not None
                assert other_range.end is not None

                if (
                    other_range.chrom < self_range.chrom
                    or other_range.end < self_range.start
                ):
                    continue
                elif (
                    self_range.chrom < other_range.chrom
                    or self_range.end < other_range.start
                ):
                    break

                if (
                    self_range.start < other_range.start
                    and self_range.end >= other_range.start
                ) or (
                    other_range.start < self_range.start
                    and other_range.end >= self_range.start
                ):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices
