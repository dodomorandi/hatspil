from typing import (Iterable, List, Optional, Sequence, Tuple, TypeVar, Union,
                    cast)

RangeType = TypeVar("RangeType", bound="Range")
RangesType = TypeVar("RangesType", bound="Ranges")


class Range:
    def __init__(self, start: Optional[int], end: Optional[int]) -> None:
        self.start = start
        self.end = end

    def __len__(self) -> int:
        assert self.start is not None
        assert self.end is not None

        return self.end - self.start

    def intersect(self: RangeType, other: RangeType) -> RangeType:
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

    def union(self: RangeType, other: RangeType)\
            -> Union[RangeType, List[RangeType]]:
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        if (self.start < other.start and self.end > other.start) or\
                (other.start < self.start and other.end > self.start):
            return cast(RangeType,
                        Range(min(self.start, other.start),
                              max(self.end, other.end)))
        else:
            return Ranges([self, other])

    def valid(self) -> bool:
        return self.start != self.end

    def __repr__(self) -> str:
        assert self.start is not None
        assert self.end is not None
        return "%d-%d" % (self.start, self.end)

    def __lt__(self: RangeType, other: RangeType) -> bool:
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        return self.start < other.start or \
            (self.start == other.start and self.end < other.end)

    def __le__(self: RangeType, other: RangeType) -> bool:
        assert self.start is not None
        assert self.end is not None
        assert other.start is not None
        assert other.end is not None

        return self.start < other.start or \
            (self.start == other.start and self.end <= other.end)


class Ranges(List[RangeType]):
    def __init__(self, ranges: Sequence[RangeType] = []) -> None:
        super().__init__(ranges)
        self.resorted = not all(self[index] <= self[index + 1]
                                for index in range(len(self) - 1))
        if self.resorted:
            super().sort()

    def overlaps(self: RangesType, other: RangesType) -> List[Tuple[int, int]]:
        overlaps_indices = []
        self_iter: Iterable[Tuple[int, RangeType]] = iter(enumerate(self))
        for self_index, self_range in self_iter:
            assert self_range.start is not None
            assert self_range.end is not None

            other_iter: Iterable[Tuple[int, RangeType]] = iter(
                enumerate(other))
            for other_index, other_range in other_iter:
                assert other_range.start is not None
                assert other_range.end is not None

                if other_range.end < self_range.start:
                    continue
                elif self_range.end < other_range.start:
                    break

                if (self_range.start < other_range.start
                    and self_range.end >= other_range.start)\
                        or (other_range.start < self_range.start
                            and other_range.end >= self_range.start):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices


class GenomicRange(Range):
    def __init__(self,
                 chrom: Optional[str],
                 start: Optional[int],
                 end: Optional[int],
                 strand: str = "*") -> None:
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

    def intersect(self, other: "GenomicRange") -> "GenomicRange":
        if self.chrom != other.chrom:
            return GenomicRange(None, None, None)

        new_range = super().intersect(other)
        if self.strand == other.strand:
            strand = self.strand
        else:
            strand = "*"

        return GenomicRange(self.chrom,
                            new_range.start,
                            new_range.end,
                            strand)

    def union(self, other: "GenomicRange")\
            -> Union["GenomicRange", "GenomicRanges"]:
        if self.chrom == other.chrom:
            union_result = super().union(other)
            if isinstance(union_result, Range):
                if self.strand == other.strand:
                    strand = self.strand
                else:
                    strand = "*"

                return GenomicRange(self.chrom, union_result.start,
                                    union_result.end, strand)
            else:
                return GenomicRanges([self, other])
        else:
            return GenomicRanges([self, other])

    def __repr__(self) -> str:
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None

        out = "%s:%d-%d" % (self.chrom, self.start, self.end)
        if self.strand != "*":
            out += "(%s)" % self.strand
        return out

    def __lt__(self, other: "GenomicRange") -> bool:
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None
        assert other.chrom is not None
        assert other.start is not None
        assert other.end is not None

        return self.chrom < other.chrom or\
            (self.chrom == other.chrom and super().__lt__(other))

    def __le__(self, other: "GenomicRange") -> bool:
        assert self.chrom is not None
        assert self.start is not None
        assert self.end is not None
        assert other.chrom is not None
        assert other.start is not None
        assert other.end is not None

        return self.chrom < other.chrom or\
            (self.chrom == other.chrom and super().__le__(other))


class GenomicRanges(Ranges):
    def __init__(self, ranges: Sequence[GenomicRange] = []) -> None:
        super().__init__(ranges)

    def overlaps(self, other: "GenomicRanges") -> List[Tuple[int, int]]:
        overlaps_indices = []
        for self_index, self_range in enumerate(self):
            assert self_range.chrom is not None
            assert self_range.start is not None
            assert self_range.end is not None

            for other_index, other_range in enumerate(other):
                assert other_range.chrom is not None
                assert other_range.start is not None
                assert other_range.end is not None

                if other_range.chrom < self_range.chrom\
                        or other_range.end < self_range.start:
                    continue
                elif self_range.chrom < other_range.chrom\
                        or self_range.end < other_range.start:
                    break

                if (self_range.start < other_range.start
                        and self_range.end >= other_range.start)\
                        or (other_range.start < self_range.start
                            and other_range.end >= self_range.start):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices
