class Range:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    def intersect(self, other):
        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)

        if new_end - new_start > max(len(self), len(other)):
            new_end = 0
            new_start = 0

        return Range(new_start, new_end)

    def union(self, other):
        if (self.start < other.start and self.end > other.start) or\
                (other.start < self.start and other.end > self.start):
            return Range(min(self.start, other.start),
                         max(self.end, other.end))
        else:
            return Ranges([self, other])

    def valid(self):
        return self.start != self.end

    def __repr__(self):
        return "%d-%d" % (self.start, self.end)

    def __lt__(self, other):
        return self.start < other.start or \
            (self.start == other.start and self.end < other.end)

    def __le__(self, other):
        return self.start < other.start or \
            (self.start == other.start and self.end <= other.end)


class Ranges(list):
    def __init__(self, ranges=[]):
        super().__init__(ranges)
        self.resorted = not all(self[index] <= self[index + 1]
                                for index in range(len(self) - 1))
        if self.resorted:
            super().sort()

    def overlaps(self, other):
        overlaps_indices = []
        for self_index, self_range in enumerate(self):
            for other_index, other_range in enumerate(other):
                if other_range.end < self_range.start:
                    continue
                elif self_range.end < other_range.start:
                    break

                if (self_range.start < other_range.start and
                        self_range.end >= other_range.start) or\
                        (other_range.start < self_range.start and
                         other_range.end >= self_range.start):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices


class GenomicRange(Range):
    def __init__(self, chrom, start, end, strand="*"):
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

    def intersect(self, other):
        if self.chrom != other.chrom:
            return GenomicRange()

        new_range = super().intersect(other)
        if self.strand == other.strand:
            strand = self.strand
        else:
            strand = "*"

        return GenomicRange(self.chrom, new_range.start, new_range.end, strand)

    def union(self, other):
        if self.chrom == other.chrom:
            union_result = super().union(other)
            if isinstance(union_result, Range):
                if self.strand == other.strand:
                    strand = self.strand
                else:
                    strand = "*"

                return GenomicRange(self.chrom,
                                    union_result.start,
                                    union_result.end,
                                    strand)
            else:
                return [self, other]
        else:
            return [self, other]

    def __repr__(self):
        out = "%s:%d-%d" % (self.chrom, self.start, self.end)
        if self.strand != "*":
            out += "(%s)" % self.strand
        return out

    def __lt__(self, other):
        return self.chrom < other.chrom or\
            (self.chrom == other.chrom and super().__lt__(other))

    def __le__(self, other):
        return self.chrom < other.chrom or\
            (self.chrom == other.chrom and super().__le__(other))


class GenomicRanges(Ranges):
    def __init__(self, ranges=[]):
        super().__init__(ranges)

    def overlaps(self, other):
        overlaps_indices = []
        for self_index, self_range in enumerate(self):
            for other_index, other_range in enumerate(other):
                if other_range.chrom < self_range.chrom or\
                        other_range.end < self_range.start:
                    continue
                elif self_range.chrom < other_range.chrom or\
                        self_range.end < other_range.start:
                    break

                if (self_range.start < other_range.start and
                        self_range.end >= other_range.start) or\
                        (other_range.start < self_range.start and
                         other_range.end >= self_range.start):
                    overlaps_indices.append((self_index, other_index))
        return overlaps_indices
