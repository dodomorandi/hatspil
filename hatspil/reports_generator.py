import datetime
import os
from itertools import accumulate
from typing import (Any, Dict, Iterable, List, Optional, Sequence, Set, TextIO,
                    Union, cast)

import colorlover as cl
import plotly.graph_objs as go
import plotly.offline as plt
import spectra
from bson import ObjectId

from . import utils
from .barcoded_filename import BarcodedFilename
from .config import Config
from .db import Db
from .db.picard_metrics import PicardMetricsType
from .report_table import ReportTable


class FigureData:
    figure: Union[go.Figure, ReportTable]
    date: str

    def __init__(self, figure: Union[go.Figure, ReportTable], date: str) -> None:
        self.figure = figure
        self.date = date


class ReportsGenerator:
    FiguresType = Dict[str, List[FigureData]]

    def __init__(
        self,
        root: str,
        config: Config,
        parameters: Dict[str, Any],
        barcoded_filenames: Optional[Iterable[BarcodedFilename]] = None,
    ) -> None:
        self.root = root
        self.config = config
        self.parameters = parameters
        self.barcoded_filenames = barcoded_filenames
        self.out_dir = os.path.join(self.root, "reports")
        self.is_plotlyjs_included = False
        self.db = Db(self.config)

    def generate_analysis_reports(self) -> None:
        if self.barcoded_filenames is None:
            return

        samples_per_directory: Dict[str, List[BarcodedFilename]] = {}
        for barcoded_filename in self.barcoded_filenames:
            report_directory = barcoded_filename.get_directory(self.out_dir)
            samples_per_directory.setdefault(report_directory, []).append(
                barcoded_filename
            )

        current_date = utils.get_overridable_current_date(self.parameters)
        for report_directory, barcoded_filenames in samples_per_directory.items():
            report_filename = os.path.join(
                report_directory, "report_{}.html".format(current_date)
            )
            if os.path.exists(report_filename):
                counter = 1
                while os.path.exists(report_filename):
                    os.path.join(
                        report_directory,
                        "report_{}_{:02}.html".format(current_date, counter),
                    )
                    counter += 1

            self.create_report_for_barcoded_filenames(
                barcoded_filenames, report_filename
            )

    def generate_global_reports(self) -> None:
        # TODO
        pass

    def create_report_for_barcoded_filenames(
        self, barcoded_filenames: Iterable[BarcodedFilename], report_filename: str
    ) -> None:
        self.is_plotlyjs_included = False
        db_data = self._get_db_data_from_barcoded_filenames(barcoded_filenames)

        figures_by_sequencing: ReportsGenerator.FiguresType = {}
        figures_by_sample: ReportsGenerator.FiguresType = {}
        figures_by_biopsy: ReportsGenerator.FiguresType = {}
        figures_by_patient: ReportsGenerator.FiguresType = {}
        figures_by_project: ReportsGenerator.FiguresType = {}

        for sample_data, sample_barcoded_filename in zip(db_data, barcoded_filenames):
            sample_barcode = sample_barcoded_filename.get_barcode()
            sequencing_figures = figures_by_sequencing.setdefault(sample_barcode, [])
            self._add_picard_data_to_figure(sample_data, sequencing_figures)

        self._add_comparative_sequencing_plots(db_data, figures_by_sequencing)
        self._add_comparative_sample_plots(db_data, figures_by_sample)
        self._add_comparative_biopsy_plots(db_data, figures_by_biopsy)
        self._add_comparative_patient_plots(db_data, figures_by_patient)
        self._add_comparative_project_plots(db_data, figures_by_project)
        self._add_variants_tables(db_data, figures_by_sequencing)

        ReportsGenerator._create_html_from_figures(
            report_filename,
            figures_by_sequencing,
            figures_by_sample,
            figures_by_biopsy,
            figures_by_patient,
            figures_by_project,
        )

    def _get_db_data_from_barcoded_filenames(
        self, barcoded_filenames: Iterable[BarcodedFilename]
    ) -> List[Dict[str, Any]]:
        db_data: List[Dict[str, Any]] = []

        for barcoded_filename in barcoded_filenames:
            from_db = self.db.from_barcoded(barcoded_filename)
            assert from_db
            db_data.append(from_db)

        return db_data

    def _add_picard_data_to_figure(
        self, sample_data: Dict[str, Any], figures: List[FigureData]
    ) -> None:
        assert "sequencing" in sample_data
        assert isinstance(sample_data["sequencing"], dict)
        sequencing_id = sample_data["sequencing"]["sequencing_id"]
        assert sequencing_id

        picard_hs_metrics = self.db.picard_metrics.find(
            {"sequencing": sequencing_id, "date": sample_data["date"], "type": "hs"}
        )
        if picard_hs_metrics is not None:
            histogram = picard_hs_metrics["histogram"]
            if "coverage" in histogram:
                histogram_params = (
                    ("baseq_count", "Base qualities", "quality", "n bases"),
                    ("count", "High quality coverages", "coverage", "n reads"),
                )
            else:
                histogram_params = (
                    ("unfiltered_baseq_count", "Base qualities", "quality", "n bases"),
                    (
                        "high_quality_coverage_count",
                        "High quality coverages",
                        "coverage",
                        "n reads",
                    ),
                )

            for coverage_param, title, x_axis_label, y_axis_label in histogram_params:
                coverage_count = histogram[coverage_param]

                before_zero = utils.rfind_if(
                    coverage_count, lambda x: cast(bool, x != 0)
                )
                assert before_zero
                coverage_count = coverage_count[:before_zero]
                max_n_reads = max(coverage_count[1:])
                max_n_reads_coverage = sum(coverage_count[1:]) * 0.99
                cum_n_reads = list(accumulate(coverage_count[1:]))
                max_coverage = None
                for index, reads in enumerate(cum_n_reads):
                    if reads > max_n_reads_coverage:
                        max_coverage = index + 3
                        break
                assert max_coverage

                figure = go.Figure(
                    data=[go.Bar(y=coverage_count)],
                    layout=go.Layout(
                        title=title,
                        xaxis={"title": x_axis_label, "range": [0, max_coverage]},
                        yaxis={"title": y_axis_label, "range": [0, max_n_reads]},
                    ),
                )

                figures.append(FigureData(figure, sample_data["date"]))

        picard_gcbias_metrics = self.db.picard_metrics.find(
            {"sequencing": sequencing_id, "date": sample_data["date"], "type": "gcbias"}
        )

        if picard_gcbias_metrics is not None:
            metrics_class = picard_gcbias_metrics["metrics class"]
            figure = go.Figure(
                data=[go.Bar(x=metrics_class["gc"], y=metrics_class["windows"])],
                layout=go.Layout(
                    title="GC percent windows",
                    xaxis={"title": "GC percentage"},
                    yaxis={"title": "n windows"},
                ),
            )
            figures.append(FigureData(figure, sample_data["date"]))

            figure = go.Figure(
                data=[go.Bar(x=metrics_class["gc"], y=metrics_class["read_starts"])],
                layout=go.Layout(
                    title="GC windows starting reads",
                    xaxis={"title": "GC percentage"},
                    yaxis={"title": "n reads"},
                ),
            )
            figures.append(FigureData(figure, sample_data["date"]))

            figure = go.Figure(
                data=[
                    go.Bar(
                        x=metrics_class["gc"],
                        y=metrics_class["normalized_coverage"],
                        error_y={
                            "type": "data",
                            "array": metrics_class["error_bar_width"],
                        },
                    )
                ],
                layout=go.Layout(
                    title="GC normalized coverage",
                    xaxis={"title": "GC percentage"},
                    yaxis={"title": "Normalized coverage"},
                ),
            )
            figures.append(FigureData(figure, sample_data["date"]))

    def _add_plots_for_multiple_samples(
        self, samples: Iterable[Dict[str, Any]], figures: List[FigureData]
    ) -> None:
        hs_metrics = [
            self.db.picard_metrics.find(
                {
                    "sequencing": sample["sequencing"]["_id"],
                    "type": PicardMetricsType.hs.name,
                    "date": sample["date"],
                }
            )
            for sample in samples
            if "sequencing" in sample
        ]
        assert all(hs_metrics)

        marked_dup_metrics = [
            self.db.picard_metrics.find(
                {
                    "sequencing": sample["sequencing"]["_id"],
                    "type": PicardMetricsType.marked_duplicates.name,
                    "date": sample["date"],
                }
            )
            for sample in samples
            if "sequencing" in sample
        ]
        assert all(marked_dup_metrics)

        barcoded_samples = [self.db.to_barcoded(sample) for sample in samples]

        raw_barcodes = [
            barcoded_sample.get_barcode()
            for barcoded_sample in barcoded_samples
            if barcoded_sample
        ]
        assert len(raw_barcodes) == len(barcoded_samples)

        high_quality_unique_aligned_filtered: List[Any] = []
        low_quality_unique_aligned_filtered: List[Any] = []
        for hs, high_quality in zip(hs_metrics, high_quality_unique_aligned_filtered):
            assert hs
            metrics_class = hs["metrics class"]
            assert metrics_class
            pf_uq_reads_aligned = metrics_class["pf_uq_reads_aligned"]
            pf_unique_reads = metrics_class["pf_unique_reads"]

            assert pf_uq_reads_aligned
            assert pf_unique_reads

            high_quality_unique_aligned_filtered.append(pf_uq_reads_aligned)
            low_quality_unique_aligned_filtered.append(pf_unique_reads - high_quality)

        paired_reads_duplicates: List[int] = []
        duplicated_unpaired_reads: List[int] = []
        unique_unpaired_reads: List[int] = []
        unaligned_reads: List[int] = []

        for marked_dup, duplicates in zip(
            marked_dup_metrics, duplicated_unpaired_reads
        ):
            assert marked_dup
            metrics_class = marked_dup["metrics class"]
            assert metrics_class
            read_pair_duplicates = metrics_class["read_pair_duplicates"]
            unpaired_read_duplicates = metrics_class["unpaired_read_duplicates"]
            unpaired_reads_examined = metrics_class["unpaired_reads_examined"]
            unmapped_reads = metrics_class["unmapped_reads"]

            assert read_pair_duplicates
            assert duplicated_unpaired_reads
            assert unpaired_reads_examined
            assert unmapped_reads

            paired_reads_duplicates.append(read_pair_duplicates * 2)
            duplicated_unpaired_reads.append(unpaired_read_duplicates)
            unique_unpaired_reads.append(unpaired_reads_examined - duplicates)
            unaligned_reads.append(unmapped_reads)

        color_scale = cl.scales["6"]["seq"]["Reds"]

        figure = go.Figure(
            data=[
                go.Bar(
                    x=raw_barcodes,
                    y=high_quality_unique_aligned_filtered,
                    name="Unique high score reads",
                    marker=dict(color=color_scale[5]),
                ),
                go.Bar(
                    x=raw_barcodes,
                    y=low_quality_unique_aligned_filtered,
                    name="Unique low score reads",
                    marker=dict(color=color_scale[4]),
                ),
                go.Bar(
                    x=raw_barcodes,
                    y=paired_reads_duplicates,
                    name="Duplicated aligned paired reads",
                    marker=dict(color=color_scale[3]),
                ),
                go.Bar(
                    x=raw_barcodes,
                    y=unique_unpaired_reads,
                    name="Unique unpaired reads",
                    marker=dict(color=color_scale[2]),
                ),
                go.Bar(
                    x=raw_barcodes,
                    y=duplicated_unpaired_reads,
                    name="Duplicated unpaired reads",
                    marker=dict(color=color_scale[1]),
                ),
                go.Bar(
                    x=raw_barcodes,
                    y=unaligned_reads,
                    name="Unaligned reads",
                    marker=dict(color=color_scale[0]),
                ),
            ],
            layout=go.Layout(
                barmode="stack",
                xaxis={"tickangle": -30, "automargin": True},
                margin={"l": 150},
                width=max(len(raw_barcodes) * 30 + 230, 510),
                legend={"xanchor": "right", "yanchor": "top"},
            ),
        )
        figures.append(
            FigureData(figure, utils.get_overridable_current_date(self.parameters))
        )

        multipliers = ("100x", "50x", "40x", "30x", "20x", "10x", "2x", "1x")
        coverage_colors = [
            color.hexcode
            for color in spectra.scale(["#003333", "#CCFFFF"]).range(len(multipliers))
        ]

        assert hs
        relative_coverages_per_sample: List[List[float]] = []
        for multiplier in multipliers:
            metrics_class = hs["metrics class"]
            assert metrics_class
            relative_coverage = metrics_class["pct_target_bases_{}".format(multiplier)]

            relative_coverages_per_sample.append(relative_coverage)

        relative_coverages_per_sample = [
            [sample_coverages[0]]
            + [
                int((sample_coverages[i] - sample_coverages[i + 1]) * 100)
                for i in range(len(sample_coverages) - 1)
            ]
            for sample_coverages in relative_coverages_per_sample
        ]

        coverage_bars = [
            go.Bar(
                x=[
                    sample_coverages[multiplier_index]
                    for sample_coverages in reversed(relative_coverages_per_sample)
                ],
                y=list(reversed(raw_barcodes)),
                name=multiplier,
                orientation="h",
                marker={"color": coverage_colors[multiplier_index]},
            )
            for (multiplier_index, multiplier) in enumerate(multipliers)
        ]
        figure = go.Figure(
            data=coverage_bars,
            layout=go.Layout(
                barmode="stack",
                yaxis={"tickangle": -30, "automargin": True, "ticksuffix": "  "},
                xaxis={"automargin": True},
                height=max(len(raw_barcodes) * 25 + 200, 360),
            ),
        )
        figures.append(
            FigureData(figure, utils.get_overridable_current_date(self.parameters))
        )

    def _add_comparative_sequencing_plots(
        self, samples: Iterable[Dict[str, Any]], sequencing_figures: FiguresType
    ) -> None:
        tissue_samples = []
        control_samples = []

        for sample in samples:
            try:
                barcoded_sample = Db.to_barcoded(sample)
            except Exception:
                continue
            assert barcoded_sample

            sequencing = sample["sequencing"]
            sequencing_id = sequencing["_id"]
            assert sequencing_id

            if not self._has_metrics_data(sequencing_id):
                continue

            if barcoded_sample.tissue.is_normal():
                control_samples.append(sample)
            else:
                tissue_samples.append(sample)

        for sample in tissue_samples:
            best_control = self._get_best_control(sample, control_samples)
            if not best_control:
                continue

            barcoded_sample = Db.to_barcoded(sample)
            assert barcoded_sample
            self._add_plots_for_multiple_samples(
                [sample, best_control],
                sequencing_figures.setdefault(barcoded_sample.get_barcode(), []),
            )

    def _add_comparative_sample_plots(
        self, samples: Iterable[Dict[str, Any]], sample_figures: FiguresType
    ) -> None:
        for current_sample in samples:
            current_samples = self._get_grouped_samples_with_controls(
                [
                    sample
                    for sample in samples
                    if sample["sample"]["_id"] == current_sample["sample"]["_id"]
                    and self._has_metrics_data(sample["sequencing"]["_id"])
                ]
            )

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                barcoded_sample.sequencing = None
                barcoded_sample.kit = None
                barcoded_sample.molecule = None
                barcoded_sample.analyte = None

                self._add_plots_for_multiple_samples(
                    current_samples,
                    sample_figures.setdefault(barcoded_sample.get_barcode(), []),
                )

    def _add_comparative_biopsy_plots(
        self, samples: Iterable[Dict[str, Any]], biopsy_figures: FiguresType
    ) -> None:
        for current_sample in samples:
            current_samples = self._get_grouped_samples_with_controls(
                [
                    sample
                    for sample in samples
                    if sample["biopsy"]["_id"] == current_sample["biopsy"]["_id"]
                    and sample["sample"]["_id"] == current_sample["sample"]["_id"]
                    and self._has_metrics_data(sample["sequencing"]["_id"])
                ]
            )

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                barcoded_sample.sequencing = None
                barcoded_sample.kit = None
                barcoded_sample.molecule = None
                barcoded_sample.analyte = None
                barcoded_sample.sample = None
                barcoded_sample.xenograft = None

                self._add_plots_for_multiple_samples(
                    current_samples,
                    biopsy_figures.setdefault(barcoded_sample.get_barcode(), []),
                )

    def _add_comparative_patient_plots(
        self, samples: Iterable[Dict[str, Any]], patient_figures: FiguresType
    ) -> None:
        for current_sample in samples:
            current_samples = self._get_grouped_samples_with_controls(
                [
                    sample
                    for sample in samples
                    if sample["patient"]["_id"] == current_sample["patient"]["_id"]
                    and sample["biopsy"]["_id"] == current_sample["biopsy"]["_id"]
                    and sample["sample"]["_id"] == current_sample["sample"]["_id"]
                    and self._has_metrics_data(sample["sequencing"]["_id"])
                ]
            )

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                barcode = "{}-{}".format(
                    barcoded_sample.project, barcoded_sample.patient
                )

                self._add_plots_for_multiple_samples(
                    current_samples, patient_figures.setdefault(barcode, [])
                )

    def _add_comparative_project_plots(
        self, samples: Iterable[Dict[str, Any]], project_figures: FiguresType
    ) -> None:
        for current_sample in samples:
            current_samples = self._get_grouped_samples_with_controls(
                [
                    sample
                    for sample in samples
                    if sample["project"]["_id"] == current_sample["project"]["_id"]
                    and sample["patient"]["_id"] == current_sample["patient"]["_id"]
                    and sample["biopsy"]["_id"] == current_sample["biopsy"]["_id"]
                    and sample["sample"]["_id"] == current_sample["sample"]["_id"]
                    and self._has_metrics_data(sample["sequencing"]["_id"])
                ]
            )

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                self._add_plots_for_multiple_samples(
                    current_samples,
                    project_figures.setdefault(barcoded_sample.get_barcode(), []),
                )

    def _has_metrics_data(self, sequencing_id: ObjectId) -> bool:
        return (
            self.db.picard_metrics.find(
                {"sequencing": sequencing_id, "type": PicardMetricsType.hs.name}
            )
            is not None
            and self.db.picard_metrics.find(
                {
                    "sequencing": sequencing_id,
                    "type": PicardMetricsType.marked_duplicates.name,
                }
            )
            is not None
        )

    def _get_best_control(
        self, sample: Dict[str, Any], control_samples: Sequence[Dict[str, Any]]
    ) -> Dict[str, Any]:
        assert sample
        assert control_samples

        if len(control_samples) == 1:
            return control_samples[0]

        matching = self._get_best_matching_samples(sample, control_samples)
        sequencing = sample["sequencing"]

        if len(matching) > 1:
            sample_date = ReportsGenerator._get_sequencing_date(sequencing)
            if sample_date is None:
                sample_date = datetime.date.today()

            normal_dates: List[datetime.date] = []
            for normal_data in control_samples:
                normal_sequencing = normal_data["sequencing"]
                if "data" in normal_sequencing:
                    splitted_date = [
                        int(value) for value in normal_sequencing["date"].split("-")
                    ]
                    if len(splitted_date) != 3:
                        normal_dates.append(datetime.date.today())
                    else:
                        normal_dates.append(
                            datetime.date(
                                splitted_date[0], splitted_date[1], splitted_date[2]
                            )
                        )
                else:
                    normal_dates.append(datetime.date.today())

            distances_from_samples = [
                abs(sample_date - normal_date) for normal_date in normal_dates
            ]
            matching_index = utils.argmin(distances_from_samples)
            assert matching_index
            return matching[matching_index]
        else:
            assert matching
            return matching[0]

    def _get_best_control_with_metrics_data(
        self, sample: Dict[str, Any], control_samples: Sequence[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        if not self._has_metrics_data(sample["sequencing"]["_id"]):
            return None

        control_samples_with_metrics = [
            control_sample
            for control_sample in control_samples
            if self._has_metrics_data(sample["sequencing"]["_id"])
        ]
        if not control_samples_with_metrics:
            return None
            return None

        return self._get_best_control(sample, control_samples_with_metrics)

    def _get_best_matching_samples(
        self, sample: Dict[str, Any], samples: Sequence[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        sequencing = sample["sequencing"]
        matching: List[Dict[str, Any]] = [
            other_sample
            for other_sample in samples
            if sequencing["sample"] == other_sample["sequencing"]["sample"]
            and sequencing["molecule"] == other_sample["sequencing"]["molecule"]
            and sequencing["analyte"] == other_sample["sequencing"]["analyte"]
        ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sequencing["sample"] == other_sample["sequencing"]["sample"]
                and sequencing["analyte"] == other_sample["sequencing"]["analyte"]
            ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sequencing["sample"] == other_sample["sequencing"]["sample"]
            ]

        if not matching:
            sample_sample = sample["sample"]
            matching = [
                other_sample
                for other_sample in samples
                if sample_sample["biopsy"] == other_sample["sample"]["biopsy"]
                and sequencing["molecule"] == other_sample["sequencing"]["molecule"]
                and sequencing["analyte"] == other_sample["sequencing"]["analyte"]
            ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sample_sample["biopsy"] == other_sample["sample"]["biopsy"]
                and sequencing["analyte"] == other_sample["sequencing"]["analyte"]
            ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sample_sample["biopsy"] == other_sample["sample"]["biopsy"]
            ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sequencing["molecule"] == other_sample["sequencing"]["molecule"]
                and sequencing["analyte"] == other_sample["sequencing"]["analyte"]
            ]

        if not matching:
            matching = [
                other_sample
                for other_sample in samples
                if sequencing["analyte"] == other_sample["sequencing"]["analyte"]
            ]

        if not matching:
            matching = list(samples)

        return matching

    def _get_grouped_samples_with_controls(
        self, _samples: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        primary_samples = []
        control_samples = []

        for _sample in _samples:
            try:
                barcoded = Db.to_barcoded(_sample)
                if not barcoded:
                    continue

                if barcoded.tissue.is_normal():
                    control_samples.append(_sample)
                else:
                    primary_samples.append((_sample, barcoded))
            except Exception:
                continue

        if not primary_samples:
            return []

        primary_samples.sort(key=lambda sample: sample[1])

        sample_best_controls = [
            self._get_best_control_with_metrics_data(sample[0], control_samples)
            for sample in primary_samples
        ]

        control_samples = [
            control_sample
            for control_sample in control_samples
            if control_sample not in sample_best_controls
        ]
        for control_sample in sample_best_controls:
            if (
                control_sample is not None
                and sample_best_controls.count(control_sample) > 1
            ):
                control_samples.append(control_sample)

        control_samples_for_primaries = [
            control_sample if control_sample not in control_samples else None
            for control_sample in sample_best_controls
        ]

        out_samples = []
        for primary_sample, control_sample in zip(
            primary_samples, control_samples_for_primaries
        ):
            out_samples.append(primary_sample[0])
            if control_sample:
                out_samples.append(control_sample)
        if len(out_samples) <= 2 and not out_samples:
            return []

        def get_sort_control_samples_key(control_sample: Dict[str, Any]) -> str:
            barcoded = Db.to_barcoded(control_sample)
            assert barcoded
            return barcoded.get_barcode()

        control_samples.sort(key=get_sort_control_samples_key)

        out_samples += control_samples
        return out_samples

    def _get_nearest_barcoded_sample(
        self, sample: Dict[str, Any], samples: Iterable[Dict[str, Any]]
    ) -> Optional[BarcodedFilename]:
        sample_barcode = Db.to_barcoded(sample)
        assert sample_barcode

        sample_date = ReportsGenerator._get_sequencing_date(sample)
        if sample_date:
            dated_samples = [
                (sample, ReportsGenerator._get_sequencing_date(sample))
                for sample in samples
            ]
            dated_samples_valid_count = sum(
                1 if dated_sample[1] is not None else 0
                for dated_sample in dated_samples
            )
            if dated_samples_valid_count:
                date_distances: List[Optional[datetime.timedelta]] = []
                for dated_sample in dated_samples:
                    if dated_sample and dated_sample[1]:
                        other_sample_date = dated_sample[1]
                    else:
                        other_sample_date = None

                    if other_sample_date:
                        date_distances.append(abs(other_sample_date - sample_date))
                    else:
                        date_distances.append(None)

                best_date_distance = min(date_distances)
                matchable_samples = [
                    sample
                    for (sample, date_distance) in zip(samples, date_distances)
                    if date_distance == best_date_distance
                ]
            else:
                matchable_samples = list(samples)

            best_sample = self._get_best_matching_samples(sample, matchable_samples)[0]
            return Db.to_barcoded(best_sample)
        else:
            best_sample = self._get_best_matching_samples(sample, list(samples))[0]
            return Db.to_barcoded(best_sample)

    @staticmethod
    def _get_sequencing_date(sequencing: Dict[str, Any]) -> Optional[datetime.date]:
        if "date" in sequencing:
            splitted_date = [int(value) for value in sequencing["date"].split("-")]
            return datetime.date(splitted_date[0], splitted_date[1], splitted_date[2])
            return datetime.date(splitted_date[0], splitted_date[1], splitted_date[2])
        else:
            return None

    def _add_variants_tables(
        self, samples: Iterable[Dict[str, Any]], figures: FiguresType
    ) -> None:
        for sample in samples:
            analyses = self.db.analyses.find_all(
                {"sequencing": sample["sequencing"]["_id"]}
            )
            if not analyses:
                continue

            annotations_list = []
            for analysis in analyses:
                assert analysis
                analysis_annotation = self.db.annotations.find(
                    {"_id": analysis["annotations"]}
                )
                annotations_list.append(analysis_annotation)

            annotations = {}
            for annotation in annotations_list:
                assert annotation
                annotation_id = annotation["id"]
                assert annotation_id
                annotations[annotation_id] = annotation

            variants = [
                variant
                for analysis in analyses
                for variant in analysis["variants"]
                if "id" in variant
                if annotations[variant["id"]].get("damaging") == "High"
            ]

            variants_annotations = [annotations[variant["id"]] for variant in variants]
            cosmic_ids = []
            for annotation in variants_annotations:
                cosmic_id = annotation.get("cosmic70")
                if cosmic_id:
                    splitted_pairs = cosmic_id.split(";")
                    cosmic_id = None

                    for raw_pairs in splitted_pairs:
                        splitted = raw_pairs.split("=")
                        if len(splitted) == 2 and splitted[0].lower() == "id":
                            cosmic_id = splitted[1]
                            break
                cosmic_ids.append(cosmic_id)

            barcoded_sample = self.db.to_barcoded(sample)
            assert barcoded_sample

            sample_barcode = barcoded_sample.get_barcode()
            assert sample_barcode

            report_table = ReportTable(
                sample_barcode,
                "Chr",
                "Start",
                "Ref",
                "Alt",
                "Gene",
                "Func",
                "snp138/cosmic70",
                "CADD13 PHRED",
                "identifiers",
                "druggable",
                "DP",
            )

            for (variant, annotation, cosmic_id) in zip(
                variants, variants_annotations, cosmic_ids
            ):
                report_table.add_row(
                    [
                        variant["Chr"],
                        variant["Start"],
                        variant["Ref"],
                        variant["Alt"],
                        annotation["Gene refGene"],
                        annotation["Func refGene"],
                        " ".join(
                            value
                            for value in [annotation.get("snp138"), cosmic_id]
                            if value
                        ),
                        annotation["CADD13_PHRED"]
                        if "CADD13_PHRED" in annotation
                        else "",
                        " ".join(annotation["hgnc_canonical_refseq"].split(":"))
                        if "hgnc_canonical_refseq" in annotation
                        else "",
                        annotation["druggable"],
                        variant["DP"],
                    ]
                )

            report_table.style = (
                "g.legend rect.bg { opacity: 0.7; }"
                "table.dataTable.display tbody tr.druggable-gene.even,"
                "table.dataTable.display tbody tr.druggable-gene.even td"
                "{ background-color: #83ff83; }"
                "table.dataTable.display tbody tr.druggable-gene.odd,"
                "table.dataTable.display tbody tr.druggable-gene.odd td"
                "{ background-color: #80e780; }"
            )
            sample_date = sample["date"]
            assert sample_date
            figures.setdefault(sample_barcode, []).append(
                FigureData(report_table, sample_date)
            )

    @staticmethod
    def _create_html_from_figures(
        output_filename: str,
        figures_by_sequencing: FiguresType,
        figures_by_sample: FiguresType,
        figures_by_biopsy: FiguresType,
        figures_by_patient: FiguresType,
        figures_by_project: FiguresType,
    ) -> None:

        styles: Set[str] = set()
        for figures_dict in (
            figures_by_sequencing,
            figures_by_sample,
            figures_by_biopsy,
            figures_by_patient,
            figures_by_project,
        ):
            for figures_data in figures_dict.values():
                for figure_data in figures_data:
                    figure = figure_data.figure
                    if not isinstance(figure, ReportTable) or not figure.style:
                        continue

                    styles.add(figure.style)

        with open(output_filename, "w") as out_fd:
            out_fd.write(
                "<!doctype html>"
                "<html>"
                "<head>"
                "<meta charset='utf-8'/>"
                "<link rel='stylesheet' type='text/css' "
                "href='https://cdn.datatables.net/1.10.19/css/jquery.dataTables.css'/>"
                "<link rel='stylesheet' type='text/css' "
                "href='https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css'/>"
                "<script src='https://code.jquery.com/jquery-3.3.1.min.js'>"
                "</script>"
                "<script src='https://code.jquery.com/ui/1.12.1/jquery-ui.min.js'>"
                "</script>"
                "<script "
                "src='https://cdn.datatables.net/1.10.19/js/jquery.dataTables.js'>"
                "</script>"
            )
            if styles:
                out_fd.write("<style>{}</style>".format("\n".join(styles)))
            out_fd.write("</head>" "<div id='top-tabs'>" "<ul>")

            if figures_by_sequencing:
                out_fd.write("<li><a href='#tab-by-sequencing'>By sequencing</a>")

            if figures_by_sample:
                out_fd.write("<li><a href='#tab-by-sample'>By sample</a>")

            if figures_by_biopsy:
                out_fd.write("<li><a href='#tab-by-biopsy'>By biopsy</a>")

            if figures_by_patient:
                out_fd.write("<li><a href='#tab-by-patient'>By patient</a>")

            if figures_by_project:
                out_fd.write("<li><a href='#tab-by-project'>By project</a>")

            out_fd.write("</ul>")

            include_plotlyjs = ReportsGenerator.write_figures_on_fd(
                out_fd, "sequencing", figures_by_sequencing, True
            )
            include_plotlyjs = ReportsGenerator.write_figures_on_fd(
                out_fd, "sample", figures_by_sample, include_plotlyjs
            )
            include_plotlyjs = ReportsGenerator.write_figures_on_fd(
                out_fd, "biopsy", figures_by_biopsy, include_plotlyjs
            )
            include_plotlyjs = ReportsGenerator.write_figures_on_fd(
                out_fd, "patient", figures_by_patient, include_plotlyjs
            )
            ReportsGenerator.write_figures_on_fd(
                out_fd, "project", figures_by_project, include_plotlyjs
            )

            out_fd.write(
                "</div>"
                "<script>"
                "$( function() {"
                "$( '#top-tabs' ).tabs();"
                "} );"
                "</script>"
                "</html>"
            )

    @staticmethod
    def write_figures_on_fd(
        fd: TextIO, figure_type: str, figures_data: FiguresType, include_plotlyjs: bool
    ) -> bool:
        if figures_data:
            fd.write(
                "<div id='tab-by-{}'>"
                "<div id='by_{}_tabs'>"
                "<ul>"
                "{}"
                "</ul>".format(
                    figure_type,
                    figure_type,
                    "".join(
                        "<li><a href='#{}_by_{}'>{}</li>".format(
                            sample_name, figure_type, sample_name
                        )
                        for sample_name in figures_data.keys()
                    ),
                )
            )

            def figure_get_date_key(figure_data: FigureData) -> datetime.date:
                if figure_data.date:
                    return datetime.datetime.strptime(figure_data.date, "%Y_%m_%d")
                else:
                    return datetime.date(datetime.MINYEAR, 1, 1)

            for sample_name, figures in figures_data.items():
                figures.sort(key=figure_get_date_key)

                last_date: Optional[str] = None
                for figure in figures:
                    if figure.date:
                        figure_date = figure.date
                    else:
                        figure_date = "Unknown date"

                    if not last_date or last_date != figure_date:
                        fd.write("<h3>{}</h3>".format(figure_date))
                        last_date = figure_date
                    fd.write("<div id='{}_by_{}'>".format(sample_name, figure_type))

                    if isinstance(figure, go.Figure):
                        fd.write(
                            plt.plot(
                                figure,
                                output_type="div",
                                include_plotlyjs=include_plotlyjs,
                            )
                        )
                        include_plotlyjs = False
                    else:
                        assert isinstance(figure, ReportTable)
                        assert figure_type == "sequencing"
                        fd.write(figure.html())

                    fd.write("</div>")

            fd.write(
                "</div>"
                "</div>"
                "<script>"
                "$(function() {{ $('#by_{}_tabs').tabs(); }} )"
                "</script>".format(figure_type)
            )
        return include_plotlyjs
