import datetime
import os
import re
from itertools import accumulate
from typing import (Any, Dict, Iterable, List, Optional, Sequence, Set, TextIO,
                    Tuple, Union, cast)

import colorlover as cl
import plotly.graph_objs as go
import plotly.offline as plt
import spectra
from bson import ObjectId

from . import utils
from .barcoded_filename import BarcodedFilename, Tissue
from .config import Config
from .db import Db
from .db.picard_metrics import PicardMetricsType
from .report_table import (ReportTable, ReportTableColumnOrdering,
                           ReportTableSortingType)


class FigureData:
    figure: Union[go.Figure, ReportTable]
    date: str

    def __init__(self, figure: Union[go.Figure, ReportTable], date: str) -> None:
        self.figure = figure
        self.date = date


class ReportsGenerator:
    FiguresType = Dict[str, List[FigureData]]
    RE_INVALID_CHAR = re.compile(r"[^A-Za-z0-9._-]")
    cached_annotations: Dict[str, Dict[str, Any]] = {}
    cached_annotations_ids: Dict[ObjectId, str] = {}
    cached_sequencing_analysis_dates: Dict[ObjectId, List[datetime.date]] = {}
    cached_sequencings_to_analyses: Dict[ObjectId, List[Tuple[ObjectId, str]]] = {}
    cached_sequencing_has_metrics_data: Dict[ObjectId, bool] = {}

    def __init__(
        self,
        root: str,
        config: Config,
        parameters: Dict[str, Any],
        barcoded_samples: Optional[Iterable[BarcodedFilename]] = None,
    ) -> None:
        self.root = root
        self.config = config
        self.parameters = parameters
        self.barcoded_samples = barcoded_samples
        self.out_dir = os.path.join(self.root, "reports")
        self.is_plotlyjs_included = False
        self.db = Db(self.config)

    def generate_analysis_reports(self) -> None:
        if self.barcoded_samples is None:
            return

        self._generate_reports("report", self.barcoded_samples)

    def generate_global_reports(self) -> None:
        barcoded_samples: List[BarcodedFilename] = []
        db_data: List[Dict[str, Any]] = []

        sequencings = self.db.sequencings.find_all()
        assert sequencings is not None
        for sequencing in sequencings:
            sample_data = self.db.from_sequencing_id(sequencing["_id"])
            assert sample_data
            try:
                barcoded_sample = Db.to_barcoded(sample_data)
            except Exception:
                continue

            if barcoded_sample:
                db_data.append(sample_data)
                barcoded_samples.append(barcoded_sample)

        self._generate_reports("global_report", barcoded_samples, db_data)

    def _generate_reports(
        self,
        filename_prefix: str,
        barcoded_samples: Iterable[BarcodedFilename],
        db_data: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        samples_per_directory: Dict[str, List[BarcodedFilename]] = {}
        for barcoded_sample in barcoded_samples:
            report_directory = barcoded_sample.get_directory(self.out_dir)
            samples_per_directory.setdefault(report_directory, []).append(
                barcoded_sample
            )

        current_date = utils.get_overridable_current_date(self.parameters)
        for report_directory, barcoded_samples in samples_per_directory.items():
            os.makedirs(report_directory, exist_ok=True)
            report_filename = os.path.join(
                report_directory, "{}_{}.html".format(filename_prefix, current_date)
            )
            if os.path.exists(report_filename):
                counter = 1
                while os.path.exists(report_filename):
                    report_filename = os.path.join(
                        report_directory,
                        "{}_{}_{:02}.html".format(
                            filename_prefix, current_date, counter
                        ),
                    )
                    counter += 1

            self._create_report_for_barcoded_samples(
                barcoded_samples, report_filename, db_data
            )

    def _create_report_for_barcoded_samples(
        self,
        barcoded_samples: Iterable[BarcodedFilename],
        report_filename: str,
        db_data: Optional[List[Dict[str, Any]]],
    ) -> None:
        self.is_plotlyjs_included = False
        if not db_data:
            db_data = self._get_db_data_from_barcoded_samples(barcoded_samples)

        figures_by_sequencing: ReportsGenerator.FiguresType = {}
        figures_by_sample: ReportsGenerator.FiguresType = {}
        figures_by_biopsy: ReportsGenerator.FiguresType = {}
        figures_by_patient: ReportsGenerator.FiguresType = {}
        figures_by_project: ReportsGenerator.FiguresType = {}

        for sample_data, sample_barcoded_sample in zip(db_data, barcoded_samples):
            sample_barcode = sample_barcoded_sample.get_barcode()
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

    def _get_db_data_from_barcoded_samples(
        self, barcoded_samples: Iterable[BarcodedFilename]
    ) -> List[Dict[str, Any]]:
        db_data: List[Dict[str, Any]] = []

        for barcoded_sample in barcoded_samples:
            from_db = self.db.from_barcoded(barcoded_sample)
            assert from_db
            db_data.append(from_db)

        return db_data

    def _add_picard_data_to_figure(
        self, sample_data: Dict[str, Any], figures: List[FigureData]
    ) -> None:
        assert "sequencing" in sample_data
        sequencing = sample_data["sequencing"]
        assert isinstance(sequencing, dict)
        assert "_id" in sequencing
        sequencing_id = sequencing["_id"]

        analysis_date_obj = self._get_best_analysis_date_for_sample(sample_data)
        if not analysis_date_obj:
            return
        analysis_date = analysis_date_obj.strftime("%Y_%m_%d")

        picard_hs_metrics = self.db.picard_metrics.find(
            {"sequencing": sequencing_id, "date": analysis_date, "type": "hs"}
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
                assert before_zero is not None
                coverage_count = coverage_count[:before_zero]
                if not coverage_count:
                    continue

                max_n_reads = max(coverage_count[1:])
                max_n_reads_coverage = sum(coverage_count[1:]) * 0.99
                cum_n_reads = list(accumulate(coverage_count[1:]))
                max_coverage = None
                for index, reads in enumerate(cum_n_reads):
                    if reads > max_n_reads_coverage:
                        max_coverage = index + 3
                        break
                assert max_coverage is not None

                figure = go.Figure(
                    data=[go.Bar(y=coverage_count)],
                    layout=go.Layout(
                        title=title,
                        xaxis={"title": x_axis_label, "range": [0, max_coverage]},
                        yaxis={"title": y_axis_label, "range": [0, max_n_reads]},
                    ),
                )

                figures.append(FigureData(figure, analysis_date))

        picard_gcbias_metrics = self.db.picard_metrics.find(
            {"sequencing": sequencing_id, "date": analysis_date, "type": "gcbias"}
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
            figures.append(FigureData(figure, analysis_date))

            figure = go.Figure(
                data=[go.Bar(x=metrics_class["gc"], y=metrics_class["read_starts"])],
                layout=go.Layout(
                    title="GC windows starting reads",
                    xaxis={"title": "GC percentage"},
                    yaxis={"title": "n reads"},
                ),
            )
            figures.append(FigureData(figure, analysis_date))

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
            figures.append(FigureData(figure, analysis_date))

    def _add_plots_for_multiple_samples(
        self, samples: Iterable[Dict[str, Any]], figures: List[FigureData]
    ) -> None:
        analyses_dates: List[Optional[str]] = []
        for sample in samples:
            if "sequencing" in sample:
                analysis_date = self._get_best_analysis_date_for_sample(sample)
                if not analysis_date:
                    analyses_dates.append(None)
                else:
                    analyses_dates.append(analysis_date.strftime("%Y_%m_%d"))
            else:
                analyses_dates.append(None)

        hs_metrics: List[Dict[str, Any]] = []
        marked_dup_metrics: List[Dict[str, Any]] = []
        barcoded_samples = []
        for sample, sample_analysis_date in zip(samples, analyses_dates):
            if "sequencing" in sample and sample_analysis_date:
                hs = self.db.picard_metrics.find(
                    {
                        "sequencing": sample["sequencing"]["_id"],
                        "type": PicardMetricsType.hs.name,
                        "date": sample_analysis_date,
                    }
                )

                if hs:
                    hs_metrics.append(hs)

                mark_dup = self.db.picard_metrics.find(
                    {
                        "sequencing": sample["sequencing"]["_id"],
                        "type": PicardMetricsType.marked_duplicates.name,
                        "date": sample_analysis_date,
                    }
                )
                if mark_dup:
                    marked_dup_metrics.append(mark_dup)

                barcoded_sample = self.db.to_barcoded(sample)
                assert barcoded_sample
                barcoded_samples.append(barcoded_sample)

        raw_barcodes = [
            barcoded_sample.get_barcode()
            for barcoded_sample in barcoded_samples
            if barcoded_sample
        ]
        assert len(raw_barcodes) == len(barcoded_samples)

        high_quality_unique_aligned_filtered: List[Any] = []
        low_quality_unique_aligned_filtered: List[Any] = []
        for hs in hs_metrics:
            assert hs
            metrics_class = hs["metrics class"]
            assert metrics_class
            pf_uq_reads_aligned = metrics_class["pf_uq_reads_aligned"]
            pf_unique_reads = metrics_class["pf_unique_reads"]

            high_quality_unique_aligned_filtered.append(pf_uq_reads_aligned)
            low_quality_unique_aligned_filtered.append(
                pf_unique_reads - pf_uq_reads_aligned
            )

        paired_reads_duplicates: List[int] = []
        duplicated_unpaired_reads: List[int] = []
        unique_unpaired_reads: List[int] = []
        unaligned_reads: List[int] = []

        for marked_dup in marked_dup_metrics:
            assert marked_dup
            metrics_class = marked_dup["metrics class"]
            assert metrics_class
            read_pair_duplicates = metrics_class["read_pair_duplicates"]
            unpaired_read_duplicates = metrics_class["unpaired_read_duplicates"]
            unpaired_reads_examined = metrics_class["unpaired_reads_examined"]
            unmapped_reads = metrics_class["unmapped_reads"]

            assert read_pair_duplicates is not None
            assert unpaired_read_duplicates is not None
            assert unpaired_reads_examined is not None
            assert unmapped_reads is not None

            paired_reads_duplicates.append(read_pair_duplicates * 2)
            duplicated_unpaired_reads.append(unpaired_read_duplicates)
            unique_unpaired_reads.append(
                unpaired_reads_examined - unpaired_read_duplicates
            )
            unaligned_reads.append(unmapped_reads)

        color_scale = cl.scales["6"]["seq"]["Reds"]

        figure_data: List[go.Bar] = []
        if high_quality_unique_aligned_filtered:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=high_quality_unique_aligned_filtered,
                    name="Unique high score reads",
                    marker=dict(color=color_scale[5]),
                )
            )

        if low_quality_unique_aligned_filtered:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=low_quality_unique_aligned_filtered,
                    name="Unique low score reads",
                    marker=dict(color=color_scale[4]),
                )
            )

        if paired_reads_duplicates:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=paired_reads_duplicates,
                    name="Duplicated aligned paired reads",
                    marker=dict(color=color_scale[3]),
                )
            )

        if unique_unpaired_reads:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=unique_unpaired_reads,
                    name="Unique unpaired reads",
                    marker=dict(color=color_scale[2]),
                )
            )

        if duplicated_unpaired_reads:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=duplicated_unpaired_reads,
                    name="Duplicated unpaired reads",
                    marker=dict(color=color_scale[1]),
                )
            )

        if unaligned_reads:
            figure_data.append(
                go.Bar(
                    x=raw_barcodes,
                    y=unaligned_reads,
                    name="Unaligned reads",
                    marker=dict(color=color_scale[0]),
                )
            )

        if figure_data:
            color_scale = cl.scales["6"]["seq"]["Reds"]

            figure = go.Figure(
                data=figure_data,
                layout=go.Layout(
                    title="Reads statistics",
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

        fractional_coverages_per_sample: List[List[float]] = []
        for hs in hs_metrics:
            assert hs
            sample_coverages: List[float] = []
            for multiplier in multipliers:
                metrics_class = hs["metrics class"]
                assert metrics_class
                relative_coverage: float = metrics_class[
                    "pct_target_bases_{}".format(multiplier)
                ]
                sample_coverages.append(relative_coverage)

            fractional_coverages_per_sample.append(sample_coverages)

        relative_coverages_per_sample = [
            [sample_coverages[0]]
            + [
                int((sample_coverages[i + 1] - sample_coverages[i]) * 100)
                for i in range(len(sample_coverages) - 1)
            ]
            for sample_coverages in fractional_coverages_per_sample
        ]

        if relative_coverages_per_sample:
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
                    title="Samples coverages",
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
            assert sequencing_id is not None

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
        added_samples: Set[str] = set()

        for current_sample in samples:
            if current_sample["sample"]["_id"] in added_samples:
                continue

            current_sample_barcoded = Db.to_barcoded(current_sample)
            assert current_sample_barcoded

            controls_to_check = []
            for sample in samples:
                if not self._has_metrics_data(sample["sequencing"]["_id"]):
                    continue

                sample_barcoded = Db.to_barcoded(sample)
                assert sample_barcoded

                if (
                    current_sample_barcoded.project == sample_barcoded.project
                    and current_sample_barcoded.patient == sample_barcoded.patient
                    and current_sample_barcoded.biopsy == sample_barcoded.biopsy
                    and current_sample_barcoded.sample == sample_barcoded.sample
                    and current_sample_barcoded.xenograft == sample_barcoded.xenograft
                ):
                    controls_to_check.append(sample)

            current_samples = self._get_grouped_samples_with_controls(controls_to_check)

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                added_samples |= set(
                    current_sample["sample"]["_id"]
                    for current_sample in current_samples
                )

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
        added_biopsies: Set[str] = set()

        for current_sample in samples:
            if current_sample["biopsy"]["_id"] in added_biopsies:
                continue

            current_sample_barcoded = Db.to_barcoded(current_sample)
            assert current_sample_barcoded

            controls_to_check = []
            for sample in samples:
                if not self._has_metrics_data(sample["sequencing"]["_id"]):
                    continue

                sample_barcoded = Db.to_barcoded(sample)
                assert sample_barcoded

                if (
                    current_sample_barcoded.project == sample_barcoded.project
                    and current_sample_barcoded.patient == sample_barcoded.patient
                    and current_sample_barcoded.biopsy == sample_barcoded.biopsy
                ):
                    controls_to_check.append(sample)

            current_samples = self._get_grouped_samples_with_controls(controls_to_check)

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                added_biopsies |= set(
                    current_sample["biopsy"]["_id"]
                    for current_sample in current_samples
                )

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
        added_patients: Set[str] = set()

        for current_sample in samples:
            if current_sample["patient"]["_id"] in added_patients:
                continue

            current_sample_barcoded = Db.to_barcoded(current_sample)
            assert current_sample_barcoded

            controls_to_check = []
            for sample in samples:
                if not self._has_metrics_data(sample["sequencing"]["_id"]):
                    continue

                sample_barcoded = Db.to_barcoded(sample)
                assert sample_barcoded

                if (
                    current_sample_barcoded.project == sample_barcoded.project
                    and current_sample_barcoded.patient == sample_barcoded.patient
                ):
                    controls_to_check.append(sample)

            current_samples = self._get_grouped_samples_with_controls(controls_to_check)

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                added_patients |= set(
                    current_sample["patient"]["_id"]
                    for current_sample in current_samples
                )
                barcode = "{}-{}".format(
                    barcoded_sample.project, barcoded_sample.patient
                )

                self._add_plots_for_multiple_samples(
                    current_samples, patient_figures.setdefault(barcode, [])
                )

    def _add_comparative_project_plots(
        self, samples: Iterable[Dict[str, Any]], project_figures: FiguresType
    ) -> None:
        added_projects: Set[str] = set()

        for current_sample in samples:
            if current_sample["project"]["_id"] in added_projects:
                continue

            current_sample_barcoded = Db.to_barcoded(current_sample)
            assert current_sample_barcoded

            controls_to_check = []
            for sample in samples:
                if not self._has_metrics_data(sample["sequencing"]["_id"]):
                    continue

                sample_barcoded = Db.to_barcoded(sample)
                assert sample_barcoded

                if current_sample_barcoded.project == sample_barcoded.project:
                    controls_to_check.append(sample)

            current_samples = self._get_grouped_samples_with_controls(controls_to_check)

            if current_samples:
                barcoded_sample = self._get_nearest_barcoded_sample(
                    current_sample, current_samples
                )
                if not barcoded_sample:
                    continue

                added_projects |= set(
                    current_sample["project"]["_id"]
                    for current_sample in current_samples
                )
                self._add_plots_for_multiple_samples(
                    current_samples,
                    project_figures.setdefault(barcoded_sample.project, []),
                )

    def _has_metrics_data(self, sequencing_id: ObjectId) -> bool:
        cached_result = ReportsGenerator.cached_sequencing_has_metrics_data.get(
            sequencing_id
        )
        if cached_result is not None:
            return cached_result

        result = (
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
        ReportsGenerator.cached_sequencing_has_metrics_data[sequencing_id] = result
        return result

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
            if not sample_date:
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
            assert matching_index is not None
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
        self, sample: Dict[str, Any], input_samples: Sequence[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        sequencing = sample["sequencing"]
        samples = [
            other_sample
            for other_sample in input_samples
            if other_sample["project"] == sample["project"]
            and other_sample["patient"] == sample["patient"]
        ]

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
            matching = list(input_samples)

        return matching

    def _get_grouped_samples_with_controls(
        self, samples: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        primary_samples = []
        control_samples = []

        for sample in samples:
            try:
                barcoded = Db.to_barcoded(sample)
                if not barcoded:
                    continue

                if barcoded.tissue.is_normal():
                    control_samples.append(sample)
                else:
                    primary_samples.append((sample, barcoded))
            except Exception:
                continue

        if not primary_samples:
            return []

        primary_samples.sort(key=lambda sample: sample[1].get_barcode())

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
                    cur_sample
                    for (cur_sample, date_distance) in zip(samples, date_distances)
                    if cur_sample["sequencing"]["_id"] != sample["sequencing"]["_id"]
                    and date_distance == best_date_distance
                    and cur_sample["project"]["_id"] == sample["project"]["_id"]
                    and cur_sample["patient"]["_id"] == sample["patient"]["_id"]
                ]
            else:
                matchable_samples = [
                    cur_sample
                    for cur_sample in samples
                    if cur_sample["sequencing"]["_id"] != sample["sequencing"]["_id"]
                    and cur_sample["project"]["_id"] == sample["project"]["_id"]
                    and cur_sample["patient"]["_id"] == sample["patient"]["_id"]
                ]

            best_sample = self._get_best_matching_samples(sample, matchable_samples)[0]
            return Db.to_barcoded(best_sample)
        else:
            matchable_samples = [
                cur_sample
                for cur_sample in samples
                if cur_sample["sequencing"]["_id"] != sample["sequencing"]["_id"]
                and cur_sample["project"]["_id"] == sample["project"]["_id"]
                and cur_sample["patient"]["_id"] == sample["patient"]["_id"]
            ]

            best_samples = self._get_best_matching_samples(sample, matchable_samples)
            if best_samples:
                sample_analysis_date = self._get_best_analysis_date_for_sample(sample)
                if sample_analysis_date:
                    matchable_samples_dates = [
                        self._get_best_analysis_date_for_sample(matchable_sample)
                        for matchable_sample in matchable_samples
                    ]

                    def get_date_distance(
                        matchable_sample_date: Optional[datetime.date]
                    ) -> datetime.timedelta:
                        assert sample_analysis_date
                        if matchable_sample_date:
                            return abs(sample_analysis_date - matchable_sample_date)
                        else:
                            return datetime.timedelta(days=999999999)

                    best_sample_index = utils.argmin(
                        matchable_samples_dates, get_date_distance
                    )
                    assert best_sample_index is not None
                    best_sample = best_samples[best_sample_index]
                else:
                    # At this point I don't have any other choice to pick one.
                    # I take the first sample in lexicographic order just to
                    # have consistent outputs
                    def get_sample_barcode(sample_data: Dict[str, Any]) -> str:
                        sample_barcoded = Db.to_barcoded(sample_data)
                        assert sample_barcoded
                        return sample_barcoded.get_barcode()

                    best_samples.sort(key=get_sample_barcode)
                    best_sample = best_samples[0]
                return Db.to_barcoded(best_sample)
            else:
                return None

    @staticmethod
    def _get_sequencing_date(sequencing: Dict[str, Any]) -> Optional[datetime.date]:
        if "date" in sequencing:
            splitted_date = [int(value) for value in sequencing["date"].split("-")]
            return datetime.date(splitted_date[0], splitted_date[1], splitted_date[2])
            return datetime.date(splitted_date[0], splitted_date[1], splitted_date[2])
        else:
            return None

    def _get_analyses_dates(self, sample_data: Dict[str, Any]) -> List[datetime.date]:
        sequencing = sample_data["sequencing"]
        sequencing_id = sequencing["_id"]
        assert sequencing_id
        cached_date = ReportsGenerator.cached_sequencing_analysis_dates.get(
            sequencing_id
        )
        if cached_date is not None:
            return list(cached_date)

        dates: List[datetime.date] = []
        self._populate_sequencing_to_analyses_cache()
        analyses_id_with_date = ReportsGenerator.cached_sequencings_to_analyses[
            sequencing_id
        ]
        assert analyses_id_with_date is not None
        if not analyses_id_with_date:
            barcoded_sample = Db.to_barcoded(sample_data)
            assert barcoded_sample
            assert barcoded_sample.tissue
            if barcoded_sample.tissue.is_normal():
                related_biopsies = self.db.biopsies.find_all(
                    {
                        "patient": sample_data["patient"]["_id"],
                        "tissue": {"$ne": sample_data["biopsy"]["tissue"]},
                    }
                )
                assert related_biopsies is not None

                related_sample_ids = set()
                for related_biopsy in related_biopsies:
                    current_samples = self.db.samples.find_all(
                        {"biopsy": related_biopsy["_id"]}
                    )
                    assert current_samples is not None
                    current_sample_ids = [
                        current_sample["_id"] for current_sample in current_samples
                    ]
                    related_sample_ids.update(current_sample_ids)

                related_sequencing_ids = set()
                for related_sample_id in related_sample_ids:
                    current_sequencings = self.db.sequencings.find_all(
                        {"sample": related_sample_id}
                    )
                    assert current_sequencings is not None
                    current_sequencing_ids = [
                        current_sequencing["_id"]
                        for current_sequencing in current_sequencings
                    ]
                    related_sequencing_ids.update(current_sequencing_ids)

                related_samples_data: List[Dict[str, Dict[str, Any]]] = []
                for sequencing_id in related_sequencing_ids:
                    current_sample_data = self.db.from_sequencing_id(sequencing_id)
                    assert current_sample_data
                    related_samples_data.append(current_sample_data)

                best_related_samples_data = self._get_best_matching_samples(
                    sample_data, related_samples_data
                )
                analyses_id_with_date = []
                for best_related_sample_data in best_related_samples_data:
                    current_analyses = ReportsGenerator.cached_sequencings_to_analyses[
                        best_related_sample_data["sequencing"]["_id"]
                    ]
                    analyses_id_with_date += current_analyses

        for analysis_id, analysis_date in analyses_id_with_date:
            splitted_date = [int(value) for value in analysis_date.split("_")]
            dates.append(
                datetime.date(splitted_date[0], splitted_date[1], splitted_date[2])
            )

        ReportsGenerator.cached_sequencing_analysis_dates[sequencing_id] = list(dates)
        return dates

    def _add_variants_tables(
        self, samples: Iterable[Dict[str, Any]], figures: FiguresType
    ) -> None:
        self._populate_sequencing_to_analyses_cache()
        for sample in samples:
            analyses = self._get_cached_analyses_from_sequencing(
                sample["sequencing"]["_id"]
            )

            if not analyses:
                continue

            for analysis in analyses:
                assert analysis
                for annotation_id in analysis["annotations"]:
                    if annotation_id not in ReportsGenerator.cached_annotations_ids:
                        annotation = self.db.annotations.find({"_id": annotation_id})
                        assert annotation
                        ReportsGenerator.cached_annotations[
                            annotation["id"]
                        ] = annotation
                        ReportsGenerator.cached_annotations_ids[
                            annotation_id
                        ] = annotation["id"]

            variants = []
            for analysis in analyses:
                for variant in analysis["variants"]:
                    variant_key = variant.get("key")
                    if not variant_key:
                        continue

                    annotation = ReportsGenerator.cached_annotations.get(variant_key)
                    if not annotation:
                        continue

                    if annotation.get("damaging") == "High":
                        variants.append(variant)

            if not variants:
                continue

            variants_annotations = [
                ReportsGenerator.cached_annotations[variant["key"]]
                for variant in variants
            ]
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
            report_table.set_order(
                ReportTableColumnOrdering(9, ReportTableSortingType.DESCENDING)
            )
            report_table.row_modifier = (
                lambda row: " class='druggable-gene'" if row["druggable"] else ""
            )

            for (variant, annotation, cosmic_id) in zip(
                variants, variants_annotations, cosmic_ids
            ):
                depth = variant.get("DP")
                report_table.add_row(
                    [
                        annotation["Chr"],
                        annotation["Start"],
                        annotation["Ref"],
                        annotation["Alt"],
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
                        depth if depth is not None else "",
                    ]
                )

            report_table.style = (
                "g.legend rect.bg { opacity: 0.7; }\n"
                "table.dataTable tbody tr.druggable-gene.even, "
                "table.dataTable tbody tr.druggable-gene.even td"
                "{ background-color: #83ff83; }\n"
                "table.dataTable tbody tr.druggable-gene.odd, "
                "table.dataTable tbody tr.druggable-gene.odd td"
                "{ background-color: #80e780; }"
            )

            analysis_date_obj = self._get_best_analysis_date_for_sample(sample)
            assert analysis_date_obj
            analysis_date = analysis_date_obj.strftime("%Y_%m_%d")
            figures.setdefault(sample_barcode, []).append(
                FigureData(report_table, analysis_date)
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
                "<title>HaTSPiL report</title>"
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
                assert all(styles)
                out_fd.write("<style>{}</style>".format("\n".join(styles)))
            out_fd.write("</head><body><div id='top-tabs'><ul>")

            if figures_by_sequencing:
                out_fd.write("<li><a href='#tab-by-sequencing'>By sequencing</a></li>")

            if figures_by_sample:
                out_fd.write("<li><a href='#tab-by-sample'>By sample</a></li>")

            if figures_by_biopsy:
                out_fd.write("<li><a href='#tab-by-biopsy'>By biopsy</a></li>")

            if figures_by_patient:
                out_fd.write("<li><a href='#tab-by-patient'>By patient</a></li>")

            if figures_by_project:
                out_fd.write("<li><a href='#tab-by-project'>By project</a></li>")

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
                "</body>"
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
                        "<li><a href='#{}_by_{}'>{}</a></li>".format(
                            ReportsGenerator.RE_INVALID_CHAR.sub("_", sample_name),
                            figure_type,
                            sample_name,
                        )
                        for sample_name, sample_figures in figures_data.items()
                        if sample_figures
                    ),
                )
            )

            def figure_get_date_key(figure_data: FigureData) -> datetime.date:
                if figure_data.date:
                    return datetime.datetime.strptime(figure_data.date, "%Y_%m_%d")
                else:
                    return datetime.date(datetime.MINYEAR, 1, 1)

            for sample_name, figures in figures_data.items():
                if not figures:
                    continue

                figures.sort(key=figure_get_date_key)

                fd.write(
                    "<div id='{}_by_{}'>".format(
                        ReportsGenerator.RE_INVALID_CHAR.sub("_", sample_name),
                        figure_type,
                    )
                )
                last_date: Optional[str] = None
                for figure in figures:
                    if figure.date:
                        figure_date = figure.date
                    else:
                        figure_date = "Unknown date"

                    if not last_date or last_date != figure_date:
                        fd.write("<h3>{}</h3>".format(figure_date))
                        last_date = figure_date

                    if isinstance(figure.figure, go.Figure):
                        fd.write(
                            plt.plot(
                                figure.figure,
                                output_type="div",
                                include_plotlyjs=include_plotlyjs,
                            )
                        )
                        include_plotlyjs = False
                    else:
                        assert isinstance(figure.figure, ReportTable)
                        assert figure_type == "sequencing"
                        fd.write(figure.figure.html())

                fd.write("</div>")

            fd.write(
                "</div>"
                "</div>"
                "<script>"
                "$(function() {{ $('#by_{}_tabs').tabs(); }} )"
                "</script>".format(figure_type)
            )
        return include_plotlyjs

    def _get_best_analysis_date_for_sample(
        self, sample: Dict[str, Any]
    ) -> Optional[datetime.date]:
        analyses_dates = self._get_analyses_dates(sample)
        if not analyses_dates:
            return None

        current_date = datetime.datetime.strptime(
            utils.get_overridable_current_date(self.parameters), "%Y_%m_%d"
        ).date()
        analyses_dates.sort(key=lambda date: abs(current_date - date))
        analysis_date = next(
            (date for date in analyses_dates if date <= current_date), analyses_dates[0]
        )

        return analysis_date

    def _populate_sequencing_to_analyses_cache(self) -> None:
        if ReportsGenerator.cached_sequencings_to_analyses:
            return

        analyses_iter = self.db.analyses.iter()
        assert analyses_iter
        for analysis in analyses_iter:
            ReportsGenerator.cached_sequencings_to_analyses.setdefault(
                analysis["sequencing"], []
            ).append((analysis["_id"], analysis["date"]))

        sequencings_iter = self.db.sequencings.iter()
        assert sequencings_iter
        for sequencing in sequencings_iter:
            ReportsGenerator.cached_sequencings_to_analyses.setdefault(
                sequencing["_id"], []
            )

    def _get_cached_analyses_from_sequencing(
        self, sequencing_id: ObjectId
    ) -> List[Dict[str, Any]]:
        self._populate_sequencing_to_analyses_cache()
        analyses_ids_and_dates = ReportsGenerator.cached_sequencings_to_analyses[
            sequencing_id
        ]
        current_analysis_ids = [
            analysis_id for analysis_id, analysis_date in analyses_ids_and_dates
        ]
        if current_analysis_ids:
            out = self.db.analyses.find_all({"_id": {"$in": current_analysis_ids}})
        else:
            out = []
        assert out is not None
        return out
