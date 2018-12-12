import os
import re
from enum import Enum
from typing import (Any, Callable, Dict, Iterable, List, Mapping, Optional,
                    Sequence, Tuple, Union, cast)

from . import utils
from .analysis import Analysis
from .barcoded_filename import BarcodedFilename
from .exceptions import DataError, PipelineError


class AnalysisType(Enum):
    Unspecified = 0
    Sample = 1
    Control = 2


class AnalysisFileData:
    def __init__(self, filename: str) -> None:
        self.filename = filename
        try:
            self.barcode = BarcodedFilename(filename)
            if self.barcode.tissue.is_normal():
                self.type = AnalysisType.Control
            elif self.barcode.tissue.is_tumor():
                self.type = AnalysisType.Sample
            else:
                self.type = AnalysisType.Unspecified
        except Exception:
            self.type = AnalysisType.Unspecified

    def __repr__(self) -> str:
        return self.filename


class SingleAnalysis(List[AnalysisFileData]):
    @property
    def sample(self) -> Optional[AnalysisFileData]:
        return next(
            (file_data for file_data in self if file_data.type == AnalysisType.Sample),
            None,
        )

    @property
    def control(self) -> Optional[AnalysisFileData]:
        return next(
            (file_data for file_data in self if file_data.type == AnalysisType.Control),
            None,
        )


AnalysesPerOrganism = Dict[str, List[SingleAnalysis]]


class ExecutorData:
    def __init__(
        self,
        command: Union[str, List[str], Callable[..., None]],
        output_format: Union[
            str, Callable[..., str], List[Union[str, Callable[..., str]]], None
        ],
        input_filenames: Optional[Sequence[str]],
        input_function: Optional[Callable[[Union[str, List[str]]], Optional[str]]],
        input_split_reads: bool,
        output_path: Optional[str],
        output_function: Optional[Callable[[str], Iterable[str]]],
        error_string: Optional[str],
        exception_string: Optional[str],
        override_last_files: bool,
        write_bam_files: bool,
        unlink_inputs: bool,
        save_only_last: bool,
        use_normals: bool,
        split_by_organism: bool,
        only_human: bool,
        split_input_files: bool,
        allow_raw_filenames: bool,
    ) -> None:

        self.command = command
        self.output_format = output_format
        self.input_filenames = input_filenames
        self.input_function = input_function
        self.input_split_reads = input_split_reads
        self.output_path = output_path
        self.output_function = output_function
        self.error_string = error_string
        self.exception_string = exception_string
        self.override_last_files = override_last_files
        self.write_bam_files = write_bam_files
        self.unlink_inputs = unlink_inputs
        self.save_only_last = save_only_last
        self.use_normals = use_normals
        self.split_by_organism = split_by_organism
        self.only_human = only_human
        self.split_input_files = split_input_files
        self.allow_raw_filenames = allow_raw_filenames


class Executor:
    RE_REPLACER = re.compile(r"\{([^}]+)\}")

    def __init__(self, analysis: Analysis) -> None:
        self.analysis = analysis
        self.data: Optional[ExecutorData] = None

    def _handle_output_filename(
        self,
        command_index: int,
        commands_len: int,
        organism: str,
        output_filename: Union[List[str], str],
        output_filenames: Dict[str, List[str]],
        output_bamfiles: Dict[str, List[str]],
    ) -> None:
        assert self.data

        if isinstance(output_filename, str):
            output_filename = [output_filename]

        new_filename: List[str] = []
        for filename in output_filename:
            dirname = os.path.dirname(filename)
            if dirname == "" or dirname == ".":
                new_filename.append(os.path.join(os.getcwd(), filename))
            else:
                new_filename.append(filename)
        output_filename = new_filename

        if not self.data.save_only_last or command_index == commands_len - 1:
            for filename in output_filename:
                if self.data.split_by_organism:
                    try:
                        output_organism = BarcodedFilename(filename).organism
                        if not output_organism:
                            output_organism = organism
                    except Exception:
                        output_organism = organism
                else:
                    output_organism = organism

                output_filenames.setdefault(output_organism, []).append(filename)
                output_extension = os.path.splitext(filename)[1]
                if output_extension == ".bam":
                    output_bamfiles.setdefault(output_organism, []).append(filename)

    def _get_output_filename(
        self, all_params: Mapping[str, Any]
    ) -> Union[List[str], str, None]:
        assert self.data

        locals().update(all_params)
        if self.data.output_format is not None:

            if isinstance(self.data.output_format, list):
                raw_output_formats = self.data.output_format[:]
            else:
                raw_output_formats = [self.data.output_format]

            output_formats: List[str] = []
            for raw_output_format in raw_output_formats:
                if isinstance(raw_output_format, str):
                    output_formats.append(raw_output_format)
                else:
                    # output_format is a Callable[..., str]
                    output_formats.append(raw_output_format(**all_params))

            if self.data.output_path is not None:
                output_formats = [
                    os.path.join(self.data.output_path, output_format)
                    for output_format in output_formats
                ]

            output_filename = []
            for s in output_formats:
                for match in Executor.RE_REPLACER.finditer(s):
                    try:
                        evaluated = eval(match.group(1))
                    except Exception:
                        raise PipelineError("cannot evaluate %s" % match.group(1))

                    if evaluated is None:
                        raise PipelineError("evaluation of %s is None" % match.group(1))
                    s = re.sub(match.group(0), str(evaluated), s)

                output_filename.append(s)

            if self.data.output_function is not None:
                output_filename = [
                    filename
                    for filenames in map(self.data.output_function, output_filename)
                    for filename in filenames
                ]

            if len(output_filename) == 1:
                return output_filename[0]
            else:
                return output_filename

        elif self.data.output_function is not None:
            input_filenames = all_params["input_filenames"]
            output_filenames = [
                filename
                for filenames in map(
                    self.data.output_function,
                    [filename.filename for filename in input_filenames],
                )
                for filename in filenames
            ]

            if len(output_filenames) == 1:
                return output_filenames[0]
            else:
                return output_filenames
        else:
            return None

    def _unlink_filename(
        self,
        input_filename: SingleAnalysis,
        real_input_filename: Optional[SingleAnalysis],
    ) -> None:
        assert self.data

        if (
            self.data.unlink_inputs
            and self.analysis.can_unlink
            and not self.analysis.run_fake
        ):
            if self.data.input_function is not None:
                assert real_input_filename
                input_filename = real_input_filename

            for file_data in input_filename:
                filename = file_data.filename
                extension = os.path.splitext(filename)[1].lower()
                os.unlink(filename)
                if extension == ".bam":
                    bai_file = filename[:-4] + ".bai"
                    if os.path.exists(bai_file):
                        os.unlink(bai_file)

    def _handle_command(
        self,
        current_command: Union[str, Callable[..., None]],
        all_params: Mapping[str, Any],
    ) -> None:
        assert self.data

        locals().update(all_params)
        if isinstance(current_command, str):
            if not self.analysis.run_fake:
                status = utils.run_and_log(current_command, self.analysis.logger)
            else:
                self.analysis.logger.info("Faking command '%s'", current_command)
                status = 0
        else:
            if not self.analysis.run_fake:
                current_command(**all_params)
            else:
                self.analysis.logger.info("Faking lambda")
            status = 0

        if status != 0:
            assert isinstance(self.data.command, str)
            arg_zero = os.path.basename(self.data.command.split(" ")[0])

            if self.data.error_string is None:
                error_string = "%s exited with status %d" % (arg_zero, status)
            else:
                error_string = self.data.error_string

            if self.data.exception_string is None:
                exception_string = "%s error" % arg_zero
            else:
                exception_string = self.data.exception_string

            for match in Executor.RE_REPLACER.finditer(error_string):
                try:
                    error_string = re.sub(
                        match.group(0), str(eval(match.group(1))), error_string
                    )
                except Exception:
                    raise PipelineError(
                        "cannot replace parameter %s" % (match.group(0))
                    )

            for match in Executor.RE_REPLACER.finditer(exception_string):
                try:
                    exception_string = re.sub(
                        match.group(0), str(eval(match.group(1))), exception_string
                    )
                except Exception:
                    raise PipelineError(
                        "cannot replace parameter %s" % (match.group(0))
                    )
            self.analysis.logger.error(error_string)
            raise PipelineError(exception_string)

    def _create_mod_input_filenames(
        self, input_filenames: AnalysesPerOrganism
    ) -> AnalysesPerOrganism:
        assert self.data
        assert self.data.input_function

        self._fix_input_filenames(input_filenames)

        mod_input_filenames: AnalysesPerOrganism = {}
        for organism, analyses in input_filenames.items():
            mod_analyses = mod_input_filenames.setdefault(organism, [])

            for analysis in analyses:
                filenames = [analysis_file.filename for analysis_file in analysis]

                if self.data.input_split_reads:
                    splitted_data: Dict[int, List[str]] = {}
                    for filename in filenames:
                        try:
                            barcoded = BarcodedFilename(filename)
                            if barcoded.read_index:
                                splitted_data.setdefault(
                                    barcoded.read_index, []
                                ).append(filename)
                            else:
                                splitted_data.setdefault(0, []).append(filename)
                        except Exception:
                            splitted_data.setdefault(0, []).append(filename)

                    for filenames in splitted_data.values():
                        param: Union[str, List[str]] = list(filenames)
                        if len(param) == 1:
                            param = param[0]

                        input_str = self.data.input_function(param)
                        if input_str:
                            mod_analyses.append(
                                SingleAnalysis([AnalysisFileData(input_str)])
                            )
                else:
                    result = cast(Callable[[List[str]], str], self.data.input_function)(
                        filenames
                    )
                    mod_analyses.append(SingleAnalysis([AnalysisFileData(result)]))

        if not mod_input_filenames:
            raise PipelineError("empty input list")

        self._fix_input_filenames(mod_input_filenames)
        return mod_input_filenames

    def _fix_input_filenames(self, input_filenames: AnalysesPerOrganism) -> None:
        assert self.data

        if not input_filenames:
            raise PipelineError("empty input list")

        if not self.data.input_split_reads:
            for organism, analyses in input_filenames.items():
                input_filenames[organism] = [
                    SingleAnalysis(
                        [file_data for analysis in analyses for file_data in analysis]
                    )
                ]

    @staticmethod
    def _get_fixed_normals_analyses(
        analyses: List[SingleAnalysis]
    ) -> List[SingleAnalysis]:

        for analysis_index, analysis in enumerate(analyses):
            if not analysis:
                continue

            file_data = next(iter(analysis))
            if file_data.type != AnalysisType.Sample:
                continue

            controls = [
                (other_analysis, other_file_data)
                for other_analysis in analyses
                for other_file_data in other_analysis
                if file_data.barcode.equals_without_tissue(other_file_data.barcode)
                and other_file_data.barcode.tissue.is_normal()
            ]

            if len(controls) == 1:
                control = controls[0]
                analysis.append(control[1])
                control[0].remove(control[1])
            else:
                sequencing_specific = [
                    control
                    for control in controls
                    if control[1].barcode.sequencing == file_data.barcode.sequencing
                ]

                if sequencing_specific:
                    controls = sequencing_specific

                analyses[analysis_index] = SingleAnalysis(
                    analysis + [control[1] for control in controls]
                )

                for control in controls:
                    control[0].remove(control[1])

        return [analysis for analysis in analyses if analysis]

    def _get_input_filenames(self) -> Tuple[AnalysesPerOrganism, AnalysesPerOrganism]:
        assert self.data

        if self.data.input_filenames is None:
            if self.analysis.last_operation_filenames is None:
                raise PipelineError(
                    "input files missing and last_operation_filenames empty"
                )

            raw_input_filenames = utils.get_sample_filenames(
                self.analysis.last_operation_filenames, self.data.split_by_organism
            )
        else:
            raw_input_filenames = utils.get_sample_filenames(
                self.data.input_filenames, self.data.split_by_organism
            )

        input_filenames: AnalysesPerOrganism = {}
        if isinstance(raw_input_filenames, dict):
            for organism, filenames in raw_input_filenames.items():
                analyses = input_filenames.setdefault(organism, [])
                if self.data.split_input_files:
                    for filename in filenames:
                        analyses.append(SingleAnalysis([AnalysisFileData(filename)]))
                else:
                    analyses.append(
                        SingleAnalysis(
                            [AnalysisFileData(filename) for filename in filenames]
                        )
                    )
        else:
            analyses = input_filenames.setdefault("", [])
            if self.data.split_input_files:
                for filename in raw_input_filenames:
                    analyses.append(SingleAnalysis([AnalysisFileData(filename)]))
            else:
                analyses.append(
                    SingleAnalysis(
                        [AnalysisFileData(filename) for filename in raw_input_filenames]
                    )
                )

        if self.analysis.parameters["use_normals"] and self.data.use_normals:
            for organism, analyses in input_filenames.items():
                input_filenames[organism] = Executor._get_fixed_normals_analyses(
                    analyses
                )

        if self.data.input_function is not None:
            return (input_filenames, self._create_mod_input_filenames(input_filenames))
        else:
            self._fix_input_filenames(input_filenames)
            return (input_filenames, {})

    def _get_additional_params(self, organism: Optional[str]) -> Dict[str, str]:
        additional_params = {}

        if not organism:
            additional_params["organism_str"] = ""
            organism = utils.get_human_annotation(self.analysis.config)
        else:
            additional_params["organism_str"] = "." + organism

        try:
            genome_ref, genome_index = utils.get_genome_ref_index_by_organism(
                self.analysis.config, organism
            )
            additional_params["genome_ref"] = genome_ref
            additional_params["genome_index"] = genome_index
        except DataError:
            pass

        try:
            additional_params["dbsnp"] = utils.get_dbsnp_by_organism(
                self.analysis.config, organism
            )
        except DataError:
            pass

        try:
            additional_params["cosmic"] = utils.get_cosmic_by_organism(
                self.analysis.config, organism
            )
        except DataError:
            pass

        return additional_params

    def _get_commands(
        self, all_params: Mapping[str, Any]
    ) -> List[Union[str, Callable[..., None]]]:
        assert self.data
        assert self.data.command

        locals().update(all_params)

        commands: List[Union[str, Callable[..., None]]] = []
        if isinstance(self.data.command, list):
            for s in self.data.command:
                if isinstance(s, str):
                    for match in Executor.RE_REPLACER.finditer(s):
                        try:
                            evaluated = eval(match.group(1))
                        except Exception:
                            raise PipelineError("cannot evaluate %s" % match.group(1))

                        if evaluated is None:
                            raise PipelineError(
                                "evaluation of %s is None" % (match.group(1))
                            )

                        s = s.replace(match.group(0), str(evaluated))

                commands.append(s)
        else:
            if isinstance(self.data.command, str):
                current_command = str(self.data.command)
                for match in Executor.RE_REPLACER.finditer(self.data.command):
                    try:
                        evaluated = eval(match.group(1))
                    except Exception:
                        raise PipelineError("cannot evaluate %s" % (match.group(1)))

                    if evaluated is None:
                        raise PipelineError(
                            "evaluation of %s is None" % (match.group(1))
                        )

                    current_command = current_command.replace(
                        match.group(0), str(evaluated)
                    )

                commands.append(current_command)
            else:
                commands.append(self.data.command)

        return commands

    def _handle_analysis(
        self,
        analysis_input: SingleAnalysis,
        mod_analysis_input: Optional[SingleAnalysis],
        output_filenames: Dict[str, List[str]],
        output_bamfiles: Dict[str, List[str]],
        analyses: AnalysesPerOrganism,
    ) -> None:
        assert self.data

        real_analysis_input: Optional[SingleAnalysis]
        if self.data.input_function is not None:
            assert mod_analysis_input

            real_analysis_input = analysis_input
            analysis_input = mod_analysis_input
        else:
            real_analysis_input = None

        input_filenames = analysis_input
        if len(analysis_input) == 1:
            input_filename = analysis_input[0]

        file_data = analysis_input.sample
        if not file_data:
            file_data = analysis_input.control

        if not self.data.allow_raw_filenames:
            assert file_data

        if file_data:
            organism = file_data.barcode.organism
        else:
            organism = None

        if not self.data.only_human or not organism or organism.startswith("hg"):

            locals().update(self._get_additional_params(organism))
            if not organism:
                organism = utils.get_human_annotation(self.analysis.config)

            local_params = {
                key: value for key, value in locals().items() if key != "self"
            }
            output_filename = self._get_output_filename(local_params)

            local_params.update(
                {"output_filename": output_filename, "config": self.analysis.config}
            )
            commands = self._get_commands(local_params)

            for command_index, current_command in enumerate(commands):
                self._handle_command(
                    current_command,
                    {key: value for key, value in locals().items() if key != "self"},
                )

                if output_filename:
                    self._handle_output_filename(
                        command_index,
                        len(commands),
                        organism,
                        output_filename,
                        output_filenames,
                        output_bamfiles,
                    )

        self._unlink_filename(analysis_input, real_analysis_input)

    def __call__(
        self,
        command: Union[str, List[str], Callable[..., None]],
        output_format: Union[
            str, Callable[..., str], List[Union[str, Callable[..., str]]], None
        ] = None,
        input_filenames: Optional[Sequence[str]] = None,
        input_function: Optional[
            Callable[[Union[str, List[str]]], Optional[str]]
        ] = None,
        input_split_reads: bool = True,
        output_path: Optional[str] = None,
        output_function: Optional[Callable[[str], Iterable[str]]] = None,
        error_string: Optional[str] = None,
        exception_string: Optional[str] = None,
        override_last_files: bool = True,
        write_bam_files: bool = True,
        unlink_inputs: bool = False,
        save_only_last: bool = True,
        use_normals: bool = False,
        split_by_organism: bool = False,
        only_human: bool = False,
        split_input_files: bool = True,
        allow_raw_filenames: bool = False,
    ) -> None:

        self.data = ExecutorData(
            command,
            output_format,
            input_filenames,
            input_function,
            input_split_reads,
            output_path,
            output_function,
            error_string,
            exception_string,
            override_last_files,
            write_bam_files,
            unlink_inputs,
            save_only_last,
            use_normals,
            split_by_organism,
            only_human,
            split_input_files,
            allow_raw_filenames,
        )

        _input_filenames, mod_input_filenames = self._get_input_filenames()

        output_filenames: Dict[str, List[str]] = {}
        output_bamfiles: Dict[str, List[str]] = {}
        for current_organism in _input_filenames.keys():
            current_input_filenames = _input_filenames[current_organism]

            iterator: Union[
                Iterable[Tuple[SingleAnalysis, SingleAnalysis]],
                Iterable[Tuple[SingleAnalysis, None]],
            ]
            if input_function is not None:
                current_mod_input_filenames = mod_input_filenames[current_organism]
                if len(current_mod_input_filenames) == len(current_input_filenames):
                    iterator = zip(current_input_filenames, current_mod_input_filenames)
                else:
                    iterator = zip(
                        current_mod_input_filenames, current_mod_input_filenames
                    )
            else:
                iterator = zip(
                    current_input_filenames, [None] * len(current_input_filenames)
                )

            for input_filename, mod_input_filename in iterator:
                self._handle_analysis(
                    input_filename,
                    mod_input_filename,
                    output_filenames,
                    output_bamfiles,
                    _input_filenames,
                )

        if override_last_files:
            self.analysis.last_operation_filenames = output_filenames
            self.analysis.can_unlink = True

        if write_bam_files and len(output_bamfiles) != 0:
            self.analysis.bamfiles = output_bamfiles

    def override_last_operation_filename(self, new_filename: str) -> None:
        if not self.analysis.last_operation_filenames:
            raise PipelineError("last operation did not leave an output file")

        if isinstance(self.analysis.last_operation_filenames, str):
            self.analysis.last_operation_filenames = new_filename
        elif isinstance(self.analysis.last_operation_filenames, list):
            if len(self.analysis.last_operation_filenames) != 1:
                raise PipelineError(
                    "last operation created a list with a number of output "
                    "files different than one"
                )

            self.analysis.last_operation_filenames[0] = new_filename
        elif isinstance(self.analysis.last_operation_filenames, dict):
            if len(self.analysis.last_operation_filenames) != 1:
                raise PipelineError(
                    "last operation created a dict using more than one organism"
                )

            organism, last_filenames = next(
                iter(self.analysis.last_operation_filenames.items())
            )
            if isinstance(last_filenames, str):
                self.analysis.last_operation_filenames[organism] = new_filename
            elif isinstance(last_filenames, list):
                if len(last_filenames) != 1:
                    raise PipelineError(
                        "last operation created a dict of lists with one "
                        "list, but the list contains a number of filenames "
                        "different than one"
                    )
                self.analysis.last_operation_filenames[organism][0] = new_filename
            else:
                raise PipelineError("last operation created an invalid dict")
        else:
            raise PipelineError("last operation created an invalid object")
