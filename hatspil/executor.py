import os
import re
from typing import Callable, Iterable, List, Optional, Sequence, Union

from . import utils
from .analysis import Analysis
from .barcoded_filename import BarcodedFilename
from .exceptions import DataError, PipelineError


class ExecutorData:
    def __init__(self,
                 command: Union[str, List[str], Callable[..., None]],
                 output_format: Optional[str],
                 input_filenames: Optional[Sequence[str]],
                 input_function: Optional[Callable[[str],
                                                   Optional[str]]],
                 input_split_reads: bool,
                 output_path: Optional[str],
                 output_function: Optional[Callable[[str],
                                                    Iterable[str]]],
                 error_string: Optional[str],
                 exception_string: Optional[str],
                 override_last_files: bool,
                 write_bam_files: bool,
                 unlink_inputs: bool,
                 save_only_last: bool,
                 use_normals: bool,
                 split_by_organism: bool,
                 only_human: bool) -> None:

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


class Executor:
    RE_REPLACER = re.compile(r"\{([^}]+)\}")

    def __init__(self, analysis: Analysis) -> None:
        self.analysis = analysis
        self.data: Optional[ExecutorData] = None

    def _handle_output_filename(self,
                                command_index,
                                commands,
                                organism,
                                output_filename,
                                output_filenames,
                                output_bamfiles):

        if not isinstance(output_filename, list):
            output_filename = [output_filename]

        new_filename = []
        for filename in output_filename:
            dirname = os.path.dirname(filename)
            if dirname == "" or dirname == ".":
                new_filename.append(
                    os.path.join(os.getcwd(), filename))
            else:
                new_filename.append(filename)
        output_filename = new_filename

        if self.data.save_only_last:
            if command_index == len(commands) - 1:
                if organism in output_filenames:
                    output_filenames[
                        organism] += output_filename
                else:
                    output_filenames[
                        organism] = output_filename

                for filename in output_filename:
                    output_extension = \
                        os.path.splitext(filename)[1]
                    if output_extension == ".bam":
                        if organism in output_bamfiles:
                            output_bamfiles[organism] += \
                                output_filename
                        else:
                            output_bamfiles[organism] = \
                                output_filename

        else:
            if organism in output_filenames:
                output_filenames[
                    organism] += output_filename
            else:
                output_filenames[
                    organism] = output_filename

            for filename in output_filename:
                output_extension = os.path.splitext(
                    filename)[1]
                if output_extension == ".bam":
                    if organism in output_bamfiles:
                        output_bamfiles[organism] += \
                            output_filename
                    else:
                        output_bamfiles[organism] = \
                            output_filename

    def _get_output_filename(self, all_params):
        locals().update(all_params)
        if self.data.output_format is not None:
            if not isinstance(self.data.output_format, list):
                output_format = [self.data.output_format]
            else:
                output_format = list(self.data.output_format)

            if self.data.output_path is not None:
                output_format = [
                    os.path.join(self.data.output_path, s)
                    for s in output_format
                ]

            output_filename = []
            for s in output_format:
                for match in Executor.RE_REPLACER.finditer(s):
                    try:
                        evaluated = eval(match.group(1))
                    except Exception:
                        raise PipelineError(
                            "cannot evaluate %s" % match.group(1))

                    if evaluated is None:
                        raise PipelineError(
                            "evaluation of %s is None" %
                            match.group(1))
                    s = re.sub(match.group(0), str(evaluated), s)

                output_filename.append(s)

            if self.data.output_function is not None:
                output_filename = [
                    filename for filename in map(
                        self.data.output_function, output_filename)
                ]

            if len(output_filename) == 1:
                output_filename = output_filename[0]

            return output_filename
        else:
            return None

    def _unlink_filename(self, input_filename, real_input_filename):
        if self.data.unlink_inputs and self.analysis.can_unlink and\
                not self.analysis.run_fake:
            if self.data.input_function is not None:
                input_filename = real_input_filename

            if not isinstance(input_filename, list):
                input_filename = [input_filename]

            new_input_filename = []
            for filename in input_filename:
                if isinstance(filename, str):
                    new_input_filename.append(filename)
                elif isinstance(filename, dict):
                    new_input_filename += filename.values()
                else:
                    new_input_filename += filename
            input_filename = new_input_filename

            for filename in input_filename:
                extension = os.path.splitext(filename)[1].lower()
                os.unlink(filename)
                if extension == ".bam":
                    bai_file = filename[:-4] + ".bai"
                    if os.path.exists(bai_file):
                        os.unlink(bai_file)

    def _handle_command(self, current_command, all_params):
        locals().update(all_params)
        if isinstance(current_command, str):
            if not self.analysis.run_fake:
                status = utils.run_and_log(
                    current_command, self.analysis.logger)
            else:
                self.analysis.logger.info(
                    "Faking command '%s'", current_command)
                status = 0
        else:
            if not self.analysis.run_fake:
                current_command(**all_params)
            else:
                self.analysis.logger.info("Faking lambda")
            status = 0

        if status != 0:
            arg_zero = os.path.basename(self.data.command.split(" ")[0])

            if self.data.error_string is None:
                error_string = "%s exited with status %d" % (
                    arg_zero, status)
            else:
                error_string = self.data.error_string

            if self.data.exception_string is None:
                exception_string = "%s error" % arg_zero
            else:
                exception_string = self.data.exception_string

            for match in Executor.RE_REPLACER.finditer(error_string):
                try:
                    error_string = re.sub(
                        match.group(0),
                        str(eval(match.group(1))),
                        error_string)
                except Exception:
                    raise PipelineError(
                        "cannot replace parameter %s" %
                        (match.group(0)))

            for match in Executor.RE_REPLACER.finditer(
                    exception_string):
                try:
                    exception_string = re.sub(
                        match.group(0),
                        str(eval(match.group(1))),
                        exception_string)
                except Exception:
                    raise PipelineError(
                        "cannot replace parameter %s" %
                        (match.group(0)))
            self.analysis.logger.error(error_string)
            raise PipelineError(exception_string)

    def _create_mod_input_filenames(self, input_filenames):
        mod_input_filenames = {}
        if self.data.input_split_reads:
            for organism, current_input_filenames in input_filenames.items(
            ):
                mod_input_filenames[organism] = \
                    [filename for filename
                     in map(self.data.input_function,
                            current_input_filenames)
                     if filename is not None]
        else:
            for organism, current_input_filenames in input_filenames.items(
            ):
                current_mod_input = self.data.input_function(
                    current_input_filenames)
                if isinstance(current_mod_input, list):
                    current_mod_input.remove(None)
                    current_mod_input.remove("")
                mod_input_filenames[organism] = current_mod_input

        if mod_input_filenames is None or \
                len(mod_input_filenames) == 0:
            raise PipelineError("empty input list")

        first_mod_input_filenames = list(mod_input_filenames.values())[0]
        if (isinstance(first_mod_input_filenames, str) and
                first_mod_input_filenames == "") or\
                (isinstance(first_mod_input_filenames, list) and
                 len(first_mod_input_filenames) == 0):
            raise PipelineError("empty input list")

        if not self.data.input_split_reads:
            for organism, current_list in input_filenames.items():
                input_filenames[organism] = [current_list]
            for organism, current_list in mod_input_filenames.items():
                mod_input_filenames[organism] = [current_list]

        return mod_input_filenames

    def _fix_input_filenames(self, input_filenames):
        if input_filenames is None or \
                len(input_filenames) == 0:
            raise PipelineError("empty input list")

        first_input_filenames = list(input_filenames.values())[0]
        if (isinstance(first_input_filenames, str) and
                first_input_filenames == "") or \
                (isinstance(first_input_filenames, list) and
                 len(first_input_filenames) == 0):
            raise PipelineError("empty input list")

        if not self.data.input_split_reads:
            for organism, current_list in input_filenames.items():
                input_filenames[organism] = [current_list]

    def _get_input_filenames(self):
        if self.data.input_filenames is None:
            if self.analysis.last_operation_filenames is None:
                raise PipelineError("input files missing and "
                                    "last_operation_filenames empty")

            input_filenames = utils.get_sample_filenames(
                self.analysis.last_operation_filenames,
                self.data.split_by_organism)
        else:
            input_filenames = utils.get_sample_filenames(
                self.data.input_filenames, self.data.split_by_organism)

        if not isinstance(input_filenames, dict):
            input_filenames = {"": input_filenames}

        if self.analysis.parameters["use_normals"] and self.data.use_normals:
            splitted_filenames = {}
            for organism, input_list in input_filenames.items():
                current_splitted = {}
                for filename in input_list:
                    barcoded_filename = BarcodedFilename(filename)

                    sample_type = barcoded_filename.tissue
                    if sample_type.is_normal():
                        current_splitted["normal"] = filename
                    else:
                        current_splitted["tumor"] = filename

                splitted_filenames[organism] = current_splitted

            input_filenames = {}
            for splitted_key, splitted_list in splitted_filenames.items():
                if len(splitted_list) == 1:
                    if list(splitted_list.keys())[0] == "tumor":
                        if splitted_key not in input_filenames:
                            input_filenames[splitted_key] = []
                        input_filenames[splitted_key].append(
                            list(splitted_list.values())[0])
                else:
                    if splitted_key not in input_filenames:
                        input_filenames[splitted_key] = []

                    input_filenames[splitted_key].append(splitted_list)

        if self.data.input_function is not None:
            return (input_filenames,
                    self._create_mod_input_filenames(input_filenames))
        else:
            self._fix_input_filenames(input_filenames)
            return (input_filenames, {})

    def _get_additional_params(self, organism):
        additional_params = {}

        if organism is None or organism == "":
            additional_params["organism_str"] = ""
            organism = "hg19"
        else:
            additional_params["organism_str"] = "." + organism

        try:
            genome_ref, genome_index = \
                utils.get_genome_ref_index_by_organism(
                    self.analysis.config, organism)
            additional_params["genome_ref"] = genome_ref
            additional_params["genome_index"] = genome_index
        except DataError:
            pass

        try:
            additional_params["dbsnp"] = utils.get_dbsnp_by_organism(
                self.analysis.config, organism)
        except DataError:
            pass

        try:
            additional_params["cosmic"] = utils.get_cosmic_by_organism(
                self.analysis.config, organism)
        except DataError:
            pass

        return additional_params

    def _get_commands(self, all_params):
        locals().update(all_params)

        commands = []
        if isinstance(self.data.command, list):
            for s in self.data.command:
                if isinstance(s, str):
                    for match in Executor.RE_REPLACER.finditer(s):
                        try:
                            evaluated = eval(match.group(1))
                        except Exception:
                            raise PipelineError(
                                "cannot evaluate %s" %
                                match.group(1))

                        if evaluated is None:
                            raise PipelineError(
                                "evaluation of %s is None" %
                                (match.group(1)))

                        s = s.replace(
                            match.group(0), str(evaluated))

                commands.append(s)
        else:
            current_command = self.data.command
            if isinstance(current_command, str):
                for match in Executor.RE_REPLACER.finditer(self.data.command):
                    try:
                        evaluated = eval(match.group(1))
                    except Exception:
                        raise PipelineError(
                            "cannot evaluate %s" %
                            (match.group(1)))

                    if evaluated is None:
                        raise PipelineError(
                            "evaluation of %s is None" %
                            (match.group(1)))

                    current_command = current_command.replace(
                        match.group(0), str(evaluated))
            commands.append(current_command)

        return commands

    def _handle_input_filename(self,
                               input_filename,
                               mod_input_filename,
                               output_filenames,
                               output_bamfiles,
                               input_filenames):

        if isinstance(input_filename, str):
            barcoded_filename = BarcodedFilename(input_filename)
            organism = barcoded_filename.organism
            read_index = barcoded_filename.read_index
        elif isinstance(input_filename, dict):
            barcoded_filename = BarcodedFilename(
                input_filename["tumor"])
            organism = barcoded_filename.organism
            read_index = barcoded_filename.read_index
        else:
            read_index = []
            for filename in input_filename:
                if isinstance(filename, str):
                    barcoded_filename = BarcodedFilename(filename)
                    organism = barcoded_filename.organism
                    index = barcoded_filename.read_index
                else:
                    barcoded_filename = BarcodedFilename(
                        filename["tumor"])
                    organism = barcoded_filename.organism
                    index = barcoded_filename.read_index
                read_index.append(index)

        if self.data.input_function is not None:
            real_input_filename = input_filename
            input_filename = mod_input_filename
        else:
            real_input_filename = None

        if not self.data.only_human \
                or organism is None \
                or organism == ""\
                or organism.startswith("hg"):

            locals().update(self._get_additional_params(organism))
            if organism is None or organism == "":
                organism = "hg19"

            local_params = {key: value for key,
                            value in locals().items() if key != "self"}
            output_filename = self._get_output_filename(local_params)

            local_params.update({"output_filename": output_filename})
            commands = self._get_commands(local_params)

            for command_index, current_command in enumerate(commands):
                self._handle_command(current_command, {
                                     key: value
                                     for key, value
                                     in locals().items()
                                     if key != "self"})

                if output_filename:
                    self._handle_output_filename(
                        command_index,
                        commands,
                        organism,
                        output_filename,
                        output_filenames,
                        output_bamfiles)

        self._unlink_filename(input_filename, real_input_filename)

    def __call__(self,
                 command: Union[str, List[str], Callable[..., None]],
                 output_format: Optional[str] = None,
                 input_filenames: Optional[Sequence[str]] = None,
                 input_function: Optional[Callable[[str],
                                                   Optional[str]]] = None,
                 input_split_reads: bool = True,
                 output_path: Optional[str] = None,
                 output_function: Optional[Callable[[str],
                                                    Iterable[str]]] = None,
                 error_string: Optional[str] = None,
                 exception_string: Optional[str] = None,
                 override_last_files: bool = True,
                 write_bam_files: bool = True,
                 unlink_inputs: bool = False,
                 save_only_last: bool = True,
                 use_normals: bool = False,
                 split_by_organism: bool = False,
                 only_human: bool = False) -> None:

        self.data = ExecutorData(command,
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
                                 only_human)

        input_filenames, mod_input_filenames = self._get_input_filenames()

        output_filenames = {}
        output_bamfiles = {}
        for current_organism in input_filenames.keys():
            current_input_filenames = input_filenames[current_organism]

            if input_function is not None:
                current_mod_input_filenames = mod_input_filenames[
                    current_organism]
                if len(current_mod_input_filenames) == len(
                        current_input_filenames):
                    iterator = zip(current_input_filenames,
                                   current_mod_input_filenames)
                else:
                    iterator = zip(current_mod_input_filenames,
                                   current_mod_input_filenames)
            else:
                iterator = zip(current_input_filenames,
                               [None] * len(current_input_filenames))

            for input_filename, mod_input_filename in iterator:
                self._handle_input_filename(
                    input_filename,
                    mod_input_filename,
                    output_filenames,
                    output_bamfiles,
                    input_filenames)

        if override_last_files:
            self.analysis.last_operation_filenames = output_filenames
            self.analysis.can_unlink = True

        if write_bam_files and len(output_bamfiles) != 0:
            self.analysis.bamfiles = output_bamfiles
