from . import utils
from .exceptions import PipelineError
from .barcoded_filename import BarcodedFilename

import os
import re


class Executor:

    def __init__(self, analysis):
        self.analysis = analysis

    def __call__(self,
                 command,
                 output_format=None,
                 input_filenames=None,
                 input_function=None,
                 input_split_reads=True,
                 output_path=None,
                 output_function=None,
                 error_string=None,
                 exception_string=None,
                 override_last_files=True,
                 write_bam_files=True,
                 unlink_inputs=False,
                 save_only_last=True,
                 use_normals=False,
                 split_by_organism=False,
                 only_human=False):

        re_replacer = re.compile(r"\{([^}]+)\}")

        if input_filenames is None:
            if self.analysis.last_operation_filenames is None:
                raise PipelineError("input files missing and "
                                    "last_operation_filenames empty")

            input_filenames = utils.get_sample_filenames(
                self.analysis.last_operation_filenames,
                split_by_organism)

        input_filenames = utils.get_sample_filenames(input_filenames,
                                                     split_by_organism)
        if not isinstance(input_filenames, dict):
            input_filenames = {"": input_filenames}

        if self.analysis.parameters["use_normals"] and use_normals:
            splitted_filenames = {}
            for organism, input_list in input_filenames.items():
                current_splitted = {}
                for filename in input_list:
                    barcoded_filename = BarcodedFilename(filename)
                    sample_id = barcoded_filename.sample
                    if sample_id not in current_splitted:
                        current_splitted[sample_id] = {}

                    sample_type = barcoded_filename.tissue
                    if sample_type.is_normal():
                        current_splitted[sample_id]["normal"] = filename
                    else:
                        current_splitted[sample_id]["tumor"] = filename

                splitted_filenames[organism] = current_splitted

            input_filenames = {}
            for splitted_key, splitted_list in splitted_filenames.items():
                for splitted_filename in splitted_list.values():
                    if len(splitted_filename) == 1:
                        if list(splitted_filename.keys())[0] == "tumor":
                            if splitted_key not in input_filenames:
                                input_filenames[splitted_key] = []
                            input_filenames[splitted_key].append(
                                list(splitted_filename.values())[0])
                    else:
                        if splitted_key not in input_filenames:
                            input_filenames[splitted_key] = []

                        input_filenames[splitted_key].append(splitted_filename)

        if input_function is not None:
            mod_input_filenames = {}

            if input_split_reads:
                for organism, current_input_filenames in input_filenames.items():
                    mod_input_filenames[organism] = \
                        [filename for filename
                         in map(input_function, current_input_filenames)
                         if filename is not None]
            else:
                for organism, current_input_filenames in input_filenames.items():
                    current_mod_input = input_function(current_input_filenames)
                    if isinstance(current_mod_input, list):
                        current_mod_input.remove(None)
                        current_mod_input.remove("")
                    mod_input_filenames[organism] = current_mod_input

            if mod_input_filenames is None or \
                    len(mod_input_filenames) == 0:
                raise PipelineError("empty input list")

            first_mod_input_filenames = list(mod_input_filenames.values())[0]
            if (isinstance(first_mod_input_filenames, str) and first_mod_input_filenames == "") or\
                    (isinstance(first_mod_input_filenames, list) and
                     len(first_mod_input_filenames) == 0):
                raise PipelineError("empty input list")

            if not input_split_reads:
                for organism, current_list in input_filenames.items():
                    input_filenames[organism] = [current_list]
                for organism, current_list in mod_input_filenames.items():
                    mod_input_filenames[organism] = [current_list]

        else:
            if input_filenames is None or \
                    len(input_filenames) == 0:
                raise PipelineError("empty input list")

            first_input_filenames = list(input_filenames.values())[0]
            if (isinstance(first_input_filenames, str) and first_input_filenames == "") or \
                    (isinstance(first_input_filenames, list) and
                     len(first_input_filenames) == 0):
                raise PipelineError("empty input list")

            if not input_split_reads:
                for organism, current_list in input_filenames.items():
                    input_filenames[organism] = [current_list]

        output_filenames = {}
        output_bamfiles = {}
        for current_organism in input_filenames.keys():
            if only_human and (current_organism != "" and
                               not current_organism.startswith("hg")):
                continue

            current_input_filenames = input_filenames[current_organism]

            if input_function is not None:
                current_mod_input_filenames = mod_input_filenames[
                    current_organism]
                if len(current_mod_input_filenames) == len(current_input_filenames):
                    iterator = zip(
                        current_input_filenames,
                        current_mod_input_filenames)
                else:
                    iterator = zip(
                        current_mod_input_filenames,
                        current_mod_input_filenames)
            else:
                iterator = zip(
                    current_input_filenames,
                    [None] *
                    len(current_input_filenames))

            for input_filename, mod_input_filename in iterator:
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
                            barcoded_filename = BarcodedFilename(
                                filename)
                            organism = barcoded_filename.organism
                            index = barcoded_filename.read_index
                        else:
                            barcoded_filename = BarcodedFilename(
                                filename["tumor"])
                            organism = barcoded_filename.organism
                            index = barcoded_filename.read_index
                        read_index.append(index)

                if only_human and (organism != "" and organism is not None and
                                   not organism.startswith("hg")):
                    continue

                if input_function is not None:
                    real_input_filename = input_filename
                    input_filename = mod_input_filename

                organism_str = None
                if organism is None or organism == "":
                    organism_str = ""
                    organism = "hg19"
                else:
                    organism_str = "." + organism

                try:
                    genome_ref, genome_index = \
                        utils.get_genome_ref_index_by_organism(
                            self.analysis.config, organism)
                except:
                    pass

                try:
                    dbsnp = utils.get_dbsnp_by_organism(
                        self.analysis.config, organism)
                except:
                    pass

                try:
                    cosmic = utils.get_cosmic_by_organism(
                        self.analysis.config, organism)
                except:
                    pass

                if output_format is not None:
                    if not isinstance(output_format, list):
                        output_format = [output_format]

                    if output_path is not None:
                        output_format = [os.path.join(output_path, s)
                                         for s in output_format]

                    output_filename = []
                    for s in output_format:
                        for match in re_replacer.finditer(s):
                            try:
                                s = re.sub(
                                    match.group(0), str(
                                        eval(
                                            match.group(1))), s)
                            except:
                                raise PipelineError(
                                    "cannot replace parameter %s"
                                    % (match.group(0)))
                        output_filename.append(s)

                    if output_function is not None:
                        output_filename = [filename for filename
                                           in map(output_function,
                                                  output_filename)]

                    if len(output_filename) == 1:
                        output_filename = output_filename[0]
                else:
                    output_filename = None

                commands = []
                if isinstance(command, list):
                    for s in command:
                        if isinstance(s, str):
                            for match in re_replacer.finditer(s):
                                try:
                                    s = s.replace(
                                        match.group(0), str(
                                            eval(
                                                match.group(1))))
                                except:
                                    raise PipelineError(
                                        "cannot replace parameter %s"
                                        % (match.group(0)))

                        commands.append(s)
                else:
                    current_command = command
                    if isinstance(current_command, str):
                        for match in re_replacer.finditer(command):
                            try:
                                current_command = current_command.replace(
                                    match.group(0), str(
                                        eval(match.group(1))))
                            except:
                                raise PipelineError(
                                    "cannot replace parameter %s"
                                    % (match.group(0)))
                    commands.append(current_command)

                for command_index, current_command in enumerate(commands):
                    if isinstance(current_command, str):
                        if not self.analysis.run_fake:
                            status = utils.run_and_log(
                                current_command,
                                self.analysis.logger)
                        else:
                            self.analysis.logger.info(
                                "Faking command '%s'", current_command)
                            status = 0
                    else:
                        if not self.analysis.run_fake:
                            local_vars = locals()
                            del local_vars["self"]
                            current_command(**local_vars)
                        else:
                            self.analysis.logger.info(
                                "Faking lambda")
                        status = 0

                    if status != 0:
                        arg_zero = os.path.basename(command.split(" ")[0])

                        if error_string is None:
                            error_string = "%s exited with status %d" % (
                                arg_zero,
                                status)

                        if exception_string is None:
                            exception_string = "%s error" % arg_zero

                        for match in re_replacer.finditer(error_string):
                            try:
                                error_string = re.sub(
                                    match.group(0),
                                    str(eval(match.group(1))),
                                    error_string)
                            except:
                                raise PipelineError(
                                    "cannot replace parameter %s"
                                    % (match.group(0)))

                        for match in re_replacer.finditer(exception_string):
                            try:
                                exception_string = re.sub(
                                    match.group(0), str(
                                        eval(
                                            match.group(1))), exception_string)
                            except:
                                raise PipelineError(
                                    "cannot replace parameter %s"
                                    % (match.group(0)))
                        self.analysis.logger.error(error_string)
                        raise PipelineError(exception_string)

                    if output_filename:
                        if not isinstance(output_filename, list):
                            output_filename = [output_filename]

                        new_filename = []
                        for filename in output_filename:
                            dirname = os.path.dirname(filename)
                            if dirname == "" or dirname == ".":
                                new_filename.append(
                                    os.path.join(
                                        os.getcwd(),
                                        filename))
                            else:
                                new_filename.append(filename)
                        output_filename = new_filename

                        if save_only_last:
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
                                output_filenames[organism] += output_filename
                            else:
                                output_filenames[organism] = output_filename

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

                if unlink_inputs:
                    if input_function is not None:
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

        if override_last_files:
            self.analysis.last_operation_filenames = output_filenames

        if write_bam_files and len(output_bamfiles) != 0:
            self.analysis.bamfiles = output_bamfiles
