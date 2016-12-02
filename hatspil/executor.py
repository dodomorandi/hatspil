from . import utils
from .exceptions import PipelineError

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
                 save_only_last=True):

        re_replacer = re.compile(r"\{([^}]+)\}")

        if input_filenames is None:
            if self.analysis.last_operation_filenames is None:
                raise PipelineError("input files missing and "
                                    "last_operation_filenames empty")

            input_filenames = utils.get_sample_filenames(
                self.analysis.last_operation_filenames)

        input_filenames = utils.get_sample_filenames(input_filenames)

        if input_function is not None:
            if input_split_reads:
                mod_input_filenames = [filename for filename
                                       in map(input_function, input_filenames)
                                       if filename is not None]
            else:
                mod_input_filenames = input_function(input_filenames)
                if isinstance(mod_input_filenames, list):
                    mod_input_filenames.remove(None)
                    mod_input_filenames.remove("")

            if mod_input_filenames is None or \
                    (isinstance(mod_input_filenames, str) and mod_input_filenames == "") or\
                    (isinstance(mod_input_filenames, list) and len(mod_input_filenames) == 0):
                raise PipelineError("empty input list")

            if not input_split_reads:
                input_filenames = [input_filenames]
                mod_input_filenames = [mod_input_filenames]

            if len(mod_input_filenames) == len(input_filenames):
                iterator = zip(input_filenames, mod_input_filenames)
            else:
                iterator = zip(mod_input_filenames, mod_input_filenames)
        else:
            if input_filenames is None or \
                    (isinstance(input_filenames, str) and input_filenames == "") or\
                    (isinstance(input_filenames, list) and len(input_filenames) == 0):
                raise PipelineError("empty input list")

            if not input_split_reads:
                input_filenames = [input_filenames]

            iterator = zip(input_filenames, [None] * len(input_filenames))

        output_filenames = {}
        output_bamfiles = {}
        for input_filename, mod_input_filename in iterator:
            if isinstance(input_filename, str):
                organism, read_index, _ = utils.get_params_from_filename(
                    input_filename, self.analysis)
            else:
                read_index = []
                for filename in input_filename:
                    organism, index, _ = utils.get_params_from_filename(
                        filename, self.analysis)
                    read_index.append(index)

            if input_function is not None:
                real_input_filename = input_filename
                input_filename = mod_input_filename

            organism_str = None
            if organism is None or organism == "":
                organism_str = ""
                organism = "hg19"
            else:
                organism_str = "_" + organism

            genome_ref, genome_index = \
                utils.get_genome_ref_index_by_organism(
                    self.analysis.config, organism)

            dbsnp = utils.get_dbsnp_by_organism(
                self.analysis.config, organism)

            cosmic = utils.get_cosmic_by_organism(
                self.analysis.config, organism)

            if output_format is not None:
                if not isinstance(output_format, list):
                    output_format = [output_format]

                if output_path is not None:
                    output_format = [os.path.join(output_path, s)
                                     for s in output_format]

                output_filename = []
                for s in output_format:
                    for match in re_replacer.finditer(s):
                        s = re.sub(match.group(0), str(eval(match.group(1))), s)
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
                            s = s.replace(match.group(0), str(eval(match.group(1))))
                    commands.append(s)
            else:
                current_command = command
                if isinstance(current_command, str):
                    for match in re_replacer.finditer(command):
                        current_command = current_command.replace(match.group(0), str(eval(match.group(1))))
                commands.append(current_command)

            for command_index, current_command in enumerate(commands):
                if isinstance(current_command, str):
                    status = utils.run_and_log(current_command, self.analysis.logger)
                else:
                    local_vars = locals()
                    del local_vars["self"]
                    current_command(**local_vars)
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
                        error_string = re.sub(match.group(0),
                                              str(eval(match.group(1))),
                                              error_string)
                    for match in re_replacer.finditer(exception_string):
                        exception_string = re.sub(match.group(0),
                                                  str(eval(match.group(1))),
                                                  exception_string)
                    self.analysis.logger.error(error_string)
                    raise PipelineError(exception_string)

                if output_filename:
                    if not isinstance(output_filename, list):
                        output_filename = [output_filename]

                    new_filename = []
                    for filename in output_filename:
                        dirname = os.path.dirname(filename)
                        if dirname == "" or dirname == ".":
                            new_filename.append(os.path.join(os.getcwd(), filename))
                        else:
                            new_filename.append(filename)
                    output_filename = new_filename

                    if save_only_last:
                        if command_index == len(commands) - 1:
                            if organism in output_filenames:
                                output_filenames[organism] += output_filename
                            else:
                                output_filenames[organism] = output_filename

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
                            output_extension = os.path.splitext(filename)[1]
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
