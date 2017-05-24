from . import utils
from .exceptions import PipelineError
from .barcoded_filename import BarcodedFilename

import logging
import os


class Analysis:
    def __init__(self, sample, root, config, parameters):
        self.sample = sample
        self.root = root
        if parameters["use_date"] is None:
            self.current = utils.get_current()
        else:
            self.current = parameters["use_date"]
        self.parameters = parameters
        self.basename = "%s.%s" % (self.sample, self.current)
        self.bam_dir = os.path.join(self.root, "BAM")
        self.out_dir = os.path.join(self.root, "Variants")
        self.bamfiles = None
        self.config = config
        self.last_operation_filenames = None
        self.run_fake = False

        if not os.path.exists(self.root):
            os.makedirs(self.root)

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        logs_dir = os.path.join(self.root, "logs")
        if not os.path.exists(logs_dir):
            os.makedirs(logs_dir)

        self.logger = logging.getLogger(self.basename)
        self.log_handler = logging.FileHandler(
            os.path.join(logs_dir, self.basename + ".steps.txt"))
        self.log_handler.setFormatter(
            logging.Formatter("%(asctime)s-15 %(message)s")
        )
        self.logger.addHandler(self.log_handler)
        self.logger.setLevel(level=logging.INFO)

    def _get_first_filename(self):
        filename = self.last_operation_filenames
        if filename is None:
            return None

        while not isinstance(filename, str):
            if isinstance(filename, list):
                if len(filename) > 0:
                    filename = filename[0]
                else:
                    return None
            elif isinstance(filename, dict):
                if len(filename) > 0:
                    filename = next(iter(filename.values()))
                else:
                    return None
            else:
                raise PipelineError(
                    "unexpected type for last_operation_filenames")

        if len(filename) > 0:
            return filename
        else:
            return None

    def _get_custom_dir(self, param):
        filename = self._get_first_filename()
        directory = getattr(self, param)
        if filename is None:
            return getattr(self, directory)

        return BarcodedFilename(filename).get_directory(directory)

    def get_bam_dir(self):
        return self._get_custom_dir("bam_dir")

    def get_out_dir(self):
        return self._get_custom_dir("out_dir")
