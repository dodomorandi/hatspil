from . import utils

import logging
import os


class Analysis:
    def __init__(self, sample, root, config, parameters):
        self.sample = sample
        self.root = root
        self.current = utils.get_current()
        self.parameters = parameters
        self.basename = "%s.%s" % (self.sample, self.current)
        self.bam_dir = os.path.join(self.root, "BAM", "Panel")
        self.out_dir = os.path.join(self.root, "Variants", "Panel")
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
