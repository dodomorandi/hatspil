import logging
import os
from typing import Any, Dict, List, Optional, Union

from . import utils
from .barcoded_filename import BarcodedFilename
from .config import Config
from .exceptions import PipelineError


class Analysis:
    def __init__(
        self, sample: str, root: str, config: Config, parameters: Dict[str, Any]
    ) -> None:
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
        self.bamfiles: Dict[str, List[str]] = {}
        self.config = config
        self.last_operation_filenames: Union[
            str, List[str], Dict[str, List[str]], None
        ] = None
        self.run_fake = False
        self.can_unlink = True

        os.makedirs(self.root, exist_ok=True)
        os.makedirs(self.out_dir, exist_ok=True)

        logs_dir = os.path.join(self.root, "logs")

        os.makedirs(logs_dir, exist_ok=True)

        self.logger = logging.getLogger(self.basename)
        self.log_handler = logging.FileHandler(
            os.path.join(logs_dir, self.basename + ".steps.txt")
        )
        self.log_handler.setFormatter(logging.Formatter("%(asctime)s-15 %(message)s"))
        self.logger.addHandler(self.log_handler)
        self.logger.setLevel(logging.INFO)

    def _get_first_filename(self) -> Optional[str]:
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
                raise PipelineError("unexpected type for last_operation_filenames")

        if len(filename) > 0:
            return filename
        else:
            return None

    def _get_custom_dir(self, param: str) -> str:
        filename = self._get_first_filename()
        directory: str = getattr(self, param)
        if filename is None:
            return directory

        return BarcodedFilename(filename).get_directory(directory)

    def get_bam_dir(self) -> str:
        return self._get_custom_dir("bam_dir")

    def get_out_dir(self) -> str:
        return self._get_custom_dir("out_dir")

    @property
    def using_normals(self) -> bool:
        return (
            self.parameters["use_normals"]
            and self.last_operation_filenames is not None
            and isinstance(self.last_operation_filenames, dict)
            and "control" in self.last_operation_filenames
            and isinstance(self.last_operation_filenames["control"], list)
            and len(self.last_operation_filenames["control"]) > 0
        )
