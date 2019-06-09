"""A small set of custom exceptions."""


class PipelineError(RuntimeError):
    """An exception generated during pipeline execution.
    
    Generally this exception is raised within `Executor` in case of
    errors.
    """


class BarcodeError(RuntimeError):
    """An exception generated in case of invalid barcoded filenames."""


class DataError(RuntimeError):
    """A generic exception for invalid data."""


class AnnotationError(RuntimeError):
    """An exception generated in case of invalid genomic annotation."""
