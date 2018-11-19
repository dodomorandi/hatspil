# Welcome to HaTSPiL!

HaTSPiL is a a Python pipeline for High-Throughput Sequencing analysis. It has been designed to be used inside our laboratory, the Salvatore Oliviero lab within the HuGef institute, TO (IT). Whoever find it useful or a good starting point to develop his own pipeline is encouraged to use and hack the code.

As for any Python package, use the *setup.py* script to build and install the package inside your system. Python >= 3.5 is required. 

One of the most important file is the config.ini. It can be obtained passing only '--configout config.ini' to the program, and it is recommended to modify the resulting file. The program will scan for a file with this name in the current directory, otherwise it must be specified explicitly with the "-c" parameter. 

## Customization of HaTSPiL
The modular structure of HaTSPiL provides a high level of customizability, but different parts of the software have different complexities. The introduction of new workflow steps is the most appealing and easy customization point, and we briefly describe some small examples. A minimum knowledge of the Python programming language is necessary to extend the features of HaTSPiL.

Suppose that a tool called *parse_bam* has to be integrated in the workflow. This tool expects a BAM file and outputs a SAM file, and it performs a filtering of the data in order to improve the information of the alignment in some way. First of all, it is necessary to understand in which phase *parse_bam* must be run, because of the input filename format. Looking at the *runner.py* module it is possible to see the following at some point:

```python
mapping = Mapping(analysis, self.fastq_dir)
mapping.run()

if (
    not self.parameters["use_normals"]
    and not self.parameters["only_mapping"]
    and barcoded_filename.analyte != Analyte.RNASEQ
):
    self._run_mutation_analysis(analysis, False)
```

Without worrying too much about details, it is possible to notice that, intuitively, a `Mapping` object is created and the mapping is run. After that, when certain conditions are met, the `_run_mutation_analysis` is called. *parse_bam* execution has to be handled after the mapping process (because it produces BAM files from FastQ files) and before the mutation analysis (the analysis has to be performed on files filtered with the custom tool). To follow the structure of HaTSPiL, it is a good practice to create a different module named *parse_bam.py* that contains the class `ParseBam`, can be initialized passing the `Analysis` object (as for the other used in the runner module) and that can be used through the `run()` method. With this idea, the code of runner can already be changed:

```python
# Start of the running.py file, with lots of imports
# ...

from parse_bam import ParseBam  # In order to import the new class

# ...

mapping = Mapping(analysis, self.fastq_dir)
mapping.run()

parse_bam = ParseBam(analysis)  # A new instance of `ParseBam` is created
parse_bam.run()                 # The `.run` method is called 

if (
    not self.parameters["use_normals"]
    and not self.parameters["only_mapping"]
    and barcoded_filename.analyte != Analyte.RNASEQ
):
    self._run_mutation_analysis(analysis, False)
```

The new module has to be written, in order to correctly perform the execution of the *parse_bam* tool. A new *parse_bam.py* is created inside the hatspil directory.

```python
from .analysis import Analysis


class ParseBam:
    def __init__(self, analysis: Analysis) -> None:  # This is the constructor of the class
        self.analysis = analysis                     # The analysis is stored inside the object

    def parse_bam(self) -> None:                     # The parse_bam method, we need to write it
        #  TODO

    def run(self) -> None:                           # The `run` helper function, in this case a
        self.parse_bam()                             # simple wrapper for `parse_bam`
```

This is the basic structure of the new module. It is worth noting that the `run` method, in this case, simply calls the `parse_bam` method. The reason behind this is consistency with respect to the other modules.
Now it is shown how the `parse_bam` method can be written using the Executor.

```python
from .analysis import Analysis
from .executor import Executor


class ParseBam:
    # ...
    def parse_bam(self) -> None:
        bam_dir = self.analysis.get_bam_dir()  # We get the BAM files directory thanks to a helper function
        executor = Executor(self.analysis)     # A new `Executor` is instantiated using the current `analysis`
        executor(
            "/usr/bin/parse_bam -i {input_filename} "
            "-o {output_filename}",  # This is the command line for the `parse_bam` tool. `{input_filename}`
                                     # and `{output_filename}` are placeholder for the `Executor`, which
                                     # replaces them with the correct filenames using all the information
                                     # it has got
            output_format=self.analysis.sample + 
            "parse_bam{organism_str}.sam",  # This is the format for the output filename. Even in this case,
                                            # the `Executor` is able to replace placeholders with the
                                            # appropriate values
            output_path=bam_dir,  # It is possible to specify where to place the output files, in this case
                                  # the BAM files directory
            unlink_inputs=True,  # We do not need the input files, they can be removed
            only_human=False,  # If supplied, this module can also handle BAM files from organisms that are
                               # not human
        )
    # ...
```

The code shows that it is pretty easy to set an Executor in order to run the desired command. It is worth noting that it is possible to use special parameters surrounded by curly braces. These are replaced internally by the Executor depending on all the arguments that are set and the input filenames. There are many parameters that can be passed to `executor`, see the documentation for further details. A critical point of this example is the fact that the *parse_bam* executable is run using a hard coded path. This can be improved customizing the config module in order to let other users specify the location of the program for their needs. The following is the code of the modified *config.py* module.

```python
# ...
class Config:
    executables = ("java", "java7", "perl", "seqtk", "fastqc", "samtools", "parse_bam")

# ...
    def __init__(self, filename: Optional[str] = None) -> None:
        # ....
        self.parse_bam = "parse_bam"
```

With this small addition, *parse_bam* becomes a mandatory executable for HaTSPiL, and its availability is checked during the initialization phase of the software. Inside the `__init__` function the default name for the executable is specified, and the path can be omitted as soon as the tool can be `run` using the PATH environment variable.
Thanks to this improvement, the *parse_bam.py* module has to be changed in order to use the new feature:

```python
# ...
        config = self.analysis.config
        executor(
            f"{config.parse_bam} -i {{input_filename}} "
            "-o {output_filename}",
        # ...
```

In this case, the Python string formatter is used in order to place the value of `config.parse_bam` inside the command. The `input_filename` now needs a double pair of curly braces, in order to avoid an immediate substitution of a local variable (which does not exist in this case) and to obtain the same behavior as before.
A small improvement that can be done is providing logging facilities. This is pretty easy, because a Python standard Logger is stored inside the analysis instance and it is already configured. With the following, the module includes a simple information for the user:

```python
# ...
    def parse_bam(self) -> None:
    self.analysis.logger.info("Starting parse_bam")
    # ...
    self.analysis.logger.info("Ended parse_bam")
```
