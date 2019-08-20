# Welcome to HaTSPiL!

HaTSPiL is a a Python pipeline for High-Throughput Sequencing analysis. It has been designed to be used inside our laboratory, the Salvatore Oliviero lab within the HuGef institute, TO (IT). Whoever find it useful or a good starting point to develop his own pipeline is encouraged to use and hack the code.

As for any Python package, use the *setup.py* script to build and install the package inside your system. Python >= 3.5 is required. 

One of the most important file is the config.ini. It can be obtained passing only '--configout config.ini' to the program, and it is recommended to modify the resulting file. The program will scan for a file with this name in the current directory, otherwise it must be specified explicitly with the "-c" parameter. 

## Required software
* JVM >= 6
* JVM 7 (yes, some software needs this version to work correctly)
* Perl 5
* Samtools
* FastQC
* SeqTK 
* Picard
* Varscan
* Genome Analysis ToolKit (GATK)
* Mutext
* Bam2TDF
* NovoAlign or BWA
* STAR for RNA-Seq analaysis (partial support)
* Xenome or Disambiguate for xenograft classification

## Installation

Download from Git repo:
```bash
git clone https://github.com/dodomorandi/hatspil.git
cd hatspil
```

Install for the current user (recommended)...

```bash
python3 ./setup.py install --user
```

... or for anyone.

```bash
sudo python3 ./setup.py install
```

## Additional dependencies

In order to use a MongoDB to store the results of the analyses, it is necessary to install the software and then install the python module:

```bash
pip3 install --user pymongo
```

## Create your configuration

HaTSPiL requires a configuration file to know how to work.
To start with a clean template you can use the following:

```bash
hatspil --configout hatspil.ini
```

This creates a configuration that can (and should) be edited to set the correct parameters.

The first section requires the path to the executables. If these are available from the current PATH environment variable, you can leave them with the default values.

The second section needs the path of the _.jar_ files. Generally they are in specific folders, therefore it is necessary to set the correct path for all of them.

The third section includes the path of the files used in the analysis. Some of them are mandatory, but others can be omitted. For instance, if you only need _hg38_ and _mm10_ annotations, you can remove the fields related to _hg19_ and _mm9_. Here a brief explanation of what is expected for each parameter:

* strelka\_basedir: where the file of Strelka can be found
* strelka\_config: the _config.ini_ file for Strelka
* {genome}\_ref: the genome fasta file
* {genome}\_index: the NovoAlign index file
* cosmic\_{genome}: the Cosmic VCF file 
* dbsnp\_{genome}: The DbSnp VCF file
* annovar\_basedir: The path of ANNOVAR

The fourth section contains the configuration for some additional parameters:

* xenome\_index: the common part of the index files for Xenome
* xenome\_threads: the number of threads for Xenome
* strelka\_threads: the number of threads for Xenome
* use\_{genome}: whether a particular annotation can be used
* mails: a comma-separated list of emails to send a notification at the end of the analysis
* use\_mongodb: whether the MongoDB can be used
* picard\_jvm\_args: the arguments to pass to the JVM when running Picard
* varscan\_jvm\_args: the arguments to pass to the JVM when running VarScan
* gatk\_jvm\_args: the arguments to pass to the JVM when running GATK
* mutect\_jvm\_args: the arguments to pass to the JVM when running MuTect
* mutect\_args: the arguments that must be passed to MuTect by default (others are appended)

The fifth section represent a kit section. Multiple kit sections can be specified in the configuration. The header of the section must be formatted using the keyword 'KIT' followed by the index of the kit and the type of the analyte. The latter can be `WHOLE_EXOME`, `GENE_PANEL`, `FUSION_PANEL` or `RNASEQ`.
The fields for each kit section are the following:

* target\_list: the list of the target genes, in case it is meaninful for the analysis. It is used by Picard.
* bait\_list: the list of baits, in case it is meaninful for the analysis. It is used by Picard.
* indels\_{genome}: a comma delimited list of VCF files used during indel realignment. Used by GATK.
* amplicons: a BED file containing the amplicons. Used during the variant calling.
* name: the name representing the current kit.
* cancer\_site: the cancer site to filter for during variant calling. Possible values are: 'esophageal', 'head\_and\_neck', 'leukemia', 'adrenocortical\_carcinoma', 'bladder', 'breast', 'prostate', 'soft\_tissue\_sarcoma', 'colorectal', 'chondroblastoma', 'angiosarcoma', 'melanoma', 'oligodendroglioma', 'pancreas', 'craniopharyngioma', 'Ewing\_sarcoma', 'glioma', 'kidney', 'glioblastoma', 'thyroid', 'ovarian', 'lung', 'endometrial', 'gastric', 'adrenocortical\_adenoma', 'cholangiocarcinoma', 'myeloma', 'meningioma', 'neuroblastoma', 'astrocytoma', 'small\_intestine', 'liver', 'multiple', 'myelodysplasia', 'lymphoma', 'medulloblastoma', 'myeloproliferative\_neoplasm', 'rhabdomyosarcoma', 'rhabdoid', 'nasopharyngeal', 'gallbladder', 'leiomyoma', 'chondromyxoid\_fibroma', 'cervix', 'thymus', 'adrenocortical\_carcinoma ', 'renal\_cell\_carcinoma', 'ameloblastoma', 'chondrosarcoma', 'intracranial\_germ\_cell', 'polymorphous\_low-grade\_adenocarcinoma'.
* adapter\_r[1|2]: the forward and reverse adapters of the analysis. Used by Cutadapt.
* [mean|sd]\_len\_library: the mean and the standard deviation of the len of the library. Used by NovoAlign.

The sixth section is related to the MongoDB configuration. The fields are self-explanatory and they are needed to specify the name of the database, the host, the post, the username and the password for the connection to the database.

Modify the configuration file according to your needs. When this step is done, you are almost ready to use HaTSPiL!

## Prepare your fastq files

HaTSPiL uses a barcode system to understand what a file represents and how to handle it. In order to run the software, it is necessary to give an appropriate name to the fastq files that are going to be analysed.

### How barcoding works

First of all, let's see how a barcoded filename is composed using an example:
```
lung-lg003-06-021-201
```
In order to understand the meaning of the barcode, it is necessary to split it into pieces:
```
+-------++------+-------+----+---+---+---+---+---+---+
| part  || lung | lg003 | 06 | 0 | 2 | 1 | 2 | 0 | 1 |
+-------++------+-------+----+---+---+---+---+---+---+
| index ||   1  |   2   | 3  | 4 | 5 | 6 | 7 | 8 | 9 |
+-------++------+-------+----+---+---+---+---+---+---+
```
* Part 1 (lung): the name of the project
* Part 2 (lg003): the pseudo-identifier for the patient/sample
* Part 3 (06): the tissue
* Part 4 (0): The molecule
* Part 5 (1): The experiment type (analyte)
* Part 6 (2): The kit
* Part 7 (2): The biopsy
* Part 8 (0): The sample
* Part 9 (1): The sequencing

The project and the patient/sample pseudo-identifier are alphanumeric identifier chosen by the user. They are mainly useful to organize the experiments.
The tissue value is a identifier, taken from [TCGA](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes). In this example, a "recurrent solid tumor" is used. Common values are _01_ for primary solid tumor, _10_ for blood derived normal and _11_ for solid tissue normal. More information can be found inside the `core.barcoded_filename` module.
The molecule identifies whether the library preparation was performed from DNA (0) or RNA (1). 
The kit, the biopsy, the sample and the sequencing are 0-based indices to identify a tree of characteristics for the sample. For instance, in this case it is the third biopsy from the same patient, it is the first sample taken from that biopsy, and it is the second sequencing (maybe the first one had a problem during the library preparation?).
The experiment type (analyte) identifies which type of analysis was performed on the sample. Possible values for now are whole exome sequencing (0), gene panel sequencing (1), fusion panel sequencing (2, unsupported for now) and RNA-seq sequencing (3, partially supported). Keep in mind that the analyte is very related to the index of the kit, because for each combination of analyte and kit it is possible to specify different behaviours using the configuration file.

### Barcode your fastq files

Let's suppose you are working on a project in which you study different types of carcinoma. This time you collected a biopsy of a lung carcinoma from a new patient. Internally your group uses a random sequence of letters to identify a patient. Moreover, this is the first biopsy from the patient. The biopsy looked quite dishomogeneous, and you decided to divide it in three samples. After some preliminary analysis, you decide to sequence only the third sample. You go for a whole exome sequencing, but unfortunately the first time you sequenced something strange happened with the indices, so you needed to sequence the sample a second time.

Here a brief overview:

* Your project is called 'carc', it stands for _carcinoma_.
* The patient is lgHRTPAF, 'lg' stands for _lung carcinoma_ and the rest is the random sequence of letters that identify the patient.
* The tissue is a primary solid tumor, coded as _01_.
* The molecule is DNA, coded as _0_.
* You have a kit to perform whole exome sequencing, it is always the same. _0_ is perfect.
* It is the first biopsy from the patient, coded as _0_ because indices are 0-based.
* The sample is the third, coded as _2_.
* It is the second sequencing, coded as _1_.

It is now possible to compose the barcode for the current analysis:
```
carc-lgHRTPAF-01-000-021
```
The sequencing was performed using paired-end technology, therefore you have got two files, one for the forward reads, on for the reverse reads. The name of the two files will be:
```
carc-lgHRTPAF-01-000-021.R1.fastq
carc-lgHRTPAF-01-000-021.R2.fastq
```

These two files can be analysed by HaTSPiL. If you have many fastq files to analyse, you can name all of them in order to perform one single run.

## Running HaTSPiL

First of all, it is necessary to create a workspace directory for HaTSPiL. For this example, the workspace will be in `/data/hatspil_workspace`. The directory containing the fastq files can be anywhere, but it is good to place it inside the hatspil workspace (other fastq files can be generated during the process). In this case, a `fastq` directory can be created inside `/data/hatspil_workspace`.

Now it is possible to run HaTSPiL in a very simple way:
```bash
hatspil --root-dir /data/hatspil_workspace --fastq-dir /data/hatspil_workspace/fastq --config hatspil.ini --scan-samples
```

If everything is fine, HaTSPiL with start analysing all the fastq files that have a valid barcode inside the fastq directory.

It is also possible to create a file with a list of partial barcodes that are intended to be analysed. For instance, it is possible to create the file 'hatspil\_carc.txt' with the following content:
```
carc
```
Then hatspil can be run using `--list-file hatspil_carc.txt` instead of `--scan-samples`.

It is possible to specify a set of whitelisted barcodes, one per line, and it is possible to use the placeholder `*` to consider valid any possible value for a specific field of the barcode. For instance:
```
carc-*-01-*21
lymp-*-10-***-**3
```
With this list, HaTSPiL only analyse the solid tumor from the project _carc_ obtained from a gene panel analysis with the third kit, and all the blood normal samples from the project _lymp_ coming from the fourth sequencing.

For more command line options, it is suggested to run `hatspil --help`.

## The results

Inside the root directory the directory 'reports' is created, containing HTML reports for the current analysis and a global report for all the analysis stored in the database.
If specific information are needed, there are many 'REPORTS' directories placed inside the fastq, the `BAM` and the `Variants` folders, which contain many different type of fine-grained information of the whole analysis.

# Customization of HaTSPiL
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
