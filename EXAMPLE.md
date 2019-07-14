# A run with HaTSPiL

Running a NGS data analysis pipeline is not trivial. There are many steps involved, and many of these do not depend on a software like HaTSPiL.
Here a brief overview that shows how to run HaTSPiL with real data.

## Read the README!

There are many software needed to run HaTSPiL. Check the readme to check which software you need to install. Obviously, install HaSTPiL as well.
For this specific example, BWA will be used as aligner. Be sure to have an available version on your system. Other software like Xenome, Disambiguate or STAR are not needed because we are going to process whole exome NGS data.
You will also need some additional files: the fasta file of the human genome, the dbSNP data in VCF format and the Cosmic DB in VCF format. For this example we are going to use the hg19/hg37 version of the genome.

## Prepare the environment

For simplicity, we are going to create a directory inside the home of the user (called 'user') and we are going to work inside this directory:
```sh
mkdir hatspil_workspace
cd hatspil_workspace
```

## Create a config for HaTSPiL
```sh
hatspil --configout config.ini
```
Now you need to edit most of the parameters in order to make everything work. Do not worry about the kits section, we will talk about it later.

## Get the data

We are going to analyze two whole exome sequencing from SRA. The two samples are:

* SRR1201243
* SRR1201245

In order to download them you need the [SRA tools](http://ncbi.github.io/sra-tools/install_config.html). One you have them installed, you can do the following:

```sh
mkdir fastq
cd fastq
fast-dump --split-files SRR1201243
fast-dump --split-files SRR1201245
cd ..
```

It will take some time to download.

## Get data from Agilent

The design of the experiment for the two dataset says that the DNA capture was performed using the 'Agilent SureSelectXT HumanAllExon 51Mb' kit.
You need to get to the [Agilent SureDesign](https://earray.chem.agilent.com/suredesign/index.htm) website and go download the 'SureSelect Human All Exon V4' zip file from the catalog (the name does not match, but this is the only 51Mb sized kit, so it should be the right one).
Let's now create an 'agilent' directory inside our workspace and unzip the file right there:

```sh
mkdir agilent
cd agilent
unzip ~/Downloads/S03723314.zip
```
Assuming that you have the file in your 'Downloads' directory.

## Transform the data from Agilent

We need to transform these files for Picard. The 'Regions' file is used as 'targets', the 'Padded' file as 'baits'.
Assuming that you have a file called `picard.jar` inside your `/opt/jars` and the sequence of the human genome stored as `/opt/db/ucsc.hg19.fasta`, the following is needed:
```sh
java -jar /opt/jars/picard.jar CreateSequenceDictionary R=/opt/db/ucsc.hg19.fasta O=/opt/db/ucsc.hg19.dict
java -jar /opt/jars/picard.jar BedToIntervalList I=S03723314_Regions.bed O=S03723314_targets.interval_list SD=/opt/db/ucsc.hg19.dict
java -jar /opt/jars/picard.jar BedToIntervalList I=S03723314_Padded.bed O=S03723314_baits.interval_list SD=/opt/db/ucsc.hg19.dict
cd ..
```

The files that have been created can be used with Picard inside HaTSPiL.

## Fix the config.ini

Now it is possible to correctly fill the section of the whole exome kit. It should look like the following:
```ini
[KIT 0 WHOLE_EXOME]
name = SRA_Kit
target_list = /home/user/hatspil_workspace/agilent/S03723314_targets.interval_list 
bait_list = /home/user/hatspil_workspace/agilent/S03723314_baits.interval_list 
indels_hg19 = /opt/db/1000G_phase1.indels.hg19.sites.vcf
cancer_site = rhabdomyosarcoma
amplicons =
adapter_r1 =
adapter_r2 =
mean_len_library =
sd_len_library =
```

The cancer site is set accordingly to the type of the original tumor. The adapter section has been left empty because they are not known. The library length section is left empty because it is only needed when using NovoAlign.

## Barcode the fastq files

As explained in the README, the fastq must be correctly barcoded in order to let HaTSPiL work correctly.

```sh
cd fastq
mv SRR1201243_1.fastq test-rm001-01-000-000.R1.fastq
mv SRR1201243_2.fastq test-rm001-01-000-000.R2.fastq

mv SRR1201245_1.fastq test-rm002-50-000-000.R1.fastq
mv SRR1201245_2.fastq test-rm002-50-000-000.R2.fastq
cd ..
```

## Run HaTSPiL

If everything is fine, HaTSPiL should work flawlessly. In this example the '--no-cutadapt' option is used because we do not know which adapters have been used, so we can skip Cutadapt.

```sh
hatspil --root-dir `pwd` --fastq-dir `pwd`/fastq --scan-samples --no-cutadapt
```
