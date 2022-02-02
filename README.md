TESSA - Transposons Exploration with Salmon's Selective Alignment
====
### Overview

TESSA runs Salmon with Selective Alignment using pre-built indices based on GENCODE sequences and RepeatMasker references for human and mouse.

### Requirements
----
TESSA requires Salmon to be installed (tested with v1.4.0) and at least 35Gb of disk space.

### Installation
----
Clone the repository to the desired directory with the following command.
```sh
git clone https://github.com/FemeniasM/TESSA
```
Click on this [link](https://mega.nz/folder/bk8U3CIB#wDBWeYMdrSO3zE3ZTF3C9Q) and download the `install` folder into the TESSA repository.

Move to the `TESSA` directory and run the installation script (`install.sh`). This script creates the indexes for the data sets and assumes that salmon is in your `$PATH`. User can change the path to the Salmon binary file with the `-s` argument as shown below. In addition, the user can define the number of threads (`-p` argument) and the k-mer size (`-k` argument) for the construction of the indexes.

```sh
cd TESSA/
bash install.sh -s path/to/salmon/binary/file -k 31 -p 12
```
The indexes building can take several minutes, if the installation was successful within the `TESSA` directory, an index sub-directory is created with the indexes for each species.

### How to use
----

The flags for the main script are described below 
```

Usage:  
bash TESSA.sh [flags]

Flags:
   -p threads [N] (1 default)     Number of therads
   -k kmer [N] (31 default)       k-mer size
   -s salmon binary path          Path to Salmon binary file (salmon default)
   -e library format ['pe']['se'] Indicates the format of the libraries: 'pe' for paired-end and 'se'
                                  for single-end reads. Supported extensions are: <.fq> or <.fastq>
                                  or <.fq.gz> or <.fastq.gz>. Paired-end file names should contain 
                                  _R1 _R2. Example: sample_R1.fq.gz, sample_R2.fq.gz (mandatory)
   -l folder with fastq files     Path to folder with fastq files (mandatory)
   -o output directory path       Path to output directory (mandatory)
   -i salmon index ['mm']['hs']   Salmon indices for mouse ('mm') or human ('hs') data respectively
```
Examples for mouse and human paired-end and single-end data are shown below assuming the reads are in `.fastq` (or `.fastq.gz`) files within the `reads` directory.

Example for single-end mouse data:
```sh
bash TESSA.sh -p <threads> -s path/to/salmon/binary/file -i mm -e se -l reads/ -o out_mm
```

Example for paired-end human data:
```sh
bash TESSA.sh -p <threads> -s path/to/salmon/binary/file -i hs -e pe -l reads/ -o out_hs
```

In the output directory a folder `quant_out` is created with the Salmon estimates and the reference file `references.csv`. These files are used to import estimates into R with functions from ExplorATE [package](https://github.com/FemeniasM/ExplorATEproject). Check the [vignette](https://femeniasm.github.io/ExplorATE_vignette/) and the [user guide](https://femeniasm.github.io/ExplorATE_user_guide/) for more information.

