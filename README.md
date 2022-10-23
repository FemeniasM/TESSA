TESSA - Transposons Exploration with Salmon's Selective Alignment
====
### Overview

TESSA runs Salmon with Selective Alignment using pre-built indices based on GENCODE sequences and RepeatMasker references for human and mouse.

### Requirements
----
TESSA requires **Salmon** to be installed (tested with v1.4.0) and at least 35Gb of disk space. 
Additionally **Bedtools** version 2.29.0 or higher is required

### Installation
----
#### Step 1 - download the repository

Clone the repository to the desired directory with the following command.
```sh
git clone https://github.com/FemeniasM/TESSA
```

#### Step 2 - edit the config.file 
Edit the config.file to set the paths to the required programs. In addition, the links to the reference files can be established to download them with the `download` function. The config.file must have the following format:
```
## Homo sapiens references
genome_hg=https://url/to/hg38.fa.gz
rm_out_hg=https://url/to/hg38.fa.out.gz

## Mus musculus references
genome_mm=https://url/to/mm39.fa.gz
rm_out_mm=https://url/to/mm39.fa.out.gz

## Drosophila melanogaster references
genome_dm=https://url/to/dm6.fa.gz
rm_out_dm=https://url/to/dm6.fa.out.gz

## Programs
salmon=/path/to/salmon-X.X.X_linux_x86_64/bin/salmon
bedtools=/path/to/bedtools2/bin/bedtools
```

#### Step 3 - make index for for the species of interest

To create indices for species of interest, users can manually download the necessary references from the web [site](https://genome.ucsc.edu/cgi-bin/hgGateway). Alternatively you can use the download function which will download the files for you. In both cases the references must be located in the `TESSA/ref` folder, as indicated below:
```
TESSA/
    |_ref/
        |_hs/
        |_mm/
        |_dm/
```
Files must contain the species identifier in the name ('hg' for human, 'mm' for mouse, and 'dm' for D. melanogaster).The supported extensions are '.fa' for genome, and '.out' for the RepeatMasker file.

Users can use the `download` argument to download the reference files. An example for human and mouse is shown below:

```
$ bash tessa download -r 'hg,mm'
```
note that when more than one reference is added they must be comma separated. 
The same command adding the *Drosophila melanogaster* references is shown below

```
$ bash tessa download -r 'hg,mm,dm'
```
The indices are created only once for each species of interest.

### How to use
----

Once the indices are built, `quant` mode is the default command for running `tessa`. The arguments are detailed below:

```
Usage:  
tessa quant [flags]

Flags:
   -t     Number of therads [N] (1 default) 
   -k     k-mer size [N] (31 default)
   -i     salmon index ['hg']['mm']['dm']
   -r     Indicates the format of the libraries: ['pe'] for paired-end and ['se']
          for single-end reads. Supported extensions are: <.fq> or <.fastq>
          or <.fq.gz> or <.fastq.gz>. Paired-end file names should contain 
          _R1 _R2. Example: sample_R1.fq.gz, sample_R2.fq.gz (mandatory) 
   -f     Path to folder with fastq files (mandatory)
   -o     Output directory path 
   -help  help
```
A typical quant mode command is:
```
bash tessa quant -p 12 -k 31 -i 'hg' -r 'pe' -f path/to/reads/ -o out/directory/
```

Input reads files must be in `.fastq` format. Supported extensions are: `<.fq>` or `<.fastq>` or `<.fq.gz>` or `<.fastq.gz>`
Paired end file names should contain `_R1` `_R2`. For example: `sample_R1.fq.gz`, `sample_R2.fq.gz`

Examples for mouse and human paired-end and single-end data are shown below assuming the reads are in `.fastq` (or `.fastq.gz`) files within the `reads` directory.

Example for single-end mouse data:
```
bash TESSA.sh -p <threads> -s path/to/salmon/binary/file -i mm -e se -l reads/ -o out_mm
```

Example for paired-end human data:

```
bash TESSA.sh -p <threads> -s path/to/salmon/binary/file -i hs -e pe -l reads/ -o out_hs
```

In the output directory a folder `quant_out` is created with the Salmon estimates and the reference file `references.csv`. These files are used to import estimates into R with functions from ExplorATE [package](https://github.com/FemeniasM/ExplorATEproject). Check the [vignette](https://femeniasm.github.io/ExplorATE_vignette/) and the [user guide](https://femeniasm.github.io/ExplorATE_user_guide/) for more information.

#### Summarize counts

TESSA includes a function to export the TE counts summarized by `repName`, `repFamily`, `repClass` or a raw table (`none`, without classification).
To use this TESSA function the following additional R packages are required: `readr`, `stringr`, `tximport`, `edgeR`, `SummarizedExperiment`, `csaw` and `DESeq2`
Users can install packages with the following command in the R environment:

```
install.packages(c('BiocManager','stringr','readr')
BiocManager::install(c('edgeR','SummarizedExperiment','DESeq2', 'tximport', 'csaw'))
```

If all R packages are available, users can use the `summarize` function with the following arguments:

```
Usage:  
tessa summarize [flags]

Flags:
   -s    Repeats classification to summarize ['repName' | 'repFamily' | 'repClass' | 'none']
         If 'none' is selected, raw count estimates for all TE locus are returned.
   -g    Groups. A string vector indicating the experimental conditions 
         of each library. Example: c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S3')
   -o    Output directory used in 'quant' step
   -help help
```

The argument -g (`group`) indicate the sample conditions. This argument must be entered as a string vector indicating the experimental conditions of each library. Example: c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3','S3')
If `none` is selected, raw counts estimates are returned. The program writes a counts table in the output directory and make a `TEs.RData` file that can be directly imported into R.

A typical command to get raw counts is:
```
bash tessa summarize -s 'repName' -g 'none' -o ./quant/out/folder
```

A typical command to get edgeR counts is:
```
bash tessa summarize -r 'repName' -g c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S3') -o ./quant/out/folder
```

