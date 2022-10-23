#!/bin/bash
threads=1
kmer=31
seqkit=$UTILDIR/seqkit
RM2Bed=$UTILDIR/RM2Bed.py
ref='hg' #'hg,mm,dm'
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
blue=`tput setaf 4`
magenta=`tput setaf 5`
reset=`tput sgr0`

print_usage () {
    echo "$0
Usage:  
install.sh ${yellow}[flags]${reset}

Flags:
   ${yellow}-t${reset} threads [N] (1 default) 
   ${yellow}-k${reset} kmer [N] (31 default)
   ${yellow}-r${reset} references comma separated ['hg,mm,dm'] ( -r 'hg' is default)
   ${yellow}-help${reset}
" 
    exit 1
}

print_help () {
    echo "$0
Usage:  
install.sh ${yellow}[flags]${reset}

Flags:
   ${yellow}-t${reset} threads [N] (1 default) 
   ${yellow}-k${reset} kmer [N] (31 default)
   ${yellow}-r${reset} references comma separated ['hg,mm,dm'] ( -r 'hg' is default)
   ${yellow}-help${reset}

Be sure to edit the corresponding paths in the config.file. This file contains the url of each reference 
(optional if you use the download function) and the paths to the salmon and bedtools programs. For example:

$ nano config.file

------------- config.file
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
----------

If Salmon is not installed on your system, follow the recommendations below:

$ wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz
$ tar –xvzf salmon-1.9.0_linux_x86_64.tar.gz

If bedtools is not installed on your system, follow the recommendations below:
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
$ tar -zxvf bedtools-2.29.1.tar.gz
$ cd bedtools2
$ make

A typical command for installation is:

$ bash tessa index -t 12 -k 31 -r 'hg'

The 'index' argument assumes that the reference files are downloaded and located in the ref folder with this 
subdirectory structure:

TESSA/
    |_ref/
        |_hs/
        |_mm/
        |_dm/

Users can use the 'download' argument to download the reference files. An example for human and mouse is shown below:

$ bash tessa download -r 'hg,mm'

note that when more than one reference is added they must be comma separated. 
The same command adding the Drosophila melanogaster references is shown below

$ bash tessa download -r 'hg,mm,dm'

Offline installation:
=====================

Make sure to download your reference files and place them inside the TESSA/ref folder with a subdirectory for each species, for example:

TESSA/
    |_ref/
        |_hs/
        |_mm/
        |_dm/

files must contain the species identifier in the name ('hg' for human, 'mm' for mouse, and 'dm' for D. melanogaster), 
The supported extensions are '.fa' for genome, and '.out' for the RepeatMasker file.

Installation times: 
==================
The installation time depends on the reference indexes that are created. For example, to index the human genome takes 
about an hour. Users can make multiple threads for this task, and run the installation script for each species individually.
" 
    exit 1
}


abort() 
{
    echo -e >&2 "
${red}
_____________________
INSTALLATION ABORTED 
=====================
${reset}"
    echo -e "${red}An error occurred. Check the arguments, and paths of the config.file file
    Exiting...${reset}" >&2
    exit 1
}

trap 'abort' 0
set -e

echo "${green} Installing TESSA with arguments:${reset}"
while getopts ":t:k:r:h:" opt; do
    case $opt in

        t)
            threads="$OPTARG"
            echo "-t <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;

        r)
            ref="$OPTARG"
            echo "-r <references> = $ref"
            ;;
        h)  
            print_help && exit 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage && exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage && exit 1
            ;;
            
    esac
done

if [ ! -s $CONFDIR/config.file ]; then
  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file config.file $fn ${reset}"
  exit 1
fi

source $CONFDIR/config.file

#checking program versions
echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Checking program versions${reset}"
bedtools_ver=$($bedtools --version | awk -F" v" '{print $2}')
salmon_ver=$($salmon --version | awk -F" " '{print $2}')

if { echo "$bedtools_ver"; echo "2.29.0"; } | sort --version-sort --check=silent; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] TESSA requires bedtools version 2.29.1 or later. 
    You are using betools version $bedtools_ver ${reset}"
    exit 1
fi

if { echo "$salmon_ver"; echo "1.4.0"; } | sort --version-sort --check=silent; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] TESSA requires bedtools version 1.4.0 or later. 
    You are using Salmon version $salmon_ver ${reset}"
    exit 1
fi

# if [[ $download = "TRUE" ]]; then
# echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Downloading reference files${reset}"
#   for i in $(echo "$ref" | sed 's/,/\t/g');do
#   
#   wget -P $CONFDIR/ref/$i $(set | grep genome_${i} | awk -F"=" '{print $2}')
#   wget -P $CONFDIR/ref/$i $(set | grep rm_out_${i} | awk -F"=" '{print $2}')
#   gunzip -d $CONFDIR/ref/$i/*.gz
# 
#   done
# fi

for i in $(echo "$ref" | sed 's/,/\t/g');do

  RM_out=$(find $CONFDIR/ref/${i}/*.out)
  genome=$(find $CONFDIR/ref/${i}/*.fa)
  
  mkdir -p $CONFDIR/pre_build_${i}
  mkdir -p $CONFDIR/index/${i}
  cd $CONFDIR/pre_build_${i}
  
   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Resolving overlaps with RM2Bed${reset}"

 python3 $RM2Bed --out_dir $CONFDIR/pre_build_${i} --out_prefix RM_ovres --sort_criterion 'size' --ovlp_resolution 'higher_score' $RM_out

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Removing Low complexity, Simple repeats and small RNAs seqs from RepeatMasker reference${reset}"

  grep -w -v "Unknown\\|rRNA\\|Satellite\\|Simple_repeat\\|Low_complexity\\|RNA\\|cRNA\\|snRNA\\|srpRNA\\|tRNA\\|Other" RM_ovres_rm.bed \
  | awk -v OFS='\t' '{print $1, $2, $3, $1":"$2":"$3":"$4"/"$7"/"$8}' | sort -k1,1 -k2,2n > RMgen.bed

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Merge RepeatMasker BED file"

  $bedtools merge -c 4 -o first -i RMgen.bed | sort -k1,1 -k2,2n > RMgen_merged.bed
  
  ### with merge# awk '{print $4}' RMgen_merged.bed | awk -F":" -v OFS=";" '{print $1":"$2":"$3":"$4, $1,$2,$3,$4}' > $CONFDIR/index/${i}/references.csv

  awk '{print $4}' RMgen.bed | awk -F":" -v OFS=";" '{print $1":"$2":"$3":"$4, $1,$2,$3,$4}' > $CONFDIR/index/${i}/references.csv
 
   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Extract target sequences"

  ### with merge# $bedtools getfasta -fi $genome -bed RMgen_merged.bed -nameOnly -fo target1.fa

  $bedtools getfasta -fi $genome -bed RMgen.bed -nameOnly -fo target1.fa

  rm RM_ovres_rm.bed RMgen.bed

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Removing duplicated and short sequences from target file"

  $seqkit rmdup -s target1.fa | $seqkit seq -m $kmer -g > target.fa

  rm target1.fa

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Extract decoys reference"

  $seqkit fx2tab --length --name -i $genome | sort -k1,1 -k2,2n > genome_len.bed

  $bedtools complement -i RMgen_merged.bed -g  genome_len.bed -L | sort -k1,1 -k2,2n | $bedtools merge | \
  awk -v OFS='\t' '{print $1, $2, $3, $1":"$2":"$3}' > decoys.bed

  $bedtools getfasta -fi $genome -bed decoys.bed -nameOnly -fo decoys1.fa

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Removing duplicated and short sequences from decoy file"

  $seqkit rmdup -s decoys1.fa | $seqkit -is replace -p "n" -r "" | $seqkit seq -g -m $kmer > decoys.fa

  rm decoys1.fa

  echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Making decoys.txt file"

  awk '/^>/ {print $0; next}' decoys.fa | sed 's/^>//g' > decoys.txt

   echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Making reference 'gentrome'"

  cat target.fa decoys.fa > reftrme.fa

  rm decoys.fa target.fa

  echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Make index"

  $salmon index -p $threads -t reftrme.fa -k $kmer -i $CONFDIR/index/${i} --decoys decoys.txt

  echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Index finished successfully"

  rm -r $CONFDIR/pre_build_${i}


 done


trap : 0
echo >&2 "${green}
The installation was completed successfully          
${magenta}===========================================
"
