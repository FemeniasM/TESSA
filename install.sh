#!/bin/bash
# Using getopt

#########################################################################################
# This script is based on Selective Alignment of sequences, 
# (Srivastava, A., Malik, L., Sarkar, H. et al. Alignment 
# and mapping methodology influence transcript abundance 
# estimation. Genome Biol 21, 239 (2020). 
# https://doi.org/10.1186/s13059-020-02151-8)
#----------------------------------------------------
# It assumes awk, bedtools and Salmon is 
# available.
# We have tested this script with awk 4.1.4, 
# bedtools v2.30.0 on an Ubuntu system. 
# 
# This script uses SeqKit (Wei Shen, Copyright ©2016-2019 Oxford Nanopore Technologies.
# W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for 
# FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962)
#
# This script uses RM2Bed.py from RepeatMaster distribution (Authors: David Ray 
# and Robert Hubley) for overlaps resolution.
########################################################################################

threads=1
kmer=31
salmon_path=="salmon"
source_dir=$(realpath $(dirname "${BASH_SOURCE[0]}"))
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`

# Argument Parsing
print_usage () {
    echo "$0
Usage:  
ExplorATE ${yellow}mo${green} [flags]${reset}

Flags:
   ${green}-p${reset} threads [N] (1 default) 
   ${green}-k${reset} kmer [N] (31 default)
   ${green}-s${reset} salmon binary path (salmon default)" 
    exit 1
}

abort() 
{
    echo -e >&2 '
===============
/// ABORTED ///
===============
'
    echo -e "${red}An error occurred. Exiting...${reset}" >&2
    exit 1
}

trap 'abort' 0
set -e

echo "${green}installing TESSA ${reset}"
while getopts ":a:p:b:c:s:f:g:r:o:k:t:u:v:h:" opt; do
    case $opt in
        s)
            salmon_path=`realpath $OPTARG`
            echo "-s <salmon path> = $salmon_path"
            ;;
        p)
            threads="$OPTARG"
            echo "-p <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;
        h)
            print_usage && exit 1
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

gzip -d $source_dir/install/*

mkdir -p $source_dir/index
cd $source_dir/index

echo "---------------------------
Genomic features processing
---------------------------"

echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] Indexing in Salmon${reset}"

echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] Indexing mouse dataset${reset}"
$salmon_path index -p $threads -t $source_dir/install/mm.fa -k $kmer -i mm --decoys $source_dir/install/mm.txt
echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] Indexing human dataset${reset}"
$salmon_path index -p $threads -t $source_dir/install/hs.fa -k $kmer -i hs --decoys $source_dir/install/hs.txt

mv $source_dir/install/mm.csv $source_dir/index/mm/references.csv
mv $source_dir/install/hs.csv $source_dir/index/hs/references.csv

rm -r $source_dir/install
trap : 0
echo >&2 "${green} Installation successfully"
