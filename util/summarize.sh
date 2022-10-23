#!/bin/bash
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
tessa ${green}[summarize]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-s${reset} Repeats classification to summarize ['repName' | 'repFamily' | 'repClass' | 'none']
                       If 'none' is selected, raw count estimates for all TE locus are returned.
   ${yellow}-g${reset} Groups. A string vector indicating the experimental conditions 
                       of each library. Example: c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S3')
   ${yellow}-o${reset} Output directory used in quant
   ${yellow}-help${reset}
" 
    exit 1
}

print_help () {
    echo "$0
Usage:  
tessa ${green}[summarize]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-s${reset} Repeats classification to summarize ['repName' | 'repFamily' | 'repClass' | 'none']
                       If 'none' is selected, raw count estimates for all TE locus are returned.
   ${yellow}-g${reset} Groups. A string vector indicating the experimental conditions 
                       of each library. Example: c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S3')
   ${yellow}-o${reset} Output directory used in quant
   ${yellow}-help${reset}

This function loads transposons counts and return a count table by repName, repFamily, repClass or raw table 
without classification. The following packages are required: readr, stringr,tximport, edgeR, SummarizedExperiment, csaw and DESeq2
Users can install packages with the following command in R enviroment:

install.packages(c('BiocManager','stringr','readr')
BiocManager::install(c('edgeR','SummarizedExperiment','DESeq2', 'tximport', 'csaw'))

The argument -g ('group') indicate the sample conditions. This argument must be entered as a string vector indicating the experimental 
conditions of each library. Example: c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3','S3')
If none is selected, raw count estimates are returned.

The program writes a counts table in the output directory and make a TEs.RData file that can be directly imported into R.


A typical command to get raw counts is:

bash tessa summarize -s 'repName' -g 'none' -o ./quant/out/folder


A typical command to get edgeR counts is:

bash tessa summarize -r 'repName' -g c('S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S3') -o ./quant/out/folder

" 
}


abort() 
{
    echo -e >&2 "
${red}
_____________________
  summarize aborted 
=====================
${reset}"
    echo -e "${red}An error occurred. Check the arguments, and paths of the config.file file
    Exiting...${reset}" >&2
    exit 1
}

trap 'abort' 0
set -e

   

echo "${green} Summarizing counts in TESSA with arguments:${reset}"
while getopts ":s:g:o:h:" opt; do
    case $opt in

        s)
            summarize="$OPTARG"
            echo "-s <Repeats classification> = $summarize"
            ;;
        g)
            groups="$OPTARG"
            echo "-g <Groups> = $groups"
            ;;
        o)
            outfolder=`realpath $OPTARG`
            echo "-o <quant output directory> = $outfolder"
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

echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} importing files and summarizing"

Rscript $UTILDIR/summarize.r $outfolder $summarize $groups

trap : 0
echo >&2 "${green}
======================================================
          TESSA summary finished successfully            
======================================================
"

