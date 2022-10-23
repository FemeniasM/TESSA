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
tessa ${green}[download]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-r${reset} references comma separated ['hg,mm,dm'] ( -r 'hg' is default)
   ${yellow}-help${reset}
" 
    exit 1
}

print_help () {
    echo "$0
Usage:  
tessa ${green}[download]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-r${reset} references comma separated ['hg,mm,dm'] ( -r 'hg' is default)
   ${yellow}-help${reset}

Be sure to edit the corresponding paths in the config.file. This file contains the url of each reference 
(optional if you use the -d argument) and the paths to the salmon and bedtools programs. For example:

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

tessa download -r 'hg'

More than one reference is allowed, for example:

tessa download -r 'hg,mm'
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
while getopts ":r:h:" opt; do
    case $opt in

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


echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Downloading reference files${reset}"
for i in $(echo "$ref" | sed 's/,/\t/g');do
  
  wget -P $CONFDIR/ref/$i $(set | grep genome_${i} | awk -F"=" '{print $2}')
  wget -P $CONFDIR/ref/$i $(set | grep rm_out_${i} | awk -F"=" '{print $2}')
  gunzip -d $CONFDIR/ref/$i/*.gz

done

trap : 0
echo >&2 "${green}
downloads completed successfully         
${magenta}===========================================
"