#!/bin/bash
# Using getopt

#################################################
# This script runs salmon with pre-built indices 
# full genomes and RepeatMasker references.
#################################################

threads=1
kmer=31
idx_path=$CONFDIR/index
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
blue=`tput setaf 4`
magenta=`tput setaf 5`
reset=`tput sgr0`

# Argument Parsing
print_usage () {
    echo "$0
Usage:  
tessa ${yellow}quant${green} [flags]${reset}

Flags:
   ${green}-t${reset} threads [N] (1 default) 
   ${green}-k${reset} kmer [N] (31 default)
   ${green}-i${reset} salmon index ['hg']['mm']['dm']
   ${green}-r${reset} library format ['pe']['se'] 
   ${green}-f${reset} folder with fastq files
   ${green}-o${reset} output directory path 
   ${green}-help${reset} help
" 
    exit 1
}

print_help () {
    echo "$0
Usage:  
tessa ${yellow}quant${green} [flags]${reset}

Flags:
   ${green}-t${reset} threads [N] (1 default) 
   ${green}-k${reset} kmer [N] (31 default)
   ${green}-i${reset} salmon index ['hg']['mm']['dm']
   ${green}-r${reset} library format ['pe']['se'] 
   ${green}-f${reset} folder with fastq files
   ${green}-o${reset} output directory path 
   ${green}-help${reset} help

After the indices are built, quant mode is the default command for running tessa. The indices are only built 
once per species and the user can use the same indices for different samples of the same species. 
A typical quant mode command is:

bash tessa quant -p 12 -k 31 -i 'hg' -r 'pe' -f path/to/reads/ -o out/directory/

Input reads files must be in .fastq format. Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>
Paired end file names should contain _R1 _R2. For example: sample_R1.fq.gz, sample_R2.fq.gz

" 
}

abort() 
{
    echo -e >&2 '
_______________________
quantification aborted
=======================
'
    echo -e "${red}An error occurred. Exiting...${reset}" >&2
   [ -f $outfolder/temp ] && rm -r $outfolder/temp
    exit 1
}

trap 'abort' 0
set -e

file_extension_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] FileExtensionError: Invalid extension${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}"
return     
}

file_name_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] Filename Error: Paired end file names should contain _R1 _R2${reset}"
echo -e "${green}Example: sample_R1.fq.gz, sample_R2.fq.gz${reset}"
return 
}

echo "${green}=============================
${yellow}|||      TESSA v0.1       |||
${green}=============================${reset}"
while getopts ":f:r:i:t:k:o:h:" opt; do
    case $opt in

        f)
            lib_folder=`realpath $OPTARG`
            echo "-f <fastq libraries path> = $lib_folder"
            ;;
        r)
            lib_format="$OPTARG"
            echo "-r <libraries format 'pe' | 'se'> = $lib_format"
            ;;
        i)
            index="$OPTARG"
            echo "-i <index 'hs' | 'mm' | 'dm'> = $index"
            ;;

        t)
            threads="$OPTARG"
            echo "-t <threads> = $threads"
            ;;
        k)
            kmer="$OPTARG"
            echo "-k <kmer> = $kmer"
            ;;

        o)
            outfolder=`realpath $OPTARG`
            echo "-o <Output files Path> = $outfolder"
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

if [ -z "$lib_format" -o -z "$index" -o -z "$lib_folder" -o -z "$threads" -o -z "$kmer" -o -z "$outfolder" ]
then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] missing required argument(s)${reset}"
    print_usage && exit 1
fi

mkdir -p $outfolder/temp
cd $outfolder/temp

source $CONFDIR/config.file

echo -e "${yellow}
Quantifying with Salmon
=======================${reset}"

if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        
        sample_name=$(basename ${fn} | sed 's/.fastq.gz\|.fq.gz\|.fastq\|.fq//g')
        
        echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} processing sample $sample_name"
        
        $salmon quant -i $idx_path/$index -l A --gcBias --useVBOpt -r ${fn} -p $threads --validateMappings 0.65 -o $outfolder/$sample_name
        
        awk 'NR>1 {print $1}' $outfolder/$sample_name/quant.sf >> $outfolder/temp/ids.txt
        
        else 
        file_extension_error
        exit 1
    fi
    done
else 
 if [[ $lib_format == 'pe' ]]; then
    for fn in $lib_folder/*R1*; do
        if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
         sample_name=$(basename ${fn} | sed 's/_R1.fastq.gz\|_R1.fq.gz\|_R1.fastq\|_R1.fq//g')
            if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
            
             ext=$(basename ${fn} | awk -F'_R1' '{print $2}')
             R1=${sample_name}_R1${ext}
             R2=${sample_name}_R2${ext}
            
             echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} processing sample $sample_name"
            
             $salmon quant -i $idx_path/$index -l A --gcBias --useVBOpt -1 $lib_folder/$R1 -2 $lib_folder/$R2 -p $threads --validateMappings 0.65 -o $outfolder/$sample_name
            
             awk 'NR>1 {print $1}' $outfolder/$sample_name/quant.sf >> $outfolder/temp/ids.txt
            else
               file_name_error
               exit 1
            fi
         else
         file_extension_error
         exit 1
         fi
     done
 else
    echo -e "\n${red} Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi
fi

cat $outfolder/temp/ids.txt | sort | uniq > $outfolder/temp/ids_sort_uniq.txt

awk -F";" 'FNR==NR { a[$1]; next } ($1 in a)' $outfolder/temp/ids_sort_uniq.txt $idx_path/$index/references.csv > $outfolder/references.csv

echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} removing temporary files"
rm -r $outfolder/temp/

trap : 0
echo >&2 "${green}
======================================================
          TESSA finished successfully            
======================================================
"
