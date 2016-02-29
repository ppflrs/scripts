#!/bin/bash
if [ "$1" == "-h" ] ; then
    echo "Help message"
    echo "Usage: `basename $0`[-h] FASTA-FILE_1.GZ BOWTIE2_INDEX THREADS_NUMBER"
    echo "-h        Show this help text"
    echo "FOLDER    Folder containing fastq.gz files"
    echo "BOWTIE2_INDEX   bowtie2 index"
    echo "THREADS_NUMBER Number of threads to use on bowtie2"
    exit 0
fi

#Argument verification code obtained from http://stackoverflow.com/questions/4341630/checking-for-the-correct-number-of-arguments
#Add option to receive files instead of folder
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 FASTA-FILE_1.GZ BOWTIE2_INDEX THREADS_NUMBER" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi
if ! [ -e "$2".1.bt2 ]; then
  echo "index $2 not found" >&2
  exit 1
fi
#re='^[0-9]+$'
#if ! [[ $3 =~ $re ]] ; then
#   echo "error: Not a valid threads number" >&2
#   exit 1
#fi

#Creates a folder to store the bam files if this folder doesn't exists
#TODO: give the user the option to choose a folder or use this as default
if ! [ -e "./bam/" ]; then
	mkdir ./bam/
fi

#name of the bowtie2_index
bowtie2_index=$2

#number of threads
threads=$3
#read files to map
dataset_1=$1
dataset_2=${dataset_1/"_1.fastq.gz"/$"_2.fastq.gz"}

_now=$(date +'%Y%m%d')
fn=$(basename $dataset_1)
dataset_name=${fn%%_1*}

log_name=./log/"$dataset_name"_mapping_info_"$_now"

if [[ -e $log_name.log ]] ; then
    i=0
    while [[ -e $log_name-$i.log ]] ; do
        let i++
    done
    name=$log_name-$i
fi

echo "Date: "$_now >  $log_name.log
echo "bowtie2 analysis." >>  $log_name.log
bowtie2 --version >> $log_name.log

echo "Analyzing "$dataset_name" at "$(date) >> $log_name.log
#echo "at "$(date)

#BAM file
sorted_bam="$dataset_name"."$bowtie2_index".bt2.sorted

bowtie2 -p $threads -x $bowtie2_index --met-stderr -a --very-sensitive-local -t -q -1 $dataset_1 -2 $dataset_2 | samtools view -Sb -F 4 - | samtools sort -@ $threads -n - ./bam/"$sorted_bam"
