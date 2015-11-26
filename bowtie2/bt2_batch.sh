#!/bin/bash

#Argument verification code obtained from http://stackoverflow.com/questions/4341630/checking-for-the-correct-number-of-arguments
#Add option to receive files instead of folder
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 DIRECTORY BOWTIE2_INDEX THREADS_NUMBER" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi
if ! [ -d "$1" ]; then
  echo "$1 not a directory" >&2
  exit 1
fi
if ! [ -e "$2".1.bt2 ]; then
  echo "index $2 not found" >&2
  exit 1
fi
re='^[0-9]+$'
if ! [[ $3 =~ $re ]] ; then
   echo "error: Not a valid threads number" >&2
   exit 1
fi

#Directory containing the fastq.gz files
dir_gz=$1

#Creates a folder to store the bam files if this folder doesn't exists
#TODO: give the user the option to choose a folder or use this as default
if ! [ -e "./bam/" ]; then
	mkdir ./bam/
fi

#name of the bowtie2_index
bowtie2_index=$2

#number of threads
threads=$3

_now=$(date -v-1d +'%Y_%m_%d')
name=mapping"$_now"
if [[ -e $name.log ]] ; then
    i=0
    while [[ -e $name-$i.log ]] ; do
        let i++
    done
    name=$name-$i
fi
touch $name.log

for file in $(find "$dir_gz" -name "*_1.fastq.gz");
do
	fn=$(basename $file)
	dataset=${fn%%_1*}
	fn1=$file
	fn2=${fn1/"_1.fastq.gz"/$"_2.fastq.gz"}
  sorted_bam="$dataset"."$index_file".bt2.sorted

  echo "Analyzing "$dataset
  echo "Analyzing "$dataset >> $name.log
  bowtie2 -p $threads -x $bowtie2_index -a --very-sensitive-local -t -q -1 $fn1 -2 $fn2 | samtools view -Sb -F 4 - | samtools sort -n - ./bam/"$sorted_bam" 2>> $name.log

  status=$?
  if test $status -ne 0
    then
      echo "Error while analyzing: '$dataset'"
  fi
done
