#!/bin/sh

#  alon.sh
#
#
#  Created by José Flores-Uribe on 8/20/15.
#

#Argument verification code obtained from http://stackoverflow.com/questions/4341630/checking-for-the-correct-number-of-arguments

if [ "$#" -ne 3 ]; then
  echo "Usage: `basename $0`[-h] ID DIR GFF">&2
  echo "-h    Show this help text"
  echo "ID    Min identity %"
  echo "DIR   Subfolder containing the bam files"
  echo "GFF   GFF file"
  exit 1
fi
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  exit 1
fi
if ! [ -d "$2" ]; then
  echo "$2 not a directory" >&2
  exit 1
fi
if ! [ -e "$3" ]; then
  echo "GFF file: $3 not found" >&2
  exit 1
fi

if [ "$1" == "-h" ] ; then
    echo "Help message"
    echo "Usage: `basename $0`[-h] ID DIR GFF"
    echo "-h    Show this help text"
    echo "ID    Min identity %"
    echo "DIR   Subfolder containing the bam files"
    echo "GFF   GFF file"
    exit 0
fi

if ! [ -e ./counts/ ]; then
	mkdir ./counts/
fi

if ! [ -e ./amb/ ]; then
        mkdir ./amb/
fi

if ! [ -e ./fasta/ ]; then
        mkdir ./fasta/
fi
if ! [ -e ./fasta_not_aligned/ ]; then
        mkdir ./fasta_not_aligned/
fi

id=$1
dir_bam=$2
gff_file=$3

echo "Files to analyze: ""$(find ./"$dir_bam" -name "*.bam" | wc -l)"

for file in $(find ./"$dir_bam" -name "*.bam")
do
file_base=$(basename $file)

echo "$file_base"
dataset=${file_base%%.sorted*}

samtools view $file | python parse_and_count_sam_alignment_atlas_final.py '-' $gff_file --fasta ./fasta/"$dataset"."$id".fasta -b ./amb/"$dataset"."$id".amb -o ./counts/"$dataset"."$id".counts --min_read_length 90 --max_clip 0.30 -p 'p' --min_id $id -u ./fasta_not_aligned/"$dataset"."$id".not_aligned.fasta

done
