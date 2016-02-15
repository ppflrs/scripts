#!/usr/bin/env python

import argparse
import sys
import HTSeq
import re
import pandas as pd

arg_parser = argparse.ArgumentParser(description='Processes a BAM file into TSV.')
arg_parser.add_argument("input_file",type=str, help='<input file>, can be a stream indicating "-"')
arg_parser.add_argument("-id","--min_id",type=float, default=95.0, help='Minimal %% of identity to reference sequence to gather the read. (Default = 95.0)')
arg_parser.add_argument("-len","--min_len",type=int, default=60, help='Minimal lenght of the read to be proccessed. (Default = 60)')
arg_parser.add_argument("-clip","--max_clip",type=float, default=0.3, help='Max clipping allowed on the alignment. (Default = 0.30)')
args = arg_parser.parse_args()

if args.input_file:
    try:
        if args.input_file == '':
            print "No input file given. exiting..."
            sys.exit(1)

        elif args.input_file == '-':

            bam_file = HTSeq.SAM_Reader(sys.stdin)

        elif args.input_file != '-':

            bam_file = HTSeq.BAM_Reader(args.input_file)
            output_filename = args.input_file.replace('.bam','.tsv')

    except Exception as e:
                    print "Failed processing SAM/BAM file"
                    raise
elif not args.input_file:
    print "No input file given. exiting..."
    sys.exit(1)

if args.min_id:
    min_id = float(args.min_id)
if args.min_len:
    min_len = int(args.min_len)
if args.max_clip:
    max_clip = float(args.max_clip)

'''Parses BAM files into TSV format using the following parameters min_len=60, max_clip=0.3, min_id=90.0.'''

def parser_cigar(cigar):
    cigar_ops_dict = {}
    cigar_ops_dict['S'] = 0
    cigar_ops_dict['M'] = 0
    for cigar_op in cigar:
        if cigar_op.type not in cigar_ops_dict.keys():
            cigar_ops_dict[cigar_op.type] = cigar_op.size
        else:
            cigar_ops_dict[cigar_op.type] += cigar_op.size
    return cigar_ops_dict

def parser_md_get_ID(md_string, read_len):
    md_list = re.split('(\D|\W)', md_string)
    md_matches = 0
    md_deletions = 0
    md_mismatches = 0
    for i in md_list:
        if i:
            if i.isdigit():
                md_matches += int(i)
            elif "^" in i:
                md_deletions += 1
            else:
                md_mismatches += 1
    return 100*float(md_matches)/read_len

def parser_aln_list(aln, aln_number, pair_pos, min_len=min_len, max_clip=max_clip, min_id=min_id):

    if aln == None:
        return None

    aln_list = list()

    query_name = aln.read.name + '.' + str(pair_pos)
    query_seq = aln.read.seq
    query_len = len(aln.read.seq)
    query_score = int(aln.optional_field('AS'))
    query_ref = aln.iv.chrom
    query_clip_pct = float(parser_cigar(aln.cigar)['S']) / query_len
    query_id = parser_md_get_ID(aln.optional_field('MD'), query_len)

    if query_len < min_len:
        return None
    elif query_clip_pct >= max_clip:
        return None
    elif query_id < min_id:
        return None
    else:
        aln_list.append(aln_number)
        aln_list.append(query_name)
        aln_list.append(query_ref)
        aln_list.append(query_seq)
        aln_list.append(query_len)
        aln_list.append(query_id)
        aln_list.append(query_score)
        aln_list.append(query_clip_pct)

        return aln_list

def bam_parser_2(bam_file):
    bam_dict = {}

    query_counter = 0

    output_list = list()

    #for aln in itertools.islice( HTSeq.pair_SAM_alignments(HTSeq.BAM_Reader(bam_file)), N ):  # printing first N reads
    for aln in HTSeq.pair_SAM_alignments(bam_file):
        query_counter += 1

        query_1, query_2 = aln

        q1_aln = parser_aln_list(query_1, aln_number = query_counter, pair_pos = 1)
        q2_aln = parser_aln_list(query_2, aln_number = query_counter, pair_pos = 2)

        alns = [q1_aln, q2_aln]

        if alns == [None, None]:
            continue
        else:
            if None in alns:
                alns.remove(None)
            output_list.append(alns)

    df_columns = ['ALN','QUERY','REF','SEQ','LEN','ID','SCORE','CLIP_PCT']
    output_list = [item for sublist in output_list for item in sublist]

    return pd.DataFrame(output_list, columns=df_columns)

t_df = bam_parser_2(bam_file)
t_df.to_csv(output_filename, sep='\t', header=False, index=False)

print 'Output file:', output_filename
