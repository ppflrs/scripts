
# coding: utf-8

# In[18]:

import sys
import HTSeq
import re
import pandas as pd

if len(sys.argv) < 2:
    sys.exit('Usage: %s BAM_file(s)' % sys.argv[0])

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


# In[3]:

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


# In[4]:

def parser_aln_list(aln, aln_number, pair_pos, min_len=60, max_clip=0.3, min_id=90.0):

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

        return '\t'.join(map(str,aln_list)) + '\n'


# In[5]:

def bam_parser_2(bam_file):
    bam_dict = {}

    query_counter = 0

    output = ''

    #for aln in itertools.islice( HTSeq.pair_SAM_alignments(HTSeq.BAM_Reader(bam_file)), 10000 ):  # printing first 5 reads
    for aln in HTSeq.pair_SAM_alignments(HTSeq.BAM_Reader(bam_file)):  # printing first 5 reads
        query_counter += 1


        query_1, query_2 = aln


        q1_aln = parser_aln_list(query_1, aln_number = query_counter, pair_pos = 1)
        q2_aln = parser_aln_list(query_2, aln_number = query_counter, pair_pos = 2)

        aln_set = set([q1_aln, q2_aln])

        if aln_set == {None}:
            continue
        else:
            aln_set.discard(None)
            output += '\n'.join(aln_set)
            #print '\n'.join(aln_set)
    return output

for bam_file in sys.argv[1:]:

    t_df2 = bam_parser_2(bam_file)
    print t_df2
