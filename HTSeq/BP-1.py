#!/usr/bin/env python

import argparse
import sys
import HTSeq
import re
import pandas as pd

'''Parses BAM files into TSV format using the following parameters min_len=60, max_clip=0.3, min_id=90.0.'''

def main():
    arg_parser = argparse.ArgumentParser(description='Processes a BAM file into TSV.')
    arg_parser.add_argument("input_file",type=str, help='<input file>, can be a stream indicating "-"')
    arg_parser.add_argument("-id","--min_id",type=float, default=95.0, help='Minimal %% of identity to reference sequence to gather the read. (Default = 95.0)')
    arg_parser.add_argument("-len","--min_len",type=int, default=60, help='Minimal lenght of the read to be proccessed. (Default = 60)')
    arg_parser.add_argument("-clip","--max_clip",type=float, default=0.3, help='Max clipping allowed on the alignment. (Default = 0.30)')
    arg_parser.add_argument("--out_dir",type=str, default='./', help='Folder where to store the output files.')
    arg_parser.add_argument("--mode",type=str, default='paired', help='Alignment type of the input files. (paired or single)')
    args = arg_parser.parse_args()

    dataset_id = ''

    if args.input_file:
        try:
            if args.input_file == '':
                print "No input file given. exiting..."
                sys.exit(1)

            elif args.input_file == '-':

                bam_file = HTSeq.SAM_Reader(sys.stdin)

            elif args.input_file != '-':

                bam_file = HTSeq.BAM_Reader(args.input_file)

                dataset_id = args.input_file.split('/')[-1].split('.')[0]

        except Exception as e:
                        print "Failed processing SAM/BAM file"
                        raise
    elif not args.input_file:
        sys.exit("No input file given. exiting...")

    if args.min_id:
        min_id = float(args.min_id)
    if args.min_len:
        min_len = int(args.min_len)
    if args.max_clip:
        max_clip = float(args.max_clip)
    if args.out_dir:
        if args.out_dir != './':
            import os

            out_dir = str(args.out_dir) + '/'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        else:
            out_dir = str(args.out_dir)
    if args.mode == 'paired':
        mode = str(args.mode)
    elif args.mode == 'single':
        mode = str(args.mode)
    else:
        sys.exit("No valid aligment type.")

    '''DF containing the raw alignments'''
    df = bam_parser_2(bam_file, min_len=min_len, max_clip=max_clip, min_id=min_id, mode=mode)

    try:
        if dataset_id == '':
            dataset_id = df.ix[0]['QUERY'].split('.')[0]
    except Exception as e:
        if args.input_file != '-':
            error_msg = 'Error: No alignments in input file.' + args.input_file
            sys.exit(error_msg)
            raise

    amb_summary = None
    aligned_aln_list = list()
    amb_list = list()

    if dataset_id == '':
        dataset_id = df.ix[0]['QUERY'].split('.')[0]

    df2 = df.sort_values(by=['ALN','SCORE'], ascending=[1,0]).drop_duplicates('ALN')

    df2['MASTER_QUERY'] = df2['QUERY'].apply(get_read_name)
    gdf2 = df2.groupby('MASTER_QUERY')

    aligned_aln_list, amb_list = dupe_remover(gdf2)

    if len(aligned_aln_list) > 0:
        unique_df = pd.concat(aligned_aln_list)
    else:
        error_msg = "Error: No relevant alignments to process in " + args.input_file
        sys.exit(error_msg)

    '''If there are ambiguous reads it will write the FASTA and TSV files'''
    if len(amb_list) > 0:
        amb_df = pd.concat(amb_list)

        g_amb_df = amb_df.groupby('MASTER_QUERY')
        amb_df = g_amb_df.apply(amb_cluster)
        amb_df = amb_df.reset_index(level=0, drop=True)

        amb_df.columns = ['ALN','QUERY','REF','SEQ','LEN','ID','SCORE','CLIP_PCT','MASTER_QUERY','AMB_STR']

        '''Counts the ambiguous reads'''
        amb_count = len(amb_df.drop_duplicates('MASTER_QUERY'))
        amb_summary = 'ambiguous\t' + str(amb_count) + '\n'
        for ref in sorted(amb_df['REF'].unique()):
            amb_count = len(amb_df.loc[amb_df['REF'] == ref])
            amb_summary += ref + '-amb\t' + str(amb_count) + '\n'

        '''FASTA file writing of ambiguously aligned reads'''
        with open(out_dir + dataset_id + '.amb.fasta','w') as fh_amb:
            ambiguous_reads = amb_df.apply(lambda x: df_2_fasta(x), axis = 1).reset_index(drop=True)
            for ambiguous_read in ambiguous_reads:
                fh_amb.write(ambiguous_read)

        output_columns = ['MASTER_QUERY','REF','SCORE','ID','AMB_STR']
        amb_df = amb_df[output_columns]

        amb_df.rename(columns={'MASTER_QUERY': 'QUERY'}, inplace=True)

        amb_df.to_csv(out_dir + dataset_id + '.amb.tsv', sep='\t', header=False, index=False)

    '''FASTA file writing'''
    with open(out_dir + dataset_id + '.fasta','w') as fh_aligned:

        aligned_reads = unique_df.apply(lambda x: df_2_fasta(x), axis = 1).reset_index(drop=True)
        for read in aligned_reads:
            fh_aligned.write(read)

        '''tsv file writing'''
        output_columns = ['QUERY','REF','SCORE','ID']
        unique_df = unique_df[output_columns]

        unique_df.to_csv(out_dir + dataset_id + '.unique_counts.tsv', sep='\t', header=False, index=False)

    '''Counts file writing'''
    with open(out_dir + dataset_id + '.counts','w') as fh_aligned_counts:
        g_unique = unique_df.groupby('REF')

        for query in sorted(unique_df['REF'].unique()):

            query_count = len(unique_df.loc[unique_df['REF'] == query])
            query_string = query + '\t' + str(query_count) + '\n'

            fh_aligned_counts.write(query_string)

        if amb_summary:
            fh_aligned_counts.write(amb_summary)

def amb_cluster(group):
    long_name_list = []
    for ref in group.REF:
        long_name_list.append(ref + '-amb')

    long_name = '-'.join(long_name_list)

    idx_to_get = group['SCORE'].idxmax()

    group['REF_AMB'] = long_name

    df_out = group.loc[idx_to_get,:]

    return df_out.reset_index(level=0, drop=True)

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

def parser_aln_list(aln, aln_number, pair_pos, min_len, max_clip, min_id):

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

def bam_parser_2(bam_file, min_len, max_clip, min_id, mode):
    bam_dict = {}

    query_counter = 0

    output_list = list()

    if mode == 'paired':
        #import itertools
        #for aln in itertools.islice( HTSeq.pair_SAM_alignments(bam_file), 1000 ):  # printing first N reads
        for aln in HTSeq.pair_SAM_alignments(bam_file):
            query_counter += 1

            query_1, query_2 = aln

            q1_aln = parser_aln_list(query_1, aln_number = query_counter, pair_pos = 1, min_len=min_len, max_clip=max_clip, min_id=min_id)
            q2_aln = parser_aln_list(query_2, aln_number = query_counter, pair_pos = 2, min_len=min_len, max_clip=max_clip, min_id=min_id)

            alns = [q1_aln, q2_aln]

            if alns == [None, None]:
                continue
            else:
                if None in alns:
                    alns.remove(None)
                output_list.append(alns)

    elif mode == 'single':
        for aln in bam_file:

            query_counter += 1

            query_1 = aln

            q1_aln = parser_aln_list(query_1, aln_number = query_counter, pair_pos = 1, min_len=min_len, max_clip=max_clip, min_id=min_id)

            alns = [q1_aln]

            if q1_aln != None:
                output_list.append(alns)

    df_columns = ['ALN','QUERY','REF','SEQ','LEN','ID','SCORE','CLIP_PCT']
    output_list = [item for sublist in output_list for item in sublist]

    return pd.DataFrame(output_list, columns=df_columns)

def get_read_name(read):
    return '.'.join(read.split('.')[:2])

def dupe_remover(grouped_df):

    aligned_aln_list = list()
    amb_list = list()

    for group in grouped_df:
        group_name, group_df = group
        if len(group_df) > 1:
            local_max_score = group_df['SCORE'].max()
            if len(group_df) != len(group_df['SCORE'].unique()):
                # It will return the alignment(s) with the top score of the group
                group_df_clean = group_df.loc[group_df['SCORE'] == local_max_score]
                if len(group_df_clean) == 1:
                    aligned_aln_list.append(group_df_clean)
                else:
                    #There are multiple top-scored alignments. These should be marked as amb
                    group_df_clean = group_df.loc[group_df['SCORE'] == local_max_score]
                    amb_list.append(group_df_clean)
            else:
                #If we have the same number of unique scores and elements in the group it means you have a top scorer.
                group_df_clean = group_df.loc[group_df['SCORE'] == local_max_score]
                aligned_aln_list.append(group_df_clean)
        else:
            #There's only one alignment
            aligned_aln_list.append(group_df)
    return aligned_aln_list, amb_list

def df_2_fasta(dataframe):

    fasta_header = '>'
    if 'AMB_STR' in set(dataframe.index):
        fasta_header += '_'.join(map(str, [dataframe['QUERY'], dataframe['AMB_STR'], dataframe['SCORE'], dataframe['ID']])) + '\n'
    else:
        fasta_header += '_'.join(map(str, [dataframe['QUERY'], dataframe['REF'], dataframe['SCORE'], dataframe['ID']])) + '\n'
    fasta_seq = dataframe['SEQ'] + '\n'

    fasta_record = ''
    fasta_record += fasta_header
    fasta_record += fasta_seq

    return fasta_record

main()
