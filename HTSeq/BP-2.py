import argparse
import sys

'''Gets a file produced from BAM2TSV and removes duplicates, and low score alignments.
produces also a fasta file.'''

def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("input_file",type=str, help='<input file> [(full path), -b/-s required], can be a stream indicating "-"')
    args = arg_parser.parse_args()

    if args.input_file:
        try:
            if args.input_file == '':
                print "No input file given. exiting..."
                sys.exit(1)

            elif args.input_file == '-':

                import pandas as pd
                df = pd.read_table(sys.stdin, names=['ALN','QUERY','REF','SEQ','LEN','ID','SCORE','CLIP_PCT'])

            elif args.input_file != '-':

                import pandas as pd
                df = pd.read_table(args.input_file, names=['ALN','QUERY','REF','SEQ','LEN','ID','SCORE','CLIP_PCT'])

        except Exception as e:
                        print "Failed processing SAM/BAM file"
                        raise
    elif not args.input_file:
        print "No input file given. exiting..."
        sys.exit(1)

    #print df

    try:
        dataset_id = df.ix[0]['QUERY'].split('.')[0]
        #print 'Analyzing ', dataset_id
    except Exception as e:
        if args.input_file != '-':
            error_msg = 'Error: No alignments in input file.' + args.input_file
            sys.exit(error_msg)
        raise

    #print 'Analyzing', args.input_file

    amb_count = None
    aligned_aln_list = list()
    amb_list = list()

    dataset_id = df.ix[0]['QUERY'].split('.')[0]

    #print 'PRE removal', len(df)
    df = df[df['ID'] >= 95.0]

    #print 'POS filtro ID > 95', len(df)

    df2 = df.sort_values(by=['ALN','SCORE'], ascending=[1,0]).drop_duplicates('ALN')

    df2['MASTER_QUERY'] = df2['QUERY'].apply(get_read_name)
    gdf2 = df2.groupby('MASTER_QUERY')

    aligned_aln_list, amb_list = dupe_remover(gdf2)
    #print 'POS dupe_remover', len(aligned_aln_list)

    if len(aligned_aln_list) > 0:
        unique_df = pd.concat(aligned_aln_list)
    else:
        error_msg = "Error: No relevant alignments to process in " + args.input_file
        sys.exit(error_msg)

    '''If there are ambiguous reads it will write the FASTA and TSV files'''
    if len(amb_list) > 0:
        amb_df = pd.concat(amb_list)

        amb_count = len(amb_df['MASTER_QUERY'].unique())

        with open('./'+ dataset_id + '.amb.fasta','w') as fh_amb:
            ambiguous_reads = amb_df.apply(lambda x: df_2_fasta(x), axis = 1).reset_index(drop=True)
            for ambiguous_read in ambiguous_reads:
                fh_amb.write(ambiguous_read)

        output_columns = ['MASTER_QUERY','REF','SCORE','ID']
        amb_df = amb_df[output_columns]
        amb_df.rename(columns={'MASTER_QUERY': 'QUERY'}, inplace=True)

        amb_df.to_csv('./'+ dataset_id + '.amb.tsv', sep='\t', header=False, index=False)

    '''FASTA file writing'''
    with open('./'+ dataset_id + '.fasta','w') as fh_aligned:

        aligned_reads = unique_df.apply(lambda x: df_2_fasta(x), axis = 1).reset_index(drop=True)
        for read in aligned_reads:
            fh_aligned.write(read)

        '''tsv file writing'''
        output_columns = ['QUERY','REF','SCORE','ID']
        unique_df = unique_df[output_columns]

        unique_df.to_csv('./'+ dataset_id + '.unique_counts.tsv', sep='\t', header=False, index=False)

    '''Counts file writing'''
    with open('./'+ dataset_id + '.counts','w') as fh_aligned_counts:
        g_unique = unique_df.groupby('REF')

        ###MAYBE USE THE .UNIQUE AND .LOC METHOD?
        for query in g_unique.groups.keys():
            query_count = query + '\t' + str(len(g_unique.get_group(query))) + '\n'
            fh_aligned_counts.write(query_count)
        if amb_count:
            fh_aligned_counts.write('ambigous\t' + str(amb_count) + '\n')


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
    fasta_header += '_'.join(map(str, [dataframe['QUERY'], dataframe['REF'], dataframe['SCORE'], dataframe['ID']])) + '\n'
    fasta_seq = dataframe['SEQ'] + '\n'

    fasta_record = ''
    fasta_record += fasta_header
    fasta_record += fasta_seq

    return fasta_record

main()
