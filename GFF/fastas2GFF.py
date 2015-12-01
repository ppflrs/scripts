#!/usr/bin/env python
#fastas2GFF.py gets a multisequence FASTA file and outputs a GFF file
import sys
import os
import datetime
from Bio import SeqIO
from Bio.Seq import Seq

if len(sys.argv) < 3 or sys.argv[1] == '-h' or sys.argv[1] not in ['-v1','-v2']:
    print('Usage:')
    print('\t%s [-h] [-v1/v2] <fasta_file> <bt2_index_name>' % sys.argv[0])
    print
    print('<fasta_file>\t\tFile containing multiple FASTA sequences.')
    print('<bt2_index_name>\tBowtie 2 index name.')
    print('\t\t\tOnly required with -v1 option.')
    print
    print('-h\tPrint this help message')
    print('-v1\tOld fastas2GFF behavior.')
    print('\tCreates a concatenated FASTA file.')
    print('-v2\tNew fastas2GFF behavior(2015.12.01).')
    print('\tCreates a GFF2 file with multiple "genomes" info on it.')
    sys.exit()

def fastas2GFF(input_file, index_name):
    handle = open(input_file, 'rU')
    dna = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    map_dict = {}
    counter = 1 #To comply with the GFF format file specifications needs to start at 1

    output_fasta_file = os.path.splitext(input_file)[0] + '.concatenated.fasta'
    with open(output_fasta_file, 'w') as cat_genes:
        cat_genes.write('>' + index_name + '\n') #FASTA header
        cat_sequence = str()

        seqid_counter = 1 #To order the output of the sequence

        for sequence in dna:
            gene_len = len(str(sequence.seq))

            description = list()
            description.append(index_name) #bowtie index name
            description.append('beja_lab') #source of the sequence
            description.append('CDS') #Feature type: This can be exon, promotor, etc for us CDS it's ok i guess
            description.append(counter) #start
            description.append(counter + gene_len) #end
            description.append('.') #score ???
            description.append('+') #strand
            description.append('.') #score ???
            description.append('GeneID '+ sequence.id + '\n')

            map_dict[seqid_counter] = description
            seqid_counter += 1
            counter += (gene_len)

            cat_genes.write(cat_sequence)

        output_gff_file = os.path.splitext(input_file)[0] + '.GFF2'
    with open(output_gff_file, 'w') as m:
        m.write('##gff-version 2\n')
        m.writelines('\t'.join(map(str,map_dict[k])) for k in map_dict.keys())
    return

def fastas2GFF_v2(input_file,):
    with open(input_file, 'r') as handle:
        dna = list(SeqIO.parse(handle, 'fasta'))

    seqs_counter = 0

    seqs_dict = dict()

    for sequence in dna:
        seqs_counter += 1

        sequence_len_counter = 1
        seq_id = sequence.id

        gene_len = len(str(sequence.seq))

        description = list()
        description.append(seq_id) #bowtie index name
        description.append(seq_id.split('.')[3]) #source of the sequence, accession, etc
        description.append('CDS') #Feature type: This can be exon, promotor, etc for us CDS it's ok i guess
        description.append(sequence_len_counter) #start
        description.append(sequence_len_counter + gene_len) #end
        description.append('.') #score ???
        description.append('+') #strand
        description.append('.') #score ???
        description.append('GeneID '+ '_'.join(sequence.id.split('.')[:3]) + '\n') #Name of the gene that will appear in the HTSeq-count report

        seqs_dict[seqs_counter] = description

    output_gff_file = os.path.splitext(input_file)[0] + '.GFF2'

    with open(output_gff_file, 'w') as gff_handle:
        gff_handle.write('##gff-version 2\n')
        gff_handle.write('#' + str(datetime.date.today()) + '\n')
        gff_handle.write('#' + 'Source: ' + input_file + '\n')

        gff_handle.writelines('\t'.join(map(str,seqs_dict[k])) for k in seqs_dict.keys())
    return

fasta_file = sys.argv[2]

if sys.argv[1] == '-v1':
    bt2_index = sys.argv[3]
    fastas2GFF(fasta_file, bt2_index)
elif sys.argv[1] == '-v2':
    fastas2GFF_v2(fasta_file)
