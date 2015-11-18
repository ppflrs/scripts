#!/usr/bin/env python
#fastas2GFF.py gets a multisequence FASTA file and outputs a GFF file
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

if len(sys.argv) != 3:
    sys.exit('Usage: %s <fasta_file> <bt2_index_name>' % sys.argv[0])

def map_bowtie(input_file, index_name):
    handle = open(input_file, 'rU')
    dna = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    map_dict = {}
    counter = 1 #To comply with the GFF format file specifications needs to start at 1

    output_fasta_file = os.path.splitext(input_file)[0] + '.concatenated.fasta'
    cat_genes = open(output_fasta_file, 'a')
    cat_genes.write('>' + index_name + '\n') #FASTA header

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
        description.append('GeneID '+ sequence.id)# + ' ' + str(seqid_counter).zfill(4)) #CDS0001, CDS0002, etc
        description.append('\n')

        map_dict[seqid_counter] = description
        seqid_counter += 1
        counter += (gene_len)
        cat_genes.write(str(sequence.seq))

    output_gff_file = os.path.splitext(input_file)[0] + '.GFF2'
    with open(output_gff_file, 'w') as m:
        m.write('##gff-version 2\n')
        m.writelines('\t'.join(map(str,map_dict[k])) for k in map_dict.keys())
    return

fasta_file = sys.argv[1]
bt2_index = sys.argv[2]

map_bowtie(fasta_file, bt2_index)
