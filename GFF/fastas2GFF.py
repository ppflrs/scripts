from Bio import SeqIO
from Bio.Seq import Seq

def map_bowtie(f, index_name):
    handle = open(f, 'rU')
    dna = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    map_dict = {}
    counter = 1 #To comply with the GFF format file specifications needs to start at 1
    cat_genes = open(f + '.seq', 'a')
    cat_genes.write('>' + index_name + '\n') #FASTA header

    seqid_counter = 1 #To order the output of the sequence

    for sequence in dna:
        gene_len = len(str(sequence.seq))

        description = list()
        description.append(sequence.id) #id
        description.append(index_name) #bowtie index name
        description.append(counter) #start
        description.append(counter + gene_len) #end
        description.append('CDS') #Feature type: This can be exon, promotor, etc for us CDS it's ok i guess
        description.append('.') #score ???
        description.append('+') #strand
        description.append('ID=CDS'+str(seqid_counter).zfill(4)) #CDS0001, CDS0002, etc
        description.append('\n')

        map_dict[seqid_counter] = description
        seqid_counter += 1
        counter += (gene_len)
        cat_genes.write(str(sequence.seq))

    with open(f + '.map', 'w') as m:
        m.write('##gff-version 3\n')
        m.writelines('\t'.join(map(str,map_dict[k])) for k in map_dict.keys())
    return
