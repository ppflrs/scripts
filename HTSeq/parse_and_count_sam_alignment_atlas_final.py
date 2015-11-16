#!/usr/bin/env python
"""
requires HTseq (http://www-huber.embl.de/users/anders/HTSeq/)
parse sam/bam alignment and extract reads/hits by quality or gene they were aligned to.
Counting is adopted from HTseq-count script and uses the "Union" mode only.
requires SAM file to have an MD string (either inserted by aligner or added using samtools calmd)
and (http://code.google.com/p/pysam/).
"""

import os
import sys
import re
import itertools
import argparse

import HTSeq  # from where, on atlas?


"""
UnknownChrom obtained from https://github.com/luispedro/htseq-copy/blob/master/HTSeq/scripts/count.py
"""


class UnknownChrom(Exception):
    pass


current_dir = os.getcwd()


def main():
    exe_parser = argparse.ArgumentParser()
    exe_parser.add_argument('infile', type=str, help='<input file> [(full path), -b/-s required]')
    exe_parser.add_argument("-u", "--not_aligned",
                            help="output reads that were not aligned, including those that were aligned multiple times(flat file).",
                            type=str)
    exe_parser.add_argument("-s", "--samout", help="output not aligned reads to [file path].", type=str)
    exe_parser.add_argument("-b", "--ambiguous_out", help="output a fasta file of ambiguous hits [file path].",
                            type=str)
    exe_parser.add_argument("-v", "--verbose", help="verbose. (default = TRUE).", action="store_true")
    exe_parser.add_argument("gff", help="<gff file> [(full path)]", type=str)
    exe_parser.add_argument("-f", "--fasta", help="output fasta file of hits (full path).", type=str)
    exe_parser.add_argument("-m", "--min_read_length", help="minimal read length to consider. (default = 60b).",
                            type=int)
    exe_parser.add_argument("-i", "--min_id", help="minimal percent id of hit to consider. (default = 80).", type=int)
    exe_parser.add_argument("-z", "--min_score", help="minimal aligner score to consider. (default = 0).", type=int)
    exe_parser.add_argument("-c", "--max_clip",
                            help="proportion of bases clipped from read for alignment. (default = 0.3).", type=float)
    exe_parser.add_argument("--stranded", help="whether the data is stranded (y, n, reverse). (default = n).", type=str,
                            choices=["y", "n", "reverse"], default="n")
    exe_parser.add_argument("--idattr", help="GFF attribute to be used as feature ID. (default = GeneID).", type=str)
    exe_parser.add_argument("--type", help="feature type (3rd column in GFF file) to be used. (default = CDS).",
                            type=str)
    exe_parser.add_argument("-a", "--minaqual", help="min. alignment quality (default = 0).", type=str)
    exe_parser.add_argument("-p", "--paired_end_mode",
                            help="input is paired end sorted by name (n) or position (p) . (default = p).", type=str,
                            choices=["p", "n"], default="p")
    exe_parser.add_argument("-o", "--out", help="name of counts output file.", type=str)
    args = exe_parser.parse_args()

    if args.paired_end_mode == 'p':
        paired_end = True
        pe_order = 'p'
    elif args.paired_end_mode == 'n':
        paired_end = True
        pe_order = 'n'

    if args.infile:
        try:
            if args.infile == '-':  # get sam on a stream
                seqfile = HTSeq.SAM_Reader(sys.stdin)
                if args.paired_end_mode:
                    # read_seq_iter = iter(seqfile)
                    # first_read = read_seq_iter.next()
                    # read_seq = itertools.chain([first_read], read_seq_iter)
                    # reader = HTSeq.pair_SAM_alignments(read_seq)
                    if pe_order == 'p':
                        reader = HTSeq.pair_SAM_alignments_with_buffer(seqfile)
                    elif pe_order == 'n':
                        reader = HTSeq.pair_SAM_alignments(seqfile)  # (read_seq)
                else:
                    reader = seqfile
            elif args.infile != '-':
                seqfile = HTSeq.SAM_Reader(args.infile)
                if args.paired_end_mode:
                    read_seq_iter = iter(seqfile)
                    first_read = read_seq_iter.next()
                    read_seq = itertools.chain([first_read], read_seq_iter)
                    reader = HTSeq.pair_SAM_alignments(read_seq)
                    if pe_order == 'p':
                        reader = HTSeq.pair_SAM_alignments_with_buffer(reader)
                    elif pe_order == 'n':
                        reader = HTSeq.pair_SAM_alignments(reader)
                else:
                    reader = seqfile
                    # fread_seq_iter = iter(reader)
                    # first_read = iter(read_seq).next()
            elif args.infile == '':
                print "no input file type given. exiting..."
                sys.exit(1)
        except:
            print "failed processing SAM/BAM file"
            raise
    elif not args.infile:
        print "no input file given. exiting..."
        sys.exit(1)

    if args.gff:
        gff_file = args.gff
    else:
        print "no gff file given. exiting..."
        sys.exit(1)

    if args.verbose:
        verbose = True
    else:
        verbose = False

    if args.min_read_length:
        min_read_len = args.min_read_length
    else:
        min_read_len = 60  # default read length

    if args.max_clip:
        max_clip_ = float(args.max_clip)
    else:
        max_clip_ = float(0.3)  # default read length

    if args.min_id:
        min_id = float(args.min_id)
    else:
        min_id = float(80)

    if args.min_score:
        min_score = int(args.min_score)
    else:
        min_score = 0

    if args.stranded == 'n':
        stranded = 'no'
    elif args.stranded == 'y':
        stranded = 'yes'
    elif args.stranded == 'reverse':
        stranded = 'reverse'

    if args.minaqual:
        minaqual = args.minaqual
    else:
        minaqual = 0

    if args.idattr:
        id_attribute = args.idattr
    else:
        id_attribute = "GeneID"
    if args.type:
        feature_type = args.type
    else:
        feature_type = 'CDS'

    # ###
    # parse GFF file
    features, counts = gff_reader(gff_file, feature_type, id_attribute, verbose, stranded)
    # ###
    if args.samout:
        samoutfile = open(args.samout, "w")
    else:
        samoutfile = None
    if args.ambiguous_out:
        ambiguousfile = open(args.ambiguous_out, "w")
    else:
        ambiguousfile = None
    if args.fasta:
        fastafile = open(args.fasta, "w")
    else:
        fastafile = None
    if args.not_aligned:
        not_aligned_file = open(args.not_aligned, "w")
    else:
        not_aligned_file = None
    if args.out:
        outfile = open(args.out, "w")
    else:
        outfile = None

        # if outfile and samoutfile and  ambiguousfile and fastafile and not_aligned_file == None:
        # print "None of the possible output file options specified. exiting..."
        # sys.exit(1)
    # #######
    # decalre counter variables
    empty = 0
    ambiguous = 0
    notaligned = 0
    lowqual = 0
    nonunique = 0
    # #######

    read_counter = 0
    for alignment in reader:  # for alignment entry (line in fact) in sam file
        # iv_seq
        # print alignment
        if not paired_end:
            if read_counter % 1000000 == 0 and verbose:
                if verbose:
                    print read_counter, 'non paired-end alignments processed'
            read_name = alignment.read.name
            # read = alignment.read  # READ. Note that def invert_strand( iv ):
            read_seq = alignment.read.seq
            read_length = len(alignment.read.seq)
            if not alignment.aligned:  # check if read is aligned to ref sequence
                if alignment is not None:
                    notaligned += 1
                    if args.samout:
                        write_to_samout(samoutfile, paired_end, alignment, "not_aligned")
                    if args.not_aligned:
                        not_aligned_file.write(read_name + '\t' + 'not_aligned' + '\n')
                        # continue
            elif alignment.aligned:

                opt_fields = alignment.optional_fields
                # flag = alignment.flag
                cigar_string = parse_cigar(alignment.original_sam_line.split('\t')[
                    5])  # just the cigar string without the fancy HTseq additions
                cigar_soft_clipped, cigar_m, cigar_insertions, cigar_deletions, cigar_insertions = parse_cigar_alignment(cigar_string)  # get alignment data from cigar string
                score, md_matches, md_deletions, md_mismatches = parse_opt_fields(
                    opt_fields)  # get alignment data from md string
                percent_id = 100.0 * (
                    float(md_matches) / (float(read_length - cigar_soft_clipped + cigar_insertions + cigar_deletions)))
                if alignment[0] is not None:  # check if read is aligned to ref sequence
                    if alignment.optional_field("NH") > 1:  # check if read is mapped more than once
                        # By default these reads are discarded. CHANGE?
                        if args.samout:
                            write_to_samout(samoutfile, paired_end, alignment, "alignment_not_unique")
                        nonunique += 1
                        if args.not_aligned:
                            not_aligned_file.write(read_name + '\t' + 'alignment_not_unique' + '\n')
                            # continue
                    if alignment.aQual < minaqual:  # check quality. default is 0
                        lowqual += 1
                        if args.samout:
                            write_to_samout(samoutfile, paired_end, alignment, "too_low_aQual")
                        if args.not_aligned:
                            not_aligned_file.write(read_name + '\t' + 'too_low_aQual' + '\n')
                            # continue
                    clipped = (float(cigar_soft_clipped) / float(read_length))
                    if read_length >= min_read_len:
                        if (float(cigar_soft_clipped) / float(read_length)) <= max_clip_:
                            if score >= args.min_score:
                                if percent_id >= float(min_id):
                                    if stranded == "reverse":
                                        iv_seq = (
                                            (invert_strand(cigar_operation.ref_iv) for cigar_operation in
                                             alignment[1].cigar
                                             if cigar_operation.type == "M" and cigar_operation.size > 0))
                                    else:
                                        iv_seq = (cigar_operation.ref_iv for cigar_operation in alignment.cigar if
                                                  cigar_operation.type == "M" and cigar_operation.size > 0)
                                    iv_seq_good = True
                                    # collects hits to chromosomes/features.
                                    """
                                    cigarOperation in HTSeq:
                                    HTSeq.parse_cigar( "20M6I10M", 1000, "chr2", "+" ) #ref_iv == genomicInterval object
                                    of htSeq
                                    [< CigarOperation: 20 base(s) matched on ref iv chr2:[1000,1020)/+,query iv[0,20)>,
                                    < CigarOperation: 6 base(s) inserted on ref iv chr2:[1020,1020)/+,query iv[20,26)>,]
                                    """
                                    # if args.fasta:
                                    # fastafile.write('>' + read_name + '\n' + read_seq + '\n')

                                else:
                                    iv_seq_good = False
                                    if args.samout:
                                        write_to_samout(samoutfile, paired_end, alignment,
                                                        "percent_id_too_low=" + str(percent_id))
                                    if args.not_aligned:
                                        not_aligned_file.write(
                                            read_name + '\t' + 'percent_id_too_low=' + str(percent_id) + '\n')
                            else:
                                iv_seq_good = False
                                if args.samout:
                                    write_to_samout(samoutfile, paired_end, alignment,
                                                    'alignment_score_too_low=' + str(score))
                                if args.not_aligned:
                                    not_aligned_file.write(
                                        read_name + '\t' + 'alignment_score_too_low=' + str(score) + '\n')
                        else:
                            iv_seq_good = False
                            if args.samout:
                                write_to_samout(samoutfile, paired_end, alignment,
                                                'too_many_bases_clipped_from_read=' + str(cigar_soft_clipped))
                            if args.not_aligned:
                                not_aligned_file.write(read_name + '\t' + 'too_many_bases_clipped_from_read=' + str(
                                    cigar_soft_clipped) + '\n')
        elif paired_end:
            # print "read counter=", read_counter
            if read_counter % 100000 == 0 and verbose:
                if verbose:
                    print read_counter, 'alignment pairs processed'
            if (alignment[0] is None) or not alignment[0].aligned:
                notaligned += 1
                try:
                    read_1_name = alignment[0].read.name
                except:
                    read_1_name = 'None'
                if args.samout:
                    write_to_samout(samoutfile, paired_end, alignment, "not_aligned")
                if args.not_aligned:
                    not_aligned_file.write(read_1_name + '\t' + 'not_aligned' + '\n')
            elif (alignment[1] is None) or not alignment[1].aligned:
                notaligned += 1
                try:
                    read_2_name = alignment[1].read.name
                except:
                    read_2_name = 'None'
                if args.samout:
                    write_to_samout(samoutfile, paired_end, alignment, "not_aligned")
                if args.not_aligned:
                    not_aligned_file.write(read_2_name + '\t' + 'not_aligned' + '\n')
            else:
                # else:
                read_1_name = alignment[0].read.name
                # read_1 = alignment[0].read  #READ.
                read_1_length = len(alignment[0].read.seq)
                read_1_seq = alignment[0].read.seq
                read_2_name = alignment[1].read.name
                # read_2 = alignment[1].read  #READ.
                # read_2_length = len(alignment[1].read.seq)
                read_2_seq = alignment[1].read.seq
                iv_seq = tuple()
                if (alignment[0] is not None) and alignment[0].aligned:  # check if read is aligned to ref sequence
                    opt_1_fields = alignment[0].optional_fields
                    # flag_1 = alignment[0].flag
                    cigar_1_string = parse_cigar(alignment[0].original_sam_line.split('\t')[
                        5])  # just the cigar string without the fancy HTseq additions
                    cigar_1_soft_clipped, cigar_1_m, cigar_1_insertions, cigar_1_deletions, cigar_1_insertions = parse_cigar_alignment(
                        cigar_1_string)
                    score_1, md_1_matches, md_1_deletions, md_1_mismatches = parse_opt_fields(
                        opt_1_fields)  # get alignment data from md string
                    percent_1_id = (100.0 * ((float(md_1_matches) / (
                        float(read_1_length - cigar_1_soft_clipped + cigar_1_insertions + cigar_1_deletions)))))
                    clipped_1 = (float(cigar_1_soft_clipped) / float(read_1_length))
                    if int(read_1_length) >= int(min_read_len):
                        if (float(cigar_1_soft_clipped) / float(read_1_length)) <= float(max_clip_):

                            # if int(score_1) >= int(args.min_score):
                            if int(score_1) >= int(min_score):
                                # if float(percent_1_id) >= float(args.min_id):
                                if float(percent_1_id) >= float(min_id):
                                    if stranded == "reverse":
                                        iv_seq = itertools.chain(iv_seq, (invert_strand(cigar_operation.ref_iv) for
                                                                          cigar_operation in alignment[0].cigar if
                                                                          cigar_operation.type == "M" and cigar_operation.size > 0))
                                    else:
                                        iv_seq = itertools.chain(iv_seq, (cigar_operation.ref_iv for cigar_operation in
                                                                          alignment[0].cigar if
                                                                          cigar_operation.type == "M" and cigar_operation.size > 0))
                                    # if args.fasta:
                                    # fastafile.write('>' + read_1_name + '\n' + read_1_seq + '\n')
                                    iv_seq_good_1 = True

                                else:
                                    iv_seq_good_1 = False
                                    if args.samout:
                                        write_to_samout(samoutfile, paired_end, alignment,
                                                        "percent_id_too_low=" + str(percent_1_id))
                                    if args.not_aligned:
                                        not_aligned_file.write(
                                            read_1_name + '\t' + 'percent_id_too_low=' + str(percent_1_id) + '\n')
                            else:
                                iv_seq_good_1 = False
                                if args.samout:
                                    write_to_samout(samoutfile, paired_end, alignment,
                                                    'alignment_score_too_low=' + str(score_1))
                                if args.not_aligned:
                                    not_aligned_file.write(
                                        read_1_name + '\t' + 'alignment_score_too_low=' + str(score_1) + '\n')
                        else:
                            iv_seq_good = False
                            if args.samout:
                                write_to_samout(samoutfile, paired_end, alignment,
                                                'too_many_bases_clipped_from_read=' + str(cigar_1_soft_clipped))
                            if args.not_aligned:
                                not_aligned_file.write(read_1_name + '\t' + 'too_many_bases_clipped_from_read=' + str(
                                    cigar_1_soft_clipped) + '\n')
                # else:
                # iv_seq = tuple()

                if (alignment[1] is not None) and alignment[1].aligned:  # check if read is aligned to ref sequence
                    opt_2_fields = alignment[1].optional_fields
                    # flag_2 = alignment[1].flag  # ',  #'bit_length', 'conjugate', 'denominator', 'imag', 'numerator', 'real']
                    cigar_2_string = parse_cigar(alignment[1].original_sam_line.split('\t')[
                        5])  # just the cigar string without the fancy HTseq additions
                    cigar_2_soft_clipped, cigar_2_m, cigar_2_insertions, cigar_2_deletions, cigar_2_insertions = parse_cigar_alignment(
                        cigar_2_string)
                    score_2, md_2_matches, md_2_deletions, md_2_mismatches = parse_opt_fields(
                        opt_2_fields)  # get alignment data from md string
                    read_2_name = alignment[1].read.name
                    read_2_length = len(alignment[1].read.seq)
                    # read_2 = alignment[1].read  # READ.
                    read_2_seq = alignment[1].read.seq
                    percent_2_id = (100.0 * (float(md_2_matches) / (
                        float(read_2_length - cigar_2_soft_clipped + cigar_2_insertions + cigar_2_deletions))))
                    clipped_2 = (float(cigar_2_soft_clipped) / float(read_2_length))
                    if int(read_2_length) >= int(min_read_len):
                        if (float(cigar_2_soft_clipped) / float(read_2_length)) <= float(max_clip_):
                            if int(score_2) >= int(min_score):
                                if float(percent_2_id) >= float(min_id):
                                    if stranded == "reverse":
                                        iv_seq = itertools.chain(iv_seq, (invert_strand(cigar_operation.ref_iv) for
                                                                          cigar_operation in alignment[1].cigar if
                                                                          cigar_operation.type == "M" and cigar_operation.size > 0))
                                    else:
                                        iv_seq = itertools.chain(iv_seq, (cigar_operation.ref_iv for cigar_operation in
                                                                          alignment[1].cigar if
                                                                          cigar_operation.type == "M" and cigar_operation.size > 0))
                                        iv_seq_good_2 = True
                                    try:
                                        if (alignment[0].optional_field("NH") > 1) or (alignment[1].optional_field(
                                                "NH") > 1):
                                            # or (alignment[1].optional_field("NH") > 1): #check if read is mapped more
                                            # than once
                                            # By default these reads are discarded. CHANGE?
                                            iv_seq_good_1 = False
                                            iv_seq_good_2 = False
                                            if args.samout:
                                                write_to_samout(samoutfile, paired_end, alignment,
                                                                "alignment_not_unique")
                                                nonunique += 1
                                            if args.not_aligned:
                                                not_aligned_file.write(read_1_name + '\t' + 'not_aligned' + '\n')
                                                not_aligned_file.write(read_2_name + '\t' + 'not_aligned' + '\n')
                                                continue
                                    except KeyError:
                                        pass
                                    if (alignment[0] and alignment[0].aQual < minaqual) or (alignment[1] and alignment[1].aQual < minaqual):
                                        # check quality. default is 0
                                        iv_seq_good_2 = False
                                        lowqual += 1
                                        if args.samout:
                                            write_to_samout(samoutfile, paired_end, alignment, "too_low_aQual")
                                        if args.not_aligned:
                                            not_aligned_file.write(read_1_name + '\t' + 'not_aligned' + '\n')
                                            not_aligned_file.write(read_2_name + '\t' + 'not_aligned' + '\n')
                                        continue
                                else:
                                    iv_seq_good_2 = False
                                    if args.samout:
                                        write_to_samout(samoutfile, paired_end, alignment,
                                                        "percent_id_too_low=" + str(percent_2_id))
                                    if args.not_aligned:
                                        not_aligned_file.write(
                                            read_2_name + '\t' + 'percent_id_too_low=' + str(percent_2_id) + '\n')
                            else:
                                iv_seq_good_2 = False
                                if args.samout:
                                    write_to_samout(samoutfile, paired_end, alignment,
                                                    'alignment_score_too_low=' + str(score_2))
                                if args.not_aligned:
                                    not_aligned_file.write(
                                        read_2_name + '\t' + 'alignment_score_too_low=' + str(score_2) + '\n')
                        else:
                            iv_seq_good_2 = False
                            if args.samout:
                                write_to_samout(samoutfile, paired_end, alignment,
                                                'too_many_bases_clipped_from_read=' + str(cigar_2_soft_clipped))
                            if args.not_aligned:
                                not_aligned_file.write(read_2_name + '\t' + 'too_many_bases_clipped_from_read=' + str(
                                    cigar_2_soft_clipped) + '\n')
        read_counter += 1

        """
        overlap_mode == "union"
        will count a hit even if read is mapped across an intron or there is an insertion.
        """
        try:
            feature_set = set()
            for iv in iv_seq:
                # print iv
                if iv.chrom not in features.chrom_vectors:  # check if alignment feaure name in features from GFF file
                    # The name of a sequence (i.e., chromosome, contig, or the like).
                    # check the gff features dictionary
                    raise UnknownChrom
                for iv2, fs2 in features[iv].steps():  # fs == feature steps.
                    """
                    from HTseq manual:
                    GenomicArray objects use by default so-called StepVectors that store the data internally in steps of
                    constant value
                    """
                    feature_set = feature_set.union(fs2)
                    # print feature_set
            if feature_set is None or len(feature_set) == 0:
                if args.samout:
                    write_to_samout(samoutfile, paired_end, alignment, "no_feature")
                if args.not_aligned:
                    not_aligned_file.write('None' + '\t' + 'no_feature' + '\n')
                empty += 1
            elif len(feature_set) > 1:
                if args.samout:
                    write_to_samout(samoutfile, paired_end, alignment, "ambiguous[" + '+'.join(feature_set) + "]")
                if ambiguousfile:
                    if paired_end:
                        if iv_seq_good_1:
                            ambiguousfile.write('>' + read_1_name + '_' + "ambiguous[" + '+'.join(
                                feature_set) + "]" + '_clipped_' + str(clipped_1) + '_score_' + str(score_2) + '_percent_id_' + str(percent_1_id) + '\n' + read_1_seq + '\n')
                        if iv_seq_good_2:
                            ambiguousfile.write('>' + read_2_name + '_' + "ambiguous[" + '+'.join(
                                feature_set) + "]" + '_clipped_' + str(clipped_2) + '_score_' + str(score_2) + '_percent_id_' + str(percent_2_id) + '\n' + read_2_seq + '\n')
                    else:
                        if iv_seq_good:
                            ambiguousfile.write('>' + alignment.read.name + '_' + "ambiguous[" + '+'.join(
                                feature_set) + "]" + '_clipped_' + str(clipped) + '_score_' + str(score) + '_percent_id_' + str(percent_id) + '\n' + read_seq + '\n')

                """
                #if args.not_aligned:
                #    if paired_end:
                #    not_aligned_file.write(alignment[0].read.name + '\t' + 'ambiguous['+'+'.join(feature_set)+']' + '\n')
                #        not_aligned_file.write(alignment[1].read.name + '\t' + 'ambiguous['+'+'.join(feature_set)+']' + '\n')
                #    else:
                #    not_aligned_file.write(alignment.read.name + '\t' + 'ambiguous['+'+'.join(feature_set)+']' + '\n')
                """
                ambiguous += 1
            elif len(feature_set) == 1:
                if args.samout:
                    write_to_samout(samoutfile, paired_end, alignment, list(feature_set)[0])
                if args.fasta:
                    if paired_end:
                        if iv_seq_good_1:
                            fastafile.write('>' + read_1_name + '_' + ''.join(list(feature_set)[0]) + '_clipped_' + str(
                                clipped_1) + '_score_' + str(score_1) + '_percent_id_' + str(percent_1_id) + '\n' + read_1_seq + '\n')
                        if iv_seq_good_2:
                            fastafile.write('>' + read_2_name + '_' + ''.join(list(feature_set)[0]) + '_clipped_' + str(
                                clipped_2) + '_score_' + str(score_2) + '_percent_id_' + str(percent_2_id) + '\n' + read_2_seq + '\n')
                    else:
                        if iv_seq_good:
                            fastafile.write('>' + read_name + '_' + ''.join(list(feature_set)[0]) + '_clipped_' + str(
                                clipped) + '_score_' + str(score) + '_percent_id_' + str(percent_id) + '\n' + read_seq + '\n')

                counts[list(feature_set)[0]] += 1
        except:
            if args.samout:
                write_to_samout(samoutfile, paired_end, alignment, "__no_feature")
            empty += 1

            # if not paired_end:
            # al = alignment
            # else:
            # al = alignment[0] if alignment[0] is not None else alignment[1]

            # if args.not_aligned:
            # not_aligned_file.write(al.read.name + '\t' + 'feature_not_in_gff_file' + '\n')
            # if not verbose:
            #    print (("Warning: Skipping read '%s', because chromosome " +
            #    "'%s', to which it has been aligned, did not appear in the GFF file.\n" ) %
            #     (al.read.name, iv.chrom) )
    print 'total', read_counter, 'alignments processed'
    if samoutfile is not None:
        samoutfile.close()
    if fastafile is not None:
        fastafile.close
    if not_aligned_file is not None:
        not_aligned_file.close()

    if outfile is not None:
        for feature in sorted(counts.keys()):
            outfile.write("%s\t%d\n" % (feature, counts[feature]))
        outfile.write("no_feature\t%d\n" % empty)
        outfile.write("ambiguous\t%d\n" % ambiguous)
        outfile.write("too_low_aQual\t%d\n" % lowqual)
        outfile.write("not_aligned\t%d\n" % notaligned)
        outfile.write("alignment_not_unique\t%d\n" % nonunique)
    if outfile is not None:
        outfile.close()


def parse_cigar_alignment(cigar_string):
    """
    parse cigar string to get alignment data: total number of insertions, clipped bases, deletions
    """
    cigar_soft_clipped = 0
    cigar_m = 0
    cigar_deletions = 0
    cigar_insertions = 0
    for c in cigar_string:  # ('5S, 12M, 3D, 4M, 1I, 45M, 3S')
        if 'S' in c:
            c = c.replace('S', '')
            cigar_soft_clipped += int(c)
        elif 'M' in c:
            c = c.replace('M', '')
            cigar_m += int(c)
        elif 'D' in c:
            c = c.replace('D', '')
            cigar_deletions += int(c)
        elif 'I' in c:
            c = c.replace('I', '')
            cigar_insertions += int(c)
    # print 'cigar_soft_clipped -->', cigar_soft_clipped, 'cigar_M -->', cigar_M, 'cigar_deletions -->',
    # cigar_deletions, 'cigar_insertion -->', cigar_insertions
    # print 'cigar total -->', cigar_soft_clipped + cigar_M + cigar_deletions + cigar_insertions
    return cigar_soft_clipped, cigar_m, cigar_insertions, cigar_deletions, cigar_insertions


def parse_opt_fields(opt_fields):
    """
    parse optional fields column in sam file.
    get MD string alignment data.
    get score --> AS
    """
    score = int()
    md_matches = 0
    md_deletions = 0
    md_mismatches = 0
    for field in opt_fields:
        if 'MD' in field[0]:
            # print 'MD -->', field[1],
            md_list = re.split('(\D|\W)', field[1])
            for i in md_list:
                if i:
                    if i.isdigit():
                        md_matches += int(i)
                    elif "^" in i:
                        md_deletions += 1
                    else:
                        md_mismatches += 1
                        # total_md = md_mismatches + md_deletions + md_matches
        if 'AS' in field[0]:
            score = field[1]
    return score, md_matches, md_deletions, md_mismatches


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError, "Illegal strand"
    return iv2


def parse_cigar(s):
    s = s  # "33S29M83S"
    l = list()
    ll = list()
    for i in s:
        if i.isdigit():
            ll.append(i)
        else:
            ll.append(i)
            l.append(''.join(ll))
            ll = list()
    return l


def write_to_samout(samoutfile, paired_end, alignment, assignment):
    if samoutfile is None:
        return
    if not paired_end:
        alignment = (alignment,)
    for read in alignment:
        if read is not None:
            samoutfile.write(read.original_sam_line.rstrip() + "\tXF:Z:" + assignment + "\n")


def gff_reader(gff_file, feature_type, id_attribute, verbose, stranded):
    """
    based on HTseq-count documentation and script
    """
    gff = HTSeq.GFF_Reader(gff_file)
    # 'gene' 'exon' 'CDS' etc. 3rd column in the gff file
    feature_type = feature_type
    # e.g 'gene_id'. the attributes (9th) column of the gff file.user inputs desired id_attribute
    id_attribute = id_attribute
    verbose = verbose
    stranded = stranded
    counts = dict()
    # features = dict()
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    counter = 0
    try:
        for feature in gff:
            if feature.type == feature_type:  # if the feature is the feature we are looking for.
                if id_attribute in feature.attr:
                    feature_id = feature.attr[id_attribute]
                else:
                    if verbose:
                        print feature, "--> the requested attribute does not exist in feature"
                        if stranded == "yes" and feature.iv.strand == ".":
                            if verbose:
                                print feature, '--> does not contain strand info.'
                # get features for iv ( genomic interval)
                features[feature.iv] += feature_id
                counts[feature.attr[id_attribute]] = 0
                counter += 1

        if counter % 1000 == 0:
            if verbose:
                print counter, "features found in gff file"
    except:
        print 'non of the desired attributes found in gff file. exiting...'
        # sys.exit(1)
        print gff
        sys.stderr.write("Error occured when processing GFF file (%s):\n" % gff.get_line_number_string())
        raise
    if verbose:
        print counter, "features found in gff file"
    return features, counts


main()
