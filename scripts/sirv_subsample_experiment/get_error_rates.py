#!/usr/bin/env python

import argparse
import sys

from collections import Counter
import re

import edlib



'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        elif type_ == 'S':
            print("error", file=sys.stderr)
            print(cigar, file=sys.stderr)
            sys.exit()
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


def get_error_profile(corr_seq, true_seq):
    res = edlib.align(corr_seq, true_seq, task="path", mode="NW")
    tot = res["editDistance"]
    cigar_string = res["cigar"]
    read_alignment, ref_alignment = cigar_to_seq(cigar_string, corr_seq, true_seq)
    insertions = ref_alignment.count("-")
    deletions = read_alignment.count("-")
    # indels =  insertions + deletions
    subs = len([1 for n1, n2 in zip(read_alignment, ref_alignment) if n1 != n2 and n1 != "-" and n2 != "-"] )
    err_rate = round(100 * float(subs + insertions + deletions) / len(true_seq), 2)
    return tot, insertions, deletions, subs, err_rate




def decide_primary_locations(sam_file, args): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    error_rates_per_read = {}
    reads_multiple_primary = set()
    reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            insertions = sum([length for type_, length in read.cigartuples if type_ == 1])
            deletions = sum([length for type_, length in read.cigartuples if type_ == 2 ])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            tot_align = insertions + deletions + subs + matches
            identity = matches/float(tot_align)
            err_rate = round(100 * float(subs + insertions + deletions) / len(true_seq), 2)
            error_rates_per_read[read.quary_name] = (insertions, deletions, subs, matches, err_rate)
    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return error_rates_per_read     



def get_error_rate_stats_per_read(reads_primary_locations, reads, annotated_splice_coordinates_pairs, args, reference = {}):
    # SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    alignments = {}
    alignments_detailed = {}
    read_index = 0
    all_called_intron_lengths = []
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # print(read.is_reverse)
            read_index += 1

            if args.align:
                read_seq = reads[read.query_name]
                ref_seq = reference[read.reference_name]
                ref_alignment, m_line, read_alignment, ins, del_, subs, matches = cigar_to_seq_mm2_local(read, ref_seq, read_seq)
                alignments_detailed[read.query_name] = (ref_alignment, m_line, read_alignment)

                # if read.query_name == 'b1d0ee62-7557-4645-8d1e-c1ccfb60c997_runid=8c239806e6f576cd17d6b7d532976b1fe830f9c6_sampleid=pcs109_sirv_mix2_LC_read=29302_ch=69_start_time=2019-04-12T22:47:30Z_strand=-':
                #     print(read.query_alignment_start)
                #     print(read.query_name, "primary", read.flag, read.reference_name) 
                #     print(ref_alignment)
                #     print(m_line)
                #     print(read_alignment)
                #     print(ins, del_, subs, matches)
            else:
                # print(read.cigartuples)

                ins = sum([length for type_, length in read.cigartuples if type_ == 1])
                del_old = sum([length for type_, length in read.cigartuples if type_ == 2 and length <  args.min_intron])
                subs = sum([length for type_, length in read.cigartuples if type_ == 8])
                matches = sum([length for type_, length in read.cigartuples if type_ == 7])
                called_intron_lengths = [length for type_, length in read.cigartuples if type_ == 3]


                # Deletion is treated separately because of potential short introns -- if in annotation
                if read.reference_name in annotated_splice_coordinates_pairs:
                    del_new = 0
                    chr_coordinate_pairs = annotated_splice_coordinates_pairs[read.reference_name]
                    read_deletion_coordinate_pairs = get_deletion_sites(read)
                    for start_del, stop_del, del_length in read_deletion_coordinate_pairs:
                        if (start_del, stop_del) not in chr_coordinate_pairs:
                            del_new += del_length
                        else:
                            print("True intron masked as deletion, length:", del_length, read.reference_name, start_del, stop_del)

                else:
                    del_new = sum([length for type_, length in read.cigartuples if type_ == 2 and length <  args.min_intron])

                # print(del_new, del_old)
                # if del_new > del_old:
                #     print("New:",del_new, "Old:", del_old) #, [length for type_, length in read.cigartuples if type_ == 2], read.cigarstring, read_deletion_coordinate_pairs)
                    # print([length for type_, length in read.cigartuples if type_ == 2 and length >=  20])

            all_called_intron_lengths.append(called_intron_lengths)
            # if [length for type_, length in read.cigartuples if type_ == 2 and length >  20]:
            #     print([length for type_, length in read.cigartuples if type_ == 2 and length >  20])
            #     print(get_deletion_coordiantes(read))
            # if tot_ins != ins or tot_subs != subs or tot_del != del_ or tot_match != matches:
            #     print(tot_ins, ins, tot_subs, subs, tot_del, del_, tot_match, matches)
            
            # if read.query_name in alignments:
            #     print("read", read_index, alignments[read.query_name], read.reference_name)
            #     print("New:", (ins, del_, subs, matches))
            #     print(read.flag, read.query_name)
            #     print("BUG")
            #     # sys.exit()
            assert read.query_name not in alignments
            alignments[read.query_name] = (ins, del_new, subs, matches, read.reference_name, read.reference_start, read.reference_end + 1, read.flag, read_index)

        else:
            pass
            # print("secondary", read.flag, read.reference_name) 
    # SAM_file.close()

    # multiple_alingment_counter = 0
    # for read in alignments:
    #     if len(alignments[read]) > 1:
    #         # print(alignments[read])
    #         multiple_alingment_counter += 1
    # print("Nr reads with multiple primary alignments:", multiple_alingment_counter)
    print( '100 smallest introns as N from the aligner:', sorted([l for lst in all_called_intron_lengths for l in lst])[:100] )
    # sys.exit()
    return alignments, alignments_detailed


def main(args):
    refs = { acc : seq for acc, (seq,_) in readfq(open(args.refs, 'r'))}
    # reads = { acc : seq for acc, (seq,qual) in readfq(open(args.reads, 'r'))}
    print("read,err,subs,ins,del,err_rate")
    decide_primary_locations(args.samfile, args)
    for acc in reads:
        transcript_id = acc.split("_")[0]

        # res = edlib.align(seq, spoa_ref, task="path", mode="NW")
        # cigar_string = res["cigar"]
        # read_alignment, ref_alignment = help_functions.cigar_to_seq(cigar_string, seq, spoa_ref)

        corr_seq = corrected[acc]
        tot, ins, del_, subs, err_rate = get_error_profile(corr_seq, true_seq)
        print("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}".format(acc, tot, ins, del_, subs, err_rate, "corrected", isoform_coverage[transcript_id], gene_coverage[gene_id], gene_fam_coverage[gene_fam_id] ))
        # print("corrected:", ed_corr, "cov:", isoform_coverage[transcript_id])
        err_rate_corr = err_rate
        orig_seq = original[acc]
        
        if err_rate_corr > 15:
            print("Over corrected!", file=sys.stderr)
            print(acc, file=sys.stderr)
            print(gene_coverage[gene_id], gene_fam_coverage[gene_fam_id], file=sys.stderr )
            print(corr_seq, tot, ins, del_, subs, err_rate, file=sys.stderr)
            print(true_seq, file=sys.stderr)
            tot, ins, del_, subs, err_rate =  get_error_profile(orig_seq, true_seq)
            print(orig_seq,tot, ins, del_, subs, err_rate, file=sys.stderr)
            print("", file=sys.stderr)
            

        tot, ins, del_, subs, err_rate = get_error_profile(orig_seq, true_seq)
        print("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}".format(acc, tot, ins, del_, subs, err_rate, "original", isoform_coverage[transcript_id], gene_coverage[gene_id], gene_fam_coverage[gene_fam_id] ))

        # res = edlib.align(orig_seq, true_seq, task="path", mode="NW")
        # ed_orig = res["editDistance"]
        # print("Original:", ed_orig, "cov:", isoform_coverage[transcript_id])





if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    parser.add_argument('refs', type=str, help='Fastq file. ')
    parser.add_argument('samfile', type=str, help='fastq file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

