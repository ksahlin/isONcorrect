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
        name, seqs, last = last[1:].split()[0], [], None #.replace(" ", "_")
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
    matches = len([1 for n1, n2 in zip(read_alignment, ref_alignment) if n1 == n2 and n1 != "-" and n2 != "-"] )
    err_rate = round(100 * float(subs + insertions + deletions) / (len(true_seq) + insertions), 2)
    read_length = len(true_seq)
    return tot, insertions, deletions, subs, matches, err_rate


def main(args):
    true = { acc.split("|")[2] : seq for acc, (seq,_) in readfq(open(args.true, 'r'))}
    corrected = { acc : seq for acc, (seq,qual) in readfq(open(args.corrected, 'r'))}
    original = { acc : seq for acc, (seq,qual) in readfq(open(args.original, 'r'))}

    isoform_coverage = Counter([acc.split("|")[2].split("_")[0] for acc in original])
    gene_coverage = Counter([acc.split("|")[1] for acc in original])
    gene_fam_coverage = Counter([acc.split("|")[0] for acc in original])
    print("read,read_length,err,ins,del,subs,matches,err_rate,type,transcript_cov,gene_cov,gene_fam_cov")
    for acc in corrected:
        transcript_id = acc.split("|")[2].split("_")[0]
        gene_id = acc.split("|")[1] 
        gene_fam_id = acc.split("|")[0] 
        true_seq = true[transcript_id]
        read_length = len(true_seq)
        # res = edlib.align(seq, spoa_ref, task="path", mode="NW")
        # cigar_string = res["cigar"]
        # read_alignment, ref_alignment = help_functions.cigar_to_seq(cigar_string, seq, spoa_ref)

        corr_seq = corrected[acc]
        tot, ins, del_, subs, matches, err_rate = get_error_profile(corr_seq, true_seq)
        print("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}".format(acc, read_length, tot, ins, del_, subs, matches, err_rate, "corrected", isoform_coverage[transcript_id], gene_coverage[gene_id], gene_fam_coverage[gene_fam_id] ))
        # print("corrected:", ed_corr, "cov:", isoform_coverage[transcript_id])
        err_rate_corr = err_rate
        orig_seq = original[acc]
        
        # if err_rate_corr > 15:
        #     print("Over corrected!", file=sys.stderr)
        #     print(acc, file=sys.stderr)
        #     print(gene_coverage[gene_id], gene_fam_coverage[gene_fam_id], file=sys.stderr )
        #     print(corr_seq, tot, ins, del_, subs, err_rate, file=sys.stderr)
        #     print(true_seq, file=sys.stderr)
        #     tot, ins, del_, subs, err_rate =  get_error_profile(orig_seq, true_seq)
        #     print(orig_seq,tot, ins, del_, subs, err_rate, file=sys.stderr)
        #     print("", file=sys.stderr)
            

        tot, ins, del_, subs, matches, err_rate = get_error_profile(orig_seq, true_seq)
        print("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}".format(acc, read_length, tot, ins, del_, subs, matches, err_rate, "original", isoform_coverage[transcript_id], gene_coverage[gene_id], gene_fam_coverage[gene_fam_id] ))

        # res = edlib.align(orig_seq, true_seq, task="path", mode="NW")
        # ed_orig = res["editDistance"]
        # print("Original:", ed_orig, "cov:", isoform_coverage[transcript_id])





if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    parser.add_argument('true', type=str, help='Fastq file. ')
    parser.add_argument('original', type=str, help='fastq file. ')
    parser.add_argument('corrected', type=str, help='Fastq file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

