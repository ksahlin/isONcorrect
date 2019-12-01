#!/usr/bin/env python

import argparse
import sys

from collections import Counter
import re

import edlib
import pysam


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
            errors = insertions + deletions + subs 
            identity = matches/float(tot_align)
            err_rate = round(100 * float(subs + insertions + deletions) / float(tot_align), 2)
            error_rates_per_read[read.query_name] = (read.reference_name, errors,insertions, deletions, subs, matches, err_rate)
    # print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return error_rates_per_read     




def main(args):
    refs = { acc : seq for acc, (seq,_) in readfq(open(args.refs, 'r'))}
    # reads = { acc : seq for acc, (seq,qual) in readfq(open(args.reads, 'r'))}
    tot_errors, tot_insertions, tot_deletions, tot_subs, tot_matches = 0,0,0,0,0
    outfile = open(args.outfile, 'w')
    outfile.write("read,ref,err,err_rate,subs,ins,del\n")
    # print(outfile.name)
    error_rates_per_read = decide_primary_locations(args.samfile, args)
    # print(error_rates_per_read)
    for acc in error_rates_per_read:
        reference_name, errors, insertions, deletions, subs, matches, err_rate = error_rates_per_read[acc]
        outfile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(acc, reference_name, errors, err_rate, subs, insertions, deletions))
        
        tot_errors +=  errors
        tot_insertions +=  insertions
        tot_deletions +=  deletions
        tot_subs +=  subs
        tot_matches += matches

    outfile.close()
    # print(tot_errors, tot_insertions, tot_deletions, tot_subs,tot_errors,tot_matches)
    print(tot_errors, tot_insertions, tot_deletions, tot_subs, round(100*tot_errors/float(tot_matches + tot_errors), 2) )


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    parser.add_argument('refs', type=str, help='Fastq file. ')
    parser.add_argument('samfile', type=str, help='fastq file. ')
    parser.add_argument('outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

