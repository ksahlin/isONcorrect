#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam

from collections import defaultdict



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
        name, seqs, last = last[1:].split()[0], [], None
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

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def main(args):
    refs = [(acc, seq) for acc, (seq, _) in readfq(open(args.fasta, 'r'))]
    sorted_refs = sorted(refs, key = lambda x: len(x[1]))
    outfile = open(args.outfasta, "w")
    for i, (acc, seq) in enumerate(sorted_refs[:-1]):
        is_redundant = False
        seq_rc = reverse_complement(seq)
        for acc2, seq2 in sorted_refs[i+1:]:
            if (seq in seq2 ) or (seq_rc in seq2):
                print(acc, "is redundant")
                is_redundant = True

        if is_redundant:
            continue
        else:
            outfile.write(">{0}\n{1}\n".format(acc,seq))
    outfile.write(">{0}\n{1}\n".format(sorted_refs[-1][0],sorted_refs[-1][1]))
    outfile.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('fasta', type=str, help='fastq file. ')
    parser.add_argument('outfasta', type=str, help='fastq file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

