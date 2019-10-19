import os
import sys

import itertools
import argparse
import errno
import math
import random
# from collections import deque

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
        name, seqs, last = last[1:].replace("_", " ").split()[0], [], None
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



def main(args):

    outfile = open(args.outfile, "w")
    for acc, (seq,qual) in readfq(open(args.fastq, 'r')):
        outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))

    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot p-minimizers shared.")
    parser.add_argument('fastq', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('outfile', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)

