from __future__ import print_function
import os,sys
import argparse
import re
# import ssw
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
from collections import Counter
import parasail
from collections import defaultdict



import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns


import edlib

'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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



def edlib_ed(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    return ed, result["cigar"]


import math
def filter_sam(sam_file, seq1, seq2, outsam, dist ):
    outfile = open(outsam, 'w')
    it, written = 0,0
    for line in open(sam_file, "r"):
        it += 1
        vals = line.split()
        read_seq = vals[9].strip()
        edit_distance1, cigar = edlib_ed(read_seq, seq1, mode="HW", task="path", k = 2*len(seq1)) 
        edit_distance2, cigar = edlib_ed(seq1, read_seq, mode="HW", task="path", k = 2*len(seq1)) 
        seq1_min = min(edit_distance1,edit_distance2)
        
        edit_distance1, cigar2 = edlib_ed(read_seq, seq2, mode="HW", task="path", k = 2*len(seq2)) 
        edit_distance2, cigar2 = edlib_ed(seq2, read_seq, mode="HW", task="path", k = 2*len(seq2)) 
        seq2_min = min(edit_distance1,edit_distance2)

        if math.fabs(seq1_min - seq2_min) < dist:
            # print("ambiguous",seq1_min, seq2_min, vals[0])
            # print(cigar)
            # print(cigar2)
            pass
        else:
            outfile.write(line)
            written += 1

    outfile.close()
    print("total iterated:", it)
    print("total written:", written)


def main(args):

    isoforms = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.isoforms, 'r'))}
    seq1, seq2 = list(isoforms.values())
    filtered_sam = filter_sam(args.sam, seq1, seq2, args.outsam, args.dist )





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('sam', type=str, help='Path to the reads SAM file')
    parser.add_argument('isoforms', type=str, help='Path to the transcript fasta file')
    parser.add_argument('dist', type=int, help='Path to the transcript fasta file')
    parser.add_argument('outsam', type=str, help='Output path of results')

    args = parser.parse_args()

    main(args)
