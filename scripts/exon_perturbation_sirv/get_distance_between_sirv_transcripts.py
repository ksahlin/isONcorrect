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

def main(args):
    transcripts = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.transcripts, 'r'))}
    for acc, seq in transcripts.items():
        for acc2, seq2 in transcripts.items():
            if acc == acc2:
                continue
            k = 2*max(len(seq), len(seq2))
            # print(k)
            edit_distance, cigar = edlib_ed(seq, seq2, mode="NW", k = k) 
            if edit_distance < 60:
                print("ED between {0} and {1} is: {2}. Cigar is: {3}".format(acc, acc2, edit_distance, cigar ))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('transcripts', type=str, help='Path to the transcript fasta file')
    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
