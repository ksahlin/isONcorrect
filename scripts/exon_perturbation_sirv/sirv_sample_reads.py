#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam
import errno
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


def get_subsamples(transcript_cov, depth, nr_isoforms):
    subsamples = {}
    already_ampled = set()
    for gene_id in transcript_cov:
        trancripts = random.sample(list(transcript_cov[gene_id]), nr_isoforms)
        for tr_acc in trancripts:
            # print(tr_acc, len(transcript_cov[tr_acc]))
            subset = random.sample(transcript_cov[gene_id][tr_acc], depth) # this is without replacement
            subsamples[tr_acc] = subset
            # print(len(subsamples[tr_acc]))
    # print(len(subsamples))
    return subsamples

    # subsamples = {}
    # already_ampled = set()
    # for tr_acc in transcript_cov:
    #     print(tr_acc, len(transcript_cov[tr_acc]))
    #     subset = random.sample(transcript_cov[tr_acc], depth) # this is without replacement
    #     subsamples[tr_acc] = subset
    # return subsamples



def get_aligned_reads(sam_file, minor_isoform, major_isoform ):
    isoform_reads = {} # defaultdict(lambda: defaultdict(set))
    for acc in [minor_isoform, major_isoform]:
        # print(acc)
        isoform_reads[acc] = set()
    

    for line in open(sam_file, "r"):
        vals = line.split()
        flag = 'minor' if vals[2] == minor_isoform else 'major'
        isoform_reads[vals[2]].add( (vals[0] +"|"+ vals[2]+"|"+ flag, vals[9], vals[10]) )

    # SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    # for read in SAM_file.fetch(until_eof=True):
    #     assert read.flag == 0
    #     isoform_reads[read.reference_name].add((read.query_name +"|"+ read.reference_name, read.query_sequence))

    return isoform_reads

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def main(args):
    nr_reads_minor = int(args.p*args.d) 
    nr_reads_major = int((1-args.p)*args.d) 
    isoforms = { acc : seq for acc, (seq, _ ) in readfq(open(args.isoforms, 'r'))}

    l = list(isoforms.keys())
    random.shuffle(l)
    minor_isoform, major_isoform = l[0], l[1]
    print("MINOR:", minor_isoform, len(isoforms[minor_isoform]), len(isoforms[major_isoform]) )
    isoform_reads = get_aligned_reads(args.alignments, minor_isoform, major_isoform)


    reads_minor = random.sample(isoform_reads[minor_isoform], nr_reads_minor)
    reads_major = random.sample(isoform_reads[major_isoform], nr_reads_major)

    all_reads = reads_minor + reads_major
    random.shuffle(all_reads)

    outfile = open(args.outfile, "w")
    for (acc, seq, qual_seq) in all_reads:
        outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual_seq))
    outfile.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('isoforms', type=str, help='fasta file. ')
    parser.add_argument('alignments', type=str, help='Sam file. ')
    parser.add_argument('p', type=float, help='fraction minor.')
    parser.add_argument('d', type=int, help='Depth per transcript.')
    parser.add_argument('outfile', type=str, help='Fastq file. ')

    args = parser.parse_args()
    dirname=os.path.dirname(args.outfile)
    mkdir_p(dirname)
    # print(dirname)
    main(args)

