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


def get_subsamples(transcript_cov):
    subsamples = {}
    already_ampled = set()
    for subsample_size in [3,5,10,20]: #,50,100,200,500]:
        large_enough_cov = [tr_id for tr_id in transcript_cov.keys() if len(transcript_cov[tr_id]) >= subsample_size]
        sampled_tr_id = random.choice(large_enough_cov)
        while sampled_tr_id in already_ampled:
            sampled_tr_id = random.choice(large_enough_cov)

        subset = random.sample(transcript_cov[sampled_tr_id], subsample_size)
        subsamples[subsample_size] = (sampled_tr_id, subset)
        already_ampled.add(sampled_tr_id)
    return subsamples


def get_abundance_aligned_reads(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    
    transcript_cov = defaultdict(set)
    amgiguous_primary = defaultdict(set)
    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16): # and read.mapping_quality > 0:
            transcript_id = read.reference_name
            transcript_cov[transcript_id].add(read.query_name)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality == 0:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)

    return transcript_cov, amgiguous_primary

def main(args):
    transcript_cov, amgiguous_primary = get_abundance_aligned_reads(args.alignments)
    for tr_id, read_aln in sorted(transcript_cov.items(), key = lambda x: len(x[1]), reverse = True):
        print(tr_id, len(read_aln) )
    sys.exit()

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('alignments', type=str, help='fastq file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

