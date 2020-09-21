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



def get_abundance_aligned_reads(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    valid_genes = set(["SIRV3", "SIRV5", "SIRV6"])
    transcript_cov = defaultdict(lambda: defaultdict(set))
    # amgiguous_primary = defaultdict(set)

    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16):
            transcript_id = read.reference_name
            if transcript_id[:5] in valid_genes:
                gene_id = transcript_id[4]
                transcript_cov[gene_id][transcript_id].add(read.query_name)
    # except:
    #     pass

    # print(transcript_cov)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality == 0:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)

    return transcript_cov

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
    transcript_cov = get_abundance_aligned_reads(args.alignments)
    print(len(transcript_cov), [len(transcript_cov[g]) for g in transcript_cov])
    subsamples = get_subsamples(transcript_cov, args.depth_per_transcript, args.nr_isoforms)
    fastq = { acc : (seq,qual) for acc, (seq,qual) in readfq(open(args.fastq, 'r'))}

    outfile = open(args.outfile, "w")
    for tr_id, set_of_reads in subsamples.items():
        for read_acc in set_of_reads:
            seq, qual  = fastq[read_acc]
            outfile.write("@{0}\n{1}\n+\n{2}\n".format(read_acc, seq, qual))

    outfile.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('fastq', type=str, help='fastq file. ')
    # parser.add_argument('ref_fastq', type=str, help='fastq file. ')
    parser.add_argument('alignments', type=str, help='fastq file. ')
    parser.add_argument('outfile', type=str, help='Fastq file. ')
    parser.add_argument('depth_per_transcript', type=int, help='Depth per transcript.')
    parser.add_argument('nr_isoforms', type=int, help='Number of nr_isoforms.')

    args = parser.parse_args()
    dirname=os.path.dirname(args.outfile)
    mkdir_p(dirname)
    print(dirname)
    main(args)

