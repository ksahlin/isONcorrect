from __future__ import print_function
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout

import itertools
from collections import defaultdict, deque
import statistics

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None
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



def kmer_counter(genome, k_size, target_kmers = {}):
	count = defaultdict(int)

	for chr_id in genome:
		# print(chr_id)
		if "_" not in chr_id:  # this assures only using kmers from primary seqs like chr1,...chr22, and chrX/Y/M
			if len(target_kmers) > 0:
				print(chr_id)
			seq = genome[chr_id]
			# seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
			# read_kmers = deque([seq[i:i+k_size] for i in range(len(seq) - k_size )])
			for i in range(0, len(seq) - k_size +1):
				kmer = seq[i:i+k_size].upper()
				if i % 10000000 == 0 and i > 1:
					print(i, chr_id)
				if len(target_kmers) == 0 or kmer in target_kmers:
					count[kmer] += 1

	return count #, position_count

def annotate_transcript_with_mean_mappability(transcripts, kmer_counts, k_size):
	ann_transcripts = {}
	for acc, seq in transcripts.items():
		sum_mappability = 0
		pos_mappabilities = []
		for i in range(0, len(seq) - k_size + 1):
			kmer = seq[i:i+k_size]
			count = kmer_counts[kmer] if kmer in kmer_counts else 1
			sum_mappability += count
			pos_mappabilities.append(count)

		avg_mappability = sum_mappability / (len(seq) - k_size + 1)
		median_mappability = statistics.median(pos_mappabilities)
		
		ann_transcripts[acc] = avg_mappability
		# ann_transcripts[acc] = median_mappability
		
	return ann_transcripts


def main(args):
    # reads = { acc.split("|")[2] : seq for acc, (seq,_) in readfq(open(args.true, 'r'))}
    transcripts = { acc : seq for acc, (seq,qual) in readfq(open(args.transcripts, 'r'))}
    transcript_kmers =  kmer_counter(transcripts, args.kmer_size)
    print("Step 1")

    genome = { acc : seq for acc, (seq,qual) in readfq(open(args.genome, 'r'))}
    print("Step 2")
    kmer_counts =  kmer_counter(genome, args.kmer_size, transcript_kmers)
    print("Step 3")

    ann_transcripts = annotate_transcript_with_mean_mappability(transcripts, kmer_counts, args.kmer_size)

    f = open(args.outfile, "w")
    for acc, mean_mapp in sorted(ann_transcripts.items(), key = lambda x: x[1], reverse = True): #sort on increasing mappability before printing?
    	# mean_mapp = ann_transcripts[acc]
    	f.write("{0},{1}\n".format(acc, mean_mapp))
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('transcripts', type=str, help='Fasta or fastq')
    parser.add_argument('genome', type=str, help='Fasta or fastq')
    parser.add_argument('outfile', type=str, help='Fasta or fastq')
    parser.add_argument('kmer_size', type=int, help='integer')


    args = parser.parse_args()

    main(args)
