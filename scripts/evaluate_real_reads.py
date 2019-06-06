
from __future__ import print_function
import os,sys
import argparse
import re
import errno

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



def cigar_to_seq_mm2(read, full_r_seq, full_q_seq):
    # print()
    r_index = 0
    q_index = 0
    r_seq = full_r_seq[read.reference_start: read.reference_end + 1]
    q_seq = full_q_seq[read.query_alignment_start: read.query_alignment_end + 1]
    r_line, m_line, q_line = [], [], []

    cigar_tuples = read.cigartuples
    # print(cigar_tuples)
    # print(full_r_seq)
    # print(full_q_seq)

    # if query is softclipped in beginning or ref seq starts earlier
    if read.query_alignment_start > 0 or read.reference_start > 0:
        r_line.append( "-"*max(0, read.query_alignment_start - read.reference_start ) + full_r_seq[ : read.reference_start ] )
        q_line.append( "-"*max(0, read.reference_start - read.query_alignment_start ) + full_q_seq[ : read.query_alignment_start ]  ) 
        m_line.append("X"*max(read.query_alignment_start, read.reference_start ) )
        if cigar_tuples[0][0] in set([2,3,4]) and read.query_alignment_start == cigar_tuples[0][1]:
            del cigar_tuples[0]

    # if query is softclipped in end or ref seq extends beyond alignment 
    ref_end_offset = len(full_r_seq) - read.reference_end
    q_end_offset = len(full_q_seq) - read.query_alignment_end

    if q_end_offset > 0 or ref_end_offset > 0:
        # print(q_end_offset, ref_end_offset, "lol")
        if ref_end_offset:
            r_end = full_r_seq[ -ref_end_offset : ] + "-"*max(0, q_end_offset - ref_end_offset ) 
        else:
            r_end =  "-" * q_end_offset 

        if q_end_offset:
            q_end = full_q_seq[ -q_end_offset : ] + "-"*max(0, ref_end_offset - q_end_offset ) 
        else:
            q_end = "-" * ref_end_offset 

        m_end = "X"*max(q_end_offset, ref_end_offset )

        if cigar_tuples[-1][0] in set([2,3,4]) and q_end_offset == cigar_tuples[-1][1]:
            # print("HAHAHAHAHAHA")
            del cigar_tuples[-1]
        # del cigar_tuples[-1]
        # print(r_end, m_end, q_end)
    else:
        r_end, m_end, q_end = "", "", ""
    # print(cigar_tuples)

    for (op_type, op_len) in cigar_tuples:
        # op_len = int(op_len)
        if op_type == 0:
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            # print(ref_piece)
            # print(query_peace)

            r_line.append(ref_piece)
            q_line.append(query_peace)
            match_seq = ''.join(['|' if r_base.upper() == q_base.upper() else '*' for (r_base, q_base) in zip(ref_piece, query_peace)]) # faster with "".join([list of str]) instead of +=
                
            m_line.append(match_seq)
            r_index += op_len
            q_index += op_len

        elif op_type == 1:
            # insertion into reference
            r_line.append('-' * op_len)
            m_line.append(' ' * op_len)
            q_line.append(q_seq[q_index: q_index + op_len])
            #  only query index change
            q_index += op_len
        elif op_type == 2 or op_type == 3 or op_type == 4:
            # deletion from reference
            r_line.append(r_seq[r_index: r_index + op_len])
            m_line.append(' ' * op_len)
            q_line.append('-' * op_len)
            #  only ref index change
            r_index += op_len
            # print("here", r_index)

    r_line.append(r_end)
    m_line.append(m_end)
    q_line.append(q_end)    

    return "".join([s for s in r_line]), "".join([s for s in m_line]), "".join([s for s in q_line])



def get_aln_stats_per_read(sam_file, reads, refs):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    alignments = {}
    for read in SAM_file.fetch(until_eof=True):
        if read.is_secondary:
            print("secondary", read.reference_name) 
        # read_index = reads_acc_to_index[read.qname]
        else:
            print("primary", read.reference_name) 
            print(read.cigartuples)
            read_seq = reads[read.query_name]
            ref_seq = refs[read.reference_name]
            ref_alignment, m_line, read_alignment = cigar_to_seq_mm2(read, ref_seq, read_seq)
            print(ref_alignment)
            print(m_line)
            print(read_alignment)
            print()

    #     if read.is_reverse:
    #         continue
    #     if read.is_unmapped:
    #         acc, read_seq, qual = reads[read_index]
    #         print(read.qname, "is umapped", len(read_seq))
    #         alignments[read_index] = (acc, "unaligned", "unaligned", [])
    #         continue


    #     acc, read_seq, qual = reads[read_index]
    #     # cigar = read.cigarstring
    #     # print(read.qname, len(read_seq))
    #     match_id = math.floor((1.0 - 0.16) * args.k)
    #     block_coverage = get_block_coverage(read_alignment, ref_alignment, k, match_id)

    #     # print(ref_alignment)
    #     # print(read_alignment)
    #     # print(block_coverage)
    #     alignments[read_index] = (acc, read_alignment, ref_alignment, block_coverage)
    # return alignments
    return alignments

def main(args):
    reads = { acc : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    refs = { acc : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}
    # print(refs)
    orig = get_aln_stats_per_read(args.orig_sam, reads, refs)
    corr = get_aln_stats_per_read(args.corr_sam, reads, refs)

    orig_stats = get_summary_stats(orig)
    corr_stats = get_summary_stats(corr)

    print( "Num_aligned_reads", "Aligned bases, tot_errors", "avg_error_rate", "median_read_error_rate", "upper_25_quant", "lower_25_quant")
    print(",".join([s for s in orig_stats]))
    print(",".join([s for s in corr_stats]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('orig_sam', type=str, help='Path to the original read file')
    parser.add_argument('corr_sam', type=str, help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)