
from __future__ import print_function
import os,sys
import argparse
import re
import errno
from collections import defaultdict

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

def cigar_to_seq_mm2_local(read, full_r_seq, full_q_seq):

    """
        Changed to parsing = and X instead of simple M.
    """
    # print()
    r_index = 0
    q_index = 0
    r_seq = full_r_seq[read.reference_start: read.reference_end + 1]
    q_seq = full_q_seq[read.query_alignment_start: read.query_alignment_end + 1]
    r_line, m_line, q_line = [], [], []

    cigar_tuples = read.cigartuples
    r_end, m_end, q_end = "", "", ""
    ins, del_, subs, matches = 0, 0, 0, 0
    for (op_type, op_len) in cigar_tuples:
        # print((op_type, op_len))
        # op_len = int(op_len)
        if op_type == 7:
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            # print(ref_piece)
            # print(query_peace)

            r_line.append(ref_piece)
            q_line.append(query_peace)
            match_seq = ''.join(['|' if r_base.upper() == q_base.upper() else '*' for (r_base, q_base) in zip(ref_piece, query_peace)]) # faster with "".join([list of str]) instead of +=
            matches += op_len 
            m_line.append(match_seq)
            r_index += op_len
            q_index += op_len

        elif op_type == 8:
            subs += op_len 


        elif op_type == 1:
            # insertion into reference
            r_line.append('-' * op_len)
            m_line.append(' ' * op_len)
            q_line.append(q_seq[q_index: q_index + op_len])
            #  only query index change
            q_index += op_len
            ins += op_len
        elif op_type == 2 or op_type == 3:
            # deletion from reference
            r_line.append(r_seq[r_index: r_index + op_len])
            m_line.append(' ' * op_len)
            q_line.append('-' * op_len)
            #  only ref index change
            r_index += op_len
            # print("here", r_index)
            if op_type == 2:
                del_ += op_len
        elif op_type == 4:
            # softclip
            pass


    r_line.append(r_end)
    m_line.append(m_end)
    q_line.append(q_end)    

    return "".join([s for s in r_line]), "".join([s for s in m_line]), "".join([s for s in q_line]), ins, del_, subs, matches



def get_abundance_aligned_reads(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    
    gene_cov_true = defaultdict(int)
    transcript_cov_true = defaultdict(int)
    gene_fam_cov_true = defaultdict(int)

    gene_cov_aligned = defaultdict(int)
    transcript_cov_aligned = defaultdict(int)
    gene_fam_cov_aligned = defaultdict(int)
    read_specific = {}
    primary_mapq_0 = set()
    optimal_cigar_str = {}
    edit_distance_to_aligned = {}
    for read in SAM_file.fetch(until_eof=True):
        
        read_acc = read.query_name
        transcript_id = read_acc.split("|")[2].split("_")[0]
        sim_read_nr = int(read_acc.split("|")[2].split("_")[1])
        
        if transcript_id not in transcript_cov_true or sim_read_nr >= transcript_cov_true[transcript_id]:
            transcript_cov_true[transcript_id] = sim_read_nr + 1

        if read.flag == 0 or read.flag == 16:
            # print(read.is_reverse)
            # print(read.cigartuples)

            optimal_cigar_str[read_acc] = read.cigarstring
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_ = sum([length for type_, length in read.cigartuples if type_ == 2])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            edit_distance_to_aligned[read_acc] = ins + del_ + subs

            if read.mapping_quality == 0:
                primary_mapq_0.add(read_acc)

            

            gene_id = read_acc.split("|")[1] 
            gene_fam_id = read_acc.split("|")[0] 
            
            gene_cov_true[gene_id] += 1
            gene_fam_cov_true[gene_fam_id] += 1

            ref_acc = read.reference_name
            
            transcript_id = ref_acc.split("|")[2]
            gene_id = ref_acc.split("|")[1] 
            gene_fam_id = ref_acc.split("|")[0] 

            transcript_cov_aligned[transcript_id] += 1
            gene_cov_aligned[gene_id] += 1
            gene_fam_cov_aligned[gene_fam_id] += 1

            assert read_acc not in read_specific
            read_specific[read_acc] = set( )
            read_specific[read_acc].add(transcript_id)
        elif read_acc not in read_specific:
            read_specific[read_acc] = set( )
            read_specific[read_acc].add("unaligned")
            edit_distance_to_aligned[read_acc] = '-'

        elif (read.flag != 0 and read.flag != 16) and read_acc in primary_mapq_0:
            if read.cigarstring == optimal_cigar_str[read_acc]:
                ref_acc = read.reference_name            
                transcript_id = ref_acc.split("|")[2]
                read_specific[read_acc].add( transcript_id )
                # print("HERE@@@", ref_acc, read_acc)
            # print("secondary", read.flag, read.reference_name)

    # read_specific_dict = {}
    # for acc in  read_specific_dict:
    #     all_ref = set([ t[0] for t in read_specific_dict[acc]])
    #     read_specific_dict[acc] = all_ref

    return transcript_cov_true, gene_cov_true, gene_fam_cov_true, transcript_cov_aligned, gene_cov_aligned, gene_fam_cov_aligned, read_specific, edit_distance_to_aligned

# def get_summary_stats(reads, quantile):
#     tot_ins, tot_del, tot_subs, tot_match = 0, 0, 0, 0
#     sorted_reads = sorted(reads.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )
#     for acc, (ins, del_, subs, matches) in sorted_reads[ : int(len(sorted_reads)*quantile)]:
#         # (ins, del_, subs, matches) = reads[acc]

#         tot_ins += ins
#         tot_del += del_
#         tot_subs += subs
#         tot_match += matches

#     sum_aln_bases = tot_ins + tot_del + tot_subs + tot_match

#     return tot_ins, tot_del, tot_subs, tot_match, sum_aln_bases

import edlib

def main(args):
    reads = { acc : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    refs = { acc.split("|")[2] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.refs, 'r')))}
    transcript_cov_true, gene_cov_true, gene_fam_cov_true, transcript_cov_aligned, gene_cov_aligned, gene_fam_cov_aligned, read_specific, edit_distance_to_aligned = get_abundance_aligned_reads(args.samfile)
    # "read_acc","aligned_to","transcript_abundance","is_tp","read_type"
    for read_acc, set_aligned_to in read_specific.items():
        true_transcript = read_acc.split("|")[2].split("_")[0]
        if true_transcript in transcript_cov_true:
            true_transcript_abundance = transcript_cov_true[true_transcript]

        is_correct = 1 if true_transcript in set_aligned_to else 0
        aligned_to = true_transcript if is_correct == 1 else set_aligned_to.pop()
        if is_correct == 1:
            ed_btw_transcripts = "-" 
            ed_read_to_true = edit_distance_to_aligned[read_acc]
            ed_read_to_aligned = edit_distance_to_aligned[read_acc]
        elif aligned_to == 'unaligned':
            res = edlib.align(reads[read_acc], refs[true_transcript], mode="NW")
            # print(args.type, true_transcript_abundance, "ed:", res, read_acc)
            continue
        else:
            res1 = edlib.align(refs[aligned_to], refs[true_transcript], mode="HW")
            res2 = edlib.align(refs[true_transcript], refs[aligned_to],  mode="HW")
            ed_btw_transcripts = min(res1["editDistance"],res2["editDistance"])

            res = edlib.align(reads[read_acc], refs[true_transcript], mode="NW")
            ed_read_to_true = res["editDistance"]
            ed_read_to_aligned = edit_distance_to_aligned[read_acc]

        print("{0},{1},{2},{3},{4},{5},{6},{7}".format(read_acc, aligned_to, true_transcript_abundance, is_correct, args.type, ed_btw_transcripts, ed_read_to_true, ed_read_to_aligned))

    # # print("id,cov_aln,cov_true,seq,type")
    # for seq_id in set(transcript_cov_true) | set(transcript_cov_aligned) :
    #     cov_aln = transcript_cov_aligned[seq_id] if seq_id in transcript_cov_aligned else 0
    #     cov_true = transcript_cov_true[seq_id] if seq_id in transcript_cov_true else 0
    #     if "," in seq_id:
    #         print("BUG", seq_id)
    #         sys.exit()
    #     if type(cov_aln) != int or type(cov_true) != int:
    #         print("BUG", type(cov_aln))
    #         sys.exit()
    #     print("{0},{1},{2},{3},{4}".format(seq_id, cov_aln, cov_true, "transcript", args.type))

    # for seq_id in set(gene_cov_true) | set(gene_cov_aligned) :
    #     cov_aln = gene_cov_aligned[seq_id] if seq_id in gene_cov_aligned else 0
    #     cov_true = gene_cov_true[seq_id] if seq_id in gene_cov_true else 0
    #     if "," in seq_id:
    #         print("BUG", seq_id)
    #         sys.exit()
    #     if type(cov_aln) != int:
    #         print("BUG", type(cov_aln))
    #         sys.exit()

    #     print("{0},{1},{2},{3},{4}".format(seq_id, cov_aln, cov_true, "gene", args.type))

    # for seq_id in set(gene_fam_cov_true) | set(gene_fam_cov_aligned) :
    #     cov_aln = gene_fam_cov_aligned[seq_id] if seq_id in gene_fam_cov_aligned else 0
    #     cov_true = gene_fam_cov_true[seq_id] if seq_id in gene_fam_cov_true else 0
    #     if "," in seq_id:
    #         print("BUG", seq_id)
    #         sys.exit()

    #     print("{0},{1},{2},{3},{4}".format(seq_id, cov_aln, cov_true, "gene_fam", args.type))      

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('samfile', type=str, help='Path to the original read file')
    parser.add_argument('refs', type=str, help='Fastq file. ')
    parser.add_argument('reads', type=str, help='fastq file. ')
    parser.add_argument('type', type=str, help='Path to the original read file')
    # parser.add_argument('reads', type=str, help='Path to the original read file')

    args = parser.parse_args()

    # outfolder = args.outfolder
    # if not os.path.exists(outfolder):
    #     os.makedirs(outfolder)
    main(args)


