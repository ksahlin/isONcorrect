
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

from collections import defaultdict

import parasail
import pysam
import gffutils

'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''
def readfq_split(fp): # this is a generator function
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
    # print(cigar_tuples)
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
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            r_line.append(ref_piece)
            q_line.append(query_peace)

            subs += op_len 
            m_line.append('*' * op_len)
            r_index += op_len
            q_index += op_len


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
        elif op_type == 0:
            sys.exit("You have to use =/X cigar operators in sam file (M was detected). In minimap2 you can to provide --eqx flag to fix this.")


    r_line.append(r_end)
    m_line.append(m_end)
    q_line.append(q_end)    

    return "".join([s for s in r_line]), "".join([s for s in m_line]), "".join([s for s in q_line]), ins, del_, subs, matches


def decide_primary_locations(sam_file, args): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    reads_primary = {}
    reads_ref_locations = {}
    reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            cigar = read.cigarstring
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_ = sum([length for type_, length in read.cigartuples if type_ == 2 ])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            tot_align = ins + del_ + subs + matches
            identity = matches/float(tot_align)
            reads_primary[read.query_name] = read
            reads_ref_locations[read.query_name] = set()
            reads_ref_locations[read.query_name].add(read.reference_name) 


            # has_large_del = [length for type_, length in read.cigartuples if type_ == 2 and length >=  args.min_intron ]
            # if has_large_del:
            #     print(has_large_del)
            # if read.query_name in reads_primary:
            #     reads_multiple_primary[read.query_name] = cigar
            #     if identity >= reads_tmp[read.query_name][0] and  matches >= reads_tmp[read.query_name][1]:
            #         reads_primary[read.query_name] = read
            #         reads_tmp[read.query_name] = (identity, matches)
            #     elif identity <= reads_tmp[read.query_name][0] and  matches <= reads_tmp[read.query_name][1]:
            #         continue
            #     else:
            #         if identity * matches > reads_tmp[read.query_name][0] * reads_tmp[read.query_name][1]:
            #             reads_primary[read.query_name] = read
            #             reads_tmp[read.query_name] = (identity, matches)
            #         else: 
            #             continue

                        
            #         # print( "Ambiguous, prefferred:", round(reads_tmp[read.query_name][0],2), reads_tmp[read.query_name][1], "over", round(identity,2), matches)

            #     #     if tot_align >= reads_tmp[read.query_name][1]:
            #     #         reads_primary[read.query_name] = read
            #     #     else:
            #     #         continue
            #     # elif tot_align >= reads_tmp[read.query_name][1]:
            #     #     if identity >= reads_tmp[read.query_name][0]:
            #     #         reads_primary[read.query_name] = read
            #     #     else:
            #     #         continue

            # else:
            #     reads_primary[read.query_name] = read
            #     reads_tmp[read.query_name] = (identity, matches)
        elif read.query_name in reads_primary and read.cigarstring == reads_primary[read.query_name].cigarstring:
            # assert read.query_name in reads_ref_locations
            reads_ref_locations[read.query_name].add(read.reference_name) 


    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_ref_locations))
    return reads_primary, reads_ref_locations    


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



def get_deletion_sites(read):
    deletion_sites = []
    ref_pos = read.reference_start
    
    for i, (t,l) in enumerate(read.cigartuples):
        if t == 2:
            deletion_start = ref_pos
            ref_pos += l
            deletion_stop = ref_pos
            deletion_sites.append( (deletion_start, deletion_stop, l) )

        elif t == 7 or t== 0 or t == 8:
            ref_pos += l
        elif t == 3:
            # splice_start = ref_pos
            ref_pos += l
            # splice_stop = ref_pos
            # splice_sites.append( (splice_start, splice_stop) )

        elif t == 1 or t == 4 or t == 5: # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return deletion_sites

def get_error_rate_stats_per_read(reads_primary_locations, reads, args, reference = {}):
    # SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    alignments = {}
    read_index = 0
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # print(read.is_reverse)
            read_index += 1
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_ = sum([length for type_, length in read.cigartuples if type_ == 2 ])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            called_intron_lengths = [length for type_, length in read.cigartuples if type_ == 3]

            assert read.query_name not in alignments
            alignments[read.query_name] = (ins, del_, subs, matches, read.reference_name, read.reference_start, read.reference_end + 1, read.flag, read_index)

        else:
            pass
    return alignments

def get_summary_stats(reads, quantile):
    tot_ins, tot_del, tot_subs, tot_match = 0, 0, 0, 0
    sorted_reads = sorted(reads.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )
    for acc, (ins, del_, subs, matches, chr_id, reference_start, reference_end, sam_flag, read_index) in sorted_reads[ : int(len(sorted_reads)*quantile)]:
        # (ins, del_, subs, matches) = reads[acc]

        tot_ins += ins
        tot_del += del_
        tot_subs += subs
        tot_match += matches

    sum_aln_bases = tot_ins + tot_del + tot_subs + tot_match

    return tot_ins, tot_del, tot_subs, tot_match, sum_aln_bases

def print_detailed_values_to_file(alignments_dict, reads, all_ref_locations,  outfile, reads_unaligned_in_other_method, reads_missing_from_clustering_correction_output, read_type):
    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )
    for (acc, (ins, del_, subs, matches, chr_id, reference_start, reference_end, sam_flag, read_index)) in alignments_sorted:
        error_rate = round( 100* (ins + del_ + subs) /float( (ins + del_ + subs + matches) ), 4 )
        read_length = len(reads[acc])
        aligned_length = int( ins + del_ + subs + matches )
        is_unaligned_in_other_method = 1 if acc in reads_unaligned_in_other_method else 0
        all_refs = ":".join([ s for s in all_ref_locations[acc]])
        # is_missing_from_clustering_or_correction = 1 if acc in reads_missing_from_clustering_correction_output else 0

        info_tuple = (acc, read_type, ins, del_, subs, matches, error_rate, read_length, aligned_length, all_refs) # 'tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'fsm', 'nic', 'ism', 'nnc', 'no_splices'  )
        # print(*info_tuple)
        # outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format(*info_tuple))
        outfile.write( ",".join( [str(item) for item in info_tuple] ) + "\n")


def print_quantile_values(alignments_dict):

    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )
    print(alignments_sorted[-5:])

    sorted_error_rates = [ sum(tup[0:3])/float(sum(tup[0:4])) for acc, tup in alignments_sorted]

    insertions = [ tup[0]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]
    deletions = [ tup[1]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]
    substitutions = [ tup[2]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]

    n = len(sorted_error_rates)
    quantiles = [0, 1.0/20, 1.0/10, 1.0/4, 1.0/2, 3.0/4, 9.0/10, 19.0/20, 1 ]
    quantile_errors = []
    quantile_insertions = []
    quantile_deletions = []
    quantile_substitutions = []

    for q in quantiles:
        if q == 1:
            quantile_errors.append(sorted_error_rates[int(q*n)-1])
            quantile_insertions.append(insertions[int(q*n)-1])
            quantile_deletions.append(deletions[int(q*n)-1])
            quantile_substitutions.append(substitutions[int(q*n)-1])
        else:
            quantile_errors.append(sorted_error_rates[int(q*n)])
            quantile_insertions.append(insertions[int(q*n)])
            quantile_deletions.append(deletions[int(q*n)])
            quantile_substitutions.append(substitutions[int(q*n)])

    return quantile_errors, quantile_insertions, quantile_deletions, quantile_substitutions


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


def parasail_alignment(read, reference, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -2, opening_penalty = 2, gap_ext = 1, ends_discrepancy_threshold = 0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(read, reference, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(read, reference, opening_penalty, gap_ext, user_matrix)
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    
    read_alignment, ref_alignment = cigar_to_seq(cigar_string, read, reference)
    return read_alignment, ref_alignment

def get_cluster_sizes(cluster_file, reads):
    cluster_sizes = defaultdict(int)
    tmp_reads = {}
    for line in open(args.cluster_file, "r"):
        cl_id, acc = line.split() 
        if acc in reads:
            cluster_sizes[cl_id] += 1  
            tmp_reads[acc] = cl_id
        else:
            cluster_sizes[cl_id] += 1  
            tmp_reads[acc.split("_")[0]] = cl_id

    reads_to_cluster_size = {}
    for acc, cl_id in tmp_reads.items():
        reads_to_cluster_size[acc] = cluster_sizes[cl_id]

    return reads_to_cluster_size




from collections import namedtuple

def pickle_dump(data, filename):
    with open(os.path.join(args.outfolder,filename), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data



def main(args):
    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}
    print(refs.keys())


    reads = { acc.split()[0] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    corr_reads = { acc.split()[0] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.corr_reads, 'r')))}
    # print(list(reads.keys()))
    # print(list(corr_reads.keys()))
    orig_primary_locations, orig_ref_locations = decide_primary_locations(args.orig_sam, args)
    corr_primary_locations, corr_ref_locations = decide_primary_locations(args.corr_sam, args)

    corr = get_error_rate_stats_per_read(corr_primary_locations, corr_reads, args)
    orig = get_error_rate_stats_per_read(orig_primary_locations, reads, args)

    print( "Reads successfully aligned (as found in sam file):", len(orig),len(corr))


    ##############
    quantile_tot_orig, quantile_insertions_orig, quantile_deletions_orig, quantile_substitutions_orig = print_quantile_values(orig)
    quantile_tot_corr, quantile_insertions_corr, quantile_deletions_corr, quantile_substitutions_corr = print_quantile_values(corr)

    orig_stats = get_summary_stats(orig, 1.0)
    corr_stats = get_summary_stats(corr, 1.0)
    
    print("Distribution of error rates (Percent)")
    print("Reads, Best, top 5%, top 10%, top 25%, Median, top 75%, top 90%, top 95%, Worst")
    print("Original,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_orig ] ))
    print("Corrected,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_corr ] ))
    ################


    
    reads_missing_from_clustering_correction_output = set(reads.keys()) - set(corr_reads.keys())
    for r_acc in reads_missing_from_clustering_correction_output:
        if r_acc in orig and r_acc not in corr:
            corr[r_acc] = orig[r_acc]
            corr_reads[r_acc] = reads[r_acc]
            corr_ref_locations[r_acc] = orig_ref_locations[r_acc]

    bug_if_not_empty = set(corr_reads.keys()) - set(reads.keys())
    reads_unaligned_in_original = set(reads.keys()) - set(orig_primary_locations.keys())
    reads_unaligned_in_correction = set(corr_reads.keys()) - set(corr_primary_locations.keys()) 

    detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read_to_transcriptome.csv"), "w")
    detailed_results_outfile.write("acc,read_type,ins,del,subs,matches,error_rate,read_length,aligned_length,chr_id\n")
    print_detailed_values_to_file(corr, corr_reads, corr_ref_locations, detailed_results_outfile, reads_unaligned_in_original, reads_missing_from_clustering_correction_output, "corrected")    
    print_detailed_values_to_file(orig, reads, orig_ref_locations, detailed_results_outfile, reads_unaligned_in_correction, reads_missing_from_clustering_correction_output, "original")
    detailed_results_outfile.close()
    for acc in orig_ref_locations:
        orig_loc = orig_ref_locations[acc]
        o = orig_loc.pop()
        corr_loc = corr_ref_locations[acc]
        c = corr_loc.pop()
        if o != c:
            print(o, c)

    print()
    print("Reads successfully aligned all reads (original/corrected):", len(orig),len(corr))
    print("Total reads (original/corrected):", len(reads),len(corr_reads))
    print("READS MISSING FROM CLUSTERING/CORRECTION INPUT:", len(reads_missing_from_clustering_correction_output))
    print("READS UNALIGNED (ORIGINAL/CORRECTED):", len(reads_unaligned_in_original), len(reads_unaligned_in_correction) )
    print("BUG IF NOT EMPTY:",len(bug_if_not_empty))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('orig_sam', type=str, help='Path to the original read file')
    parser.add_argument('corr_sam', type=str, help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('corr_reads', type=str, help='Path to the corrected read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

