
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

from collections import defaultdict

import parasail
import pysam
import gffutils

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


def get_aln_stats_per_read(sam_file, reads, refs, args):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    alignments = {}
    alignments_detailed = {}
    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            # print(read.is_reverse)
            # print(read.cigartuples)
            read_seq = reads[read.query_name]
            ref_seq = refs[read.reference_name]
            ref_alignment, m_line, read_alignment, ins, del_, subs, matches = cigar_to_seq_mm2_local(read, ref_seq, read_seq)
            # if read.query_name == 'b1d0ee62-7557-4645-8d1e-c1ccfb60c997_runid=8c239806e6f576cd17d6b7d532976b1fe830f9c6_sampleid=pcs109_sirv_mix2_LC_read=29302_ch=69_start_time=2019-04-12T22:47:30Z_strand=-':
            #     print(read.query_alignment_start)
            #     print(read.query_name, "primary", read.flag, read.reference_name) 
            #     print(ref_alignment)
            #     print(m_line)
            #     print(read_alignment)
            #     print(ins, del_, subs, matches)

            if read.query_name in alignments:
                print("BUG")
                sys.exit()
            
            alignments[read.query_name] = (ins, del_, subs, matches)
            if args.align:
                alignments_detailed[read.query_name] = (ref_alignment, m_line, read_alignment)
            # print()
            # return
        else:
            pass
            # print("secondary", read.flag, read.reference_name) 
    SAM_file.close()
    return alignments, alignments_detailed

def get_summary_stats(reads, quantile):
    tot_ins, tot_del, tot_subs, tot_match = 0, 0, 0, 0
    sorted_reads = sorted(reads.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )
    for acc, (ins, del_, subs, matches) in sorted_reads[ : int(len(sorted_reads)*quantile)]:
        # (ins, del_, subs, matches) = reads[acc]

        tot_ins += ins
        tot_del += del_
        tot_subs += subs
        tot_match += matches

    sum_aln_bases = tot_ins + tot_del + tot_subs + tot_match

    return tot_ins, tot_del, tot_subs, tot_match, sum_aln_bases

def print_detailed_values_to_file(alignments_dict, annotations_dict, reads_to_cluster_size, reads, outfile, read_type):
    # read_calss is FSM, NIC, NNC, ISM
    # donwnload human gtf file to compare against, check how sqanti does it.
    # also sent isONclust tsv file to this script to get cluster size 
    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )

    for (acc, (ins, del_, subs, matches)) in alignments_sorted:
        error_rate = (ins + del_ + subs) /float( (ins + del_ + subs + matches) ) 
        read_class = "NA" # annotations_dict[acc]
        nr_perfect_splice_sites = "NA" #annotations_dict[acc]
        cluster_size = reads_to_cluster_size[acc]
        read_length = len(reads[acc])
        info_tuple = (acc, read_type, ins, del_, subs, matches, error_rate, read_length, cluster_size, read_class, nr_perfect_splice_sites )
        outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format(*info_tuple))


def print_quantile_values(alignments_dict):

    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )
    print(alignments_sorted[-5:])

    sorted_error_rates = [ sum(tup[0:3])/float(sum(tup)) for acc, tup in alignments_sorted]

    insertions = [ tup[0]/float(sum(tup)) for acc, tup in alignments_sorted]
    deletions = [ tup[1]/float(sum(tup)) for acc, tup in alignments_sorted]
    substitutions = [ tup[2]/float(sum(tup)) for acc, tup in alignments_sorted]

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

    # quantile_errors = [sorted_error_rates[0], sorted_error_rates[int(n/20)], sorted_error_rates[int(n/10)], 
    #             sorted_error_rates[int(n/4)], sorted_error_rates[int(n/2)], 
    #             sorted_error_rates[int(3*n/4)], sorted_error_rates[int(9*n/10)], sorted_error_rates[int(19*n/20)], 
    #             sorted_error_rates[-1]]

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

def get_cluster_sizes(cluster_file):
    cluster_sizes = defaultdict(int)
    tmp_reads = {}
    for line in open(args.cluster_file, "r"):
        cl_id, acc = line.split() 
        cluster_sizes[cl_id] += 1  
        tmp_reads[acc] = cl_id

    reads_to_cluster_size = {}
    for acc, cl_id in tmp_reads.items():
        reads_to_cluster_size[acc] = cluster_sizes[cl_id]

    return reads_to_cluster_size


def get_splice_sites(cigar_tuples, first_exon_start, minimum_annotated_intron):
    splice_sites = []
    ref_pos = first_exon_start
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
            if l >= minimum_annotated_intron -1:
                # print("long del", l)
                splice_start = ref_pos
                ref_pos += l
                splice_stop = ref_pos
                splice_sites.append( (splice_start, splice_stop) )
            else:
                ref_pos += l

                # if l > 15:
                #     print("Large deletion!!!", l)


        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            splice_sites.append( (splice_start, splice_stop) )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return splice_sites

def get_read_splice_sites(sam_file, minimum_annotated_intron):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    read_splice_sites = {}
    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            # compare cs tag at intron sites
            q_cigar = read.cigarstring
            q_start = read.reference_start
            q_end = read.reference_end
            read_cigar_tuples = []
            result = re.split(r'[=DXSMINH]+', q_cigar)
            i = 0
            for length in result[:-1]:
                i += len(length)
                type_ = q_cigar[i]
                i += 1
                read_cigar_tuples.append((int(length), type_ ))  
            read_splice_sites[read.query_name] = {}  
            read_splice_sites[read.query_name][read.reference_name] = get_splice_sites(read_cigar_tuples, q_start, minimum_annotated_intron)
            # print("read", read_splice_sites[read.query_name][read.reference_name])
    return read_splice_sites


def get_annotated_splicesites(ref_gff_file):
    fn = gffutils.example_filename(ref_gff_file)
    db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    db = gffutils.FeatureDB('test.db', keep_order=True)
    
    # for tr in db.children(gene, featuretype='transcript', order_by='start'):
    #     # print(tr.id, dir(tr)) 
    #     splice_sites = []
    #     for j in db.children(tr, featuretype='exon', order_by='start'):
    #         # print(j, j.start, j.end)
    #         splice_sites.append(j.start -1)
    #         splice_sites.append(j.end)
    #     splice_sites_tmp = splice_sites[1:-1]
    #     splice_sites = []
    #     for i in range(0, len(splice_sites_tmp),2):
    #         splice_sites.append( (splice_sites_tmp[i], splice_sites_tmp[i+1]) )
    #     # splice_sites = [item for item in zip(splice_sites[:-1], splice_sites[1:])]
    #     ref_isoforms[tr.id] = splice_sites


    splice_coordinates = {} # to calc individual fraction of correct sites and NIC
    ref_isoforms = {} # To calculate Full splice matches

    # gene_graphs = {} # gene_id : { (exon_start, exon_stop) : set() }
    # collapsed_exon_to_transcript = {}
    minimum_annotated_intron = 1000000000
    for gene in db.features_of_type('gene'):
        chromosome = str(gene.seqid)
        if chromosome not in ref_isoforms:
            ref_isoforms[chromosome] = {}
        if chromosome not in splice_coordinates:
            splice_coordinates[chromosome] = set()

        
        #add splice sites
        for exon in db.children(gene, featuretype='exon', order_by='start'):
            exon_start, exon_end = exon.start - 1, exon.stop # double check if gff is 1-indexed here!
            splice_coordinates[chromosome].add(exon_start)
            splice_coordinates[chromosome].add(exon_end)
            
            # collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])
            # if (exon.start, exon.stop) in already_parsed_exons: 
            # gene_graph.add_node( (exon.start, exon.stop), weight=1  )
            # print(gene_graph.nodes[(exon.start, exon.stop)])

        #add annotated transcripts
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            # print('transcript', transcript.id, transcript.start, transcript.stop, [ (exon.start, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')])
            tmp_splice_sites = []
            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                tmp_splice_sites.append( (e1.stop, e2.start -1 ))

                if e2.start -1 - e1.stop < minimum_annotated_intron:
                    minimum_annotated_intron = e2.start -1 - e1.stop
                # print('exon', exon.id, exon.start, exon.stop)
            
            ref_isoforms[chromosome][tuple(tmp_splice_sites)] = transcript.id
      
    return ref_isoforms, splice_coordinates, minimum_annotated_intron


def get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, all_reads_splice_sites):
    total_true = 0
    total = 0
    total_reads = 0
    total_fsm = 0
    for read_acc in all_reads_splice_sites:
        total_reads += 1
        assert len(all_reads_splice_sites[read_acc]) == 1
        for chr_id in all_reads_splice_sites[read_acc]:
            # print(chr_id)
            # print(annotated_splice_coordinates)
            annotated_sites = annotated_splice_coordinates[chr_id]
            # print(chr_id)
            for read_splice_sites in all_reads_splice_sites[read_acc][chr_id]:
                # print(read_splice_sites)
                for sp in read_splice_sites:
                    # print(sp)
                    # print(annotated_sites)
                    total += 1
                    if sp in annotated_sites:
                        total_true += 1

            # check set intersection between read splice sites and annotated splice sites
            if tuple(all_reads_splice_sites[read_acc][chr_id]) in annotated_ref_isoforms[chr_id]:
                total_fsm += 1
                # print("FSM!!")
    # print(annotated_ref_isoforms[chr_id])
    # print( tuple(all_reads_splice_sites[read_acc][chr_id]))
    print("Total splice sizes found in cigar in reads:", total, "total matching annotations:", total_true, "total reads aligned:", total_reads)
    print("total FSM:", total_fsm)




def main(args):
    reads = { acc : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    corr_reads = { acc : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.corr_reads, 'r')))}
    refs = { acc : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}
    # print(refs)
    orig, orig_detailed = get_aln_stats_per_read(args.orig_sam, reads, refs, args)
    corr, corr_detailed = get_aln_stats_per_read(args.corr_sam, corr_reads, refs, args)

    if args.align:
        tot = 0
        for acc in reads:
            if acc not in corr_reads:
                print(acc, "Not in corrected reads")
            else:
                orig_seq = reads[acc]
                corr_seq = corr_reads[acc]
                read_alignment, ref_alignment = parasail_alignment(orig_seq, corr_seq)
                orig_dlen = [d_len for ch, d_len in  [(ch, len(list(_))) for ch, _ in itertools.groupby(read_alignment)] if ch == "-" ]
                corr_dlen = [d_len for ch, d_len in  [(ch, len(list(_))) for ch, _ in itertools.groupby(ref_alignment)] if ch == "-" ] 
                if corr_dlen:
                    max_orig_dlen = max( orig_dlen )
                else:
                    max_orig_dlen = 0
                if corr_dlen:
                    max_corr_dlen = max( corr_dlen )
                else:
                    max_corr_dlen = 0

                if max_corr_dlen > 10 or max_orig_dlen > 10:
                    print("orig", acc, read_alignment)
                    print(orig[acc])
                    print("corr", acc, ref_alignment)
                    print(corr[acc])
                    print(corr_detailed[acc][0])
                    print(corr_detailed[acc][1])
                    print(corr_detailed[acc][2])

                    tot += 1
                    print()
        print("TOT structural diffs:", tot)
    print( "Reads successfully aligned:", len(orig),len(corr))

    quantile_tot_orig, quantile_insertions_orig, quantile_deletions_orig, quantile_substitutions_orig = print_quantile_values(orig)
    quantile_tot_corr, quantile_insertions_corr, quantile_deletions_corr, quantile_substitutions_corr = print_quantile_values(corr)

    
    ## Splice site analysis
    annotated_ref_isoforms, annotated_splice_coordinates, minimum_annotated_intron = get_annotated_splicesites(args.gff_file)
    # print(annotated_ref_isoforms)
    # print(annotated_splice_coordinates)
    print("SHORTEST INTRON:", minimum_annotated_intron)
    corrected_splice_sites = get_read_splice_sites(args.corr_sam, minimum_annotated_intron)
    original_splice_sites = get_read_splice_sites(args.orig_sam, minimum_annotated_intron)
    corr_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, corrected_splice_sites)
    orig_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, original_splice_sites)

    annotations_dict_corrected = {} # get_annotation_of_reads(gff_file, corr_reads)
    annotations_dict_original = {} # get_annotation_of_reads(gff_file, reads)
    reads_to_cluster_size = get_cluster_sizes(args.cluster_file)

    detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read.csv"), "w")
    print_detailed_values_to_file(corr, annotations_dict_corrected, reads_to_cluster_size, reads, detailed_results_outfile, "corrected")    
    print_detailed_values_to_file(orig, annotations_dict_original, reads_to_cluster_size, corr_reads, detailed_results_outfile, "original")
    detailed_results_outfile.close()

    # alignments_stats = get_summary_stats(alignments_dict, 1.0)
    # print("Original reads (total): ins:{0}, del:{1}, subs:{2}, match:{3}".format(*alignments_stats[:-1]), "tot aligned region (ins+del+subs+match):", alignments_stats[-1] )
    # print("Original reads percent:{0}, del:{1}, subs:{2}, match:{3}".format(*[round(100*float(s)/alignments_stats[-1] , 1) for s in alignments_stats[:-1]]))
    # print( "Num_aligned_reads", "Aligned bases", "tot_errors", "avg_error_rate", "median_read_error_rate", "upper_25_quant", "lower_25_quant", "min", "max")


    orig_stats = get_summary_stats(orig, 1.0)
    # print("Original reads (total): ins:{0}, del:{1}, subs:{2}, match:{3}".format(*orig_stats[:-1]), "tot aligned region (ins+del+subs+match):", orig_stats[-1] )
    # print("Original reads percent:{0}, del:{1}, subs:{2}, match:{3}".format(*[round(100*float(s)/orig_stats[-1] , 1) for s in orig_stats[:-1]]))

    corr_stats = get_summary_stats(corr, 1.0)
    # print("Corrected reads (total): ins:{0}, del:{1}, subs:{2}, match:{3}".format(*corr_stats[:-1]), "tot aligned region (ins+del+subs+match):", corr_stats[-1])
    # print("Corrected reads percent:{0}, del:{1}, subs:{2}, match:{3}".format(*[round(100*float(s)/corr_stats[-1] , 1) for s in corr_stats[:-1]]))


    # print( "Num_aligned_reads", "Aligned bases", "tot_errors", "avg_error_rate", "median_read_error_rate", "upper_25_quant", "lower_25_quant", "min", "max")
    # orig_sorted = sorted(orig.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )
    # print(orig_sorted[-5:])

    # orig_sorted_error_rates = [ sum(tup[0:3])/float(sum(tup)) for acc, tup in orig_sorted]
    # orig_vals = [orig_sorted_error_rates[0], orig_sorted_error_rates[int(len(orig_sorted_error_rates)/20)], orig_sorted_error_rates[int(len(orig_sorted_error_rates)/10)], 
    #             orig_sorted_error_rates[int(len(orig_sorted_error_rates)/4)], orig_sorted_error_rates[int(len(orig_sorted_error_rates)/2)], 
    #             orig_sorted_error_rates[int(3*len(orig_sorted_error_rates)/4)], orig_sorted_error_rates[int(9*len(orig_sorted_error_rates)/10)], orig_sorted_error_rates[int(19*len(orig_sorted_error_rates)/20)], 
    #             orig_sorted_error_rates[-1]]


    # corr_sorted = sorted(corr.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1])) )
    # print(corr_sorted[-5:])
    # corr_sorted_error_rates = [ sum(tup[0:3])/float(sum(tup)) for acc, tup in corr_sorted]
    # corr_vals = [corr_sorted_error_rates[0],corr_sorted_error_rates[int(len(corr_sorted_error_rates)/20)], corr_sorted_error_rates[int(len(corr_sorted_error_rates)/10)],
    #              corr_sorted_error_rates[int(len(corr_sorted_error_rates)/4)], corr_sorted_error_rates[int(len(corr_sorted_error_rates)/2)],
    #               corr_sorted_error_rates[int(3*len(corr_sorted_error_rates)/4)], corr_sorted_error_rates[int(9*len(corr_sorted_error_rates)/10)], corr_sorted_error_rates[int(19*len(corr_sorted_error_rates)/20)],
    #               corr_sorted_error_rates[-1]]
    
    print("Distribution of error rates (Percent)")
    print("Reads, Best, top 5%, top 10%, top 25%, Median, top 75%, top 90%, top 95%, Worst")
    print("Original,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_orig ] ))
    print("Corrected,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_corr ] ))

    outfile = open(os.path.join(args.outfolder, "results.csv"), "w")
    outfile.write("Original,tot,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_tot_orig ], *orig_stats, len(orig)))
    outfile.write("Corrected,tot,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_tot_corr ], *corr_stats, len(corr)))
    
    outfile.write("Original,ins,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_insertions_orig ], *orig_stats, len(orig)))
    outfile.write("Corrected,ins,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_insertions_corr ], *corr_stats, len(corr)))

    outfile.write("Original,del,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_deletions_orig ], *orig_stats, len(orig)))
    outfile.write("Corrected,del,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_deletions_corr ], *corr_stats, len(corr)))

    outfile.write("Original,subs,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_substitutions_orig ], *orig_stats, len(orig)))
    outfile.write("Corrected,subs,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_substitutions_corr ], *corr_stats, len(corr)))

    outfile.close()

    print("Reads successfully aligned:", len(orig),len(corr))

    # print(orig_sorted)
    # print(",".join([s for s in orig_stats]))
    # print(",".join([s for s in corr_stats]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('orig_sam', type=str, help='Path to the original read file')
    parser.add_argument('corr_sam', type=str, help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('corr_reads', type=str, help='Path to the corrected read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('cluster_file', type=str, help='Path to the refs file')
    parser.add_argument('gff_file', type=str, help='Path to the refs file')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('--align', action= "store_true", help='Include pairwise alignment of original and corrected read.')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)