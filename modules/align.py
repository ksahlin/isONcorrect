import sys
import os
import tempfile
import subprocess
import re
import math
from collections import deque

import pysam
import parasail


from modules import help_functions

phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.5)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
phred_char_to_val = {chr(i) : max( (ord(chr(i)) - 33), 3)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing


def minimap2(reference, reads, outfile):
    print('Aligning with minimap2.')
    sys.stdout.flush()
    # work_dir = "/tmp/" #tempfile.mkdtemp() 
    # print(work_dir)
    stderr_file = open(reference + ".minimap.stderr", 'w')
    # print('Stderr file: ', stderr_file)
    sys.stdout.flush()
    # print(type(reference))
    processes = 1 # mp.cpu_count()
    print("minimap with {0} processes.".format(processes))
    # print("minimap -f 0.00000001 -Sw5 -L40 -m0 -t {0} {1} {2}".format(processes, reference, reads))
    print("minimap2 -ax  splice -uf -k14 --cs=long -t {0} {1} {2}".format(processes, reference, reads))

    with open(outfile, "w") as minimap_file:
        sys.stdout.flush()
        subprocess.check_call([ "minimap2", "-ax", "splice", "-un", "-k13", "-w 1", "-n1", "-m 80", "-r 1000", "-s40", "-O4,24", "-E2,1",
                                "-t", str(processes),
                               reference, reads ],
                                stdout=minimap_file,
                                stderr=stderr_file)
        sys.stdout.flush()

    return outfile


def align_mm2(reads, ref, args):

    ref_file = os.path.join(args.outfolder, "augmented_ref_{0}.fa".format(args.iteration))
    r = open(ref_file, "w")
    r.write(">{0}\n{1}".format("ref", ref))
    r.close()
    minimap2_out = os.path.join(args.outfolder, "mm2_aligned_{0}.sam".format(args.iteration))
    sam_file = minimap2(ref_file, reads,  minimap2_out)
    
    reads = { i : (acc, seq, qual) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(reads,"r")))}
    reads_acc_to_index = { acc : i for i, (acc, seq, qual) in reads.items()}
    k = args.k
    best_matches = {}
    current_min_ed = {}
    print()
    print()

    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    alignments = {}
    for read in SAM_file.fetch(until_eof=True):
        read_index = reads_acc_to_index[read.qname]
        if read.is_reverse:
            # print("read {0} is reverse complemented aligned to {1}.".format(read.query_name, references[read.reference_id]))
            # print(read.tostring(SAM_file))
            # sys.exit()
            continue
        if read.is_unmapped:
            acc, read_seq, qual = reads[read_index]
            print(read.qname, "is umapped", len(read_seq))
            alignments[read_index] = (acc, "unaligned", "unaligned", [])
            continue


        acc, read_seq, qual = reads[read_index]
        # cigar = read.cigarstring
        # print(read.qname, len(read_seq))
        ref_alignment, m_line, read_alignment = cigar_to_seq_mm2(read, ref, read_seq)
        match_id = math.floor((1.0 - 0.16) * args.k)
        block_coverage = get_block_coverage(read_alignment, ref_alignment, k, match_id)

        # print(ref_alignment)
        # print(read_alignment)
        # print(block_coverage)
        alignments[read_index] = (acc, read_alignment, ref_alignment, block_coverage)
    return alignments


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


def get_block_vector(match_vector, k, match_id):
    match_window = deque(match_vector[:k]) # initialization
    current_match_count = sum(match_window)
    aligned_region = []
    if current_match_count >= match_id:
        aligned_region.append(1)
    else:
        aligned_region.append(0)


    for new_m_state in match_vector[k:]:
        prev_m_state = match_window.popleft()
        current_match_count = current_match_count - prev_m_state + new_m_state 
        match_window.append(new_m_state)
        
        if current_match_count >= match_id:
            aligned_region.append(1)
        else:        
            aligned_region.append(0)

    return aligned_region

def get_block_coverage(read_alignment, ref_alignment, k, match_id):
    match_vector = [ 1 if n1 == n2 else 0 for n1, n2 in zip(read_alignment, ref_alignment) ]
    aligned_region = get_block_vector(match_vector, k, match_id)
    aligned_region_reverse = get_block_vector(match_vector[::-1], k, match_id)
    block_coverage_f = [m for m in aligned_region + [aligned_region[-1]]*(k-1) ]
    block_coverage_reverse = [m for m in aligned_region_reverse + [aligned_region_reverse[-1]]*(k-1) ]
    assert len(block_coverage_f) == len(block_coverage_reverse)
    block_coverage = [ 1 if n1 +n2 >0 else 0 for n1, n2 in zip(block_coverage_f, block_coverage_reverse[::-1])]
    return block_coverage

import re
def get_block_coverage2(read_alignment, ref_alignment, k, match_id):
    match_vector = "".join([ "-" if n1 == "-" or n2 == "-" else "|" for n1, n2 in zip(read_alignment, ref_alignment) ])
    p=re.compile("[-]+")
    # all_indels = re.findall(p, match_vector)
    all_indels = []
    for m in p.finditer(match_vector):
        if len(m.group()) >= 10:
            all_indels.append((m.start(), len(m.group())))
            print(m.start(), m.group())


    block_coverage = [ 1 for i in range(len(match_vector))]
    prev_stop = 0
    for pos, i_len in all_indels:
        if pos - prev_stop < 7:
            s_pos = prev_stop
        else:
            s_pos = pos

        for i in range(i_len):
            block_coverage[s_pos+i] = 0
        
        prev_stop = s_pos + i_len

    # aligned_region = get_block_vector(match_vector, k, match_id)
    # aligned_region_reverse = get_block_vector(match_vector[::-1], k, match_id)
    # block_coverage_f = [m for m in aligned_region + [aligned_region[-1]]*(k-1) ]
    # block_coverage_reverse = [m for m in aligned_region_reverse + [aligned_region_reverse[-1]]*(k-1) ]
    # assert len(block_coverage_f) == len(block_coverage_reverse)
    # block_coverage = [ 1 if n1 +n2 >0 else 0 for n1, n2 in zip(block_coverage_f, block_coverage_reverse[::-1])]
    return block_coverage


def parasail_block_alignment(read, reference, k, match_id, x_acc = "", y_acc = "", match_score = 4, mismatch_penalty = -4, opening_penalty = 5, gap_ext = 1, ends_discrepancy_threshold = 0):
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
    block_coverage = get_block_coverage(read_alignment, ref_alignment, k, match_id)
    block_coverage2 = get_block_coverage2(read_alignment, ref_alignment, k, match_id)
    print(block_coverage)
    print(block_coverage2)
    print()
    return read, reference, read_alignment, ref_alignment, block_coverage2


def block_align_parasail(reference, read, reference_qual, read_qual, args):
    poisson_mean = sum([ read_qual.count(char_) * phred_char_to_p[char_] for char_ in set(read_qual)])
    poisson_mean2 = sum([ reference_qual.count(char_) * phred_char_to_p[char_] for char_ in set(reference_qual)])
    error_rate_sum = poisson_mean/float(len(read)) + poisson_mean2/float(len(reference))  # k = max(int(mean_plus_two_stdvs_q2 + mean_plus_two_stdvs_q1) + 1 + int(len(read)*args.variant_rate) , 40)
    # if error_rate_sum <= 0.01:
    #     gap_opening_penalty = 5
    # elif  0.01 < error_rate_sum <= 0.04:
    #     gap_opening_penalty = 4
    # elif  0.04 < error_rate_sum <= 0.1:
    #     gap_opening_penalty = 3
    # elif  0.1 < error_rate_sum:
    #     gap_opening_penalty = 5
    # elif  0.2 < error_rate_sum:
    #     gap_opening_penalty = 1
    
    gap_opening_penalty = 5
    gap_ext = 0
    # print(error_rate_sum)
    # gap_opening_penalty = 2
    match_id_tailored = math.floor((1.0 - error_rate_sum) * args.k)
    match_id_tailored = math.floor((1.0 - 0.16) * args.k)
    read, ref, read_alignment, ref_alignment, block_coverage  = parasail_block_alignment(read, reference, args.k, \
                                                                                        match_id_tailored, opening_penalty = gap_opening_penalty, gap_ext = gap_ext,\
                                                                                        match_score = 5, mismatch_penalty = -4 )
    # print()
    # print("Expected errors:", poisson_mean, poisson_mean2)
    # print(read_alignment)
    # print(block_coverage)
    # print(ref_alignment)
    # print()
    return error_rate_sum, read_alignment, ref_alignment, block_coverage



def align_parasail(reads, ref, args):
    reads = { i : (acc, seq, qual) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(reads,"r")))}
    ref_qual = "".join(["+" for i in range(len(ref))])
    alignments = {}
    for j in range(0, len(reads)):
        acc, read_seq, read_qual = reads[j]
        error_rate_sum, read_alignment, ref_alignment, block_coverage = block_align_parasail(ref, read_seq, ref_qual, read_qual, args)
        # print(acc)
        # print(ref_alignment)
        # print(read_alignment)
        # print("".join([str(m) for m in block_coverage]))

        alignments[j] = (acc, read_alignment, ref_alignment, block_coverage)
            
    return alignments






