import edlib
import re
import math

from collections import deque
import itertools

def annotate_with_quality_values(alignment_matrix, seq_to_acc, qual_dict):
    alignment_matrix_of_qualities = {}
    alignment_matrix_of_max_qualities = {}

    for s in alignment_matrix:
        s_accessions = seq_to_acc[s]
        all_quals = [ qual_dict[s_acc] for s_acc in s_accessions ]
        sum_quals_vector = [sum(t) for t in zip(*all_quals)]
        max_quals_vector = [max(t) for t in zip(*all_quals)]
        # sum all probabilities from all reads equal to s here using all accessions in s to acc s_to_acc 
        
        list_sum_quals = []
        list_max_quals = []
        current_pos_in_s = 0

        for j in range(len(alignment_matrix[s])):
            # print(current_pos_in_s, len(sum_quals_vector) )
            current_quality = sum_quals_vector[current_pos_in_s]
            current_max_quality = max_quals_vector[current_pos_in_s]

            list_sum_quals.append(current_quality)
            list_max_quals.append(current_max_quality)

            char_at_pos = alignment_matrix[s][j]
            if char_at_pos != "-":
                if current_pos_in_s < len(sum_quals_vector) - 1:
                    current_pos_in_s += 1 

        alignment_matrix_of_qualities[s] = list_sum_quals
        alignment_matrix_of_max_qualities[s] = list_max_quals

    PFM_qualities = []
    PFM_max_qualities = []
    for j in range(len(alignment_matrix[s])): # for each column
        PFM_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "U": 0, "-": 0})
        PFM_max_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "U": 0, "-": 0})
        for s in alignment_matrix_of_qualities:
            nucl = alignment_matrix[s][j]
            sum_quality_at_position = alignment_matrix_of_qualities[s][j]
            PFM_qualities[j][nucl] += sum_quality_at_position

            max_quality_at_position = alignment_matrix_of_max_qualities[s][j]
            PFM_max_qualities[j][nucl] += max_quality_at_position



    # get all the differences to majority here  
    list_of_majority_nucleotides = []
    for j in range(len(PFM_qualities)):
        max_v_j = max(PFM_qualities[j], key = lambda x: PFM_qualities[j][x] )
        majority_count =  PFM_qualities[j][max_v_j]
        max_v_j_set = set([v for v in PFM_qualities[j] if PFM_qualities[j][v] == majority_count ])
        all_major = "".join(max_v_j_set)
        list_of_majority_nucleotides.append(all_major)

    assert len(list_of_majority_nucleotides) == len(PFM_qualities)
    global_all_difference_qualities = []    
    for s in alignment_matrix:
        s_aligned_vector = alignment_matrix[s]
        for j in range(len(s_aligned_vector)):
            if s_aligned_vector[j] not in list_of_majority_nucleotides[j] and len(list_of_majority_nucleotides[j]) == 1:
                global_all_difference_qualities.append(alignment_matrix_of_qualities[s][j])

    global_all_difference_qualities.sort()
    if len(global_all_difference_qualities) > 0:
        global_correction_threshold = global_all_difference_qualities[ int( math.ceil( len(global_all_difference_qualities)/2.0) ) - 1 ]
    else:
        global_correction_threshold = -1
        print("nothing to correct!")
    print("GLOBAL QUAL THRESH:", global_correction_threshold)
    return alignment_matrix_of_qualities, PFM_qualities, PFM_max_qualities, global_correction_threshold



def create_position_frequency_matrix(alignment_matrix, partition):
    nr_columns = len( alignment_matrix[ list(alignment_matrix)[0] ]) # just pick a key
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "U" : 0, "-": 0} for j in range(nr_columns)]
    for s in alignment_matrix:
        s_aln = alignment_matrix[s]
        indegree = partition[s][3]
        for j in range(nr_columns): # for each column
            nucl = s_aln[j]
            PFM[j][nucl] += indegree
    return PFM


def min_ed(max_insertion, q_ins):
    result = edlib.align(max_insertion, q_ins, task="path", mode= "NW")
    cigar = result["cigar"]
    tuples = []
    # do not allow deletions in max_insertion: because we do not want to alter this sequence
    if "D" in cigar:
        return ""
    matches = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in matches[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        tuples.append((int(length), type_ ))

    q_insertion_modified = ""
    q_ins_pos = 0
    for length, type_ in tuples:
        # if we reach here we are guaranteed no deletions in alignment of max_insertion
        # we therefore simply thread in the matching or mismatching characters (if type '=' or 'X')
        # or we put a "-" (if type is 'I')
        if type_ == "I":
            q_insertion_modified += "-"*length
        else:
            q_insertion_modified += q_ins[q_ins_pos : q_ins_pos + length]
            q_ins_pos += length
    return q_insertion_modified



def get_best_solution(max_insertion, q_ins):

    if q_ins == "-":
        return ["-" for i in range(len(max_insertion))]

    else:
        pos = max_insertion.find(q_ins) 
        if pos >= 0:    
            q_insertion_modified = "-"*pos + max_insertion[ pos : pos + len(q_ins) ] + "-"* len(max_insertion[ pos + len(q_ins) : ])
            return [character for character in q_insertion_modified]        

        # else, check if smaller deletion can be aligned from left to write, e.g. say max deletion is GACG
        # then an insertion AG may be aligned as -A-G. Take this alignment instead
        q_insertion_modified = min_ed(max_insertion, q_ins)
        if q_insertion_modified:
            return [character for character in q_insertion_modified]

        # otherwise just shift left
        # check if there is at least one matching character we could align to
        max_p = 0
        max_matches = 0
        for p in range(0, len(max_insertion) - len(q_ins) + 1 ):
            nr_matches = len([1 for c1, c2 in zip(q_ins, max_insertion[p: p + len(q_ins) ] ) if c1 == c2])
            if nr_matches > max_matches:
                max_p = p
                max_matches = nr_matches

        if max_p > 0:
            q_insertion_modified = "-"*max_p + q_ins + "-"*len(max_insertion[max_p + len(q_ins) : ])
            # print("specially solved: q:{0} max:{1} ".format(q_insertion_modified, max_insertion) )
            return [character for character in q_insertion_modified]


        q_insertion_modified = []
        for p in range(len(max_insertion)):
            # all shorter insertions are left shifted -- identical indels are guaranteed to be aligned
            # however, no multialignment is performed among the indels
            if p < len(q_ins):
                q_insertion_modified.append(q_ins[p])
            else:
                q_insertion_modified.append("-")
        return [character for character in q_insertion_modified]




def create_multialignment_format_NEW(query_to_target_positioned_dict, start, stop):
    """
            1. From segments datastructure, get a list (coordinate in MAM) of lists insertions that contains all the insertions in that position
            2. Make data structure unique_indels =  set(insertions) to get all unique insertions in a given position
            3. Make a list max_insertions containing the lengths of each max_indel in (max_indels >1bp should be padded). From this list we can get the total_matrix_length
            4. In a data structure position_solutions : {max_indel : {q_ins : solution_list, q_ins2: }, }, 
                find oplimal solution for each unique indel in unique_indels to max_indel in datastructure  (a list of dicts in each position) 
            5. for each pos in range(total_matrix_length):
                1. for q_acc in segments:
                    1. [ ] if pos odd:
                        1. [ ] trivial
                    2. [ ] elif pos even and max_indel == 1:
                        1. [ ] trivial
                    3. [ ] else:
                        1. [ ] [alignmet[q_acc].append(char) for char in solution]
    """
    assert len(query_to_target_positioned_dict) > 0
    target_vector_length = len( list(query_to_target_positioned_dict.values())[0][0])
    assert stop < target_vector_length # vector coordinates are 0-indexed

    # get allreads alignments covering interesting segment
    # row_accessions = []
    segments = {}
    segment_lists = []
    for q_acc, (q_to_t_pos, t_vector_start, t_vector_end) in query_to_target_positioned_dict.items():
        if t_vector_start <= start and t_vector_end >= stop: # read cover region
            segment = q_to_t_pos[start - t_vector_start : stop - t_vector_start +1 ]
            # row_accessions.append(q_acc)
            segments[q_acc] = segment
            segment_lists.append(segment)

    # 1
    unique_insertions = [set(i) for i in zip(*segment_lists)]
    # unique_insertions = [set() for i in range(len(segment))]
    nr_pos = len(segments[q_acc])
    # for q_acc in segments:
    #     segment = segments[q_acc]
    #     [ unique_insertions[p].add(segment[p]) for p in range(nr_pos) ]


    # 2
    # unique_insertions = []
    # for p in range(nr_pos):
    #     unique_insertions.append(set(position_insertions[p]))

    # 3
    max_insertions = []
    for p in range(nr_pos):
        max_ins_len = len(max(unique_insertions[p], key=len))
        if max_ins_len > 1:
            max_ins = sorted([ins for ins in unique_insertions[p] if len(ins) == max_ins_len ])[0]
            assert p % 2 == 0 # has to be between positions in reference
            max_ins =  "-" + max_ins + "-" 
            max_insertions.append(max_ins)    
        else:
            max_insertions.append("-")    

    # total_matrix_length = sum([len(max_ins) for max_ins in max_insertions]) # we now have all info to get total_matrix_length  

    # 4
    position_solutions = {max_ins : {} for max_ins in max_insertions}

    for nucl in ["A", "G", "C", "T", "U", "-"]:
        position_solutions[nucl] = {"A": ["A"], "G": ["G"], "C": ["C"], "T": ["T"], "U": ["U"], "-": ["-"], "N": ["N"]}

    for p in range(nr_pos):
        max_ins = max_insertions[p]
        if len(max_ins) > 1:
            for ins in unique_insertions[p]:
                if ins not in position_solutions[max_ins]:
                    position_solutions[max_ins][ins] = get_best_solution(max_ins, ins)

    # 5 create alignment matrix
    alignment_matrix = {}
    for q_acc in segments:
        alignment_matrix[q_acc] = []
        # tmp_list = []
        seqs = segments[q_acc]
        tmp_list = [ position_solutions[max_insertions[p]][seqs[p]] for p in range(nr_pos)]

        # for p in range(nr_pos):
        #     max_ins = max_insertions[p]
        #     seq = seqs[p]
        #     tmp_list.append(position_solutions[max_ins][seq])

        alignment_matrix[q_acc] = [character for sublist in tmp_list for character in sublist]

    # assert total_matrix_length == len(alignment_matrix[q_acc])
    return alignment_matrix


def create_multialignment_matrix(partition):
    """
        a partition is a dictionary of pairwise alignments for a given center repr_seq. "partition has the following
        structure:  partition = {s : (edit_distance, m_alignment, s_alignment, degree_of_s)}
        s can only have a degree larger than 1 if s=repr_seq, otherwise s has a degree of 1.

        This function does the following

        query_to_target_positioned_dict = {}
        for each seq in partition we call function
            position_query_to_alignment(query_aligned, target_aligned, target_start=0)
            returns the following data
            query_to_target_positioned, target_vector_start_position = 0, target_vector_end_position 
            we store all these pairwise alignments in query_to_target_positioned_dict

        where this function retuns a dictionary with query seq as key in the following form:
        query_to_target_positioned_dict[query_accession] = (query_to_target_positioned, target_vector_start_position, target_vector_end_position)

        then it calls create_multialignment_format(query_to_target_positioned_dict, start, stop)
        This function returns an alignment matric in the following smaple format:
        alignment_matrix = {"q1" : ["-", "A","-", "C", "-", "G","A","C","C","G", "G", "-", "A", "T","T","T"],
                            "q2" : ["-", "A","-", "C", "-", "G","A","G","-","-", "G", "-", "A", "T","T","T"],
                            "q3" : ["-", "A","-", "C", "-", "G","A","-","-","-", "G", "-", "A", "T","T","T"],
                            "q4" : ["-", "A","-", "C", "-", "G","C","C","-","-", "G", "-", "A", "-","-","-"],
                            "q5" : ["-", "A","-", "C", "-", "G","-","-","-","-", "G", "-", "A", "T","-","-"],
                            "q6" : ["G", "A","-", "C", "-", "G","C","-","-","-", "G", "-", "A", "-","-","-"]
                            }

        finally, we transform the alignment_matrix into a PFM by multiplying each query sequnece with the correct degree.

        PFM is a list of dicts, where the inner dicts are one per column in the PFM matrix, its keys by characters A, C, G, T, -
        alignment_matrix is a representation of all alignments in a partition. this is a dictionary where sequences s_i belonging to the 
        partition as keys and the alignment of s_i with respect to the alignment matix.
    """
    repr_seq = partition["ref"][1]
    query_to_target_positioned_dict = {}
    for q_acc in partition:
        (edit_distance, m_alignment, s_alignment, degree_of_s) = partition[q_acc]
        # s_positioned_new, target_vector_start_position, target_vector_end_position = position_query_to_alignment_NEW(s_alignment, m_alignment, 0)
        # print(s_alignment)
        # print(m_alignment)
        # print(repr_seq)
        s_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(s_alignment, m_alignment, 0)
        # if s_alignment == m_alignment:
        #     print("POS:",s_positioned)
        #     sys.exit()

        # assert s_positioned_new == s_positioned
        # print(s_positioned)
        # print()
        assert target_vector_start_position == 0
        # print(target_vector_end_position + 1, len(repr_seq), 2*len(repr_seq)+1)
        assert target_vector_end_position + 1 == 2*len(repr_seq) + 1 # vector positions are 0-indexed
        query_to_target_positioned_dict[q_acc] = (s_positioned, target_vector_start_position, target_vector_end_position)

    alignment_matrix = create_multialignment_format_NEW(query_to_target_positioned_dict, 0, 2*len(repr_seq))
    # print(alignment_matrix)
    # sys.exit()
    return alignment_matrix


def position_query_to_alignment(query_aligned, target_aligned, target_alignment_start_position):
    """
        input:      0-indexed target positions
        returns:    target vector is a list of 2*t_len +1 positions to keep insertions between base pairs    
                list of strings of nucleotides for each position, start position in target vector list, end position in target vector list
    """

    query_positioned = [] 
    target_position = target_alignment_start_position # 0-indexed
    temp_ins = ""
    # iterating over alignment positions
    for p in range(len(target_aligned)):
        if target_aligned[p] == "-":
            temp_ins += query_aligned[p]
        else:
            if not temp_ins:
                query_positioned.append("-") #[2*target_position] = "-"
            else:
                query_positioned.append(temp_ins)
                temp_ins = ""

            query_positioned.append(query_aligned[p])

            target_position += 1

    if not temp_ins:
        query_positioned.append("-")
    else:
        query_positioned.append(temp_ins)

    target_vector_start_position = 2*target_alignment_start_position
    target_vector_end_position = 2*(target_position-1) + 2

    return query_positioned, target_vector_start_position, target_vector_end_position



def get_block_vector(match_vector, k, match_id):
    # initialization
    pos = 0
    valid_positions = 0
    current_match_count = 0
    while valid_positions <=k:
        if match_vector[pos] == -1:
            pass
        else:
            current_match_count += match_vector[pos]
            valid_positions +=1
        pos +=1
    match_window = deque([ val for val in match_vector[:pos] if val != -1])

    # match_window = deque(match_vector[:k]) # initialization
    # current_match_count = sum([m for m in match_window)
    aligned_region = []
    if current_match_count >= match_id:
        aligned_region.append(1)
    else:
        aligned_region.append(0)
    curr_window_state = 1 if current_match_count >= match_id else 0
    for new_m_state in match_vector[pos:]:
        if new_m_state == -1:
            aligned_region.append(curr_window_state)
        else:
            prev_m_state = match_window.popleft()
            current_match_count = current_match_count - prev_m_state + new_m_state 
            curr_window_state = 1 if current_match_count >= match_id else 0
            match_window.append(new_m_state)
            aligned_region.append(curr_window_state)

    # print(len(aligned_region))
    aligned_region = aligned_region + [m for m in [aligned_region[-1]]*(pos-1) ]
    # print(match_vector)
    # print(aligned_region)
    # print(len(aligned_region),len(match_vector))
    assert len(aligned_region) == len(match_vector)
    return aligned_region


def get_block_coverage(read_alignment, ref_alignment, k, match_id):
    match_vector = [ 0 if n1 != n2 else -1 if n2 == "-" else 1 for n1, n2 in zip(read_alignment, ref_alignment) ]
    block_coverage_f = get_block_vector(match_vector, k, match_id)
    block_coverage_reverse = get_block_vector(match_vector[::-1], k, match_id)
    assert len(block_coverage_f) == len(block_coverage_reverse)
    block_coverage = [ 1 if n1 + n2 > 0 else 0 for n1, n2 in zip(block_coverage_f, block_coverage_reverse[::-1])]
    return block_coverage


def get_block_coverage2(read_alignment, ref_alignment, k, match_id):
    ref_nucl_since_gap_open = 0
    deteted_blocks = []
    start_coord = 0
    for i, (n1, n2) in enumerate(zip(read_alignment, ref_alignment)):
        if n2 == "-":
            continue
        elif n1 == "-":
            ref_nucl_since_gap_open +=1
        else:
            if ref_nucl_since_gap_open >= 10:
                deteted_blocks.append((start_coord, i))
            ref_nucl_since_gap_open = 0
            start_coord = i+1

    block_coverage = [ 1 for i in range(len(read_alignment))]
    prev_stop = 0
    for start, stop in deteted_blocks:
        if start - prev_stop < 14:
            start_pos = prev_stop
        else:
            start_pos = start

        for i in range(start_pos, stop):
            block_coverage[i] = 0
        prev_stop = stop

    return block_coverage



def get_homopolymer_factor(ref_aln):
    homopol_correction_vector = []
    ref_aln_to_str = "".join([n for n in ref_aln if n != "-"])
    str_to_list_index = [i for i, n in enumerate(ref_aln) if n != "-" ]
    # print(ref_aln_to_str)
    # print(str_to_list_index)
    # print([(ch, sum(1 for _ in g)) for ch, g in itertools.groupby(ref_aln_to_str)])
    curr_pos = 0
    start_coord = 0
    for ch, g in itertools.groupby(ref_aln_to_str):
        h_len =  sum(1 for _ in g)
        stop_coord = str_to_list_index[curr_pos + h_len -1]
        # print(start_coord, stop_coord)
        for i in range(stop_coord - start_coord +1):
            homopol_correction_vector.append(h_len)
        curr_pos += h_len
        start_coord = stop_coord +1

    # last positions corner case
    for i in range(len(ref_aln) - len(homopol_correction_vector)):
        homopol_correction_vector.append(h_len)

    # print(ref_aln)
    # print(homopol_correction_vector)
    # print(len(homopol_correction_vector), len(ref_aln))
    return homopol_correction_vector


def blocks_from_msa(partition, k, match_id):
    block_matrix = {}
    ref_aln = partition["ref"][0]
    for acc in partition:
        if acc == "ref":
            homopol_correction_vector = get_homopolymer_factor(ref_aln)
        else:
            seq_aln, cnt = partition[acc]
            # block_vector = get_block_coverage(seq_aln, ref_aln, k, match_id)
            block_vector = get_block_coverage2(seq_aln, ref_aln, k, match_id)
            # print("".join([str(b) for b in block_vector]))
            # print("".join([str(b) for b in block_vector2]))
            # print("".join([str(b) for b in ref_aln]))
            # print("".join([str(b) for b in seq_aln]))
            # print()
            block_matrix[acc] = (block_vector, cnt)
    return block_matrix, homopol_correction_vector


def get_block_freqs_and_majority(block_matrix, partition):
    nr_columns = len( list(block_matrix.values())[0][0])
    block_frequency_matrix = [{"A": 0, "C": 0, "G": 0, "T": 0, "U" : 0, "-": 0} for j in range(nr_columns)]
    BLOCK_FREQ = [0 for j in range(nr_columns)]
    for s in block_matrix:
        read_aln_vector, _  = partition[s]
        block_v, cnt = block_matrix[s]
        for j in range(nr_columns):
            BLOCK_FREQ[j] += cnt*block_v[j]
            if block_v[j] == 1:
                n = read_aln_vector[j]
                block_frequency_matrix[j][n] += cnt

    block_majority_vector = [max(d.items(), key=lambda x: x[1])[0] for d in block_frequency_matrix]
    # for j in range(nr_columns):
    #     print(j, BLOCK_FREQ[j], block_majority_vector[j], block_frequency_matrix[j])

    # print("".join([n for n in block_majority_vector if n != "-"]))
    # print([(n, BLOCK_FREQ[i]) for i, n in enumerate(block_majority_vector)])
    # sys.exit()
    return BLOCK_FREQ, block_majority_vector, block_frequency_matrix

def get_global_probs(read_errors):
    tot_ins = sum([ read_errors[seq][0] for seq in read_errors])
    tot_del = sum([ read_errors[seq][1] for seq in read_errors])
    tot_subs = sum([ read_errors[seq][2] for seq in read_errors])
    tot_aln_length = sum([ read_errors[seq][3] for seq in read_errors])
    print(tot_ins,tot_del, tot_subs, tot_aln_length)
    p_ins, p_del, p_subs = tot_ins/float(tot_aln_length), tot_del/float(tot_aln_length), tot_subs/float(tot_aln_length)
    print(p_ins, p_del, p_subs)
    # sys.exit()
    return p_ins, p_del, p_subs

def correct_to_consensus(partition, read_errors, args):
    """
         partition[acc] = (edit_dist, aln_rep, aln_s, depth_of_string)
    """
    N_t = sum([container_tuple[3] for acc, container_tuple in partition.items()]) # total number of sequences in partition
    if len(partition) > 1:
        # all strings has not converged
        alignment_matrix = create_multialignment_matrix(partition) 
        partition = { acc : (alignment_matrix[acc], partition[acc][3]) for acc in alignment_matrix}
        PFM = PFM_from_msa(partition)
        # PFM = create_position_frequency_matrix(alignment_matrix, partition)
        # for i,pos_dict in enumerate(PFM):
        #     print(i, pos_dict)
        p_ins, p_del, p_subs = get_global_probs(read_errors)
        tot_prob = p_ins + p_del + p_subs 
        print(tot_prob,args.k, int( (1- tot_prob)*args.k))
        match_id = int( (1- tot_prob)*args.k)
        block_matrix, homopol_correction_vector = blocks_from_msa(partition, args.k, match_id)
        block_freq_vector, block_majority_vector, BFM = get_block_freqs_and_majority(block_matrix, partition)
        pos_cutoffs = [{} for p in range(len(block_freq_vector))]
        # print(block_freq_vector)
        # print(homopol_correction_vector)
        for i, c in enumerate(block_freq_vector): # poisson rates
            p_del_pos = 1.0 - (1.0 - p_del)**homopol_correction_vector[i]
            
            if homopol_correction_vector[i] > 1:
                p_del_lambda = 1.0 - (1.0 - p_del)**homopol_correction_vector[i]
                p_corr_lambda = (1.0 - p_del)**(homopol_correction_vector[i] - 1)
                pos_cutoffs[i]["h"] = (c*p_del_lambda, c*p_corr_lambda)

            
            # if (c * p_del_pos) + 3*math.sqrt((c * p_del_pos))> 50:
            #     print(i,c, homopol_correction_vector[i], (c * p_del_pos) + 3*math.sqrt((c * p_del_pos)))
            # print(p_del_pos)
            # if i == 628:
            #     print( c,  max(1, (c * p_del) + 3*math.sqrt((c * p_del)) ),  max(1, (c * p_subs) + 3*math.sqrt((c * p_subs)) ),max(1, (c * p_ins) + 3*math.sqrt((c * p_ins)) ) )
            pos_cutoffs[i]["d"] = max(1, (c * p_del_pos) + 3*math.sqrt((c * p_del_pos)) )
            pos_cutoffs[i]["mm"] = max(1, (c * p_subs) + 3*math.sqrt((c * p_subs)) )
            pos_cutoffs[i]["i"] = max(1, (c * p_ins) + 3*math.sqrt((c * p_ins)) )

        # print(BFM)
        # print(block_freq_vector)
        # print(pos_cutoffs)
        # sys.exit()
        
        # print(pos_cutoffs)
        S_prime_partition = correct_from_msa(partition, BFM, pos_cutoffs = pos_cutoffs, block_matrix = block_matrix, block_majority_vector= block_majority_vector, homopol_correction_vector = homopol_correction_vector)
        # sys.exit()
    else:
        print("Partition converged: Partition size(unique strings):{0}, partition support: {1}.".format(len(partition), N_t))

    
    # only_deletions_pos = set([p for p, d in enumerate(PFM) if d["A"] == d["C"] == d["G"] == d["T"] == 0] )
    # print("only deleted pos:", only_deletions_pos)
    # for seq, almnt in alignment_matrix.items():
    #     if seq not in seq_to_acc:
    #         print(">augmented_ref")
    #     else:
    #         print(">{0}".format(seq_to_acc[seq][0]))
    #     print("".join([n for p, n in enumerate(almnt) if p not in only_deletions_pos]))
    # sys.exit()  

    return S_prime_partition



def PFM_from_msa(partition):
    nr_columns = len( list(partition.values())[0][0]) # just pick a key
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "U" : 0, "-": 0} for j in range(nr_columns)]
    for s in partition:
        seq_aln, cnt = partition[s]
        for j, n in enumerate(seq_aln):
            PFM[j][n] += cnt

    # cov_tot = len(partition)
    # for j in range(nr_columns):
    #     cov = sum([ PFM[j][n] for n in PFM[j] if n != "-"])
    #     print(j,cov, cov_tot, [ PFM[j][n] for n in PFM[j] if n != "-"])
    # print("A", "".join([str(min(p["A"], 9)) for p in PFM]))
    # print("C", "".join([str(min(p["C"], 9)) for p in PFM]))
    # print("G", "".join([str(min(p["G"], 9)) for p in PFM]))
    # print("T", "".join([str(min(p["T"], 9)) for p in PFM]))
    # sys.exit()

    return PFM


# from collections import defaultdict
# def get_block_freqs_and_majority_window(block_matrix, partition, k):
#     nr_columns = len( list(block_matrix.values())[0][0]) - k + 1
#     block_frequency_matrix = defaultdict(int) 
#     BLOCK_FREQ = [0 for j in range(nr_columns)]
#     for s in block_matrix:
#         read_aln_vector, _  = partition[s]
#         block_v, cnt = block_matrix[s]
#         for j in range(nr_columns):
#             BLOCK_FREQ[j] += cnt*block_v[j]
#             if block_v[j] == 1:
#                 segment = "".join([n for n in read_aln_vector[j:j+k]])
#                 block_frequency_matrix[j][segment] += cnt

#     block_majority_vector = [max(d.items(), key=lambda x: x[1])[0] for d in block_frequency_matrix]
#     # for j in range(nr_columns):
#     #     print(j, BLOCK_FREQ[j], block_majority_vector[j], block_frequency_matrix[j])

#     # print("".join([segment for segment in block_majority_vector if segment != "-"]))
#     # print([(segment, BLOCK_FREQ[i]) for i, segment in enumerate(block_majority_vector)])
#     # sys.exit()
#     return BLOCK_FREQ, block_majority_vector, block_frequency_matrix


# def correct_from_msa_window(ref_seq, partition, BFM, seq_to_acc, k, pos_cutoffs = [], block_matrix = {}, block_majority_vector = []):
#     nr_columns = len(BFM)
#     S_prime_partition = {}
#     ref_alignment = partition[ref_seq][0]
#     ref_alignment_kmers = ["".join([n for n in ref_alignment[i:i+k]]) for i in range(len(ref_alignment) -k)]
#     subs_tot = 0
#     ins_tot = 0
#     del_tot = 0
#     # print(sum([ block_matrix[s][0][628] for s in block_matrix ]))
#     # sys.exit()
#     for s in partition:
#         if s == ref_seq:
#             continue
#         subs_pos_corrected = 0
#         ins_pos_corrected = 0
#         del_pos_corrected = 0
#         seq_aln, cnt = partition[s]
#         seq_alignment_kmers = ["".join( n for n in seq_aln[i:i+k]) for i in range(len(seq_aln) -k)]
#         s_new_kmers = [ "".join( n for n in seq_aln[i:i+k]) for i in range(len(seq_aln) -k)]
#         block_v_kmers = [ block_v[i:i+k] for i in range(len(block_v) -k)]
#         # s_new = [n for n in seq_aln]
#         # block_v, _ = block_matrix[s]
#         assert len(block_v_kmers) == len(seq_alignment_kmers)

#         if cnt == 1:
#             for j, kmer in enumerate(seq_alignment_kmers):
#                 if block_v_kmers[j] == 0:
#                     # print(j,"not present")
#                     continue

#                 ref_kmer = ref_alignment_kmers[j]
#                 # print(j, pos_cutoffs[j]["mm"],pos_cutoffs[j]["i"], pos_cutoffs[j]["d"] )
#                 # print(BFM[j])
#                 # print()

#                 # kmer_count = [ c if in block:] implement code for finding elegible reads to count kmers over. Elegible if has at least one 1 in block vector
#                 if kmer_count < 2:
#                     tmp_dict = {segment:cnt for segment, cnt in BFM[j].items()}
#                     majority_correct_to = block_majority_vector[j] 
#                     s_new_kmers[j] = ref_kmer
        
#         for kmer1, kmer2 in zip(s_new_kmers[:-1], s_new_kmers[1:]):
#             if kmer1[1:] == kmer2[:-1]:
#                 pass
#             else:
#                 print(kmer1, kmer2, "not matching")

#         # accessions_of_s = seq_to_acc[s] 
#         # for acc in accessions_of_s:
#         #     S_prime_partition[acc] = "".join([kmer for kmer in s_new if kmer != "-"])
#     sys.exit()
#     print("Corrected {0} subs pos, {1} ins pos, and {2} del pos corrected in partition".format(subs_tot, ins_tot, del_tot))

#     return S_prime_partition


# from scipy.stats import poisson
def correct_from_msa(partition, BFM, pos_cutoffs = [], block_matrix = {}, block_majority_vector = [], homopol_correction_vector = []):
    nr_columns = len(BFM)
    S_prime_partition = {}
    ref_alignment = partition['ref'][0]
    subs_tot = 0
    ins_tot = 0
    del_tot = 0
    # print(sum([ block_matrix[s][0][628] for s in block_matrix ]))
    # sys.exit()
    print(homopol_correction_vector)
    for read_acc in partition:
        if read_acc == "ref":
            continue
        subs_pos_corrected = 0
        ins_pos_corrected = 0
        del_pos_corrected = 0
        seq_aln, cnt = partition[read_acc]
        s_new = [n for n in seq_aln]
        block_v, _ = block_matrix[read_acc]
        assert len(block_v) == len(seq_aln)

        if cnt == 1:
            for j, n in enumerate(seq_aln):
                if block_v[j] == 0:
                    # print(j,"not present")
                    continue

                ref_nucl = ref_alignment[j]
                # print(j, pos_cutoffs[j]["mm"],pos_cutoffs[j]["i"], pos_cutoffs[j]["d"] )
                # print(BFM[j])
                # print()

                if n != ref_nucl:
                    if n == "-":  # deletion
                        # if BFM[j][n] <= pos_cutoffs[j]["d"]:
                        #     s_new[j] = ref_nucl
                        #     del_pos_corrected += 1
                        if homopol_correction_vector[j] > 1:
                            lambda_del = pos_cutoffs[j]["h"][0]
                            lambda_corr = pos_cutoffs[j]["h"][1]
                            p_del = poisson.pmf(BFM[j][n], lambda_del)
                            p_corr = poisson.pmf(BFM[j][n], lambda_corr)
                            # if homopol_correction_vector[j] == 4:
                            #     print(BFM[j][n], pos_cutoffs[j]["h"], p_del, p_corr)
                            if p_corr < p_del:
                                s_new[j] = ref_nucl
                                del_pos_corrected += 1

                            # if BFM[j][n] <= pos_cutoffs[j]["d"] and p_corr > p_del:
                            #     print(j, homopol_correction_vector[j],)
                            #     print(n, BFM[j][n], pos_cutoffs[j]["h"], p_del, p_corr)
                            #     print(pos_cutoffs[j]["d"])
                        else:
                            if BFM[j][n] <= pos_cutoffs[j]["d"]:
                                s_new[j] = ref_nucl
                                del_pos_corrected += 1
                    else:
                        if ref_nucl == "-": # insertion
                            if BFM[j][n] <= pos_cutoffs[j]["i"]:
                                s_new[j] = ref_nucl
                                ins_pos_corrected += 1
                        else: # substitution                                
                            if BFM[j][n] <= pos_cutoffs[j]["mm"]:
                                s_new[j] = ref_nucl
                                subs_pos_corrected +=1

                            # tmp_dict = {n:cnt for n, cnt in BFM[j].items()}
                            # base_correct_to = max(tmp_dict, key = lambda x: tmp_dict[x] )
                            # s_new[j] = base_correct_to
                            # if base_correct_to != "-":
                            #     subs_pos_corrected +=1
                            # else:
                            #     ins_pos_corrected += 1

                else:
                    tmp_dict = {n:cnt for n, cnt in BFM[j].items()}
                    majority_correct_to = block_majority_vector[j] 
                    if majority_correct_to != n:
                        if majority_correct_to == "-":
                            if BFM[j][n] <= pos_cutoffs[j]["i"]:
                                s_new[j] = "-"
                                ins_pos_corrected += 1
                        elif n != "-": # substitution: both majority and n and ref are not "-"
                            if BFM[j][n] <= pos_cutoffs[j]["mm"]:
                                s_new[j] = majority_correct_to
                                subs_pos_corrected += 1                       
                        else: # deletion: majority has nucleotide but this is not present in reference, this should be rare. Unsure if this should be classified even as a bug
                            print(j, pos_cutoffs[j],BFM[j][n])
                            print("Unexpected! Majority: {0}, base in read and reference:{1}. Total counts over position: {2}".format(majority_correct_to, n, tmp_dict) )
                            s_new[j] = majority_correct_to

        subs_tot += subs_pos_corrected
        ins_tot += ins_pos_corrected
        del_tot += del_pos_corrected
        print("Corrected {0} subs pos, {1} ins pos, and {2} del pos corrected in seq of length {3}".format(subs_pos_corrected, ins_pos_corrected, del_pos_corrected, len(read_acc)))

        # accessions_of_s = seq_to_acc[read_acc] 
        # for acc in accessions_of_s:
        S_prime_partition[read_acc] = "".join([n for n in s_new if n != "-"])

    print("Corrected {0} subs pos, {1} ins pos, and {2} del pos corrected in partition".format(subs_tot, ins_tot, del_tot))
    # sys.exit()
    return S_prime_partition


def msa(reference_seq, partition, seq_to_acc):

    """
         partition[seq] = (edit_dist, aln_rep, aln_s, depth_of_string)
    """
    N_t = sum([container_tuple[1] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    if len(partition) > 1:
        # all strings has not converged
        alignment_matrix = { seq : [n for n in partition[seq][0]] for seq in partition }
        lengths = [len(v) for v in alignment_matrix.values()]
        assert len(set(lengths)) == 1

        PFM = PFM_from_msa(partition)
        S_prime_partition = correct_from_msa(partition, PFM, seq_to_acc)

    else:
        print("Partition converged: Partition size(unique strings):{0}, partition support: {1}.".format(len(partition), N_t))

    
    return S_prime_partition


