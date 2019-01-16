import edlib
import re
import math

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
        PFM_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
        PFM_max_qualities.append({"A": 0, "C": 0, "G": 0, "T": 0, "-": 0})
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
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "-": 0} for j in range(nr_columns)]
    for s in alignment_matrix:
        s_aln = alignment_matrix[s]
        indegree = partition[s][3]
        for j in range(nr_columns): # for each column
            nucl = s_aln[j]
            PFM[j][nucl] += indegree

    # As = [PFM[p]["A"] for p in range(len(PFM)) ]
    # Cs = [PFM[p]["C"] for p in range(len(PFM)) ]
    # Gs = [PFM[p]["G"] for p in range(len(PFM)) ]
    # Ts = [PFM[p]["T"] for p in range(len(PFM)) ]
    # dels = [PFM[p]["-"] for p in range(len(PFM))]
    # for i in range(len(As)):
    #     if i % 2 == 1:
    #         print(As[i], Cs[i], Gs[i], Ts[i], dels[i])
    # print(As)
    # print(Cs)
    # print(Gs)
    # print(Ts)
    # print(dels)

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

    for nucl in ["A", "G", "C", "T", "-"]:
        position_solutions[nucl] = {"A": ["A"], "G": ["G"], "C": ["C"], "T": ["T"], "-": ["-"]}

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


def create_multialignment_matrix(repr_seq, partition):
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
    query_to_target_positioned_dict = {}
    for q_acc in partition:
        (edit_distance, m_alignment, s_alignment, degree_of_s) = partition[q_acc]
        # s_positioned_new, target_vector_start_position, target_vector_end_position = position_query_to_alignment_NEW(s_alignment, m_alignment, 0)
        # print(s_alignment)
        # print(m_alignment)
        # print(repr_seq)
        s_positioned, target_vector_start_position, target_vector_end_position = position_query_to_alignment(s_alignment, m_alignment, 0)
        # assert s_positioned_new == s_positioned
        # print(s_positioned)
        # print()
        assert target_vector_start_position == 0
        # print(target_vector_end_position + 1, len(repr_seq), 2*len(repr_seq)+1)
        assert target_vector_end_position + 1 == 2*len(repr_seq) + 1 # vector positions are 0-indexed
        query_to_target_positioned_dict[q_acc] = (s_positioned, target_vector_start_position, target_vector_end_position)

    alignment_matrix = create_multialignment_format_NEW(query_to_target_positioned_dict, 0, 2*len(repr_seq))
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

def correct_to_consensus(repr_seq, partition, seq_to_acc, qual_dict):
    """
         partition[seq] = (edit_dist, aln_rep, aln_s, depth_of_string)
    """
    S_prime_partition = {}
    S_prime_quality_vector = {}
    N_t = sum([container_tuple[3] for s, container_tuple in partition.items()]) # total number of sequences in partition
    
    if len(partition) > 1:
        # all strings has not converged
        alignment_matrix = create_multialignment_matrix(repr_seq, partition) 
        PFM = create_position_frequency_matrix(alignment_matrix, partition)

        for s_before in partition:
            s_after = "".join([n for n in alignment_matrix[s_before] if n != "-"])
            assert s_before == s_after
        # print(len(partition), N_t)
        alignment_matrix_of_qualities, PFM_qualities, PFM_max_qualities, global_correction_threshold = annotate_with_quality_values(alignment_matrix, seq_to_acc, qual_dict)

        assert len(alignment_matrix_of_qualities) == len(alignment_matrix)
        if global_correction_threshold < 0:
            return S_prime_partition, S_prime_quality_vector

        majority_vector = []
        for j in range(len(PFM_qualities)):
            max_v_j = max(PFM_qualities[j], key = lambda x: PFM_qualities[j][x] )

            majority_count =  PFM_qualities[j][max_v_j]
            max_v_j_set = set([v for v in PFM_qualities[j] if PFM_qualities[j][v] == majority_count ])
            all_major = "".join(max_v_j_set)
            majority_vector.append( all_major )
        assert len(majority_vector) == len(PFM_qualities)
        ############################
        ############################


        for s in sorted(partition):
            if partition[s][3] > 1: # at least 2 identical sequences --> its a nearest_neighbor of the partition, has converged, and should not be corrected
                print("not correcting converged sequence!")
                continue

            s_alignment_in_matrix = alignment_matrix[s]
            # s_min = 0
            for i, p in enumerate(s_alignment_in_matrix):
                if p != "-":
                    s_min = i
                    break
            for i, p in enumerate(s_alignment_in_matrix[::-1]):
                if p != "-":
                    s_max = len(s_alignment_in_matrix) - i
                    break

            # print("S", s_min, s_max)
            s_quals = alignment_matrix_of_qualities[s]
            # ALL POSITIONS with sum of probabilities lower than the largest probability
            minority_positions = [ (j,majority_vector[j], s_alignment_in_matrix[j]) for j in range(len(majority_vector)) if s_alignment_in_matrix[j] not in majority_vector[j] and s_min <= j <= s_max ]
            minority_positions_correctable = [ j  for j in range(len(majority_vector)) if (len(majority_vector[j]) == 1 and majority_vector[j] != s_alignment_in_matrix[j] ) ]
            minority_positions_correctable = [ (j, PFM_qualities[j][ s_alignment_in_matrix[j] ])  for j in minority_positions_correctable ]
            nr_pos_to_correct = int(math.ceil( len(minority_positions_correctable) * 0.5)) # (step/ float(step +1)) ))
            # print("positions to correct:", nr_pos_to_correct) 

            if nr_pos_to_correct  == 0:
                print("Edit distance to nearest_neighbor:", partition[s][0], "is nearest_neighbor:", s ==repr_seq, "Minority positions:", minority_positions)
                continue
            if len(minority_positions_correctable) == 0:
                print("no unambiguous majority positions")
                continue

            minority_positions_correctable.sort(key=lambda x: x[1])
            # print(len(minority_positions_correctable) ,minority_positions_correctable)
            _, quality_threshold_to_correct = minority_positions_correctable[ nr_pos_to_correct - 1 ]
            minority_positions_to_correct = [ (j, qual_j) for j, qual_j in minority_positions_correctable if qual_j <= quality_threshold_to_correct ]
            print(quality_threshold_to_correct, len(minority_positions_to_correct))
            # minority_positions_to_correct = [ (j, qual_j) for j, qual_j in minority_positions_correctable if qual_j <= global_correction_threshold ]
            # print("actual:", len(minority_positions_to_correct))
            # minority_positions_to_correct = sorted(minority_positions_correctable, key=lambda x: x[1])[:nr_pos_to_correct]  # sorted list with the smallest probabilities first

            # print(minority_positions_to_correct)
            s_new = alignment_matrix[s]
            s_qual_new = alignment_matrix_of_qualities[s]
            for j, qual_j in minority_positions_to_correct:
                highest_prob_character_at_j = majority_vector[j]
                assert len(majority_vector[j]) == 1
                s_new[j] = highest_prob_character_at_j
                s_qual_new[j] = PFM_max_qualities[j][highest_prob_character_at_j] 
            
            s_modified = "".join([nucl for nucl in s_new if nucl != "-" ])
            s_qual_modified = [s_qual_new[j] for j in range(len(s_new)) if s_new[j] != "-" ]

            # only unique strings can change in this step
            accessions_of_s = seq_to_acc[s] 
            for acc in accessions_of_s:
                S_prime_partition[acc] = s_modified
                S_prime_quality_vector[acc] = s_qual_modified
    else:
        print("Partition converged: Partition size(unique strings):{0}, partition support: {1}.".format(len(partition), N_t))

    
    return S_prime_partition, S_prime_quality_vector

