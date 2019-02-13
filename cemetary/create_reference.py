
def create_position_frequency_matrix_modified_count_vector(repr_seq, alignment_matrix, partition):
    nr_columns = len( alignment_matrix[ list(alignment_matrix)[0] ]) # just pick a key
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "U" : 0, "-": 0} for j in range(nr_columns)]
    for s in alignment_matrix:
        s_aln = alignment_matrix[s]
        if s != repr_seq:
            indegree = partition[s][3]
            assert indegree == 1
            for j in range(nr_columns): # for each column
                nucl = s_aln[j]
                PFM[j][nucl] += indegree
        else:
            ref_indegrees = partition[s][3]
            ref_pos = 0
            for j in range(nr_columns): # for each column
                # if j > 24000:
                #     print(ref_pos, len(ref_indegrees), len(repr_seq), j, nr_columns, s_aln[j:])
                nucl = s_aln[j]
                PFM[j][nucl] += ref_indegrees[ref_pos]
                if nucl != "-" and ref_pos < len(ref_indegrees) -1:
                    ref_pos += 1



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


def update_reference(repr_seq, partition, seq_to_acc, qual_dict, ref_count_vector):
    """
         partition[seq] = (edit_dist, aln_rep, aln_s, depth_of_string)
         Update to majority position at each position except if the combination: majority is "-" and current_repr has a non "-" character in that position
    """
    # for s in partition:
    #     print("len:",len([n for n in partition[s][1] if n != "-" ]))
    alignment_matrix = create_multialignment_matrix(repr_seq, partition) 
    
    repr_vector = []
    for n in alignment_matrix[repr_seq]:
        repr_vector.append(n)
    # print( "Old ref count vector:",partition[repr_seq][3])
    print(sorted(ref_count_vector)[-50:], sum(ref_count_vector))
    # alignment_matrix_of_qualities, PFM_qualities, PFM_max_qualities, global_correction_threshold = annotate_with_quality_values(alignment_matrix, seq_to_acc, qual_dict)
    # assert len(alignment_matrix_of_qualities) == len(alignment_matrix)

    PFM = create_position_frequency_matrix_modified_count_vector(repr_seq, alignment_matrix, partition)

    majority_vector = []
    new_ref_count_vector = []




    # for j in range(len(PFM)): 
    #     max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
    #     majority_count = PFM[j][max_v_j]
    #     max_v_j_set = set([v for v in PFM[j] if PFM[j][v] == majority_count ])
    #     all_major = "".join(max_v_j_set)
    #     if len(all_major) > 1:
    #         print("OMG!!", PFM[j])
    #         all_major = max_v_j_set.pop()
    #     majority_vector.append( all_major )

    #     if all_major != "-": 
    #         new_ref_count_vector.append(majority_count)
    # new_ref = "".join([n for n in majority_vector if n !="-"])





    for j in range(len(PFM)): 
        if repr_vector[j] == "-":
            # all_poss_insertions = [(v,PFM[j][v]) for v in PFM[j] if (PFM[j][v] > 8 and v != "-")]

            # if all_poss_insertions:
            #     print(j, all_poss_insertions, PFM[j]["-"], len(alignment_matrix) , PFM[j])
            #     max_v_j, majority_count = max(all_poss_insertions, key = lambda x: x[1] )
            #     max_v_j_set = set([v for v,c in all_poss_insertions if c == majority_count ])
            #     all_major = "".join(max_v_j_set)
            #     # print(all_major)
            # else:
            max_v_j = max(PFM[j], key = lambda x: PFM[j][x] )
            majority_count = PFM[j][max_v_j]
            max_v_j_set = set([v for v in PFM[j] if PFM[j][v] == majority_count ])
            all_major = "".join(max_v_j_set)
            if len(all_major) > 1:
                print("OMG!!", PFM[j])
                all_major = "".join([n for n in max_v_j_set if n != "-"])
        else:
            tmp_dict = PFM[j]
            del tmp_dict["-"]
            max_v_j = max(tmp_dict, key = lambda x: tmp_dict[x] )
            majority_count = tmp_dict[max_v_j]
            max_v_j_set = set([v for v in tmp_dict if tmp_dict[v] == majority_count ])
            all_major = "".join(max_v_j_set)
            if len(all_major) > 1:
                print("OK", j, PFM[j])
                # max_v_j_set.remove(repr_vector[j])
                # all_major = repr_vector[j] # take the representative nucleotide
                all_major = max_v_j_set.pop()
        majority_vector.append( all_major )
        if all_major != "-": 
            new_ref_count_vector.append(majority_count)
        # else:
        #     majority_vector.append( all_major )
        # if len(all_major) >1:
        #     print(PFM[j])

    # print(majority_vector)
    new_ref = "".join([n for n in majority_vector if n !="-"])

    assert len(majority_vector) == len(PFM)
    print(len(new_ref), len(new_ref_count_vector))
    assert len(new_ref) == len(new_ref_count_vector)
    return new_ref, new_ref_count_vector






def main(args):

    print("started sorting seqs")
    start = time()
    reads = { i : (acc, seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.fastq, 'r')))}
    # read_array = [ (i, 0, acc, seq, qual, float(acc.split("_")[-1])) for i, (acc, (seq, qual)) in enumerate(readfq(open(sorted_reads_fastq_file, 'r')))]
    print("Correcting {0} reads.".format(len(reads)))
    start = time()

    ref_file = os.path.join(args.outfolder, "reference.fa")
    reference_seq = create_augmented_reference.run_spoa(args.fastq, ref_file, "spoa")

    origin_acc, origin_seq, origin_qual = reads[0]

    qual_dict = {origin_acc : [ phred_char_to_val[q] for q in origin_qual] }    
    block_dict = {origin_acc : [ 1 for q in origin_qual] }    
    seq_to_acc = get_seq_to_acc(reads)
    partition = {}

    ### Augment reference starting with origin as basis
    reference_seq = origin_seq
    ref_count_vector = [1]*len(reference_seq)
    reference_qual = origin_qual
    partition[reference_seq] = (0, reference_seq, reference_seq, [1]*len(reference_seq))

    for i in range(1, len(reads)):
        print(i)
        acc, seq, qual = reads[i]
        qual_dict[acc] =  [ phred_char_to_val[q] for q in qual]
        error_rate_sum, read_alignment, reference_alignment, alignment_ratio, block_coverage = block_align(reference_seq, seq, reference_qual, qual, args)

        if segment_should_be_added(reference_alignment, read_alignment, qual):
            # correct current ref
            print("Using {0} to update consensus".format(len(partition)))
            reference_seq, ref_count_vector = correct_seqs.update_reference(reference_seq, partition, seq_to_acc, qual_dict, ref_count_vector)
            print(reference_seq)

            reference_qual = "".join(["+" for i in range(len(reference_seq))]) # TODO: Maybe update to actually use quals here

            # realign read to new ref and add segment to updated reference
            error_rate_sum, read_alignment, reference_alignment, alignment_ratio, block_coverage = block_align(reference_seq, seq, reference_qual, qual, args)
            reference_seq, reference_qual, ref_count_vector = add_segment(reference_alignment, read_alignment, qual, reference_qual, ref_count_vector)
            
            # Reinitiate data structure partition
            partition = {} 
            partition[reference_seq] = (0, reference_seq, reference_seq, ref_count_vector)

            # realign read that cause to update reference to the new reference
            error_rate_sum, read_alignment, reference_alignment, alignment_ratio, block_coverage = block_align(reference_seq, seq, reference_qual, qual, args)

        
        mismatches = len([ 1 for n1, n2 in zip(reference_alignment,read_alignment) if n1 != n2 and n1 != "-" and n2 != "-" ])
        matches = len([ 1 for n1, n2 in zip(reference_alignment,read_alignment) if n1 == n2 and n1 != "-"])
        indels = len(reference_alignment) - mismatches - matches
        partition[seq] = (matches+indels, reference_alignment, read_alignment, len(seq_to_acc[seq]))

    # last iteration

    print("Using {0} to update consensus".format(len(partition)))
    reference_seq, ref_count_vector = correct_seqs.update_reference(reference_seq, partition, seq_to_acc, qual_dict, ref_count_vector)
    reference_qual = "".join(["+" for i in range(len(reference_seq))]) # TODO: Maybe update to actually use quals here
    partition = {} 

    print(sorted(ref_count_vector)[-50:], sum(ref_count_vector))
    print(ref_count_vector)
    print(reference_seq)
    r2 = []
    for i, c in enumerate(ref_count_vector):
        if c > 10:
            r2.append(reference_seq[i])
    print("r2", "".join(r2))

    sys.exit()




def add_segment(reference_alignment, read_alignment, qual, reference_qual, ref_count_vector):
    curr_ref_pos = 0
    curr_aln_pos = 0
    new_ref = []
    new_quals = []
    new_ref_count_vector = []
    for ch, g in itertools.groupby(reference_alignment):
        l = list(g)
        if l[0] == "-" and len(l) > 10 : # add quality values of read check here
            # print(l)
            read_pos = len([read_char for read_char in read_alignment[:curr_aln_pos] if read_char != "-"])
            q_part = qual[read_pos:read_pos+ len(l)]
            # print([round(phred_char_to_p[q_v], 4) for q_v in q_part])
            print("AVG:", sum([phred_char_to_p[q_v] for q_v in q_part])/float(len(q_part)))
            print(read_alignment[curr_aln_pos: curr_aln_pos + len(l)] )
            if sum([phred_char_to_p[q_v] for q_v in q_part])/float(len(q_part)) < 0.2:
                print("Adding")
                new_ref.append(read_alignment[curr_aln_pos: curr_aln_pos + len(l)] )
                new_quals.append(qual[read_pos:read_pos+ len(l)])
                new_ref_count_vector.append([1]*len(l))
        if l[0] != "-":
            new_ref.append(l)
            new_quals.append(reference_qual[curr_ref_pos : curr_ref_pos + len(l)])
            new_ref_count_vector.append(ref_count_vector[curr_ref_pos : curr_ref_pos + len(l)])
            curr_ref_pos += len(l)

        curr_aln_pos += len(l)

    reference_seq = "".join([c for sublist in new_ref for c in sublist])
    reference_qual = "".join([c for sublist in new_quals for c in sublist])
    reference_count_vector = [c for sublist in new_ref_count_vector for c in sublist]
    print(len(reference_seq), len(reference_qual), len(reference_count_vector))  
    assert len(reference_seq) == len(reference_qual) == len(reference_count_vector)

    return reference_seq, reference_qual, reference_count_vector


# def segment_should_be_added(reference_alignment, read_alignment, qual):
#     curr_ref_pos = 0
#     curr_aln_pos = 0
#     for ch, g in itertools.groupby(reference_alignment):
#         l = list(g)
#         if l[0] == "-" and len(l) > 10 : # add quality values of read check here
#             read_pos = len([read_char for read_char in read_alignment[:curr_aln_pos] if read_char != "-"])
#             q_part = qual[read_pos:read_pos+ len(l)]
#             if sum([phred_char_to_p[q_v] for q_v in q_part])/float(len(q_part)) < 0.2:
#                 return True
#         if l[0] != "-":
#             curr_ref_pos += len(l)
#         curr_aln_pos += len(l)
#     return False

