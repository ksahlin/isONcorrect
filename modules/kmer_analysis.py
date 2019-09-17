
from collections import defaultdict

import edlib

from modules import help_functions


def get_kmers_over_pos(read_aln, ref_aln, k_size):
    k_list = []
    ref_pos = 0

    ct = 0
    for i, n in enumerate(ref_aln[::-1]):
        if n != "-":
            ct += 1
        if ct >= k_size:
            break
    last_kmer_pos = len(ref_aln) - i - 1
    # print(last_kmer_pos,ref_aln[last_kmer_pos:] )

    for i in range(last_kmer_pos):
        if ref_aln[i] != "-":
            it = i
            read_kmer = []
            while len(read_kmer) < k_size:
                if read_aln[it] != "-":
                    read_kmer.append(read_aln[it])
                it += 1
                if it >= len(read_aln):
                    break

            if len(read_kmer) == k_size:
                read_kmer = "".join([b for b in read_kmer])
                k_list.append(read_kmer) 
            else:
                break
            ref_pos += 1
        else:
            pass
    # print(len(k_list))
    return k_list

def kmers(seq,k):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

def get_kmer_read_pos_on_ref(reads, reference_seq, k_size):
    DBG = defaultdict(int)

    # r_pos = { i: {} for i in range(len(reference_seq) - k_size+1)}
    # print(len(r_pos), "lol")
    for i, (acc, (seq, qual, pos1, pos2)) in enumerate(reads.items()):
        for km in kmers(seq,k_size):
            DBG[km] +=1
        # res = edlib.align(seq, reference_seq, task="path", mode="NW")
        # cigar_string = res["cigar"]
        # read_aln, ref_aln = help_functions.cigar_to_seq(cigar_string, seq, reference_seq)

        # # read_aln, ref_aln = parasail_aln(seq, reference_seq)
        # k_list = get_kmers_over_pos(read_aln, ref_aln, k_size)
        # # print(len(k_list), len(reference_seq) - k_size+1)
        # # assert len(k_list) == len(reference_seq) - k_size+1
        # for i in range(len(k_list)):
        #     kmer = k_list[i] 
        #     if kmer in r_pos[i]:
        #         r_pos[i][kmer] +=1
        #     else:
        #         r_pos[i][kmer] = 1
    return DBG #r_pos 



def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def traverse(d, start_anchor, stop_anchor, max_depth = 40):
    curr_kmer = start_anchor
    max_path_kmers = [curr_kmer]
    for i in range(max_depth):
        kmer_candidates = [(d[x],x) for x in fw(curr_kmer) if x in d]
        if kmer_candidates:
            k_count, curr_kmer = max(kmer_candidates) 
            if curr_kmer == stop_anchor:
                break
            else:    
                max_path_kmers.append(curr_kmer)
        else:
            break
    max_path_kmers.append(stop_anchor)

    is_valid_path = all([k[1:] == k_fw[:-1] for k,k_fw in zip(max_path_kmers[:-1], max_path_kmers[1:]) ])
    if is_valid_path:
        print("SUCCESS", [d[k] for k in max_path_kmers])
        corr_spoa = max_path_kmers[0] + "".join([km[-1] for km in max_path_kmers[1:] ])
        print(corr_spoa)
        return corr_spoa
    else:
        print("NOT VALID PATH",max_path_kmers)
        return False


