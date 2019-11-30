#! /usr/bin/env python

from __future__ import print_function
import os,sys
import argparse

import errno
from time import time
import itertools
import tempfile
import shutil

import math
import re
from collections import deque
from collections import defaultdict

import edlib

from modules import create_augmented_reference, help_functions, correct_seqs #,align

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers

def get_kmer_maximizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = max(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = max(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer > curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers




def get_minimizers_and_positions(reads, w, k, hash_fcn):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        if hash_fcn == "lex":
            minimizers = get_kmer_minimizers(seq, k, w)
        elif hash_fcn == "rev_lex":
            minimizers = get_kmer_maximizers(seq, k, w)

        M[r_id] = minimizers

    return M


def get_minimizer_combinations_database(reads, M, k, x_low, x_high):
    M2 = defaultdict(lambda: defaultdict(list))
    tmp_cnt = 0
    forbidden = 'A'*k
    for r_id in M:
        minimizers = M[r_id]
        for (m1,p1), m1_curr_spans in  minimizers_comb_iterator(minimizers, k, x_low, x_high):
            for (m2, p2) in m1_curr_spans:
                if m2 == m1 == forbidden:
                    continue

                tmp_cnt +=1
                M2[m1][m2].append((r_id, p1, p2))

    print(tmp_cnt, "MINIMIZER COMBINATIONS GENERATED")

    avg_bundance = 0
    singleton_minimzer = 0
    cnt = 1
    abundants=[]
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            if len(M2[m1][m2]) > 1:
                avg_bundance += len(M2[m1][m2])
                cnt +=1
            else:
                del M2[m1][m2]
                singleton_minimzer += 1

            # if len(M2[m1][m2]) > len(reads):
            #     abundants.append((m1,m2, len(M2[m1][m2])))
            #     if m2 == forbidden: # poly A tail
            #         del M2[m1][m2]
    for m1,m2,ab in sorted(abundants, key=lambda x: x[2], reverse=True):
        print("Too abundant:", m1, m2, ab, len(reads))

    print("Average abundance for non-unique minimizer-combs:", avg_bundance/float(cnt))
    print("Number of singleton minimizer combinations filtered out:", singleton_minimzer)
    return M2



def minimizers_comb_iterator(minimizers, k, x_low, x_high):
    # print("read")
    for i, (m1, p1) in enumerate(minimizers[:-1]):
        m1_curr_spans = []
        for j, (m2, p2) in enumerate(minimizers[i+1:]):
            if x_low < p2 - p1 and p2 - p1 <= x_high:
                m1_curr_spans.append( (m2, p2) )
                # yield (m1,p1), (m2, p2) 
            elif p2 - p1 > x_high:
                break
        yield (m1, p1), m1_curr_spans[::-1]


def edlib_alignment(x, y, k):
    # if i == 100 and j % 1000 == 0:
    #     print("Edlib processed alignments: {0}".format(j+1))

    result = edlib.align(x,y, "NW", 'dist', k) # , task="path")
    ed = result["editDistance"]
    # locations = result["locations"]
    return ed #, locations




def fill_p2(p, all_intervals_sorted_by_finish):

    stop_to_max_j = {stop : j for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish)}
    all_choord_to_max_j = []
    j_max = 0
    for i in range(0, all_intervals_sorted_by_finish[-1][1]):
        if i in stop_to_max_j:
            j_max = stop_to_max_j[i]
        
        all_choord_to_max_j.append(j_max)

    for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish):
        j_max = all_choord_to_max_j[start]
        p.append(j_max)
    return p

def solve_WIS2(all_intervals_sorted_by_finish):
    # print("instance size", len(all_intervals_sorted_by_finish))
    # p = [None]
    # fill_p(p, all_intervals_sorted_by_finish)
    p = [None]    
    fill_p2(p, all_intervals_sorted_by_finish)
    # if p != p2:
    #     print(p)
    #     print(p2)
    # assert p == p2

    v = [None] + [w*(stop-start) for (start, stop, w, _) in all_intervals_sorted_by_finish]
    OPT = [0]
    for j in range(1, len(all_intervals_sorted_by_finish) +1):
        OPT.append( max(v[j] + OPT[ p[j] ], OPT[j-1] ) )

    # assert len(p) == len(all_intervals_sorted_by_finish) + 1 == len(v) == len(OPT)

    # Find solution
    opt_indicies = []
    j = len(all_intervals_sorted_by_finish)
    while j >= 0:
        if j == 0:
            break
        if v[j] + OPT[p[j]] > OPT[j-1]:
            opt_indicies.append(j - 1) # we have shifted all indices forward by one so we neew to reduce to j -1 because of indexing in python works
            j = p[j]
        else:
            j -= 1
    return opt_indicies





def find_most_supported_span(r_id, m1, p1, m1_curr_spans, minimizer_combinations_database, reads, all_intervals, k_size, tmp_cnt, read_complexity_cnt, already_computed):

    # curr_m_pos = read_min_comb[0][0][1]
    # curr_m_pos2 = read_min_comb[0][1][1]
    acc, seq, qual = reads[r_id]
    curr_best_seqs = {}
    curr_best_seqs_to_id = {}
    # cnt = 0
    for (m2,p2) in m1_curr_spans:
        # print(p1,p2)
        relevant_reads = minimizer_combinations_database[m1][m2]
        seqs = {} #defaultdict(list)
        added_strings = {} 
        # not_added_strings = set() 
        if len(relevant_reads) >= max(3,len(curr_best_seqs)): #max(3, previously_calculated_regions_read[p2]): # 
            # cnt += 1
            ref_seq = seq[p1  : p2 + k_size]
            ref_qual = qual[p1 : p2 + k_size]            

            seqs["curr_read"] = (ref_seq, ref_qual, p1, p2)
            added_strings[ref_seq] = 0
            reads_visited = {}
            for relevant_read_id, pos1, pos2 in relevant_reads:
                if r_id  == relevant_read_id:
                    continue
                
                read_seq = reads[relevant_read_id][1][pos1: pos2 + k_size]
                read_qual = reads[relevant_read_id][2][pos1: pos2 + k_size]

                if read_seq == ref_seq:
                    seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
                    reads_visited[relevant_read_id] = 0
                    already_computed[relevant_read_id] = (p1,p2,pos1,pos2, 0)
                    continue
                elif relevant_read_id in reads_visited:
                    # print("Prev:", reads_visited[relevant_read_id])
                    # print("Act:", edlib_alignment(ref_seq, read_seq, p_error_sum_thresh*len(ref_seq)) )
                    pass
        # Implement if we see this to recompute all the aligments exact ed here instead!! Thats the only way to guarantee exactly the same
        # or maybe use this traceback to get exact: https://github.com/Martinsos/edlib/pull/132#issuecomment-522258271
                elif read_seq in added_strings:  #== ref_seq:
                    # seqs[relevant_read_id] = []
                    seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
                    reads_visited[relevant_read_id] = added_strings[read_seq]
                    already_computed[relevant_read_id] = (p1,p2,pos1,pos2, added_strings[read_seq])
                    # print("Saved!!")
                    continue
                # elif read_seq in not_added_strings:
                #     # print("save")
                #     continue
                elif relevant_read_id in already_computed:
                    curr_ref_start, curr_ref_end, curr_read_start, curr_read_end, curr_ed = already_computed[relevant_read_id]
                    if (curr_read_start <= pos1 and pos2 <= curr_read_end) and (curr_ref_start <= p1 and p2 <=  curr_ref_end):
                        read_beg_diff = pos1 - curr_read_start
                        read_end_diff = pos2 - curr_read_end
                        ref_beg_diff = p1 - curr_ref_start
                        ref_end_diff = p2 - curr_ref_end

                        ed_est = curr_ed + math.fabs(ref_end_diff - read_end_diff) + math.fabs(read_beg_diff - ref_beg_diff) 
                        if 0 <= ed_est <= 0.1*len(ref_seq): # < curr_p_error_sum_thresh*len(ref_seq):
                            # print("saved:",pos1, pos2, p1,p2, curr_ref_start, curr_ref_end, curr_read_start, curr_read_end, curr_ed, curr_p_error_sum_thresh)
                            # print("saved", (curr_ref_end - curr_ref_start) - (p2-p1), (curr_read_end - curr_read_start) - (pos2 - pos1))

                            # read_seq = reads[relevant_read_id][1][pos1: pos2 + k_size]
                            # read_qual = reads[relevant_read_id][2][pos1: pos2 + k_size]
                            seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
                            added_strings[read_seq] = ed_est
                            reads_visited[relevant_read_id] = ed_est

                            # est = curr_ed + math.fabs(ref_end_diff - read_end_diff) + math.fabs(read_beg_diff - ref_beg_diff)
                            # act = edlib_alignment(ref_seq, read_seq, p_error_sum_thresh*len(ref_seq))
                            # if est != act:
                            #     print("estimated", est)
                            #     print("Actual", act)

                            continue

                    else:
                        pass
                
                editdist = edlib_alignment(ref_seq, read_seq, 0.1*len(ref_seq))

                # if editdist == 0:
                #     print("Here!")
                tmp_cnt += 1
                if editdist >= 0:    # passing second edit distance check
                    if relevant_read_id in reads_visited: # we have already seen the minimizer combination
                        prev_read_seq, read_read_qual, prev_pos1, prev_pos2 = seqs[relevant_read_id]
                        editdist_prev = edlib_alignment(ref_seq, prev_read_seq, len(ref_seq))
                        tmp_cnt += 1
                        read_complexity_cnt += 1
                        # print("Read already visited lool", "prev:",reads_visited[relevant_read_id], "this:", editdist)
                        # print(seqs[relevant_read_id])
                        # print(p_error_sum_thresh)
                        # print(ref_seq,p1,p2)
                        # print(read_seq, pos1, pos2)
                        # print("current:", editdist, "prev", reads_visited[relevant_read_id])
                        if editdist < editdist_prev:
                            # seqs[relevant_read_id] = []
                            seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
                            added_strings[read_seq] = editdist
                            reads_visited[relevant_read_id] = editdist
                            already_computed[relevant_read_id] = (p1,p2,pos1,pos2, editdist)
                            # print("REPLACED OLD MATCH")
                        else:
                            # seqs[relevant_read_id] = []
                            seqs[relevant_read_id] = (prev_read_seq, read_read_qual, prev_pos1, prev_pos2)
                            added_strings[prev_read_seq] = editdist_prev
                            reads_visited[relevant_read_id] = editdist_prev
                            already_computed[relevant_read_id] = (p1,p2,prev_pos1, prev_pos2, editdist_prev)
                    else:
                        seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
                        added_strings[read_seq] = editdist
                        reads_visited[relevant_read_id] = editdist
                        already_computed[relevant_read_id] = (p1,p2,pos1,pos2, editdist)
                # else:
                #     not_added_strings.add(read_seq)
                #     # print(read_seq, "not considered, store?")


            # zz_cnt = 0

            # # make sure all identical strings are added to he same instance
            # for relevant_read_id, pos1, pos2 in relevant_reads:
            #     if r_id  == relevant_read_id:
            #         continue
            #     read_seq = reads[relevant_read_id][1][pos1: pos2 + k_size]
            #     read_qual = reads[relevant_read_id][2][pos1: pos2 + k_size]
            #     if read_seq in added_strings and relevant_read_id not in seqs:
            #         seqs[relevant_read_id] = (read_seq, read_qual, pos1, pos2)
            #         # zz_cnt += 1
            #         # print("Not added but others are!!!", zz_cnt, read_seq, ref_seq, read_seq == ref_seq)


            all_intervals.append( (p1 + k_size, p2,  len(seqs), seqs) )

    return tmp_cnt, read_complexity_cnt


D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}



def batch(dictionary, size):
    batches = []
    sub_dict = {}
    for i, (acc, seq) in enumerate(dictionary.items()):
        if i > 0 and i % size == 0:
            batches.append(sub_dict)
            sub_dict = {}
            sub_dict[acc] = seq
        else:
            sub_dict[acc] = seq

    if i/size != 0:
        sub_dict[acc] = seq
        batches.append(sub_dict)
    
    return batches

def main(args):
    # start = time()
    all_reads = { i : (acc, seq, qual) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(args.fastq, 'r')))}
    eprint("Total cluster of {0} reads.".format(len(all_reads)))
    if len(all_reads) <= args.exact_instance_limit:
        args.exact = True
    if args.set_w_dynamically:
        args.w = args.k + min(7, int( len(all_reads)/500))

    eprint("ARGUMENT SETTINGS:")
    for key, value in args.__dict__.items():
        eprint("{0}: {1}".format(key, value))
        # setattr(self, key, value)
    eprint()

    work_dir = tempfile.mkdtemp()
    print("Temporary workdirektory:", work_dir)
    anchors_to_read_acc = {}
    k_size = args.k
    for batch_id, reads in enumerate(batch(all_reads, args.max_seqs)):
        print("correcting {0} reads in a batch".format(len(reads)))
        batch_start_time = time()
        w = args.w
        x_high = args.xmax
        x_low = args.xmin
        hash_fcn = "lex"

        minimizer_database  = get_minimizers_and_positions(reads, w, k_size, hash_fcn)
        minimizer_combinations_database = get_minimizer_combinations_database(reads, minimizer_database, k_size, x_low, x_high)
        # print(minimizer_database)
        if args.verbose:
            eprint("done creating minimizer combinations")

        tot_corr = 0
        previously_corrected_regions = defaultdict(list)
        tmp_cnt = 0

        for r_id in sorted(reads): #, reverse=True):
            read_min_comb = [ ((m1,p1), m1_curr_spans) for (m1,p1), m1_curr_spans in  minimizers_comb_iterator(minimizer_database[r_id], k_size, x_low, x_high)]
            # print(read_min_comb)
            # sys.exit()
            if args.exact:
                previously_corrected_regions = defaultdict(list)
            # stored_calculated_regions = defaultdict(list)
    
            #  = stored_calculated_regions[r_id]
            corr_pos = []
            (acc, seq, qual) = reads[r_id]
            # print("starting correcting:", seq)

            # print(r_id, sorted(previously_corrected_regions[r_id], key=lambda x:x[1]))
            read_previously_considered_positions = set([tmp_pos for tmp_p1, tmp_p2, w_tmp, _ in previously_corrected_regions[r_id] for tmp_pos in range(tmp_p1, tmp_p2)])
            
            if args.verbose:
                if read_previously_considered_positions:
                    eprint("not corrected:", [ (p1_, p2_) for p1_, p2_ in zip(sorted(read_previously_considered_positions)[:-1], sorted(read_previously_considered_positions)[1:]) if p2_ > p1_ + 1 ] )
                else:
                    eprint("not corrected: entire read", )

            if previously_corrected_regions[r_id]:
                read_previously_considered_positions = set([tmp_pos for tmp_p1, tmp_p2, w_tmp, _ in previously_corrected_regions[r_id] for tmp_pos in range(tmp_p1, tmp_p2)])
                group_id = 0
                pos_group = {}
                sorted_corr_pos = sorted(read_previously_considered_positions)
                for p1, p2 in zip(sorted_corr_pos[:-1], sorted_corr_pos[1:]):
                    if p2 > p1 + 1:
                       pos_group[p1] = group_id 
                       group_id += 1
                       pos_group[p2] = group_id 
                    else:
                       pos_group[p1] = group_id 
                if p2 == p1 + 1:
                    pos_group[p2] = group_id 
            else:
                read_previously_considered_positions= set()
                pos_group = {}

            already_computed = {}
            read_complexity_cnt = 0
            # test_cnt = 0
            # old_cnt = 0
            # test_cnt2 = 0
            all_intervals = []
            # prev_visited_intervals = []

            for (m1,p1), m1_curr_spans in read_min_comb: 
                # If any position is not in range of current corrections: then correct, not just start and stop
                not_prev_corrected_spans = [(m2,p2) for (m2,p2) in m1_curr_spans if not (p1 + k_size in read_previously_considered_positions and p2 - 1 in read_previously_considered_positions) ] 
                set_not_prev = set(not_prev_corrected_spans)
                not_prev_corrected_spans2 = [(m2,p2) for (m2,p2) in m1_curr_spans if (m2,p2) not in set_not_prev and (p1 + k_size in pos_group and p2 - 1 in pos_group and pos_group[p1 + k_size] != pos_group[p2 - 1]) ] 
                not_prev_corrected_spans += not_prev_corrected_spans2


                if not_prev_corrected_spans: # p1 + k_size not in read_previously_considered_positions:
                    tmp_cnt, read_complexity_cnt = find_most_supported_span(r_id, m1, p1, not_prev_corrected_spans, minimizer_combinations_database, reads, all_intervals, k_size, tmp_cnt, read_complexity_cnt, already_computed)

            # sys.exit()
            if args.verbose:
                print("{0} edlib invoked due to repeated anchors for this read.".format(read_complexity_cnt))
                print(tmp_cnt, "total computed editdist.")
                eprint("Correcting read", r_id)

            # add prev_visited_intervals to intervals to consider
            # all_intervals.extend(prev_visited_intervals)

            if previously_corrected_regions[r_id]: # add previously corrected regions in to the solver
                all_intervals.extend(previously_corrected_regions[r_id])
                del previously_corrected_regions[r_id]

            if not all_intervals:
                # eprint("Found nothing to correct")
                corrected_seq = seq
            else:
                all_intervals_sorted_by_finish = sorted(all_intervals, key = lambda x: x[1])
                opt_indicies = solve_WIS2(all_intervals_sorted_by_finish) # solve Weighted Interval Scheduling here to find set of best non overlapping intervals to correct over
                sol = []
                prev_stop = 0
                for j in opt_indicies:
                    start, stop, weights, instance = all_intervals_sorted_by_finish[j]
                    sol.append( ( instance["curr_read"][0][:k_size], instance["curr_read"][0][-k_size:] ) )
                    
                    if start - k_size > prev_stop and prev_stop > 0:
                        if verbose:
                            eprint("Gap in correction:", start-k_size - prev_stop, "between positions:", prev_stop, start, )
                    prev_stop = stop + k_size

                hashable_sol = tuple(sol)
                if hashable_sol not in anchors_to_read_acc:
                    anchors_to_read_acc[hashable_sol] = set()
                    anchors_to_read_acc[hashable_sol].add(acc) 
                else:
                    anchors_to_read_acc[hashable_sol].add(acc) 
            print("processed:", r_id, len(anchors_to_read_acc), len(hashable_sol), hashable_sol)
                
    eprint( "Number of unique transcripts (based on anchor solution):", len(anchors_to_read_acc))
    for anchors, read_set in anchors_to_read_acc.items():
        print(anchors, read_set)        
    eprint( "Number of unique transcripts (based on anchor solution):", len(anchors_to_read_acc))

    # To final spoa generation consensus
    
    # outfile = open(os.path.join(args.outfolder, "corrected_reads.fastq"), "w")
    # for r_id, (acc, seq, qual) in corrected_reads.items():
    #     outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))
    # outfile.close()
    # print("removing temporary workdir")
    # shutil.rmtree(work_dir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Collapsing long-read transcriptome reads into transcripts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

    parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq file with reads')
    # parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores allocated for clustering')

    parser.add_argument('--k', type=int, default=9, help='Kmer size')
    parser.add_argument('--w', type=int, default=10, help='Window size')
    parser.add_argument('--xmin', type=int, default=14, help='Upper interval length')
    parser.add_argument('--xmax', type=int, default=80, help='Lower interval length')
    parser.add_argument('--max_seqs', type=int, default=1000,  help='Maximum number of seqs to correct at a time (in case of large clusters).')
    parser.add_argument('--T', type=float, default=0.1, help='Minimum fraction keeping substitution')
    # parser.add_argument('--C', type=float, default=0.05, help='Minimum fraction of keeping alternative refernece contexts')
    parser.add_argument('--exact', action="store_true", help='Get exact solution for WIS for evary read (recalculating weights for each read (much slower but slightly more accuracy,\
                                                                 not to be used for clusters with over ~500 reads)')
    parser.add_argument('--exact_instance_limit', type=int, default=0,  help='Activates slower exact mode for instance smaller than this limit')
    parser.add_argument('--set_w_dynamically', action="store_true", help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')
    parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()


    if args.xmin < 2*args.k:
        args.xmin = 2*args.k
        eprint("xmin set to {0}".format(args.xmin))

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if not args.fastq and not args.flnc and not  args.ccs:
        parser.print_help()
        sys.exit()




    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    if 100 < args.w or args.w < args.k:
        eprint('Please specify a window of size larger or equal to k, and smaller than 100.')
        sys.exit(1)

    main(args)

