#! /usr/bin/env python

from __future__ import print_function
import os,sys
import argparse

import errno
from time import time, sleep
import itertools
import tempfile
import shutil

import math
import re
from collections import deque
from collections import defaultdict

import edlib
import parasail

from isoncorrect import create_augmented_reference, help_functions, correct_seqs #,align

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def rindex(lst, value):
    return len(lst) - operator.indexOf(reversed(lst), value) - 1

def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([hash(seq[i:i+k_size]) for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizer_pos = rindex(list(window_kmers), curr_min)
    minimizers = [ (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ] # get the last element if ties in window

    for i in range(w+1,len(seq) - k_size):
        new_kmer = hash(seq[i:i+k_size])
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous window's minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer and minimizer_pos < i - w: 
            curr_min = min(window_kmers)
            minimizer_pos = rindex(list(window_kmers), curr_min) + i - w  
            minimizers.append( (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ) # get the last element if ties in window

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (seq[i: i+k_size], i) )

    return minimizers

def get_kmer_minimizers_lex(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizer_pos = rindex(list(window_kmers), curr_min)
    minimizers = [ (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ] # get the last element if ties in window

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous window's minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer and minimizer_pos < i - w: 
            curr_min = min(window_kmers)
            minimizer_pos = rindex(list(window_kmers), curr_min) + i - w  
            minimizers.append( (seq[minimizer_pos: minimizer_pos+k_size], minimizer_pos) ) # get the last element if ties in window

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (seq[i: i+k_size], i) )

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


def get_minimizers_and_positions_compressed(reads, w, k, hash_fcn):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]

        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))

        if hash_fcn == "random":
            minimizers = get_kmer_minimizers(seq_hpol_comp, k, w)
        elif hash_fcn == "lex":
            minimizers = get_kmer_minimizers_lex(seq_hpol_comp, k, w)
        elif hash_fcn == "rev_lex":
            minimizers = get_kmer_maximizers(seq_hpol_comp, k, w)

        indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2] # indicies we want to take quality values from to get quality string of homopolymer compressed read 
        indices.append(len(seq) - 1)
        positions_in_non_compressed_sring = [(m, indices[p]) for m, p in minimizers ]
        M[r_id] = positions_in_non_compressed_sring

    return M


def get_minimizers_and_positions(reads, w, k, hash_fcn):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        if hash_fcn == "random":
            minimizers = get_kmer_minimizers(seq, k, w)
        elif hash_fcn == "lex":
            minimizers = get_kmer_minimizers_lex(seq, k, w)
        elif hash_fcn == "rev_lex":
            minimizers = get_kmer_maximizers(seq, k, w)

        M[r_id] = minimizers

    return M



from array import array
def get_minimizer_combinations_database(reads, M, k, x_low, x_high):
    # M2 = defaultdict(lambda: defaultdict(list))
    M2 = defaultdict(lambda: defaultdict(lambda :array("I")))
    tmp_cnt = 0
    forbidden = 'A'*k
    for r_id in M:
        minimizers = M[r_id]
        for (m1,p1), m1_curr_spans in  minimizers_comb_iterator(minimizers, k, x_low, x_high):
            for (m2, p2) in m1_curr_spans:
                if m2 == m1 == forbidden:
                    continue

                tmp_cnt +=1
                # t = array('I', [r_id, p1, p2])
                # M2[m1][m2].append( t )
                # M2[m1][m2].append((r_id, p1, p2))

                M2[m1][m2].append(r_id)
                M2[m1][m2].append(p1)
                M2[m1][m2].append(p2)

    print(tmp_cnt, "MINIMIZER COMBINATIONS GENERATED")
    # import time
    # time.sleep(10)
    # sys.exit()

    avg_bundance = 0
    singleton_minimzer = 0
    high_abundance = 0
    cnt = 1
    abundants=[]
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            if len(M2[m1][m2]) > 3:
                avg_bundance += len(M2[m1][m2])//3
                cnt +=1
            else:
                del M2[m1][m2]
                singleton_minimzer += 1

            if len(M2[m1][m2])// 3 > len(reads):
                abundants.append((m1,m2, len(M2[m1][m2])//3 ))
                if m2 == forbidden or len(M2[m1][m2])// 3 > 10*len(reads): # poly A tail or highly abundant
                    del M2[m1][m2]
                    high_abundance += 1
    for m1,m2,ab in sorted(abundants, key=lambda x: x[2], reverse=True):
        print("Not unique within reads:", m1, m2, ab, len(reads))

    print("Average abundance for non-unique minimizer-combs:", avg_bundance/float(cnt))
    print("Number of singleton minimizer combinations filtered out:", singleton_minimzer)
    print("Number of highly abundant minimizer combinations (10x more frequent than nr reads) or poly-A anchors filtered out :", high_abundance)

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


import operator
def argmin(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value


def randstrobe_order2(subseq, m1, k_size, prime):
    min_index, min_value = argmin([ hash(m1+ subseq[i:i+k_size]) % prime for i in range(len(subseq) - k_size + 1)])
    min_m2 = subseq[min_index:min_index+k_size]
    # print(len(m1 + min_m2))
    return  min_m2, min_index

def randstrobe_order2_list(hash_seq_list, start, stop, hash_m1, k_size, prime):
    min_index, min_value = argmin([ (hash_m1+ hash_seq_list[i]) % prime for i in range(start, stop)])
    # min_m2 = hash_subseq_list[min_index]
    # print(len(m1 + min_m2))
    return min_index

# def seq_to_strobes(r_id, seq, M2, k_size, w_min, w_max, prime, forbidden, reverse = False):
#     for p1 in range(len(seq) - 2*k_size + 1):
#         m1 = seq[p1:p1+k_size]
#         window_p_start = p1 + k_size + w_min if p1+k_size+w_max <= len(seq) else max( (p1 + k_size + w_min) -  (p1+k_size+w_max - len(seq)), p1 + k_size )
#         window_p_end = min(p1 + k_size + w_max, len(seq))
#         # assert window_p_end - window_p_start >= k_size, window_p_end - window_p_start  
#         # assert window_p_start >= p1+k_size, (window_p_start, p1+k_size)
#         # print(window_p_start, window_p_end)
#         window_seq = seq[window_p_start : window_p_end]
#         m2, min_index2 = randstrobe_order2(window_seq, m1, k_size, prime)
#         # print("min_index2", min_index2, m1,m2)
#         p2 = window_p_start + min_index2
#         # assert p2 >= p1 + k_size, (p2, p1 + k_size)
#         if m2 == m1 == forbidden:
#             continue
#         if reverse:
#             M2[m2[::-1]][m1[::-1]].append(r_id)
#             M2[m2[::-1]][m1[::-1]].append(len(seq) - p2 - k_size)
#             M2[m2[::-1]][m1[::-1]].append(len(seq) - p1 - k_size)

#         else:   
#             M2[m1][m2].append(r_id)
#             M2[m1][m2].append(p1)
#             M2[m1][m2].append(p2)


def get_randstrobes_with_positions_database_2way(reads, k_size, w_min, w_max, primes):
    
    M2 = defaultdict(lambda: defaultdict(lambda :array("I")))
    tmp_cnt = 0
    forbidden = 'A'*k_size
    # read_to_randstrobes = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        rand_strobes = defaultdict(set)
        for i, prime in enumerate(primes):
            if i % 2 == 0:
                # fw = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter(seq, k_size, w_min, w_max, prime, forbidden)]
                # fw2 = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter2(seq, k_size, w_min, w_max, prime, forbidden)]
                # seq_to_strobes_iter_list2(seq, k_size, w_min, w_max, prime, forbidden)
                # print( [ (z1, z2) for z1, z2 in zip(fw, fw2) if z1 != z2])
                # assert fw == fw2
                for m1, p1, m2, p2 in seq_to_strobes_iter2(seq, k_size, w_min, w_max, prime, forbidden):
                # for (m1, p1, (m2, p2)) in fw:
                    rand_strobes[(m1, p1)].add((m2,p2))
            else:
                # rv = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter(seq[::-1], k_size, w_min, w_max, prime, forbidden, reverse = True)]
                # rv2 = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter2(seq[::-1], k_size, w_min, w_max, prime, forbidden, reverse = True)]
                # assert rv == rv2
                for m1, p1, m2, p2 in seq_to_strobes_iter2(seq[::-1], k_size, w_min, w_max, prime, forbidden, reverse = True):
                # for (m1, p1, (m2, p2)) in rv:
                    rand_strobes[(m1, p1)].add((m2,p2))

        for (m1, p1), s in rand_strobes.items():
            for (m2, p2) in s:
                M2[m1][m2].append(r_id)
                M2[m1][m2].append(p1)
                M2[m1][m2].append(p2)
                tmp_cnt += 1

        # read_to_randstrobes[r_id] = [((m1, p2), list(s))  for (m1, p2), s in rand_strobes.items()]

    print(tmp_cnt, "RANDSTROBES 2WAY COMBINATIONS GENERATED")

    avg_bundance = 0
    singleton_minimzer = 0
    cnt = 1
    abundants=[]
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            if len(M2[m1][m2]) > 3:
                avg_bundance += len(M2[m1][m2])//3
                cnt +=1
            else:
                del M2[m1][m2]
                singleton_minimzer += 1
            # if len(M2[m1][m2])// 3 > 5:
            #     print(len(M2[m1][m2])// 3)
            if len(M2[m1][m2])// 3 > len(reads):
                abundants.append((m1,m2, len(M2[m1][m2])//3 ))
                if m2 == forbidden: # poly A tail
                    del M2[m1][m2]
    for m1,m2,ab in sorted(abundants, key=lambda x: x[2], reverse=True):
        print("Too abundant:", m1, m2, ab, len(reads))

    print("Average abundance for non-unique randstrobes:", avg_bundance/float(cnt))
    print("Number of singleton randstrobes filtered out:", singleton_minimzer)
    # sys.exit()
    return M2 #, read_to_randstrobes


def seq_to_strobes_iter(seq, k_size, w_min, w_max, prime, forbidden, reverse = False):
    for p1 in range(len(seq) - 2*k_size + 1):
        m1 = seq[p1:p1+k_size]
        window_p_start = p1 + k_size + w_min if p1+k_size+w_max <= len(seq) else max( (p1 + k_size + w_min) -  (p1+k_size+w_max - len(seq)), p1 + k_size )
        window_p_end = min(p1 + k_size + w_max, len(seq))
        # assert window_p_end - window_p_start >= k_size, window_p_end - window_p_start  
        # assert window_p_start >= p1+k_size, (window_p_start, p1+k_size)
        # print(window_p_start, window_p_end)
        window_seq = seq[window_p_start : window_p_end]
        m2, min_index2 = randstrobe_order2(window_seq, m1, k_size, prime)
        # print("min_index2", min_index2, m1,m2)
        p2 = window_p_start + min_index2
        # assert p2 >= p1 + k_size, (p2, p1 + k_size)
        if m2 == m1 == forbidden:
            continue
        if reverse:
            yield m2[::-1], len(seq) - p2 - k_size, m1[::-1], len(seq) - p1 - k_size
        else:   
            yield m1, p1, m2, p2


def seq_to_strobes_iter2(seq, k_size, w_min, w_max, prime, forbidden, reverse = False):
    hash_seq_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
    for p1 in range(len(seq) - 2*k_size + 1):
        hash_m1 = hash_seq_list[p1]
        window_p_start = p1 + k_size + w_min if p1 + w_max <= len(hash_seq_list) else max( (p1 + k_size + w_min) -  (p1+k_size+w_max - len(hash_seq_list)), p1+ k_size )
        window_p_end = min(p1 + w_max, len(hash_seq_list))
        # assert window_p_end - window_p_start >= k_size, window_p_end - window_p_start  
        # assert window_p_start >= p1+k_size, (window_p_start, p1+k_size)
        # print(window_p_start, window_p_end)
        # window_seq_list = hash_seq_list[window_p_start : window_p_end+1]
        # print(len(hash_seq_list), window_p_start, window_p_end)
        min_index2 = randstrobe_order2_list(hash_seq_list, window_p_start, window_p_end, hash_m1, k_size, prime)
        # print("min_index2", min_index2, m1,m2)
        p2 = window_p_start + min_index2
        assert p2 >= p1 + k_size, (p2, p1 + k_size)
        assert p2 <= len(seq) -k_size, (p2, len(seq) -k_size)
        m1 = seq[p1:p1+k_size]
        m2 = seq[p2:p2+k_size]
        if m2 == m1 == forbidden:
            continue
        if reverse:
            yield m2[::-1], len(seq) - p2 - k_size, m1[::-1], len(seq) - p1 - k_size
        else:   
            yield m1, p1, m2, p2

# def seq_to_strobes_iter_list2(seq, k_size, w_min, w_max, prime, forbidden, reverse = False):
#     hash_seq_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
#     return [ (i, argmin([ (hash_seq_list[i]+ hash_seq_list[j]) % prime for j in range( min(i + k_size + w_min, len(hash_seq_list)- k_size) , min(i + w_max, len(hash_seq_list) - 1))] ) ) for i in range(len(hash_seq_list) - k_size) ]
#     # return  min_index, min_value = argmin([ (hash_m1+ hash_seq_list[i]) % prime for i in range(start, stop)])


def randstrobes_read_2way(seq, k_size, w_min, w_max, primes):
    forbidden = 'A'*k_size
    rand_strobes = defaultdict(set)
    for i, prime in enumerate(primes):
        if i % 2 == 0:
            # fw = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter(seq, k_size, w_min, w_max, prime, forbidden)]
            for m1, p1, m2, p2 in seq_to_strobes_iter2(seq, k_size, w_min, w_max, prime, forbidden):
                rand_strobes[(m1, p1)].add((m2,p2))
        else:
            # rv = [(m1, p1, (m2, p2)) for m1, p1, m2, p2 in seq_to_strobes_iter(seq[::-1], k_size, w_min, w_max, prime, forbidden, reverse = True)]
            for m1, p1, m2, p2 in seq_to_strobes_iter2(seq[::-1], k_size, w_min, w_max, prime, forbidden, reverse = True):
                rand_strobes[(m1, p1)].add((m2,p2))

    return [((m1, p2), list(s))  for (m1, p2), s in rand_strobes.items()]


def get_randstrobes_with_positions_database(reads, k_size, w_min, w_max, primes):
    
    M2 = defaultdict(lambda: defaultdict(lambda :array("I")))
    tmp_cnt = 0
    forbidden = 'A'*k_size
    read_to_randstrobes = {}

    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        rand_strobes = defaultdict(set)
        for p1 in range(len(seq) - 2*k_size + 1):
            m1 = seq[p1:p1+k_size]
            window_p_start = p1 + k_size + w_min if p1+k_size+w_max <= len(seq) else max( (p1 + k_size + w_min) -  (p1+k_size+w_max - len(seq)), p1 + k_size )
            window_p_end = min(p1 + k_size + w_max, len(seq))
            # assert window_p_end - window_p_start >= k_size, window_p_end - window_p_start  
            # assert window_p_start >= p1+k_size, (window_p_start, p1+k_size)
            # print(window_p_start, window_p_end)
            window_seq = seq[window_p_start : window_p_end]
            added_m2 = set()
            added_p2 = set()
            for prime in primes:
                m2, min_index2 = randstrobe_order2(window_seq, m1, k_size, prime)
                # print("min_index2", min_index2, m1,m2)
                p2 = window_p_start + min_index2

                if m2 == m1 == forbidden:
                    continue
                if m2 in added_m2 and p2 in added_p2:
                    # print("LOOOL")
                    continue
                tmp_cnt +=1

                # assert p2 >= p1 + k_size, (p2, p1 + k_size)
                M2[m1][m2].append(r_id)
                M2[m1][m2].append(p1)
                M2[m1][m2].append(p2)
                added_m2.add(m2)
                added_p2.add(p2)
                rand_strobes[(m1, p1)].add((m2,p2))

        read_to_randstrobes[r_id] = [((m1, p2), list(s))  for (m1, p2), s in rand_strobes.items()]

    print(tmp_cnt, "RANDSTROBES COMBINATIONS GENERATED")

    avg_bundance = 0
    singleton_minimzer = 0
    cnt = 1
    abundants=[]
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            if len(M2[m1][m2]) > 3:
                avg_bundance += len(M2[m1][m2])//3
                cnt +=1
            else:
                del M2[m1][m2]
                singleton_minimzer += 1

            if len(M2[m1][m2])// 3 > len(reads):
                abundants.append((m1,m2, len(M2[m1][m2])//3 ))
                if m2 == forbidden: # poly A tail
                    del M2[m1][m2]
    for m1,m2,ab in sorted(abundants, key=lambda x: x[2], reverse=True):
        print("Too abundant:", m1, m2, ab, len(reads))

    print("Average abundance for non-unique randstrobes:", avg_bundance/float(cnt))
    print("Number of singleton randstrobes filtered out:", singleton_minimzer)
    # sys.exit()
    return M2, read_to_randstrobes



def randstrobe_iterator(seq, k_size, w_min, w_max, primes):
    forbidden = 'A'*k_size
    for p1 in range(len(seq) - 2*k_size + 1):
        m1 = seq[p1:p1+k_size]
        window_p_start = p1 + k_size + w_min if p1+k_size+w_max <= len(seq) else max( (p1 + k_size + w_min) -  (p1+k_size+w_max - len(seq)), p1 + k_size )
        window_p_end = min(p1 + k_size + w_max, len(seq))
        # assert window_p_end - window_p_start >= k_size, window_p_end - window_p_start  
        # assert window_p_start >= p1+k_size, (window_p_start, p1+k_size)
        window_seq = seq[window_p_start : window_p_end]
        m1_curr_spans = []
        added_m2 = set()
        added_p2 = set()
        for prime in primes:
            m2, min_index2 = randstrobe_order2(window_seq, m1, k_size, prime)
            p2 = window_p_start + min_index2

            if m2 == m1 == forbidden:
                continue

            if m2 in added_m2 and p2 in added_p2:
                # print("LEWL")
                continue

            # assert p2 >= p1 + k_size, (p2, p1 + k_size)
            m1_curr_spans.append((m2, p2))
            added_m2.add(m2)
            added_p2.add(p2)

        yield m1, p1, m1_curr_spans




def edlib_alignment(x, y, k):
    # if i == 100 and j % 1000 == 0:
    #     print("Edlib processed alignments: {0}".format(j+1))

    result = edlib.align(x,y, "NW", 'dist', k) # , task="path")
    ed = result["editDistance"]
    # locations = result["locations"]
    return ed #, locations


def get_context_offset(vector, k):
    nuc_obs = 0
    if not vector:
        return 0
    for offset, n in enumerate(vector):
        if n != '-':
            nuc_obs += 1

        if nuc_obs > k:
            break
    return offset

# from modules import kmer_analysis

def get_contexts(alignment_matrix, k_size):
    ref_aln = alignment_matrix["ref"]
    contexts_per_pos = []
    for i_tmp in range(len(ref_aln)):
        nuc_obs = 0
        # j_tmp = i_tmp
        p1 = get_context_offset(ref_aln[:i_tmp][::-1], k_size)
        p2 = get_context_offset(ref_aln[i_tmp + 1:], k_size) + 1
        # print(p1,p2)
        # print(ref_aln)
        contexts_per_pos.append((i_tmp - p1, i_tmp + p2))
    return contexts_per_pos

def is_substring(seq, set_of_seqs):
    for s in set_of_seqs:
        if seq in s or s in seq:
            return True
    return False

import numpy as np
# import numba
# from numba import jit
# @jit(nopython=True)
def test_numba(A, contexts_per_pos, n, context_threshold):
    FCM = [] # frequency context matrix

    prev_context_start, prev_context_stop = -1,-1
    prev_context = []

    for i in range(n):
        eligible_contexts = []
        context_start, context_stop = contexts_per_pos[i][0], contexts_per_pos[i][1]
        variant_pos = i - context_start
        if context_start == prev_context_start and context_stop == prev_context_stop:
            # use previous context and depths but shift variant position one forward
            for v, cn, dep in FCM[i-1]:
                # print(v, cn, dep, variant_pos, i, context_start)
                # print(len(cn), variant_pos, i, "window",context_start, context_stop, n)
                # eligible_contexts.append((cn[min(variant_pos, len(cn)-1)] ,cn, dep))
                eligible_contexts.append((cn[variant_pos] ,cn, dep))

            FCM.append( eligible_contexts )
            continue
        else:
            prev_context_start, prev_context_stop = context_start, context_stop

        context = A[ :, context_start :context_stop] 
        b = np.ascontiguousarray(context).view(np.dtype((np.void, context.dtype.itemsize * context.shape[1])))
        unq_a, unq_cnt = np.unique(b, return_counts=True)
        unq_a = unq_a.view(context.dtype).reshape(-1, context.shape[1])
        # print( unq_a)
        # print('OK', unq_cnt)         
        # print('Unique Values along with occurrence Count')
        # Iterate over the zip object
        for (j, count) in sorted(enumerate(unq_cnt), key = lambda x: x[1], reverse=True):
            # print(count)
            if count > max(context_threshold/10, 5):
                # print("im here")
                context = tuple(unq_a[j].tolist())
                # try:
                # print('lol', context, variant_pos, count)
                # variant = context[min(variant_pos, len(context)-1)] #context[variant_pos]
                variant = context[variant_pos] #context[variant_pos]
                eligible_contexts.append( (variant, context, count ) )
                # FCM.[i][ (variant, context ) ] = count
                # except:
                #     pass
            else: 
                # print('Breaking', count)
                break
        FCM.append(eligible_contexts)
    return FCM


def sep_function_test(alignment_matrix, FCM, contexts_per_pos):
    for acc, aln_tuple in alignment_matrix.items():
        # aln_tuple = tuple(aln_list)
        if acc == "ref":
            continue
        subtuples = [ (j, aln_tuple[contexts_per_pos[j][0] :contexts_per_pos[j][1]]) for j in range(1, len(aln_tuple)) if contexts_per_pos[j][0] != contexts_per_pos[j-1][0] or contexts_per_pos[j][1] != contexts_per_pos[j-1][1] ]
        for j in range(len(aln_tuple)):
            FCM[j][ (aln_tuple[j], subtuples[j]) ] += 1

    #     prev_context_start, prev_context_stop = -1,-1
    #     prev_context = tuple()
    #     for j in range(len(aln_tuple)):
    #         context_start, context_stop = contexts_per_pos[j][0], contexts_per_pos[j][1]
    #         if context_start == prev_context_start and context_stop == prev_context_stop:
    #             context = prev_context
    #         else:
    #             context =  aln_tuple[context_start :context_stop] #prev_context
    #             # context =  aln_tuple[context_start :context_stop] #prev_context
    #             # context = "".join([s for s in aln_tuple[context_start :context_stop]]) #tuple( aln_tuple[context_start :context_stop] )
            
    #         variant = aln_tuple[j]
    #         FCM[j][ (variant, context) ] += 1
    #         prev_context_start, prev_context_stop = context_start, context_stop
    #         prev_context = context
    # # return FCM

def sep_function(alignment_matrix, FCM, contexts_per_pos):
    for acc, aln_list in alignment_matrix.items():
        aln_tuple = tuple(aln_list)
        # if acc == "ref":
        #     continue
        prev_context_start, prev_context_stop = -1,-1
        prev_context = []
        for j in range(len(aln_tuple)):
            context_start, context_stop = contexts_per_pos[j][0], contexts_per_pos[j][1]
            if context_start == prev_context_start and context_stop == prev_context_stop:
                context = prev_context
            else:
                context = tuple( aln_tuple[context_start :context_stop] )
            
            variant = aln_tuple[j]
            FCM[j][ (variant, context) ] += 1
            prev_context_start, prev_context_stop = context_start, context_stop
            prev_context = context
    # return FCM

# def sep_function_only_str(partition):
#     FCM = [defaultdict(int) for j in range(nr_columns)] # frequency context matrix

#     for (q_id, pos1, pos2) in partition = 
#     (res["editDistance"], ref_alignment, read_alignment, 1)
#     return FCM

def get_alternative_ref_contexts(alignment_matrix, contexts_per_pos, context_threshold, disable_numpy):

    ref_aln = alignment_matrix["ref"]
    nr_columns = len(ref_aln)
    if disable_numpy:
        FCM = [defaultdict(int) for j in range(nr_columns)] # frequency context matrix
        sep_function(alignment_matrix, FCM, contexts_per_pos)
    else:
        A = np.array([l for l in alignment_matrix.values()])
        n = len(A[0])
        FCM = test_numba(A, contexts_per_pos, n, context_threshold)
    # print(len(FCM), len(FCM3))
    # assert len(FCM) == len(FCM3)

    alternative_contexts = [set() for j in range(nr_columns)]
    for j in range(len(FCM)):
        # eligible_contexts2 = [ (variant, tuple([c for c in cn]) , dep) for (variant,cn), dep in FCM[j].items() if dep > max(context_threshold/10, 5)] # remove very low abundant noise before calculation
        eligible_contexts = FCM[j]

        # print(eligible_contexts)
        # if len(eligible_contexts2) != len(eligible_contexts):
        #     print(eligible_contexts)
        #     print(eligible_contexts2)
        #     print(max(context_threshold/10, 5))
        #     print( ref_aln[j-9:j+9])
        # assert len(eligible_contexts2) == len(eligible_contexts)
        # if set( [(v,cn) for v,cn,d in eligible_contexts]) !=  set( [(v,cn) for v,cn,d in eligible_contexts2]):
        #     print("Different" )
        #     print(j,eligible_contexts)
        #     print(j,eligible_contexts2)
        #     print(len(FCM), ref_aln[j], ref_aln)
        #     print()
        # else:
        #     pass
        #     # print('Same')

        # print(eligible_contexts)
        # print(FCM3[j])
        # print()
        # filter FCM to only positions with more than one variants
        if len(eligible_contexts) == 0 or len( {v for (v, ref_tmp, dep) in eligible_contexts}) == 1:
            # print(j)
            pass
        # elif len(eligible_contexts) == 1:
        #     # print(eligible_contexts)
        #     # print(eligible_contexts[0])
        #     alternative_contexts[j].add(eligible_contexts[0])
        else:
            eligible_contexts_hcomp = {}
            context_start, context_stop = contexts_per_pos[j][0], contexts_per_pos[j][1]
            consensus_context = ''.join( [ c for c in ref_aln[context_start :context_stop] if c != '-'] )
            consensus_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(consensus_context))
            eligible_contexts_hcomp[consensus_hpol_comp] = set()
            eligible_contexts_hcomp[consensus_context] = set()
            # for variant, aln_seq, depth in eligible_contexts:
            #     print(variant, aln_seq, depth, j)
            # if j ==199:
            #     for acc_z, aln_list_z in alignment_matrix.items():
            #         print(aln_list_z[context_start-15:context_stop+15])
            #     print(ref_aln[context_start-15:context_stop+15], "consensus")
            for variant, aln_seq, depth in sorted(eligible_contexts, key = lambda x: x[2], reverse=True): # highest depth first
                seq = "".join([c for c in aln_seq if c != '-'])
                seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
                if seq_hpol_comp in eligible_contexts_hcomp:
                    continue
                else:
                    min_ed = 1000
                    new_alt_ref  = ''
                    for prev_seq in eligible_contexts_hcomp:
                        ed = edlib_alignment( prev_seq, seq, 1000)
                        # print(ed)
                        if ed < min_ed:
                            min_ed = ed
                            new_alt_ref = seq
                        ed = edlib_alignment( prev_seq, seq_hpol_comp, 1000)
                        # print(ed)
                        if ed < min_ed:
                            min_ed = ed
                            new_alt_ref = seq_hpol_comp

                    threshold = context_threshold/max(1.0, min_ed)
                    if depth >= threshold:
                        # print(depth, context_threshold, threshold, min_ed)
                        eligible_contexts_hcomp[new_alt_ref] = (variant, depth, aln_seq, threshold, seq)

                # print(variant,seq, seq_hpol_comp, depth)
                # if seq_hpol_comp not in eligible_contexts_hcomp:
                #     # if not is_substring(seq_hpol_comp,eligible_contexts_hcomp):
                #     eligible_contexts_hcomp[seq_hpol_comp] = (variant, depth, aln_seq)
                # elif depth > eligible_contexts_hcomp[seq_hpol_comp][1]:
                #     eligible_contexts_hcomp[seq_hpol_comp] = (variant, depth, aln_seq)
            # print(j, eligible_contexts_hcomp)
            del eligible_contexts_hcomp[consensus_hpol_comp]
            if consensus_hpol_comp != consensus_context:
                del eligible_contexts_hcomp[consensus_context]

            if len(eligible_contexts_hcomp) > 0:
                for seq_tmp_hcomp, (variant, depth, aln_seq, threshold, seq) in eligible_contexts_hcomp.items():
                    alternative_contexts[j].add( (variant, aln_seq, depth, threshold) )

            # if len( {v for (v, ref_tmp, dep, thresh_) in alternative_contexts[j]}) > 0:
            #     print("MAAAAAAAAAADE IT", j, len(alignment_matrix), context_threshold, [(variant,depth,threshold, "".join([c for c in aln_seq if c != '-']))  for (variant, aln_seq, depth, threshold)  in alternative_contexts[j]] )
            #     print(consensus_context)
                # if j ==199:
                #     sys.exit()
        # BUUUG no alternative refs
    # print([ ("".join([c for c in cn if c != '-']) , dep) for cn, dep in FCM[50].items() if ])
    # sys.exit()
    # for i_tmp in range(len(ref_aln)):
    #     depths_per_pos[i_tmp] = defaultdict(int)
    #     for 
    return alternative_contexts


def get_best_corrections(curr_best_seqs, reads, k_size, work_dir,  v_depth_ratio_threshold = 0.1, max_seqs_to_spoa = 200, disable_numpy=False, use_racon = False):
    """
        curr_best_seqs is an array with q_id, pos1, pos2
        the current read is on index 0 in curr_best_seqs array
    """
    weight = len(curr_best_seqs)//3
    curr_read_id = curr_best_seqs[0]

    # print()
    # print()
    # print(weight)
    # print()
    # print()

    reads_path = open(os.path.join(work_dir, "reads_tmp.fa"), "w")
    for i, (q_id, pos1, pos2) in  enumerate(grouper(curr_best_seqs, 3)):
        seq = reads[q_id][1][pos1: pos2 + k_size]
        if i > max_seqs_to_spoa:
            break
        reads_path.write(">{0}\n{1}\n".format(str(q_id)+str(pos1)+str(pos2), seq))
    reads_path.close()
    # print(reads_path.name)
    # sys.exit()
    spoa_ref = create_augmented_reference.run_spoa(reads_path.name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
    # print(spoa_ref)
    # spoa_ref_m = create_augmented_reference.run_spoa_m(reads_path.name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
    # spoa_ref_m2 = create_augmented_reference.run_spoa_m2(reads_path.name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
    if use_racon and weight > 2:
        read_alignments_paf = open(os.path.join(work_dir, "reads_tmp.paf"), "w")
        for i, (q_id, pos1, pos2) in  enumerate(grouper(curr_best_seqs, 3)):
            seq = reads[q_id][1][pos1: pos2 + k_size]
            if i > max_seqs_to_spoa:
                break
            read_alignments_paf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(str(q_id)+str(pos1)+str(pos2), len(seq), 0, len(seq) - 1, "+", "consensus", len(spoa_ref), 0, len(spoa_ref) - 1, len(spoa_ref) - k_size, len(spoa_ref), 255 ))
        read_alignments_paf.close()
        consensus_file = open(os.path.join(work_dir,"spoa.fa"), 'w')
        consensus_file.write(">{0}\n{1}\n".format('consensus', spoa_ref))
        consensus_file.close()
        racon_ref = create_augmented_reference.run_racon(reads_path.name, read_alignments_paf.name, consensus_file.name, work_dir, 1, 1)  #run_spoa_m2(reads_path.name, os.path.join(work_dir,"spoa_tmp.fa"), "spoa")
        if spoa_ref != racon_ref:
            print(weight)
            print(spoa_ref)
            print(racon_ref)
            print()
        if racon_ref:
            spoa_ref = racon_ref

    partition = {"ref" : (0, spoa_ref, spoa_ref, 1)}
    for q_id, pos1, pos2 in  grouper(curr_best_seqs, 3):
        seq = reads[q_id][1][pos1: pos2 + k_size]
    # for q_id, (seq, pos1, pos2) in curr_best_seqs.items():
        res = edlib.align(seq, spoa_ref, task="path", mode="NW")
        cigar_string = res["cigar"]
        read_alignment, ref_alignment = help_functions.cigar_to_seq(cigar_string, seq, spoa_ref)
        partition[(q_id, pos1, pos2)] = (res["editDistance"], ref_alignment, read_alignment, 1)

    alignment_matrix = correct_seqs.create_multialignment_matrix(partition)

    # if spoa_ref != spoa_ref_m or spoa_ref != spoa_ref_m2:
    #     print(weight)
    #     print(spoa_ref)
    #     print(spoa_ref_m)
    #     print(spoa_ref_m2)
    #     # print(racon_ref)
    # if spoa_ref == "GGTCGGCGACCGGAGTCACAGCGCGACCAACGGGCAAAGGCCCATAGGCTTTTCCATCGGCA":
    #     for acc_z, aln_list_z in alignment_matrix.items():
    #         print("".join([c for c in aln_list_z]) )
    #     print("".join([c for c in alignment_matrix["ref"]]), "consensus")      
    #     print()

    nr_columns = len(alignment_matrix["ref"])
    PFM = [{"A": 0, "C": 0, "G": 0, "T": 0, "U" : 0, "-": 0, "N": 0} for j in range(nr_columns)]
    for r_tmp, aln_list in alignment_matrix.items():
        if r_tmp == "ref":
            continue
        for j, n in enumerate(aln_list):
            PFM[j][n] += 1

    other_corrections = defaultdict(list)

    # TODO: potentially add context around substitutions here
    variant_threshold = max(3, weight * v_depth_ratio_threshold)
    context_threshold = max(3, weight * v_depth_ratio_threshold) # max(5, len(curr_best_seqs) * context_depth_ratio_threshold)
    # print(variant_threshold, context_threshold, len(curr_best_seqs), spoa_ref)

    contexts_per_pos = get_contexts(alignment_matrix, int(k_size/2))
    alternative_refs = get_alternative_ref_contexts(alignment_matrix, contexts_per_pos, context_threshold, disable_numpy)
    ref_aln = alignment_matrix["ref"]
    # if spoa_ref == "GGTCGGCGACCGGAGTCACAGCGCGACCAACGGGCAAAGGCCCATAGGCTTTTCCATCGGCA":
    #     for zz, r in enumerate(alternative_refs):
    #         print(zz, r)

    for i, d in enumerate(PFM):
        # calculate reference context here the first thing we do from ref_aln
        if len(alternative_refs[i]) >= 1 and len( {v for (v, ref_tmp, dep, thresh_) in alternative_refs[i]}) >= 1:
            # if spoa_ref == "GGTCGGCGACCGGAGTCACAGCGCGACCAACGGGCAAAGGCCCATAGGCTTTTCCATCGGCA":
            #     print(alternative_refs[i])
            for q_id_tuple in alignment_matrix:    
                if q_id_tuple == "ref":
                    continue

                read_aln = alignment_matrix[q_id_tuple]
                read_nucl = read_aln[i]
                ref_nucl = alignment_matrix["ref"][i]

                if read_nucl == ref_nucl:
                    other_corrections[q_id_tuple].append(ref_nucl)
                    # print(read_nucl, ref_nucl, [ (v, "".join([c for c in ref_tmp if c != '-']), dep) for (v, ref_tmp,dep) in  alternative_refs[i]])
                    continue

                
                read_context =  tuple(read_aln[contexts_per_pos[i][0] : contexts_per_pos[i][1]])
                consensus_context =  tuple(ref_aln[contexts_per_pos[i][0] : contexts_per_pos[i][1]])
                min_ed = edlib_alignment( "".join([c for c in consensus_context if c != '-']), "".join([c for c in read_context if c != '-']), 1000)
                ref_to_correct_to = ref_nucl
                # ref_to_correct_to = read_nucl
                corr_d = 0
                for (variant, ref_context, depth, thresh_) in alternative_refs[i]:
                    if read_context == ref_context:
                        min_ed = 0
                        ref_to_correct_to = variant
                        corr_d = depth
                        # other_corrections[q_id_tuple].append(read_nucl)
                        # print('Corr (identical) to', (variant, "".join([c for c in ref_context if c != '-']), depth), [ (v, "".join([c for c in ref_tmp if c != '-']), dep) for (v, ref_tmp,dep) in  alternative_refs[i]])
                        break

                    ed = edlib_alignment( "".join([c for c in ref_context if c != '-']), "".join([c for c in read_context if c != '-']), 1000)
                    if ed < min_ed:
                        min_ed = ed
                        ref_to_correct_to = variant
                        corr_d = depth
                
                other_corrections[q_id_tuple].append(ref_to_correct_to)
                # print()
                # print( [ (v, "".join([c for c in ref_tmp if c != '-']), dep) for (v, ref_tmp,dep) in  alternative_refs[i]] )
                # print('Correcting to:',ref_to_correct_to, min_ed, corr_d, "".join([c for c in read_context if c != '-']))
                # print()
        else:

            ref_nucl = alignment_matrix["ref"][i]
            if ref_nucl != "-":
                for q_id_tuple in alignment_matrix:
                    if q_id_tuple == "ref":
                        continue
                    read_aln = alignment_matrix[q_id_tuple]
                    read_nucl = read_aln[i]
                    if read_nucl != "-" and read_nucl != ref_nucl:
                        if variant_threshold <= d[read_nucl]:
                            if ref_aln[contexts_per_pos[i][0]: i] == read_aln[contexts_per_pos[i][0]: i] and ref_aln[i+1: contexts_per_pos[i][1]] == read_aln[i+1: contexts_per_pos[i][1]]: # potential variant position
                                other_corrections[q_id_tuple].append(read_nucl)
                            else:
                                other_corrections[q_id_tuple].append(ref_nucl)

                        else:
                            other_corrections[q_id_tuple].append(ref_nucl)
                    else:
                        other_corrections[q_id_tuple].append(ref_nucl)



    other_corrections_final = defaultdict(list)
    for q_id_tuple in other_corrections:
        other_read_corr = other_corrections[q_id_tuple]
        other_read_corr = "".join([n for n in other_read_corr if n != "-"]) #spoa_ref #
        # print(q_id_tuple)
        q_id, q_p1, q_p2 = q_id_tuple
        # print()
        other_corrections_final[q_id].append( (q_p1 + k_size, q_p2, weight, other_read_corr[k_size:-k_size] ))

        if q_id == curr_read_id:
            curr_read_corr = "".join([n for n in other_read_corr if n != "-"]) #spoa_ref #

    return curr_read_corr[k_size:-k_size], other_corrections_final




def fill_p2(p, all_intervals_sorted_by_finish):

    stop_to_max_j = {stop : j for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish) if start < stop }
    # print(stop_to_max_j)
    all_choord_to_max_j = []
    j_max = 0
    # print("L", all_intervals_sorted_by_finish[-1][1])
    for i in range(0, all_intervals_sorted_by_finish[-1][1] +1):
        if i in stop_to_max_j:
            j_max = stop_to_max_j[i]
        
        all_choord_to_max_j.append(j_max)
    # print("AAAAAAAAA", len(all_choord_to_max_j))
    # print(all_choord_to_max_j)
    for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish):
        # print(start)
        j_max = all_choord_to_max_j[start]
        p.append(j_max)
    return p

def solve_WIS(all_intervals_sorted_by_finish):
    # Using notation from https://courses.cs.washington.edu/courses/cse521/13wi/slides/06dp-sched.pdf
    # print("instance size", len(all_intervals_sorted_by_finish))
    # p = [None]
    # fill_p(p, all_intervals_sorted_by_finish)
    p = [None]    
    fill_p2(p, all_intervals_sorted_by_finish)
    # if p != p2:
    #     print(p)
    #     print(p2)
    # assert p == p2
    epsilon = 0.0001
    # w - 1 since the read interval isself is included in the instance
    v = [None] + [(w - 1)*(stop-start + epsilon) for (start, stop, w, _) in all_intervals_sorted_by_finish]
    OPT = [0]
    # print(len(v), len(p), len(all_intervals_sorted_by_finish) +1)
    # print(p)
    for j in range(1, len(all_intervals_sorted_by_finish) +1):
        # print(v[j])
        # print(p[j])
        # print( len(p), j, len(OPT),p[j] )
        # print(OPT[p[j]])
        # print(OPT[j-1])
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


from itertools import zip_longest
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def add_items(seqs, r_id, p1, p2):
    seqs.append(r_id)
    seqs.append(p1)
    seqs.append(p2)

def find_most_supported_span(r_id, m1, p1, m1_curr_spans, minimizer_combinations_database, reads, all_intervals, k_size, tmp_cnt, read_complexity_cnt, quality_values_database, already_computed):

    acc, seq, qual = reads[r_id]
    for (m2,p2) in m1_curr_spans:
        # print('here', p1,p2)
        relevant_reads = minimizer_combinations_database[m1][m2]
        # seqs = array("I") #{} #defaultdict(list)
        to_add = {}
        added_strings = {} 
        # locations = {}
        # not_added_strings = set() 
        if len(relevant_reads)//3 >= 3: 
            # cnt += 1
            ref_seq = seq[p1  : p2 + k_size]
            # ref_qual = qual[p1 : p2 + k_size]            
            p_error_ref = (quality_values_database[r_id][p2 + k_size] - quality_values_database[r_id][p1])/(p2 + k_size - p1)

            # seqs["curr_read"] = (p1, p2)
            # add_items(seqs, "curr_read", p1, p2)
            # add_items(seqs, r_id, p1, p2) 
            to_add[r_id] = (r_id, p1, p2, 0)
            # locations[0] = len(seqs) - 3
            added_strings[ref_seq] = 0
            reads_visited = {}
            for relevant_read_id, pos1, pos2 in grouper(relevant_reads, 3): #relevant_reads:
                if r_id  == relevant_read_id:
                    continue
                
                read_seq = reads[relevant_read_id][1][pos1: pos2 + k_size]
                #read_qual = reads[relevant_read_id][2][pos1: pos2 + k_size]

                if read_seq == ref_seq:
                    # seqs[relevant_read_id] = (pos1, pos2)
                    # add_items(seqs, relevant_read_id, pos1, pos2)
                    to_add[relevant_read_id] = (relevant_read_id, pos1, pos2, 0)
                    # locations[relevant_read_id] = len(seqs) - 3
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
                    # seqs[relevant_read_id] = (pos1, pos2)
                    # add_items(seqs, relevant_read_id, pos1, pos2)
                    if relevant_read_id in to_add and added_strings[read_seq] >= to_add[relevant_read_id][3]:
                        continue
                    else:
                        to_add[relevant_read_id] = (relevant_read_id, pos1, pos2, added_strings[read_seq])

                    # locations[relevant_read_id] = len(seqs) - 3
                    reads_visited[relevant_read_id] = added_strings[read_seq]
                    already_computed[relevant_read_id] = (p1,p2,pos1,pos2, added_strings[read_seq])
                    continue

                elif relevant_read_id in already_computed:
                    curr_ref_start, curr_ref_end, curr_read_start, curr_read_end, curr_ed = already_computed[relevant_read_id]
                    if (curr_read_start <= pos1 and pos2 <= curr_read_end) and (curr_ref_start <= p1 and p2 <=  curr_ref_end):
                        p_error_read = (quality_values_database[relevant_read_id][pos2 + k_size] - quality_values_database[relevant_read_id][pos1])/(pos2 + k_size - pos1)
                        p_error_sum_thresh = (p_error_ref + p_error_read)*len(ref_seq) #max(8, (p_error_ref + p_error_read)*(1/3)*len(ref_seq)) # roughly a 1/3 of the errors are indels # curr_p_error_sum_thresh*len(ref_seq)
                        read_beg_diff = pos1 - curr_read_start
                        read_end_diff = pos2 - curr_read_end
                        ref_beg_diff = p1 - curr_ref_start
                        ref_end_diff = p2 - curr_ref_end
                        ed_est = curr_ed + math.fabs(ref_end_diff - read_end_diff) + math.fabs(read_beg_diff - ref_beg_diff) 
                        if 0 <= ed_est <= p_error_sum_thresh: # max(8, p_error_sum_thresh*len(ref_seq)):
                            # seqs[relevant_read_id] = (pos1, pos2)
                            # add_items(seqs, relevant_read_id, pos1, pos2)
                            if relevant_read_id in to_add and ed_est >= to_add[relevant_read_id][3]:
                                continue
                            else:
                                to_add[relevant_read_id] = (relevant_read_id, pos1, pos2, ed_est)
                            # locations[relevant_read_id] = len(seqs) - 3
                            added_strings[read_seq] = ed_est
                            reads_visited[relevant_read_id] = ed_est

                            continue

                    else:
                        pass
                


                p_error_read = (quality_values_database[relevant_read_id][pos2 + k_size] - quality_values_database[relevant_read_id][pos1])/(pos2 + k_size - pos1)
                # p_error_sum_thresh = (p_error_ref + p_error_read)*len(ref_seq)  #max(8, (p_error_ref + p_error_read)*(1/3)*len(ref_seq)) # max(8,(p_error_ref + p_error_read)*(1/3)*len(ref_seq)) # roughly a 1/3 of the errors are indels  #sum([D[char_] for char_ in read_qual])/len(read_qual) #+ 0.1
                editdist = edlib_alignment(ref_seq, read_seq, len(ref_seq))
                # p_error_sum_thresh = p_error_ref + p_error_read #sum([D[char_] for char_ in read_qual])/len(read_qual) #+ 0.1
                # if p_error_sum_thresh*len(ref_seq) < 5:
                #     print(p_error_sum_thresh*len(ref_seq), p_error_sum_thresh,len(ref_seq))
                # editdist = edlib_alignment(ref_seq, read_seq, p_error_sum_thresh*len(ref_seq))

                tmp_cnt += 1
                if editdist >= 0:    # passing second edit distance check
                    if False: #relevant_read_id in reads_visited: # we have already seen the minimizer combination
                        # prev_pos1, prev_pos2 = seqs[relevant_read_id]
                        prev_pos1, prev_pos2 = seqs[ locations[relevant_read_id] + 1], seqs[ locations[relevant_read_id] + 2]
                        prev_read_seq = reads[relevant_read_id][1][prev_pos1: prev_pos2 + k_size]
                        editdist_prev = edlib_alignment(ref_seq, prev_read_seq, len(ref_seq))
                        tmp_cnt += 1
                        read_complexity_cnt += 1

                        if editdist < editdist_prev:
                            # seqs[relevant_read_id] = (pos1, pos2)
                            seqs[locations[relevant_read_id] + 1] = pos1
                            seqs[locations[relevant_read_id] + 2] = pos2
                            added_strings[read_seq] = editdist
                            reads_visited[relevant_read_id] = editdist
                            already_computed[relevant_read_id] = (p1, p2, pos1, pos2, editdist)
                            # print("REPLACED OLD MATCH")
                        # else:
                        #     # seqs[relevant_read_id] = (prev_pos1, prev_pos2)
                        #     added_strings[prev_read_seq] = editdist_prev
                        #     reads_visited[relevant_read_id] = editdist_prev
                        #     already_computed[relevant_read_id] = (p1,p2,prev_pos1, prev_pos2, editdist_prev)
                    else:
                        # seqs[relevant_read_id] = (pos1, pos2)
                        # add_items(seqs, relevant_read_id, pos1, pos2)
                        if relevant_read_id in to_add and editdist >= to_add[relevant_read_id][3]:
                            continue
                        else:
                            to_add[relevant_read_id] = (relevant_read_id, pos1, pos2, editdist)

                        # locations[relevant_read_id] = len(seqs) - 3
                        added_strings[read_seq] = editdist
                        reads_visited[relevant_read_id] = editdist
                        already_computed[relevant_read_id] = (p1,p2,pos1,pos2, editdist)

            seqs = array("I")
            # print(to_add)
            for relev_r_id in to_add:
                add_items(seqs, to_add[relev_r_id][0], to_add[relev_r_id][1], to_add[relev_r_id][2])

            all_intervals.append( (p1 + k_size, p2,  len(seqs)//3, seqs) )
    # del seqs
    return tmp_cnt, read_complexity_cnt

def get_intervals_to_correct(opt_indicies, all_intervals_sorted_by_finish):
    intervals_to_correct =[]
    for j in opt_indicies:
        start, stop, weights, instance = all_intervals_sorted_by_finish[j]
        intervals_to_correct.append( (start, stop, weights, instance))

    return intervals_to_correct



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

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples


def parasail_alignment(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 24, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1, s2)
    # print()
    # print(s1_alignment)
    # print(s2_alignment)
    # print(cigar_string)
    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score


def fix_correction(orig, corr):
    seq = []
    o_segm = []
    c_segm = []
    l = 0
    for o, c in zip(orig,corr):
        if o != '-' and c != '-':
            if l > 10: # take original read segment
                seq.append( ''.join([x for x in o_segm if x != '-']) )
            elif l > 0: # take corrected read segment
                seq.append( ''.join([x for x in c_segm if x != '-']) )
            seq.append(c)
            l=0
            o_segm = []
            c_segm = []
        elif o == '-':
            c_segm.append(c)
            l += 1
        elif c == '-':
            o_segm.append(o)
            l += 1
        else:
            raise("Parsing alignment error of parasail's alignment")

    # if ending in an indel
    if l > 10: # take original read segment
        seq.append( ''.join([x for x in o_segm if x != '-']) )
    elif l > 0: # take corrected read segment
        seq.append( ''.join([x for x in c_segm if x != '-']) )
    l=0
    o_segm = []
    c_segm = []
    return ''.join([s for s in seq])


def correct_read(seq, reads, intervals_to_correct, k_size, work_dir, v_depth_ratio_threshold,  max_seqs_to_spoa, disable_numpy, verbose, use_racon):
    corr_seq = []
    # print(opt_indicies)
    other_reads_corrected_regions = defaultdict(list)
    prev_stop = 0
    for start, stop, weights, instance in intervals_to_correct:
    # for j in opt_indicies:
    #     start, stop, weights, instance = all_intervals_sorted_by_finish[j]
        if start - k_size > prev_stop and prev_stop > 0:
            # print()
            if verbose:
                eprint("Gap in correction:", start-k_size - prev_stop, "between positions:", prev_stop, start, )
            # print()
            # sys.exit()
        prev_stop = stop + k_size

        if isinstance(instance, str): # already corrected
            best_corr = instance
        else:
            best_corr, other_corrections = get_best_corrections(instance, reads, k_size, work_dir, v_depth_ratio_threshold, max_seqs_to_spoa, disable_numpy, use_racon) # store all corrected regions within all reads in large container and keep track when correcting new read to not re-compute these regions     
            for other_r_id, other_corr_regions in other_corrections.items():
                for region in other_corr_regions:
                    other_reads_corrected_regions[other_r_id].append(region)
        # print(seq[start: stop],  best_corr)
        corr_seq.append((start,stop, best_corr))
    # print(corr_seq)
    tmp = [seq[0 : corr_seq[0][0]] ]
    for cnt, (start_, stop_, seq_segment) in enumerate(corr_seq):
        tmp.append(seq_segment)
        # print(tmp)
        # if math.fabs( (stop_ - start_) - len(seq_segment)) > 20:
        #     print("Structural correction", stop_ - start_, len(seq_segment))
        #     print(seq_segment)
        #     print(seq[start_: stop_])

        if cnt == len(corr_seq) - 1:
            tmp.append( seq[ stop_ : ] )
        else:
            # print(cnt)
            tmp.append( seq[ stop_ : corr_seq[cnt+1][0]] )
                    
    corr = "".join([s for s in tmp])

    # check for structural overcorrections
    start = time()
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(seq, corr, match_score=4, mismatch_penalty=-8, opening_penalty=12, gap_ext=1)
    # print('Alignment took: ', time() - start )
    adjusted_corr = fix_correction(s1_alignment, s2_alignment)
    s1_alignment, s2_alignment, cigar_string, cigar_tuples, score = parasail_alignment(seq, adjusted_corr, match_score=4, mismatch_penalty=-8, opening_penalty=12, gap_ext=1)

    return adjusted_corr, other_reads_corrected_regions


D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}

def get_qvs(reads):
    quality_values_database = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        quality_values_database[r_id] = [0]
        tmp_tot_sum = 0
        for char_ in qual:
            qv = D[char_]
            quality_values_database[r_id].append( tmp_tot_sum + qv )  #= [D[char_] for char_ in qual]
            tmp_tot_sum += qv
    return quality_values_database



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
    elif len(dictionary) == 1:
        batches.append(sub_dict)        
    
    return batches

def get_primes(n, nprimes):
    primes = []
    for num in range(n, 2, -1):
        if all(num%i!=0 for i in range(2,int(math.sqrt(num))+1)):
            primes.append(num)
            if len(primes) >= nprimes:
                return primes


def isoncorrect_main(args):
    # start = time()
    all_reads = { i + 1 : (acc, seq, qual) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(args.fastq, 'r')))}
    eprint("Total cluster of {0} reads.".format(len(all_reads)))
    max_seqs_to_spoa = args.max_seqs_to_spoa
    if len(all_reads) <= args.exact_instance_limit:
        args.exact = True

    eprint("ARGUMENT SETTINGS:")
    for key, value in args.__dict__.items():
        eprint("{0}: {1}".format(key, value))
        # setattr(self, key, value)
    eprint()

    work_dir = tempfile.mkdtemp()
    print("Temporary workdirektory:", work_dir)

    # start = time()
    corrected_reads = {}
    v_depth_ratio_threshold = args.T
    # context_depth_ratio_threshold = args.C
    k_size = args.k
    for batch_id, reads in enumerate(batch(all_reads, args.max_seqs)):
        print("correcting {0} reads in a batch".format(len(reads)))
        batch_start_time = time()
        if args.set_w_dynamically:
            # Activates for 'small' clusters with less than 700 reads
            if len(reads) >= 100:
                w = min(args.w, args.k + (len(reads)//100 + 4))  # min(args.w,)  #args.w = args.k + min(7, int( len(all_reads)/500))
            else:
                w = args.k + 1 + len(reads)//30  # min(args.w,)  #args.w = args.k + min(7, int( len(all_reads)/500))
        else:
            w = args.w
        print("Window used for batch:", w)
        x_high = args.xmax
        x_low = args.xmin
        hash_fcn = "lex"
        # for hash_fcn in ["lex"]: # ["lex"]: #  add "rev_lex" # add four others
        if args.randstrobes:
            if args.set_layers_manually:
                primes = get_primes(1000, args.layers)
            else:
                layers = 1 if len(reads) >= 1000 else 2
                print("Using {0} layers.".format(layers))
                primes = get_primes(1000, layers)
            # print(primes)
            # primes = [97,73]
            # minimizer_combinations_database, read_to_randstrobes = get_randstrobes_with_positions_database(reads, k_size, x_low, x_high, primes)
            minimizer_combinations_database = get_randstrobes_with_positions_database_2way(reads, k_size, x_low, x_high, primes)
        else:
            if args.compression:
                minimizer_database  = get_minimizers_and_positions_compressed(reads, w, k_size, hash_fcn)
            else:
                minimizer_database  = get_minimizers_and_positions(reads, w, k_size, hash_fcn)

            minimizer_combinations_database = get_minimizer_combinations_database(reads, minimizer_database, k_size, x_low, x_high)

        quality_values_database = get_qvs(reads)
        # print(minimizer_database)
        if args.verbose:
            eprint("done creating minimizer combinations")

        # print( [ (xx, len(reads_to_M2[xx])) for xx in reads_to_M2 ])
        # sys.exit()
        # corrected_reads = {}
        # tot_errors_before = {"subs" : 0, "del": 0, "ins": 0}
        # tot_errors_after = {"subs" : 0, "del": 0, "ins": 0}
        tot_corr = 0
        previously_corrected_regions = defaultdict(list)
        # stored_calculated_regions = defaultdict(lambda: defaultdict(int))
        tmp_cnt = 0

        for r_id in sorted(reads): #, reverse=True):
            # print()
            # print(reads[r_id][0])
            if args.randstrobes:
                seq = reads[r_id][1]
                # print("seq length:", len(seq))
                # read_min_comb = [ ((m1,p1), m1_curr_spans) for m1, p1, m1_curr_spans in randstrobe_iterator(seq, k_size, x_low, x_high, primes)] 
                read_min_comb = randstrobes_read_2way(seq, k_size, x_low, x_high, primes)
                # read_min_comb = read_to_randstrobes[r_id]
            else:
                # seq = reads[r_id][1]
                # print("seq length:", len(seq))
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
                if len(read_previously_considered_positions) > 1:
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

            # if r_id in interval_database:
            #     print(len(interval_database[r_id]))
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
                    tmp_cnt, read_complexity_cnt = find_most_supported_span(r_id, m1, p1, not_prev_corrected_spans, minimizer_combinations_database, reads, all_intervals, k_size, tmp_cnt, read_complexity_cnt, quality_values_database, already_computed)


            # from pympler import asizeof
            # print("reads", asizeof.asizeof(reads)/1000000)
            # print("not_prev_corrected_spans", asizeof.asizeof(not_prev_corrected_spans)/1000000)
            # # print("other_reads_corrected_regions", asizeof.asizeof(other_reads_corrected_regions)/1000000)
            # print("previously_corrected_regions", asizeof.asizeof(previously_corrected_regions)/1000000)
            # print("all_intervals", asizeof.asizeof(all_intervals)/1000000)
            # print("read_min_comb", asizeof.asizeof(read_min_comb)/1000000)
            # print("quality_values_database", asizeof.asizeof(quality_values_database)/1000000)
            # print("already_computed", asizeof.asizeof(already_computed)/1000000)   
            # print("minimizer_database", asizeof.asizeof(minimizer_database)/1000000)            
            # print("minimizer_combinations_database", asizeof.asizeof(minimizer_combinations_database)/1000000)
            # sleep(100)
            # sys.exit()

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
                all_intervals.sort(key = lambda x: (x[1],x[0]))
                # print(r_id)
                # print([(xxx, yyy, www, zzz) for (xxx, yyy, www, zzz)  in all_intervals])
                # print()
                # sys.exit()
                # if r_id == 1:
                #     all_intervals = [(3191, 3191, 3, all_intervals[0][3])]
                #     print("changed:", all_intervals )
                opt_indicies = solve_WIS(all_intervals) # solve Weighted Interval Scheduling here to find set of best non overlapping intervals to correct over
                # print(opt_indicies)
                # assert opt_indicies == opt_indicies2
                # print(opt_indicies)
                if not opt_indicies: # all intervals had only the read to be corrected itself as support (i.e., 1 in support)
                    corrected_seq = seq
                else:
                    intervals_to_correct = get_intervals_to_correct(opt_indicies[::-1], all_intervals)
                    del all_intervals
                    all_intervals = []
                    corrected_seq, other_reads_corrected_regions = correct_read(seq, reads, intervals_to_correct, k_size, work_dir, v_depth_ratio_threshold, max_seqs_to_spoa, args.disable_numpy, args.verbose, args.use_racon)
                    del intervals_to_correct
                    for other_r_id, corrected_regions in other_reads_corrected_regions.items():
                        for corr_region in corrected_regions:
                            previously_corrected_regions[other_r_id].append(corr_region)

            # from pympler import asizeof
            # print("reads", asizeof.asizeof(reads)/1000000)
            # # print("read_min_comb", asizeof.asizeof(read_min_comb)/1000000)
            # print("not_prev_corrected_spans", asizeof.asizeof(not_prev_corrected_spans)/1000000)
            # print("other_reads_corrected_regions", asizeof.asizeof(other_reads_corrected_regions)/1000000)
            # print("previously_corrected_regions", asizeof.asizeof(previously_corrected_regions)/1000000)
            # print("all_intervals", asizeof.asizeof(all_intervals)/1000000)
            # print("read_min_comb", asizeof.asizeof(read_min_comb)/1000000)
            # print("quality_values_database", asizeof.asizeof(quality_values_database)/1000000)
            # print("already_computed", asizeof.asizeof(already_computed)/1000000)   
            # print("minimizer_database", asizeof.asizeof(minimizer_database)/1000000)            
            # print("minimizer_combinations_database", asizeof.asizeof(minimizer_combinations_database)/1000000)

            corrected_reads[r_id] = (acc, corrected_seq, "+"*len(corrected_seq))
            if args.verbose:
                print("@{0}\n{1}\n+\n{2}".format(acc, corrected_seq, "+"*len(corrected_seq) ))
                eprint("{0},{1}".format(r_id,corrected_seq))
        print()
        print("Done with batch_id:", batch_id)
        print("Took {0} seconds.".format(time()- batch_start_time))
        # from pympler import asizeof
        # print("reads", asizeof.asizeof(reads)/1000000)
        # # print("read_min_comb", asizeof.asizeof(read_min_comb)/1000000)
        # print("not_prev_corrected_spans", asizeof.asizeof(not_prev_corrected_spans)/1000000)
        # print("other_reads_corrected_regions", asizeof.asizeof(other_reads_corrected_regions)/1000000)
        # print("previously_corrected_regions", asizeof.asizeof(previously_corrected_regions)/1000000)
        # print("all_intervals", asizeof.asizeof(all_intervals)/1000000)
        # print("read_min_comb", asizeof.asizeof(read_min_comb)/1000000)
        # print("quality_values_database", asizeof.asizeof(quality_values_database)/1000000)
        # print("already_computed", asizeof.asizeof(already_computed)/1000000)   
        # # print("minimizer_database", asizeof.asizeof(minimizer_database)/1000000)            
        # print("minimizer_combinations_database", asizeof.asizeof(minimizer_combinations_database)/1000000)
        # minimizer_combinations_database.clear()
        # quality_values_database.clear()
        # previously_corrected_regions.clear()
                # eval_sim2(corrected_seq, seq, qual, tot_errors_before, tot_errors_after)
                # if r_id == 10:
                #     sys.exit()

    # eprint("tot_before:", tot_errors_before)
    # eprint("tot_after:", sum(tot_errors_after.values()), tot_errors_after)
    eprint( len(corrected_reads))
    outfile = open(os.path.join(args.outfolder, "corrected_reads.fastq"), "w")

    for r_id, (acc, seq, qual) in corrected_reads.items():
        outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))
    outfile.close()

    print("removing temporary workdir")
    shutil.rmtree(work_dir)


def main():
    parser = argparse.ArgumentParser(description="De novo error correction of long-read transcriptome reads", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.3.4')

    parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq file with reads')
    # parser.add_argument('--t', dest="nr_cores", type=int, default=8, help='Number of cores allocated for clustering')

    parser.add_argument('--k', type=int, default=9, help='Kmer size')
    parser.add_argument('--w', type=int, default=20, help='Window size')
    parser.add_argument('--xmin', type=int, default=18, help='Lower interval length')
    parser.add_argument('--xmax', type=int, default=80, help='Upper interval length')
    parser.add_argument('--T', type=float, default=0.1, help='Minimum fraction keeping substitution')
    # parser.add_argument('--C', type=float, default=0.05, help='Minimum fraction of keeping alternative refernece contexts')
    parser.add_argument('--exact', action="store_true", help='Get exact solution for WIS for evary read (recalculating weights for each read (much slower but slightly more accuracy,\
                                                                 not to be used for clusters with over ~500 reads)')
    parser.add_argument('--disable_numpy', action="store_true", help='Do not require numpy to be installed, but this version is about 1.5x slower than with numpy.')

    parser.add_argument('--max_seqs_to_spoa', type=int, default=200,  help='Maximum number of seqs to spoa')
    parser.add_argument('--max_seqs', type=int, default=2000,  help='Maximum number of seqs to correct at a time (in case of large clusters).')
    parser.add_argument('--use_racon', action="store_true", help='Use racon to polish consensus after spoa (more time consuming but higher accuracy).')

    parser.add_argument('--exact_instance_limit', type=int, default=0,  help='Activates slower exact mode for instance smaller than this limit')
    # parser.add_argument('--w_equal_k_limit', type=int, default=0,  help='Sets w=k which is slower and more memory consuming but more accurate and useful for smalled clusters.')
    parser.add_argument('--set_w_dynamically', action="store_true", help='Set w = k + max(2*k, floor(cluster_size/1000)).')
    parser.add_argument('--verbose', action="store_true", help='Print various developer stats.')
    parser.add_argument('--randstrobes', action="store_true", help='EXPERIMENTAL PARAMETER: IsONcorrect uses paired minimizers (described in isONcorrect paper). This experimental option\
                                                                 uses randstrobes instead of paired minimizers to find shared regions. Randstrobes \
                                                                 reduces memory footprint substantially (and runtime) with only slight increase in post correction quality.')


    parser.add_argument('--layers', type=int, default=argparse.SUPPRESS, help='EXPERIMENTAL PARAMETER: Active when --randstrobes specified.\
                                                                How many "layers" with randstrobes we want per sequence to sample.\
                                                               More layers gives more accureate results but is more memory consuming and slower.\
                                                               It is not reccomended to specify more than 5. ')
    parser.add_argument('--set_layers_manually', action="store_true", help='EXPERIMENTAL PARAMETER: By default isONcorrect sets layers = 1 if nr seqs in batch to be corrected is >= 1000, else layers = 2.\
                                                                            This command will manually pick the number of layers specified with the --layers parameter.')


    parser.add_argument('--compression', action="store_true", help='Use homopolymenr compressed reads. (Deprecated, because we will have fewer \
                                                                        minmimizer combinations to span regions in homopolymenr dense regions. Solution \
                                                                        could be to adjust upper interval legnth dynamically to guarantee a certain number of spanning intervals.')
    parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='isoncorrect_main')
    args = parser.parse_args()


    if args.xmin < 2*args.k:
        args.xmin = 2*args.k
        eprint("xmin set to {0}".format(args.xmin))

    # if not args.paired_minimizers and 'max_seqs' not in args:
    #     print("max_seqs was not specified and paired_minimizer setting not used. Setting max_seqs to 2000")
    #     args.max_seqs = 2000
    # elif args.paired_minimizers and 'max_seqs' not in args:
    #     print("max_seqs was not specified. Setting max_seqs to 1000")
    #     args.max_seqs = 1000

    if args.set_layers_manually and 'layers' not in args:
        args.layers = 2

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if not args.fastq and not args.flnc and not  args.ccs:
        parser.print_help()
        sys.exit()




    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    # edlib_module = 'edlib'
    # parasail_module = 'parasail'
    # if edlib_module not in sys.modules:
    #     print('You have not imported the {0} module. Only performing clustering with mapping, i.e., no alignment.'.format(edlib_module))
    # if parasail_module not in sys.modules:
    #     eprint('You have not imported the {0} module. Only performing clustering with mapping, i.e., no alignment!'.format(parasail_module))
    #     sys.exit(1)
    if 100 < args.w or args.w < args.k:
        eprint('Please specify a window of size larger or equal to k, and smaller than 100.')
        sys.exit(1)

    isoncorrect_main(args)


if __name__ == '__main__':
    main()


