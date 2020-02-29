from __future__ import print_function
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout

# import networkx as nx

from collections import defaultdict
# import numpy as np

# Botonds DP
def cutoff(x):
    return x

def longest_path_botond(dot_graph_path):
    G = nx.drawing.nx_pydot.read_dot(dot_graph_path)   
    topo = list(nx.topological_sort(G))   
    data = nx.get_node_attributes(G, 'label')
    # data = {x: chr(int(y.split(" - ")[1].split("\"")[0])) for x, y in data.items()}
    data = {x: y.split(" - ")[1].split("\"")[0] for x, y in data.items()}
    scores = defaultdict(int)
    bt = defaultdict(lambda: -1)
     
    for n in topo:
        sc = []
        for p in G.predecessors(n):
            tmp = G.get_edge_data(p, n, key=0)
            if 'label' in tmp:
                weight = int(tmp['label'].strip("\""))
                sc.append((p, scores[p] + cutoff(weight)))
                # sc.append((p, scores[p] + weight))
        sc = sorted(sc, key=lambda x: x[1], reverse=True)
        if len(sc) < 1:
            continue
        bt[n] = sc[0][0]
        scores[n] = sc[0][1]
     
    max_score = max(scores.values())
    lookup = {y: x for x, y in scores.items()}
     
    prev = lookup[max_score]
    nodes = []
     
    while prev != -1:
        nodes.append(prev)
        prev = bt[prev]
     
    nodes = list(reversed(nodes))
     
    print(">cons_" + str(len(nodes)))
     
    bases = [data[n] for n in nodes]
    cons = ''.join(bases)
    print(cons)
    return cons


# def longest_path(dot_graph_path):
#     G = nx.DiGraph()
#     for line in open(dot_graph_path, "r"):
#         if len(line.split("->")) == 2: # edge
#             n1, info = line.split("->")
#             n1 = n1.strip()
#             n2, info = info.strip().split("[")
#             n2 = n2.strip()
#             if "dotted" in info: # edge should not be created for substitutions 
#                 continue  
#             elif "," not in info: # not in consensus path
#                 weight = info.split("=")[1].strip().strip("]").strip('"')
#                 # print("ok", weight)
#             else:
#                 weight = info.split(",")[0].split("=")[1].strip().strip('"')
#                 # print("bla", weight)
#             if int(weight) > 4:
#                 # print(weight)
#                 G.add_edge(n1, n2, weight=int(weight))

#         elif "label" in line: # at node
#             node_id, info = line.split("[")
#             node_id = node_id.strip()
#             nucl = info.split(" - ")[1][0]
#             # print(node_id, nucl)
#             G.add_node(node_id, nucleotide=nucl)

#     # G = nx.drawing.nx_pydot.read_dot(dot_graph_path)
#     # print(G.edges(data=True))
#     # print(G.nodes(data=True))
#     order = list(nx.algorithms.dag.topological_sort(G))
#     # print([ [G[n][nbr]["weight"] for nbr in G.neighbors(n)] for n in order ])
#     # print(order)
#     longest_path = nx.algorithms.dag.dag_longest_path(G)
#     # print("longest path", [(i,n) for i,n in enumerate(longest_path)])
#     augmented_ref = "".join([ G.node[n]["nucleotide"] for n in longest_path])
#     # augmented_ref = "".join([ info["label"][-2] for n, info in G.nodes(data=True)])
#     print(augmented_ref)
#     print(len(augmented_ref))
#     return augmented_ref


def run_spoa_affine(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call([ spoa_path, "-q", reads, "-l", "2", "-r", "2", "-x", "-8", "-m", "10", "-o","-8", "-e", "-1", "-d", dot_graph_path], stderr=output_file, stdout=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print("affine heaviest path:", consensus, len(consensus))
    print()
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()
    return consensus, msa

def run_spoa_convex(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa convex...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        # print(spoa_path, "-l", "2", "-r", "2", "-g", "-4", "-e", "0", "-d", dot_graph_path, reads)
        subprocess.check_call([ spoa_path, "-l", "2", "-r", "2", "-m", "10", "-g", "-8", "-e", "-2", "-q", "-24", "-c", "-1" , "-d", dot_graph_path, reads], stdout=output_file, stderr=null)
        # subprocess.check_call([ spoa_path, "-q", reads, "-l", "2", "-r", "2", "-o","-4", "-e", "0", "-d", dot_graph_path], stderr=output_file, stdout=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print("convex heaviest path:", consensus, len(consensus))
    print()
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()
    return consensus, msa

def run_spoa_affine_v2_0_3(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        # print(spoa_path, "-l", "2", "-r", "2", "-g", "-4", "-e", "0", "-d", dot_graph_path, reads)
        subprocess.check_call([ spoa_path, "-l", "2", "-r", "2", "-x", "-4", "-m", "10", "-g", "-8", "-e", "-1", "-d", dot_graph_path, reads], stdout=output_file, stderr=null)
        # subprocess.check_call([ spoa_path, "-q", reads, "-l", "2", "-r", "2", "-o","-4", "-e", "0", "-d", dot_graph_path], stderr=output_file, stdout=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print("affine_v2_0_2 heaviest path:", consensus, len(consensus))
    print()
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()
    return consensus, msa

def run_spoa_with_msa(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        # subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "2", "-g", "-2", "-d", dot_graph_path], stdout=output_file, stderr=null)
        subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "2", "-g", "-2"], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    # print("regular spoa:", consensus)
    # print(len(consensus))
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()

    return consensus, msa

def run_spoa(reads, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "0", "-g", "-2"], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    # output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    del l
    return consensus


# For eventual De Bruijn graph approach
import itertools
from collections import defaultdict, deque
def kmer_counter(reads, k_size):
    count = defaultdict(int)
    position_count = defaultdict(list)

    for r_i in reads:
        (acc, seq, qual) = reads[r_i]
        # seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        read_kmers = deque([seq[i:i+k_size] for i in range(len(seq) - k_size )])
        for i, kmer in enumerate(read_kmers):
            count[kmer] += 1

            position_count[kmer].append( (r_i, i))

    # cnt_sorted = sorted(count.items(), key = lambda x: x[1], reverse=True)
    # print(len(cnt_sorted),cnt_sorted)
    return count, position_count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('reads', type=str, help='Fasta or fastq')
    parser.add_argument('outfile', type=str, help='Fasta or fastq')
    parser.add_argument('--spoa_path', type=str, default='spoa', required=False, help='Path to spoa binary with bwa binary name at the end.')


    args = parser.parse_args()

    run_spoa(args.reads, args.outfile, args.spoa_path)
