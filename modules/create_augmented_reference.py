from __future__ import print_function
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout

import networkx as nx

def longest_path(dot_graph_path):
    G = nx.DiGraph()
    for line in open(dot_graph_path, "r"):
        if len(line.split("->")) == 2: # edge
            n1, info = line.split("->")
            n1 = n1.strip()
            n2, info = info.strip().split("[")
            n2 = n2.strip()
            if "dotted" in info: # edge should not be created for substitutions 
                continue  
            elif "," not in info: # not in consensus path
                weight = info.split("=")[1].strip().strip("]").strip('"')
                # print("ok", weight)
            else:
                weight = info.split(",")[0].split("=")[1].strip().strip('"')
                # print("bla", weight)
            if int(weight) > 4:
                # print(weight)
                G.add_edge(n1, n2, weight=int(weight))

        elif "label" in line: # at node
            node_id, info = line.split("[")
            node_id = node_id.strip()
            nucl = info.split(" - ")[1][0]
            # print(node_id, nucl)
            G.add_node(node_id, nucleotide=nucl)

    # G = nx.drawing.nx_pydot.read_dot(dot_graph_path)
    # print(G.edges(data=True))
    # print(G.nodes(data=True))
    order = list(nx.algorithms.dag.topological_sort(G))
    # print([ [G[n][nbr]["weight"] for nbr in G.neighbors(n)] for n in order ])
    print(order)
    longest_path = nx.algorithms.dag.dag_longest_path(G)
    print("longest path", [(i,n) for i,n in enumerate(longest_path)])
    augmented_ref = "".join([ G.node[n]["nucleotide"] for n in longest_path])
    # augmented_ref = "".join([ info["label"][-2] for n, info in G.nodes(data=True)])
    print(augmented_ref)
    print(len(augmented_ref))
    # sys.exit()
    return augmented_ref


def run_spoa_affine(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call([ spoa_path, "-q", reads, "-l", "0", "-r", "2", "-o","-4", "-e", "0", "-d", dot_graph_path], stderr=output_file, stdout=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print("affine:", consensus)
    print(len(consensus))
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()

    return consensus, msa

def run_spoa(reads, ref_out_file, spoa_out_file, spoa_path, dot_graph_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "2", "--gap", "-2", "-d", dot_graph_path], stdout=output_file, stderr=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print("regular spoa:", consensus)
    print(len(consensus))
    # print(msa)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()

    return consensus, msa

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('reads', type=str, help='Fasta or fastq')
    parser.add_argument('outfile', type=str, help='Fasta or fastq')
    parser.add_argument('--spoa_path', type=str, default='spoa', required=False, help='Path to spoa binary with bwa binary name at the end.')


    args = parser.parse_args()

    run_spoa(args.reads, args.outfile, args.spoa_path)
