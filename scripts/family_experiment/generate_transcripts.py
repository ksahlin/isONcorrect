
from __future__ import print_function
import os,sys
import argparse
import re
import misc_functions
from collections import deque
import random
import itertools
# import math



def mutate_member(exons, config, params, mut = True):
    leaf_node = {}
    total_mutations = 0
    total_length = 0
    log_string = ""
    for ex_nr, seq in misc_functions.iteritems(exons):
        if mut:
            new_seq, mutation_log, exon_mutations, exon_indels, total_del_length = misc_functions.mutate_sequence(seq, config["mut"], config["ins"], config["del"])
            # print(new_seq, mutation_log, exon_mutations)
            leaf_node[ex_nr] = new_seq 
        else:
            new_seq = seq
            leaf_node[ex_nr] = seq 
            mutation_log, exon_mutations = "-", 0
        log_string += "mutations in exon: {0}, mutation places: {1}\n".format(exon_mutations, mutation_log)
        # params.logfile.write("mutations in exon: {0}, mutation places: {1}\n".format(exon_mutations, mutation_log))
        total_mutations += exon_mutations
        total_length += len(new_seq)

    if total_mutations > 0 and mut:
        params.logfile.write(log_string)
        params.logfile.write("Total mutations: {0}\n".format(total_mutations))
        params.logfile.write("mutation rate all exons: {0}\n".format(total_mutations/ float(total_length)))
        return leaf_node
    elif not mut:
        params.logfile.write(log_string)
        params.logfile.write("Total mutations: {0}\n".format(total_mutations))
        params.logfile.write("mutation rate all exons: {0}\n".format(total_mutations/ float(total_length)))
        return leaf_node        
    else:
        print("NO mutations!")
        return False

def create_family(exons, config):
    tree_evolution = deque() # implemented as a FIFO queue

    variants_in_family = {}  # variants labelled as interger number. each varinat is a list of (mutated) exon strings
    for acc, seq in exons:
        ex_nr = int(acc.split("|")[0])
        config["original_family"] = acc.split("|")[1]
        variants_in_family[ex_nr] = seq
    tree_evolution.append({0 : variants_in_family}) 

    node_counter = 0
    while len(tree_evolution) < config["family_size"]:
        node = tree_evolution.popleft()
        print("popping node", list(node.keys()))
        exons = list(node.values())[0]

        node_counter += 1
        params.logfile.write("Node_" + str(node_counter) + "\n")
        while True:
            left_leaf_node = mutate_member(exons, config, params, mut=False)
            if left_leaf_node:
                tree_evolution.append({node_counter : left_leaf_node}) 
                break

        node_counter += 1
        params.logfile.write("Node_" + str(node_counter) + "\n")
        while True:
            right_leaf_node = mutate_member(exons, config, params)
            if right_leaf_node:
                tree_evolution.append({node_counter : right_leaf_node}) 
                break


    # print(tree_evolution)
    print(len(tree_evolution))
    return tree_evolution

def get_isoforms(gene_copy, nr_isoforms, node, config, prev_simulated_isoforms):
        ### new drop method #####
        isoforms = {}
        nr_simulated = 0
        while True:
            for i in range(len(gene_copy), 0, -1):
                for isoform in itertools.combinations(list(gene_copy.keys()), i):
                    print(isoform)
                    isoform_copy_accession = config["original_family"] + ";" + "member:{0};exons:".format(list(node.keys())[0]) + ",".join([str(ex_nr) for ex_nr in isoform])
                    isoform_copy_sequence = "".join([gene_copy[ex_nr] for ex_nr in isoform])
                    if isoform_copy_sequence in prev_simulated_isoforms:
                        continue
                    else:
                        isoforms["{0}".format(isoform_copy_accession)] = isoform_copy_sequence
                        nr_simulated += 1
                        prev_simulated_isoforms.add(isoform_copy_sequence)

                    if nr_simulated >= nr_isoforms:
                        return isoforms

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def drop_exons(tree_evolution, config):
    final_isoforms = {}
    isoform_counts = config["isoform_counts"]
    current_isoforms_seqs = set() 
    for i, node in enumerate( tree_evolution ):
        current_isoforms_seqs = set([seq for seq in final_isoforms.values()])
        # print(list(node.keys())[0])
        gene_copy = list(node.values())[0]
        # print(len(gene_copy), gene_copy)
        already_generated_isoforms = set()
        nr_isoforms = isoform_counts[ i % len(isoform_counts) ]  # random.choice(isoform_counts) #
        print("nr_isoforms", nr_isoforms)

        assert nr_isoforms < 2**(len(gene_copy))
        isoforms_ = get_isoforms(gene_copy, nr_isoforms, node, config, current_isoforms_seqs)

        final_isoforms = merge_two_dicts(final_isoforms, isoforms_)


        # print(isoforms)

        # otherwise we will try forever...
        # print(nr_isoforms, 2**(len(gene_copy)-1) )
        # assert nr_isoforms < 2**(len(gene_copy)-1)

    #     for i in range(nr_isoforms):
    #         sample_again = True
    #         while sample_again:
    #             nr_exons = random.randint(2, len(gene_copy))
    #             isoform = random.sample(list(gene_copy.keys()), nr_exons)
    #             # for acc, seq in gene_copy.iteritems():
    #             # print(isoform.sort())
    #             isoform.sort()
    #             # print(isoform)
    #             current_isoform = "".join([str(i) for i in isoform])
    #             if current_isoform not in already_generated_isoforms:
    #                 already_generated_isoforms.add(current_isoform)
    #                 isoform_copy_accession = config["original_family"] + ";" + "member:{0};exons:".format(list(node.keys())[0]) + ",".join([str(ex_nr) for ex_nr in isoform])
    #                 # print(isoform_copy_accession)
    #                 isoform_copy_sequence = "".join([gene_copy[ex_nr] for ex_nr in isoform])
    #                 # print(isoform_copy_sequence)
    #                 isoforms[">{0}".format(isoform_copy_accession)] = isoform_copy_sequence
    #                 sample_again = False
    #             else:
    #                 pass
    #                 print(current_isoform, "already sampled")
    # print(isoforms)

    # print(len(isoforms2), len(isoforms))

    return final_isoforms


def set_parameters(params):
    config_parameters = {}
    type_dict = {}
    isoform_counts = []

    if params.isoform_distribution == "exponential":
        type_dict["exponential"] = [2**i for i in range(1, params.family_size +1) ]
        config_parameters["isoform_counts"] = type_dict[params.isoform_distribution]
    elif  params.isoform_distribution == "constant":
        type_dict["constant"] =[params.nr_isoforms for i in range( params.family_size) ]
        config_parameters["isoform_counts"] = type_dict["constant"]

    if params.mutation_rate:
        config_parameters["mut"] = params.mutation_rate / float(3)
        config_parameters["ins"] = params.mutation_rate / float(3)
        config_parameters["del"] = params.mutation_rate / float(3)
    config_parameters["original_family"] = params.gene_member
    config_parameters["family_size"] = params.family_size
    return config_parameters


def main(params):
    # config = misc_functions.read_config(params.config)
    config = set_parameters(params)
    print(params.exon_file, params.transcript_file)
    if params.exon_file:
        exons = misc_functions.read_fasta(open(params.exon_file,"r"))
        tree_evolution = create_family(exons, config)
        isoforms = drop_exons(tree_evolution, config)
    elif params.transcript_file:
        isoforms = {acc : seq for acc, seq in misc_functions.read_fasta(open(params.transcript_file,"r"))} 
    out_file = open(params.outfile, "w")
    for acc, seq in misc_functions.iteritems(isoforms):
        out_file.write(">{0}\n{1}\n".format(acc,seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a gene family and isoforms from a set of original exons.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--exon_file', type=str, help='The fasta file with original exons.')
    group.add_argument('--transcript_file', type=str, help='The fasta file with a full transcript family.')
    # parser.add_argument('exon_file', type=str, help='The fasta file with original exons.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')
    # parser.add_argument('config', type=str, help='config file')
    parser.add_argument('--gene_member', type=str, help='gene member')
    parser.add_argument('--family_size', type=int, help='Family size')
    parser.add_argument('--isoform_distribution', type=str, help='isoform distribution, either constant or exponential')
    parser.add_argument('--nr_isoforms', type=int, default = 1, help='Family size')
    parser.add_argument('--mutation_rate', type=float, help='Mutation rate')

    params = parser.parse_args()
    path_, file_prefix = os.path.split(params.outfile)
    misc_functions.mkdir_p(path_)
    params.logfile = open(os.path.join(path_, file_prefix[:-3] + ".log"), "w")
    main(params)