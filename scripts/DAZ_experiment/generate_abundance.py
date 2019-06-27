
from __future__ import print_function
import os,sys
import argparse
import re
import misc_functions
from collections import deque
import random
import itertools
# import math


def generate_abundance(isoforms, config, params):
    sequence_pool = {}
    params.logfile.write("ABUNDANCES\n")

    for i, (acc, seq) in enumerate(misc_functions.iteritems(isoforms)):
        # abundance = random.choice(config["abundance"])
        abundance = config["abundance"][ i % len(config["abundance"]) ]
        params.logfile.write("{0}, abundance: {1}\n".format(acc, abundance))

        for i in range(1, abundance +1):
            new_acc = acc + "_" + str(i)
            sequence_pool[new_acc] = seq
    return sequence_pool

def set_parameters(nr_transcripts):
    config_parameters = {}
    type_dict = {}
    isoform_counts = []

    if params.abundance == "constant":
        config_parameters["abundance"] = [1 for i in range(nr_transcripts) ]
    elif params.abundance == "exponential":
        config_parameters["abundance"] = [2**i for i in range(8) ]

    return config_parameters


def main(params):
    isoforms = {acc : seq for acc, seq in misc_functions.read_fasta(open(params.transcript_file,"r"))} 
    config = set_parameters(len(isoforms))
    all_transcripts_for_sequencing = generate_abundance(isoforms, config, params)
    out_file = open(params.outfile, "w")
    for acc, seq in misc_functions.iteritems(all_transcripts_for_sequencing):
        out_file.write(">{0}\n{1}\n".format(acc,seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a gene family and isoforms from a set of original exons.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--transcript_file', type=str, help='The fasta file with a full transcript family.')
    parser.add_argument('--abundance', type=str, help='either constant or exponential.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')

    params = parser.parse_args()
    path_, file_prefix = os.path.split(params.outfile)
    misc_functions.mkdir_p(path_)
    params.logfile = open(os.path.join(path_, file_prefix[:-3] + ".log"), "w")
    main(params)