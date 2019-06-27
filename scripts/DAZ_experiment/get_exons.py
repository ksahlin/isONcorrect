from __future__ import print_function
import os,sys
import argparse
import re
from misc_functions import read_fasta
from collections import defaultdict

def main(params):
    targeted = set(["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    exons = read_fasta(open(params.database, 'r'))
    # create output files here..
    out_files = defaultdict(str) #{}
    # for fam in targeted:
    #     out_files[fam] = open(os.path.join(params.outfolder, fam + "_exons.fa"), 'w')

    exon_ordering = defaultdict(str) #{}
    # for fam in targeted:
    #     exon_ordering[fam] = {}

    #store exons
    for acc, seq in exons:
        gene_name, chr_name, exon_rank = acc.split("|")
        try:
            exon_ranks = [int(exon_rank)]
        except:
            exon_ranks = [int(rank) for rank in exon_rank.split(";")]
        for family in targeted:
            if re.search(family, gene_name):
                for rank in exon_ranks:
                    if exon_ordering[gene_name]:
                        print("lol")
                        # print(gene_name,rank, seq)
                        if rank in exon_ordering[gene_name]:
                            print(rank, "already exists", gene_name)
                            if seq != exon_ordering[gene_name][rank]:
                                print("also, its a different sequence:")
                                print(seq)
                                print(exon_ordering[gene_name][rank])
                                # sys.exit()
                        exon_ordering[gene_name][rank] = seq
                    else:
                        print("hmm")
                        exon_ordering[gene_name] = {}
                        exon_ordering[gene_name][rank] = seq
                        out_files[gene_name] = open(os.path.join(params.outfolder, gene_name + "_exons.fa"), 'w')
                        # print(gene_name,rank, seq)

                    # print(exon_ordering[gene_name])



    for gene_name in exon_ordering:
        gene_length = 0
        tot_exons = 0
        for rank in sorted(exon_ordering[gene_name].keys()):
            seq = exon_ordering[gene_name][rank]
            out_files[gene_name].write(">{0}\n{1}\n".format(str(rank) + "|" + gene_name , seq))
            gene_length += len(seq)
            tot_exons += 1
        print(gene_name, "nr exons: ", tot_exons, "tot length: ", gene_length)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a gene family and isoforms from a set of original exons.")
    parser.add_argument('database', type=str, help='Exon database.')
    parser.add_argument('outfolder', type=str, help='Output path to fasta file')


    params = parser.parse_args()

    outfolder = params.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(params)