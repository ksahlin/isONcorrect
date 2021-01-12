from __future__ import print_function
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import numpy as np
import seaborn as sns
import pandas as pd


def plot_err_per_mappability(file, args):

    sns.set(style="whitegrid")
    data = pd.read_csv(file)

    bins = [1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000, 10000]
    labels = [ '{0}-{1}'.format(i,j) for i,j in zip(bins[:-1], bins[1:]) ]
    print(bins)
    print(labels)
    mappability_bins = pd.cut(data['mappability'], bins, labels = labels )
    data["mappability_bins"] = pd.Series(mappability_bins, index=data.index)
    print(mappability_bins)

    # error rateh per mappability bins
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=10)
    fig, ax = plt.subplots()
    ax = sns.boxplot(x="mappability_bins", y="err_rate", data=data, showfliers=False)
    plt.xlabel('Average k-mer copy number in transcript (binned)', fontsize=11)
    plt.ylabel('Error rate after correction',fontsize=14)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.tight_layout()
    # plt.ylim(1, 100)
    plt.savefig( args.outfile, dpi = 350)
    # plt.clf()
    plt.close()

    # y=sys.argv[3]
    # g = sns.catplot(x="Depth", y=y, col="w", row = "k",
    #             data=df, hue="type", hue_order= ["exact", "approx", "original"],
    #             kind="box", aspect=1)

    # g.set(ylim=(0,12))
    # g.set_ylabels("Error rate (%)")
    # g.set_xlabels("Read depth")
    # g.set_titles("{row_var} = {row_name}, {col_var} = {row_var} + {col_name}")
    # plt.savefig(sys.argv[2])
    # plt.close()

def main(args):
    reads = {}
    for line in open(args.results_file, 'r'):
        vals = line.strip().split(',')
        acc, err_rate, read_type = vals[0], vals[7], vals[8]
        if read_type == 'corrected':
            reads[acc] = [err_rate]

    transcripts = {}
    for line in open(args.mappability_file, 'r'):
        acc, mappability = line.strip().split(',')
        transcripts[acc] = mappability

    # example read_header ID: SYNE1|ENSG00000131018|ENST00000367255_1
    for read_acc in list(reads.keys()):
        tr_id = read_acc.split("_")[0]
        mappability = transcripts[tr_id]
        reads[read_acc].append(mappability)


    f = open(args.outfile[:-3] + "csv", "w")
    f.write("read_acc,err_rate,mappability\n")
    for read_acc, (err_rate, mappability) in reads.items(): 
    	f.write("{0},{1},{2}\n".format(read_acc, err_rate, mappability))
    f.close()

    plot_err_per_mappability(f.name, args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('results_file', type=str, help='csv')
    parser.add_argument('mappability_file', type=str, help='csv')
    parser.add_argument('outfile', type=str, help='pdf')


    args = parser.parse_args()

    main(args)
