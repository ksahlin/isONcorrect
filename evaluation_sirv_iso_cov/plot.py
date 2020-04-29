
## Various plots from large table

import sys
import argparse
import os
import random
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
from matplotlib import pyplot


def sirv_error_rate_per_transcript(input_csv, outfolder):
    # get_error_rates(input_csv)
    # get_error_rates(input_csv, read_to_infer = 'original' )
    # sys.exit()
    pd.set_option("display.precision", 8)
    indata = pd.read_csv(input_csv)

    # new controlled experiment:
    g = sns.FacetGrid(indata, col="nr_isoforms", hue="read_type", col_wrap=3)
    g.map(sns.lineplot, "nr_reads","error_rate",ci = 'sd', estimator= 'median').add_legend()
    g.set(ylim=(0, 12))
    ticks = [1,2,3,4,5,10,20]
    labels = ["1","","","","5","10","20"]
    g.set_axis_labels("Reads per transcript", "Error rate (%)")
    g.set(xscale ='log')
    g.set(xticks = ticks,xticklabels = labels)

    ticks = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    labels = ["","1", "","3","","5","","7","","9","","11",""]   
    g.set(yticks = ticks, yticklabels = labels)
    

    # ax = sns.lineplot(x="nr_reads", y="error_rate",  hue="read_type",
    #                   ci = None, estimator= 'median',  data=indata)
    # ax.set_ylim(0,12)
    # ax.set_xscale('log')
    # ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
    # ax.set_ylabel("Error rate (%)")
    # ax.set_xlabel("Reads per transcript")

    plt.savefig(os.path.join(outfolder, "sirv_iso_cov.eps"))
    plt.savefig(os.path.join(outfolder, "sirv_iso_cov.pdf"))
    plt.close()



def main(args):
    sns.set(style="whitegrid")

    sirv_error_rate_per_transcript(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sirv subsampling reads.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

