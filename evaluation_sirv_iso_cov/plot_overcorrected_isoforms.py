
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


def sirv_overcorrected_isoforms(input_csv, outfolder):
    # get_error_rates(input_csv)
    # get_error_rates(input_csv, read_to_infer = 'original' )
    # sys.exit()
    pd.set_option("display.precision", 8)
    indata = pd.read_csv(input_csv)
    indata["nr_isoforms"] = ["$%s$" % x for x in indata["nr_isoforms"]]

    ax = sns.lineplot(x="nr_reads", y="is_overcorrected",  hue="nr_isoforms",  data=indata)
    ax.set_ylim(0,0.1)
    ax.set_xscale('log')
    ax.set_xticks( [1,2,3,4,5,10,20]) #([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
    x_labels = ["1","","","","5","10","20"]
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Fraction overcorreced reads")
    ax.set_xlabel("Reads per transcript")

    plt.savefig(os.path.join(outfolder, "sirv_iso_overcorrected.eps"))
    plt.savefig(os.path.join(outfolder, "sirv_iso_overcorrected.pdf"))
    plt.close()



def main(args):
    sns.set(style="whitegrid")

    sirv_overcorrected_isoforms(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sirv subsampling reads.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

