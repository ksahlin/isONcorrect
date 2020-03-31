
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
    ax = sns.lineplot(x="nr_reads", y="error_rate",  hue="read_type",
                      ci = 'sd', estimator= 'median',  data=indata)
    ax.set_ylim(0,12)
    ax.set_xscale('log')
    ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
    ax.set_ylabel("Error rate (%)")
    ax.set_xlabel("Reads per transcript")

    plt.savefig(os.path.join(outfolder, "sirv_subsampling_error_rates.eps"))
    plt.savefig(os.path.join(outfolder, "sirv_subsampling_error_rates.pdf"))
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

