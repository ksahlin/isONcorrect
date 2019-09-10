
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


def total_error_rate(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))

    g = sns.catplot(x="read_type", y="error_rate", #col="Depth",
                data=indata,  #hue="read_type", hue_order= ["corrected", "original"],
                kind="violin", aspect=1)

    g.set(ylim=(0,15))
    g.set_ylabels("Error rate (%)")
    g.set_xlabels("Method")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "total_error_rate.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate.pdf"))
    plt.close()

def splice_site_classification(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))

    g = sns.catplot(x="read_type", y="error_rate", #col="Depth",
                data=indata,  #hue="read_type", hue_order= ["corrected", "original"],
                kind="violin", aspect=1)

    g.set(ylim=(0,15))
    g.set_ylabels("Error rate (%)")
    g.set_xlabels("Method")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "total_error_rate.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate.pdf"))
    plt.close()



def total_error_rate2(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))
    data =indata[indata.read_type == 'corrected']
    sns.distplot(data['error_rate'], kde=True, label='Corrected', bins=200, hist_kws=dict(alpha=0.5))
    data = indata[indata.read_type == 'original']
    sns.distplot(data['error_rate'], kde=True, label='Original', bins=200, hist_kws=dict(alpha=0.5))


    plt.xticks(np.arange(0, 20, step=1))
    plt.xlim(0,20)
    plt.xlabel("Difference to HG38 (%)")
    plt.ylabel("Frequency")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "total_error_rate2.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate2.pdf"))
    plt.close()

def main(args):
    
    sns.set(style="whitegrid")

    # total_error_rate(args.input_csv, args.outfolder)
    total_error_rate2(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

