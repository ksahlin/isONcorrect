
## Various plots from large table
from time import time
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


# def get_accuracies(f):
#     #  ~/Documents/data/ont/drosophila/more_dros/dros_acc/A_*   # 7 datasets
#     dataset = {}
#     for i, line in enumarate(f):
#         if i == 0:
#             continue
#         else:
#             read_id, accuracy = line.split()

#     return dataset


def plot_histogram(input_csv, d_index, args):
    # from https://stackoverflow.com/questions/51417483/mean-median-mode-lines-showing-only-in-last-graph-in-seaborn

    # df = pd.DataFrame({"Acc": [1, 2, 3, 4, 6, 7, 9, 9, 9, 10], "dummy": range(10)})
    df = pd.read_csv(input_csv, sep='\t')

    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw= {"height_ratios": (0.2, 1)})
    mean=df['Acc'].mean()
    median=df['Acc'].median()
    mode=df['Acc'].mode().get_values()[0]
    q1=df['Acc'].quantile(0.25)
    q3=df['Acc'].quantile(0.75)


    sns.boxplot(df["Acc"], ax=ax_box)
    ax_box.axvline(mean, color='r', linestyle='--')
    ax_box.axvline(median, color='g', linestyle='-')
    ax_box.axvline(mode, color='b', linestyle='-')

    sns.distplot(df["Acc"], ax=ax_hist)
    ax_hist.axvline(mean, color='r', linestyle='--')
    ax_hist.axvline(median, color='g', linestyle='-')
    ax_hist.axvline(mode, color='b', linestyle='-')

    plt.legend({'Mean':mean,'Median':median,'Mode':mode})

    ax_box.set(xlabel='')
    plt.savefig(os.path.join(args.outfolder, str(d_index)+".pdf"))
    plt.clf()
    return mean, median, mode, df['Acc'].values.tolist(), q1, q3

def box_plot(indata):
    boxplot = indata.boxplot(grid=False, rot=30, fontsize=12)
    plt.tight_layout()
    boxplot.set_ylabel("Accuracy (%)")
    plt.savefig(os.path.join(args.outfolder, "boxplot.pdf"))


def main(args):
    sns.set(style="whitegrid")
    dfs = {}
    for i, file in enumerate(args.input_csvs):
        basename = os.path.basename(file)
        run_id = basename.split("_")[2]
        mean, median, mode, accuracies, q1, q3 = plot_histogram(file, basename, args)
        dfs[run_id] = accuracies 
        print(q1, median, q3, mean, mode, basename)

    # start = time()
    # data = pd.DataFrame.from_dict(dfs, orient='index')
    # indata2 = data.transpose()
    # print("method 1:", time() - start)

    start = time()
    indata = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dfs.items() ]))
    print("method 2:", time() - start)
    print("Converting done")
    box_plot(indata)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate full datasets.")
    parser.add_argument('input_csvs', type=str, nargs='*', help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

