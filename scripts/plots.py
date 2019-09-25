
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


def label_transcript(row):
   if row['fsm'] == 1 :
      return 'FSM'
   if row['nic'] == 1 :
      return 'NIC'
   if row['ism'] == 1:
      return 'ISM'
   if row['nnc']  == 1:
      return 'NNC'
   if row['no_splices']  == 1:
      return 'NO_SPLICE'

def splice_site_classification_plot(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    indata['transcript_type'] = indata.apply (lambda row: label_transcript(row), axis=1)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))

    g = sns.catplot(x="transcript_type", #col="Depth",
                data=indata,  hue="read_type", hue_order= ["corrected", "original"],
                order= ["FSM", "ISM", "NIC", "NNC", 'NO_SPLICE'], kind="count", aspect=1)

    # g.set(ylim=(0,15))
    g.set_ylabels("Count")
    g.set_xlabels("Transcript type")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "splice_site_classification.eps"))
    plt.savefig(os.path.join(outfolder, "splice_site_classification.pdf"))
    plt.close()

def number_splices_fsm(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    indata['transcript_type'] = indata.apply (lambda row: label_transcript(row), axis=1)
    indata = indata[indata['transcript_type']=='FSM']
    g = sns.catplot(x="tot_splices", #col="Depth",
                data=indata,  hue="read_type", hue_order= ["corrected", "original"], kind="count", aspect=1)
    # axes = g.axes
    g.set_ylabels("Count")
    g.set_xlabels("Number of splice sites")
    # axes.set_xticks(np.arange(0, 70, step=5) )
    # axes.set_xlim(xlim=(0, 70))
    # g.set_xlim(0,70)
    # g.set_xticks(np.arange(0, 70, step=5))
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "nr_splices.eps"))
    plt.savefig(os.path.join(outfolder, "nr_splices.pdf"))
    plt.close()


def unique_fsm(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    orig = indata[indata['read_type']=='original']
    corr = indata[indata['read_type']=='corrected']

    print('orig:', orig['transcript_fsm_id'].nunique())
    print('corr', corr['transcript_fsm_id'].nunique())
    # print(set(orig['transcript_fsm_id'].unique()))
    # print(set(corr['transcript_fsm_id'].unique()))

    orig_reads = pd.Series(orig.transcript_fsm_id.values,index=orig.acc).to_dict()
    corr_reads = pd.Series(corr.transcript_fsm_id.values,index=corr.acc).to_dict()

    all_fsm_orig = set(orig['transcript_fsm_id'].unique())
    all_fsm_corr = set(corr['transcript_fsm_id'].unique())
    # print(orig_reads)
    fsm_read_absent = set()
    fsm_read_overcorrected = set()
    for read in orig_reads:
        if read not in corr_reads:
            # print('here')
            if orig_reads[read] not in all_fsm_corr:
                # print('orig read not in corr and FSM not found in corr')
                fsm_read_absent.add(orig_reads[read])
        else:
            if orig_reads[read] not in all_fsm_corr:
                # print('orig read not in corr and FSM not found in corr')
                fsm_read_overcorrected.add(orig_reads[read])


    print('Original reads was a FSM but was not in the input data of the corrected reads:',len(fsm_read_absent))
    print('Original reads was a FSM but probably overcorrected/modified in the corrected reads (BAD):', len(fsm_read_overcorrected))

    fsm_read_absent = set()
    fsm_read_overcorrected = set()
    for read in corr_reads:
        if read not in orig_reads:
            # print('here')
            if corr_reads[read] not in all_fsm_orig:
                # print('orig read not in corr and FSM not found in corr')
                fsm_read_absent.add(corr_reads[read])
        else:
            if corr_reads[read] not in all_fsm_orig:
                # print('orig read not in corr and FSM not found in corr')
                fsm_read_overcorrected.add(corr_reads[read])


    print('Corrected reads was a FSM and was not aligned at all in the original data (Good):', len(fsm_read_absent))
    print('Corrected reads was a FSM but did not align as FSM in original data (Good):', len(fsm_read_overcorrected))



def total_error_rate2(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))
    data =indata[indata.read_type == 'corrected']
    sns.distplot(data['error_rate'], norm_hist=False, kde=False, label='Corrected', bins=500, hist_kws=dict(alpha=0.5))
    data = indata[indata.read_type == 'original']
    sns.distplot(data['error_rate'], norm_hist=False,  kde=False, label='Original', bins=500, hist_kws=dict(alpha=0.5))


    plt.xticks(np.arange(0, 10, step=1))
    plt.xlim(0,10)
    plt.xlabel("Error rate (%)")
    # plt.xlabel("Difference to HG38 (%)")
    plt.ylabel("Frequency")
    plt.legend(prop={'size': 12})

    orig = indata[indata['read_type']=='original']
    print(orig.median(axis = 0))
    corr = indata[indata['read_type']=='corrected']
    print(corr.median(axis = 0))

    plt.savefig(os.path.join(outfolder, "total_error_rate2.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate2.pdf"))
    plt.close()


def error_rate_per_cluster_size(input_csv, outfolder):
    indata = pd.read_csv(input_csv)
    # indata = indata[indata['read_type']=='corrected']
    ax = sns.pointplot(x="cluster_size", y="error_rate", hue="read_type", data=indata);
    
    # ax = sns.lineplot(x="cluster_size", y="error_rate", err_style="bars",
    #                 hue="read_type", markers=True, dashes=False, data=indata)

    # g = sns.catplot(x="cluster_size", y="error_rate", #col="Depth",
    #             data=indata,  #hue="read_type", hue_order= ["corrected", "original"],
    #             kind="violin", aspect=1)

    # g.set(ylim=(0,15))
    # g.set_ylabels("Error rate (%)")
    # g.set_xlabels("Cluster size")

    # ax.set_ylim(0,15)
    ax.set_ylabel("Error rate %")
    ax.set_xlabel("Cluster size")
    # ax.set_xscale('log')
    # print([x for x in ax.get_xticklabels()])
    print([x.set_text('{0}'.format(int(int(x.get_text())/1000)) + 'K') for x in ax.get_xticklabels() if int(x.get_text()) > 1000])
    # print([x for x in ax.get_xticklabels()])
    print([x.set_text('') for i, x in enumerate(ax.get_xticklabels()) if i > 3 and i % 2 == 0 ])

    # xlabels = ['{:,.2f}'.format(x/1000) + 'K' if x > 1000 else x for x in ax.get_xticks()]
    # print(xlabels)
    ax.set_xticklabels(ax.get_xticklabels())

    plt.savefig(os.path.join(outfolder, "error_rate_per_cluster_size.eps"))
    plt.savefig(os.path.join(outfolder, "error_rate_per_cluster_size.pdf"))
    plt.close()


def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)

    # total_error_rate2(args.input_csv, args.outfolder)
    error_rate_per_cluster_size(args.input_csv, args.outfolder)
    # total_error_rate(args.input_csv, args.outfolder)
    # splice_site_classification_plot(args.input_csv, args.outfolder)

    # unique_fsm(args.input_csv, args.outfolder)
    # number_splices_fsm(args.input_csv, args.outfolder)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

