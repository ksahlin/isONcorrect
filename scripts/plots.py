
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

    # g = sns.catplot(x="transcript_type", #col="Depth",
    #             data=indata,  hue="read_type", hue_order= ["corrected", "original"],
    #             order= ["FSM", "ISM", "NIC", "NNC", 'NO_SPLICE'], kind="count", aspect=1)

    g = sns.catplot(x="transcript_type", #col="Depth",
                data=indata,  hue="read_type", hue_order= ["corrected", "original"],
                order= ["FSM", "ISM", "NIC", "NNC"], kind="count", aspect=1)
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


from collections import defaultdict

def transpose_count(dct):
    d = defaultdict(int)
    for key1, value in dct.items():
        d[value]  += 1
    return d

def unique_fsm(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    orig = indata[indata['read_type']=='original']
    corr = indata[indata['read_type']=='corrected']

    # print('orig:', orig['transcript_fsm_id'].nunique())
    # print('corr', corr['transcript_fsm_id'].nunique())
    # print(set(orig['transcript_fsm_id'].unique()))
    # print(set(corr['transcript_fsm_id'].unique()))

    orig_reads = pd.Series(orig.transcript_fsm_id.values,index=orig.acc).to_dict()
    corr_reads = pd.Series(corr.transcript_fsm_id.values,index=corr.acc).to_dict()

    all_fsm_orig = set(orig['transcript_fsm_id'].unique())
    all_fsm_corr = set(corr['transcript_fsm_id'].unique())
    print('In orig but not in corr', len(all_fsm_orig - all_fsm_corr))
    print('In corr but not in orig', len(all_fsm_corr - all_fsm_orig))

    all_fsm_orig = {x for x in all_fsm_orig if x==x}
    all_fsm_corr = {x for x in all_fsm_corr if x==x}
    print('Total unique in orig', len(all_fsm_orig))
    print('Total unique in corr', len(all_fsm_corr))

    # Categories: (0. not present/aligned in Corr, 0'. not present/aligned in Orig)  1. Both nAn, 2. Both FSM same, 2. Both FSM different, 3. Corr_FSM, Orig nAn, 3. Orig_FSM, Corr nAn,
    # print(orig_reads)
    not_present_corr = set(orig_reads) - set(corr_reads) 
    not_present_orig = set(corr_reads) - set(orig_reads) 
    in_both = set(corr_reads) & set(orig_reads)
    both_nan = [ acc for acc in in_both if orig_reads[acc] not in all_fsm_orig and corr_reads[acc] not in all_fsm_corr ]
    both_fsm_same = [ acc for acc in in_both if orig_reads[acc] == corr_reads[acc] and  corr_reads[acc] in all_fsm_corr ]
    both_fsm_different = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] in all_fsm_corr and  orig_reads[acc] in all_fsm_orig ]
    corr_fsm_orig_nan = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] in all_fsm_corr and  orig_reads[acc] not in all_fsm_orig ]
    orig_fsm_corr_nan = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] not in all_fsm_corr and  orig_reads[acc] in all_fsm_orig ]

    print()
    print('not_present_corr', len(not_present_corr))
    print('not_present_orig', len(not_present_orig))
    print('in_both', len(in_both))
    print('both_nan', len(both_nan))
    print('both_fsm_same', len(both_fsm_same))
    print('both_fsm_different', len(both_fsm_different))
    print('corr_fsm_orig_nan', len (corr_fsm_orig_nan))
    print('orig_fsm_corr_nan', len(orig_fsm_corr_nan))
    print()
    reads_per_transcript_orig = transpose_count(orig_reads)
    reads_per_transcript_corr = transpose_count(corr_reads)
    print("transcript_id,FSM_orig,FSM_corr")
    fsm_depth_orig = []
    fsm_depth_corr = []
    occurr_in_orig_read_depth = defaultdict(int)
    for transcript_id in set(all_fsm_orig | all_fsm_corr):
        if transcript_id not in reads_per_transcript_corr:
            print(transcript_id, "not in corrected, depth:{0}.".format(reads_per_transcript_orig[transcript_id]))
            occurr_in_orig_read_depth[reads_per_transcript_orig[transcript_id]] += 1
        elif transcript_id not in reads_per_transcript_orig:
            print(transcript_id, "not in ORIGINAL, depth in corr:{0}.".format(reads_per_transcript_corr[transcript_id]))
        fsm_depth_orig.append( reads_per_transcript_orig[transcript_id] + 1 )
        fsm_depth_corr.append( reads_per_transcript_corr[transcript_id] + 1 )
            # print(transcript_id,reads_per_transcript_orig[transcript_id], reads_per_transcript_corr[transcript_id])
    print(occurr_in_orig_read_depth)
    d = {"fsm_orig" : fsm_depth_orig, "fsm_corr" : fsm_depth_corr}
    df = pd.DataFrame(d)
    # ax.set(xscale="log", yscale="log")
    plt.xscale('log')
    plt.yscale('log')
    ax = sns.scatterplot(x='fsm_orig', y='fsm_corr', alpha= 0.5, data = df)
    # ax.set_title("FSM before and after correction")
    # ax.set_xlim(1,10)
    # ax.set_ylim(1,10)
    ax.set_ylabel("Read depth corrected (log(1+x) scale)")
    ax.set_xlabel("Read depth original (log(1+x) scale)")
    plt.savefig(os.path.join(outfolder, "fsm_breakdown.eps"))
    plt.savefig(os.path.join(outfolder, "fsm_breakdown.pdf"))
    plt.close()

    print("sum",sum([reads_per_transcript_orig[tr_id] for tr_id in set(all_fsm_orig | all_fsm_corr)]), sum([reads_per_transcript_corr[tr_id] for tr_id in set(all_fsm_orig | all_fsm_corr)]) )



def total_error_rate2(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))
    data = indata[indata.read_type == 'original']
    sns.distplot(data['error_rate'], norm_hist=False,  kde=False, label='Original', bins=100, hist_kws=dict(alpha=0.5))
    data =indata[indata.read_type == 'corrected']
    sns.distplot(data['error_rate'], norm_hist=False, kde=False, label='Corrected', bins=100, hist_kws=dict(alpha=0.5))

    plt.xticks(np.arange(0, 10, step=1))
    plt.xlim(0,10)
    plt.xlabel("Error rate (%)")
    # plt.xlabel("Difference to HG38 (%)")
    plt.ylabel("Frequency")
    plt.legend(prop={'size': 12})

    orig = indata[indata['read_type']=='original']
    print(orig.median(axis = 0))
    print(orig.sum(axis = 0))

    corr = indata[indata['read_type']=='corrected']
    print(corr.median(axis = 0))
    print(corr.sum(axis = 0))

    plt.savefig(os.path.join(outfolder, "total_error_rate2.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate2.pdf"))
    plt.close()

def total_error_rate(input_csv, outfolder):
    pd.set_option("display.precision", 8)
    df = pd.read_csv(input_csv)
    df_corr = df.loc[df['read_type'] == 'corrected']
    
    median_error =  df_corr['error_rate'].median()
    df_corr['subs_rate'] = df_corr['subs']/df_corr['aligned_length']
    df_corr['ins_rate'] = df_corr['ins']/df_corr['aligned_length']
    df_corr['del_rate'] = df_corr['del']/df_corr['aligned_length']

    print("median error rate/subs/ins/del corrected:", median_error, 100*df_corr['subs_rate'].median(), 100*df_corr['ins_rate'].median(), 100*df_corr['del_rate'].median())
    # sys.exit()
    
    df_orig = df.loc[df['read_type'] == 'original']

    median_error =  df_orig['error_rate'].median()
    df_orig['subs_rate'] = df_orig['subs']/df_orig['aligned_length']
    df_orig['ins_rate'] = df_orig['ins']/df_orig['aligned_length']
    df_orig['del_rate'] = df_orig['del']/df_orig['aligned_length']
    # print(df_orig['del_rate'])

    print("median error rate/subs/ins/del original:",median_error, 100*df_orig['subs_rate'].median(), 100*df_orig['ins_rate'].median(), 100*df_orig['del_rate'].median())

    error_rate_orig = df_orig['error_rate'].tolist()
    error_rate_corr = df_corr['error_rate'].tolist()
    print("total number of original reads aligned:", len(error_rate_orig))
    print("total number of corrected reads aligned:", len(error_rate_corr))

    # bins = [0.1*i for i in range(300)]
    pyplot.hist(error_rate_corr, 100, range=[0, 20], alpha=0.5, label='Corrected')
    pyplot.hist(error_rate_orig, 100, range=[0, 20], alpha=0.5, label='Original')
    pyplot.legend(loc='upper right')
    # pyplot.xlabel("Difference to genome (%)")
    pyplot.xlabel("Error rate (%)")
    pyplot.ylabel("Read count")
    plt.savefig(os.path.join(outfolder, "error_rate.eps"))
    plt.savefig(os.path.join(outfolder, "error_rate.pdf"))
    

# def error_rate_per_cluster_size(input_csv, outfolder):
#     indata = pd.read_csv(input_csv)
#     # indata = indata[indata['read_type']=='corrected']
#     ax = sns.pointplot(x="cluster_size", y="error_rate", hue="read_type", data=indata);
    
#     # ax = sns.lineplot(x="cluster_size", y="error_rate", err_style="bars",
#     #                 hue="read_type", markers=True, dashes=False, data=indata)

#     # g = sns.catplot(x="cluster_size", y="error_rate", #col="Depth",
#     #             data=indata,  #hue="read_type", hue_order= ["corrected", "original"],
#     #             kind="violin", aspect=1)

#     # g.set(ylim=(0,15))
#     # g.set_ylabels("Error rate (%)")
#     # g.set_xlabels("Cluster size")

#     # ax.set_ylim(0,15)
#     ax.set_ylabel("Error rate %")
#     ax.set_xlabel("Cluster size")
#     # ax.set_xscale('log')
#     # print([x for x in ax.get_xticklabels()])
#     print([x.set_text('{0}'.format(int(int(x.get_text())/1000)) + 'K') for x in ax.get_xticklabels() if int(x.get_text()) > 1000])
#     # print([x for x in ax.get_xticklabels()])
#     print([x.set_text('') for i, x in enumerate(ax.get_xticklabels()) if i > 3 and i % 2 == 0 ])

#     # xlabels = ['{:,.2f}'.format(x/1000) + 'K' if x > 1000 else x for x in ax.get_xticks()]
#     # print(xlabels)
#     ax.set_xticklabels(ax.get_xticklabels())

#     plt.savefig(os.path.join(outfolder, "error_rate_per_cluster_size.eps"))
#     plt.savefig(os.path.join(outfolder, "error_rate_per_cluster_size.pdf"))
#     plt.close()

def get_error_rates(input_csv, read_to_infer = "corrected"):
    dels =[]
    inses =[]
    subses =[]
    matcheses =[]
    # aln_lengths =[]
    read_lengths =[]
    for line in open(input_csv, 'r'):

        # if SIRV
        # experiment_id,acc,read_type,ins,del_,subs,matches,error_rate,read_length,aligned_length,chr_id = line.split(',')

        # if SIRV full        
        # acc,read_type,ins,del_,subs,matches,error_rate,read_length,aligned_length,chr_id = line.split(',')
        
        # if Drosophila
        # acc,read_type,ins,del_,subs,matches,error_rate,read_length,aligned_length,cluster_size, is_unaligned_in_other_method,tot_splices,read_sm_junctions,read_nic_junctions,fsm,nic,ism,nnc,no_splices,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag = line.split(',')
        
        # if SIM
        read,read_length,err,ins,del_,subs,matches,error_rate,read_type,transcript_cov,gene_cov,gene_fam_cov = line.split(',')

        if read_type == read_to_infer:
            # read_length =  (int(subs)+ int(ins) + int(del_))/ (float(err_rate)/100.0)
            # matches = read_length - (int(subs)+ int(ins) + int(del_))

            dels.append(int(del_))
            inses.append(int(ins))
            subses.append(int(subs))
            matcheses.append(int(matches))
            # aln_lengths.append(int(aligned_length))
            read_lengths.append(int(read_length))

    # if SIM:
    tot_bases = sum(read_lengths)
    # tot_bases = sum(read_lengths)


    # err_rates1 = [ (dels[i] + inses[i] + subses[i])/ float(aln_lengths[i]) for i in range(len(matcheses))]
    err_rates2 = [ (dels[i] + inses[i] + subses[i])/ (dels[i] + inses[i] + subses[i] + matcheses[i]) for i in range(len(matcheses))]
    dels_rates = [ (dels[i])/ (dels[i] + inses[i] + subses[i] + matcheses[i]) for i in range(len(matcheses))]
    inses_rates = [ (inses[i])/ (dels[i] + inses[i] + subses[i] + matcheses[i]) for i in range(len(matcheses))]
    subses_rates = [ (subses[i])/ (dels[i] + inses[i] + subses[i] + matcheses[i]) for i in range(len(matcheses))]
    # print(sorted(err_rates1)[int(len(err_rates1)/2) ])
    print(sorted(err_rates2)[int(len(err_rates2)/2) ])
    print("median dels_rates", sorted(dels_rates)[int(len(dels_rates)/2) ])
    print("median inses_rates",sorted(inses_rates)[int(len(inses_rates)/2) ])
    print("median subses_rates", sorted(subses_rates)[int(len(subses_rates)/2) ])
    print("type,ins,del,subs,matches,tot_bases")
    print( "{0},{1},{2},{3},{4},{5}".format(read_to_infer,sum(inses),sum(dels), sum(subses), sum(matcheses), tot_bases))
    # print(sorted(subses_rates))

def sirv_error_rate_per_transcript(input_csv, outfolder):
    # get_error_rates(input_csv)
    # get_error_rates(input_csv, read_to_infer = 'original' )
    # sys.exit()
    pd.set_option("display.precision", 8)
    indata = pd.read_csv(input_csv)

    # df_corr = indata.loc[indata['read_type'] == 'corrected']
    # median_error =  df_corr['error_rate'].median()
    # df_corr['subs_rate'] = df_corr['subs']/df_corr['aligned_length']
    # df_corr['ins_rate'] = df_corr['ins']/df_corr['aligned_length']
    # df_corr['del_rate'] = df_corr['del']/df_corr['aligned_length']
    # df_corr['error_rate_new'] = (df_corr['subs'] + df_corr['ins'] + df_corr['del'] ) /df_corr['aligned_length']

    # df_orig = indata.loc[indata['read_type'] == 'original']
    # median_error =  df_orig['error_rate'].median()
    # df_orig['subs_rate'] = df_orig['subs']/df_orig['aligned_length']
    # df_orig['ins_rate'] = df_orig['ins']/df_orig['aligned_length']
    # df_orig['del_rate'] = df_orig['del']/df_orig['aligned_length']
    # df_orig['error_rate_new'] = (df_orig['del'] + df_orig['ins'] + df_orig['del'] ) /df_orig['aligned_length']


    # df2 = indata.groupby(['chr_id', "experiment_id"])['chr_id', "experiment_id"].transform('count')
    # df2 = df2.rename(columns={"chr_id": "transcript_cov", "experiment_id": "transcript_cov2"})
    # indata = pd.concat([indata, df2], axis=1)
    # ax = sns.lineplot(x="transcript_cov", y="error_rate",  hue="read_type",
    #                   ci = 'sd', estimator= 'median',  data=indata)

    # new controlled experiment:
    ax = sns.lineplot(x="nr_reads", y="error_rate",  hue="read_type",
                      ci = 'sd', estimator= 'median',  data=indata)
    ax.set_ylim(0,12)
    ax.set_xscale('log')
    ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
    ax.set_ylabel("Error rate (%)")
    ax.set_xlabel("Reads per transcript")

    plt.savefig(os.path.join(outfolder, "error_rates.eps"))
    plt.savefig(os.path.join(outfolder, "error_rates.pdf"))
    plt.close()



def main(args):
    sns.set(style="whitegrid")
    # sns.set()
    # flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    # sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)

    get_error_rates(args.input_csv, read_to_infer = 'corrected')
    get_error_rates(args.input_csv, read_to_infer = 'original' )

    # splice_site_classification_plot(args.input_csv, args.outfolder)
    # unique_fsm(args.input_csv, args.outfolder)
    # total_error_rate(args.input_csv, args.outfolder)
    # sirv_error_rate_per_transcript(args.input_csv, args.outfolder)


    # total_error_rate2(args.input_csv, args.outfolder)
    # error_rate_per_cluster_size(args.input_csv, args.outfolder)
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

