
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




def total_error_rate(input_csv, outfolder, dataset):

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

    plt.savefig(os.path.join(outfolder, dataset+ "_full.eps"))
    plt.savefig(os.path.join(outfolder, dataset+ "_full.pdf"))
    plt.close()



def get_error_rates(input_csv, dataset, read_to_infer = "corrected"):
    dels =[]
    inses =[]
    subses =[]
    matcheses =[]
    # aln_lengths =[]
    read_lengths =[]
    for line in open(input_csv, 'r'):
        if dataset == 'sirv':
            # if SIRV full        
            acc,read_type,ins,del_,subs,matches,error_rate,read_length,aligned_length,chr_id = line.split(',')
        
        if dataset == 'drosophila':
            # if Drosophila
            acc,read_type,ins,del_,subs,matches,error_rate,read_length,aligned_length,cluster_size, is_unaligned_in_other_method,tot_splices,read_sm_junctions,read_nic_junctions,fsm,nic,ism,nnc,no_splices,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag = line.split(',')

        
        if dataset == 'sim':    
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

    tot_bases = sum(read_lengths)
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

    errors = [sum(inses), sum(dels), sum(subses), sum(matcheses), tot_bases]
    error_types = ["ins","del","subs","matches","tot_bases"]
    for err, err_type in zip(errors, error_types) :
        print("{0},{1},{2},{3}".format(dataset, read_to_infer, err_type, err/float(tot_bases)))

        # print("type,ins,del,subs,matches,tot_bases")
    # print( "{0},{1},{2},{3},{4},{5}".format(read_to_infer,sum(inses),sum(dels), sum(subses), sum(matcheses), tot_bases))
    # print(sorted(subses_rates))


def main(args):
    sns.set(style="whitegrid")
    # sns.set()
    # flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    # sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)

    get_error_rates(args.input_csv, args.dataset, read_to_infer = 'corrected')
    get_error_rates(args.input_csv, args.dataset, read_to_infer = 'original')

    total_error_rate(args.input_csv, args.outfolder, args.dataset)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate full datasets.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    parser.add_argument('dataset', type=str, help='which dataset')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

