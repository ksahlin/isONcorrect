
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

sns.set(style="whitegrid")

data=sys.argv[1]
# seq=sys.argv[2]

# f = open(data, "r")
indata = pd.read_csv(data) #, dtype={'id': str, "cov_true": int, "cov_aln": int, "type": str})
list_of_list = indata.values.tolist()
# print(list_of_list[-1])
read_annot = {}
orig_good = 0
corr_good = 0
both_good = 0
both_bad = 0
for row in list_of_list:
    acc, annot, ab, is_corr, read_type, ed_to_correct = row
    if acc not in read_annot:
        read_annot[acc] = [None, None, ab, None, None]

    if read_type == 'corrected':
        read_annot[acc][1] =  is_corr
        read_annot[acc][4] =  ed_to_correct
    else:
        read_annot[acc][0] =  is_corr
        read_annot[acc][3] =  ed_to_correct

wrong_in_both = []
wrong_in_orig = []
wrong_in_corr = []

for acc in read_annot:
    orig, corr, ab, ed_orig, ed_corr = read_annot[acc]
    if orig == 1 and corr == 1:
        both_good += 1
    elif orig == 1 and corr == 0:
        orig_good += 1
        wrong_in_corr.append(ed_corr)
    elif  orig == 0 and corr == 1:
        corr_good += 1
        wrong_in_orig.append(ed_orig)
    else:
        both_bad += 1
        wrong_in_both.append(ed_orig)

# print(wrong_in_both)
print(wrong_in_orig)
# print(sum(wrong_in_both)/float(len(wrong_in_both)))
# print(sum(wrong_in_orig)/float(len(wrong_in_orig)))
# print(sum(wrong_in_corr)/float(len(wrong_in_corr)))
print(wrong_in_corr)
print(both_good, orig_good, corr_good, both_bad)
pyplot.hist(wrong_in_corr, 200, alpha=0.5, label='Corrected')
pyplot.hist(wrong_in_orig, 200, alpha=0.5, label='Original')
# plt.xscale('log')
pyplot.legend(loc='upper right')
pyplot.xlabel("Local edit distance to true ref")
pyplot.ylabel("Count")
plt.savefig(sys.argv[2] + '.eps')
plt.savefig(sys.argv[2])




# ax = sns.lineplot(x="transcript_abundance", y="is_tp",  hue="read_type", 
#                   ci = 'sd', data=indata)
# ax.set_ylabel("Fraction correctly aligned back to true transcript")
# ax.set_xlabel("True transcript abundance")
# plt.savefig(sys.argv[2])
# plt.close()

