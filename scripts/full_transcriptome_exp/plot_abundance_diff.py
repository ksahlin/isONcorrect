
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
from collections import defaultdict

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
x = []
y = []
overcorrection_amount = []
abundance = []
overcorr_and_ab = defaultdict(lambda: defaultdict(int))
# overcorr_and_ab = []
# for i in range(17):
#     overcorr_and_ab.append([])
#     for j in range(19):
#         overcorr_and_ab[-1].append(0)
for row in list_of_list:
    acc, annot, ab, is_corr, read_type, ed_btw_transcripts, ed_read_to_true, ed_read_to_aligned = row
    if acc not in read_annot:
        read_annot[acc] = [None, None, ab, None, None]

    if read_type == 'corrected':
        read_annot[acc][1] =  is_corr
        read_annot[acc][4] =  ed_read_to_true
    else:
        read_annot[acc][0] =  is_corr
        read_annot[acc][3] =  ed_read_to_true

    if is_corr == 0 and read_type == 'corrected' and ed_read_to_aligned < ed_read_to_true:
        x.append(ed_read_to_aligned)
        y.append(ed_read_to_true)
        overcorrection_amount.append(ed_read_to_true - ed_read_to_aligned)
        abundance.append(ab)

        if ed_read_to_true - ed_read_to_aligned > 15:
            overcorr_and_ab[16][ab] += 1
        else:
            overcorr_and_ab[ed_read_to_true - ed_read_to_aligned][ab] += 1
        # if ab > 10:
        #     index = 10 + int(ab/10) - 2
        # else:
        #     index = ab - 1
        # if ed_read_to_true - ed_read_to_aligned > 15:
        #     overcorr_and_ab[-1][index] += 1
        # else:
        #     overcorr_and_ab[ed_read_to_true - ed_read_to_aligned][index] += 1

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

print("over corrected:", len(overcorrection_amount) )
# print(len(wrong_in_both), len(wrong_in_orig), len(wrong_in_corr))
print(abundance)
# # print(wrong_in_both)
# print(wrong_in_orig)
# # print(sum(wrong_in_both)/float(len(wrong_in_both)))
# # print(sum(wrong_in_orig)/float(len(wrong_in_orig)))
# # print(sum(wrong_in_corr)/float(len(wrong_in_corr)))
# print(wrong_in_corr)
# print(both_good, orig_good, corr_good, both_bad)
# pyplot.hist(wrong_in_corr, 200, alpha=0.5, label='Corrected')
# pyplot.hist(wrong_in_orig, 200, alpha=0.5, label='Original')
# # plt.xscale('log')
# pyplot.legend(loc='upper right')
# pyplot.xlabel("Local edit distance to true ref")
# pyplot.ylabel("Count")
# plt.savefig(sys.argv[2] + '.eps')
# plt.savefig(sys.argv[2])

print(len([ed for ed in overcorrection_amount if ed <=15]))
print(len([ed for ed in overcorrection_amount if ed <=5]))

print(len([ab for ab in abundance if ab <=15]))
print(len([ab for ab in abundance if ab <=5]))

pyplot.hist(overcorrection_amount, 40, alpha=0.5)
# plt.xscale('log')
# pyplot.legend(loc='upper right')
pyplot.xlabel("Overcorrected (edit distance)")
pyplot.ylabel("Count")
plt.savefig(sys.argv[2] + '.eps')
plt.savefig(sys.argv[2])

print(overcorr_and_ab)

# import numpy as np; np.random.seed(0)
# uniform_data = np.random.rand(10, 12)
# print(uniform_data)
# flights = sns.load_dataset("flights")
# print(flights)
# d = {"overcorrection_amount" : overcorrection_amount, "abundance" : abundance}
df = pd.DataFrame(overcorr_and_ab)
df = df.reindex(sorted(df.columns), axis=1)
df = df.reindex(sorted(df.index, reverse=True), axis=0)
# df.sort_index(level=1, ascending=True, inplace=True)

print(df)
plt.clf()
xticks= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,">15"]
ax = sns.heatmap(df, cmap='coolwarm', annot=True, xticklabels = xticks)
ax.set_ylabel("Abundance")
ax.set_xlabel("Overcorrection (edit distance)")

plt.savefig(sys.argv[2]+ '_heatmap.pdf')
# ax = sns.scatterplot(x="overcorrection_amount", y="abundance", data=df)
# plt.savefig(sys.argv[2]+ '_scatter.pdf')
# plt.savefig(sys.argv[2]+ '_scatter.eps')

# d = {"ed_read_to_aligned" : x, "ed_read_to_true" : y}
# df = pd.DataFrame(d)
# ax = sns.scatterplot(x='ed_read_to_true', y='ed_read_to_aligned', alpha= 0.5, data = df)
# plt.savefig(sys.argv[2] + '.eps')
# plt.savefig(sys.argv[2])

# ax = sns.lineplot(x="transcript_abundance", y="is_tp",  hue="read_type", 
#                   ci = 'sd', data=indata)
# ax.set_ylabel("Fraction correctly aligned back to true transcript")
# ax.set_xlabel("True transcript abundance")
# plt.savefig(sys.argv[2])
# plt.close()

