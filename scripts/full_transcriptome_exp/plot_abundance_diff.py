
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
    acc, annot, ab, is_corr, read_type = row
    if acc not in read_annot:
        read_annot[acc] = [None,None, ab]

    if read_type == 'corrected':
        read_annot[acc][1] =  is_corr
    else:
        read_annot[acc][0] =  is_corr

wrong_in_both = []
wrong_in_orig = []
wrong_in_corr = []

for acc in read_annot:
    orig, corr, ab = read_annot[acc]
    if orig == 1 and corr == 1:
        both_good += 1
    elif  orig == 1 and corr == 0:
        orig_good += 1
        wrong_in_corr.append(acc)
    elif  orig == 0 and corr == 1:
        corr_good += 1
        wrong_in_orig.append(acc)
    else:
        both_bad += 1
        wrong_in_both.append(acc)

print(wrong_in_both)
print(wrong_in_orig)
# print(sum(wrong_in_both)/float(len(wrong_in_both)))
# print(sum(wrong_in_orig)/float(len(wrong_in_orig)))
# print(sum(wrong_in_corr)/float(len(wrong_in_corr)))
print(wrong_in_corr)
print(both_good, orig_good, corr_good, both_bad)



# indata = df.loc[df['seq'] == 'transcript']
# print(indata.info())
# indata[["id"]] = indata.apply(pd.to_str)
# indata[["cov_true", "cov_aln"]] = indata.apply(pd.to_numeric)

# indata["cov_aln"] = indata.to_numeric(indata["cov_aln"])
# y=sys.argv[3]

# g = sns.catplot(x="abundance_original", y="abundance_corrected", col="Depth", col_wrap=3,
#             data=indata, kind="point", aspect=1)
# ax = sns.scatterplot(x="abundance_original", y="abundance_corrected", #col="Depth", # col_wrap=3, #hue="time",
#                       data=indata)
# g = sns.lmplot(x="cov_true", y="cov_aln",  hue="type", scatter = False, x_ci= 'sd', #x_bins=50, #x_jitter = 0.3, #fit_reg= True, x_jitter = 0.2, y_jitter = 0.2, #col="Depth", 
#                data=indata, scatter_kws={'alpha':0.3}) #, col_wrap=2, height=3)

ax = sns.lineplot(x="transcript_abundance", y="is_tp",  hue="read_type", 
                  ci = 'sd', data=indata)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylabel("Fraction correctly aligned back to true transcript")
ax.set_xlabel("True transcript abundance")
# g.set(xscale="log", yscale="log")
# g.set(ylim=(0,100))
# g.set_ylabels("Aligned transcript abundance")
# g.set_xlabels("True transcript abundance")


# g = sns.catplot(x="cov_aln", y="cov_true", col="Depth", col_wrap=3,
#             data=indata, hue="type", hue_order= ["isoncorrect", "original"],
#             kind="point", aspect=1)
# g.set_ylabels("% Reads switched transcript")
# g.set_xlabels("True transcript abundance in simulation")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()

# g.set(ylim=(0,100))
# ax.set_ylabel("Abundance after correction")
# ax.set_xlabel("Abundance before correction")

# g.set_ylabels("Abundance inferred from alignment")
# g.set_xlabels("True abundance per transcript")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

# plt.savefig(sys.argv[3])
# plt.close()

