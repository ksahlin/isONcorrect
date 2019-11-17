
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
df = pd.read_csv(data) #, dtype={'id': str, "cov_true": int, "cov_aln": int, "type": str})
indata = df.loc[df['seq'] == 'transcript']
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

ax = sns.lineplot(x="cov_true", y="cov_aln",  hue="type", 
                  ci = 'sd', data=indata)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylabel("Aligned transcript abundance")
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

