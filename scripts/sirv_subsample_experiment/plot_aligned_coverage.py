
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
# indata = df.loc[df['seq'] == str(seq)]
# print(indata.info())
# indata[["id"]] = indata.apply(pd.to_str)
# indata[["cov_true", "cov_aln"]] = indata.apply(pd.to_numeric)

# indata["cov_aln"] = indata.to_numeric(indata["cov_aln"])
# y=sys.argv[3]

# g = sns.catplot(x="abundance_original", y="abundance_corrected", col="Depth", col_wrap=3,
#             data=indata, kind="point", aspect=1)
# ax = sns.scatterplot(x="abundance_original", y="abundance_corrected", #col="Depth", # col_wrap=3, #hue="time",
#                       data=indata)
# g = sns.lmplot(x="cov_true", y="cov_aln",  hue="type", #x_jitter = 0.3, #fit_reg= True, x_jitter = 0.2, y_jitter = 0.2, #col="Depth", 
#                 data=indata) #, col_wrap=2, height=3)

ax = sns.barplot(x="transcript",  y="depth", hue="type", data=indata)
# g = sns.catplot(x="transcript", y="depth", col="Depth", col_wrap=3,
#             data=indata, hue="type", hue_order= ["isoncorrect", "original"],
#             kind="point", aspect=1)

# # g.set(ylim=(0,100))
# g.set_ylabels("Reads aligned")
# g.set_xlabels("Transcript ID")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
ax.set_ylabel("Reads aligned")
ax.set_xlabel("Transcript ID")

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

