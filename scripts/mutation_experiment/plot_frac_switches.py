
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
# f = open(data, "r")
df = pd.read_csv(data)
minor_indata = df.loc[df['minor'] == 1]
# y=sys.argv[3]

g = sns.catplot(x="p", y="mutation_present", col="Depth", col_wrap=3,
            data=minor_indata, hue="type", hue_order= ["exact", "approx", "original"],
            kind="bar", aspect=1)

# g.set(ylim=(0,100))
g.set_ylabels("Minor mutation retained %")
g.set_xlabels("Fraction mutation present in data")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=minor_indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()


# g = sns.catplot(x="p", y="Major_Retained", col="Depth", col_wrap=3,
#             data=indata, hue="type", hue_order= ["exact", "approx", "original"],
#             kind="violin", aspect=1)

# g.set(ylim=(0,100))
# g.set_ylabels("Major mutation retained %")
# g.set_xlabels("Fraction mutation present in data")

# plt.savefig(sys.argv[3])
# plt.close()