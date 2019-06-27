
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
indata = pd.read_csv(data)
y=sys.argv[3]
g = sns.catplot(x="p", y=y, col="Depth", col_wrap=3,
            data=indata, hue="type", hue_order= ["exact", "approx", "original"],
            kind="violin", aspect=1)

g.set(ylim=(0,15))
g.set_ylabels("Error rate %")
g.set_xlabels("Mutation rate")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()
