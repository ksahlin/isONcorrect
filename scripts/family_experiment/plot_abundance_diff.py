
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
indata = pd.read_csv(data, sep='\t')
# y=sys.argv[3]

g = sns.catplot(x="abundance_original", y="abundance_corrected", col="Depth", col_wrap=3,
            data=indata, kind="strip", aspect=1)

# g.set(ylim=(0,100))
g.set_ylabels("Abundance after correction")
g.set_xlabels("Abundance before correction")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()

