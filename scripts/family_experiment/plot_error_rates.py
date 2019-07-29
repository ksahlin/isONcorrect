
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
df = pd.read_csv(data)
print(len(df))
df['q_acc'] = df['q_acc'].apply(lambda x: x.split("_")[0])
indata = df.loc[df['q_acc'] == df['r_acc']]
print(len(indata))

y=sys.argv[3]
g = sns.catplot(x="Depth", y=y, #col="Depth", col_wrap=3,
            data=indata, hue="type", hue_order= ["exact", "approx", "original"],
            kind="violin", aspect=1)

g.set(ylim=(0,15))
# g.ax.set_yscale("log")
# g.set(ylim=(0.01,15))
# g.ax.set_yticks([1e-1,1e0,1e1,1e2])
# g.set_yticklabels(("0.1%","1%","10%","100%"), visible = True)
g.set_ylabels("Error rate")
g.set_xlabels("Read depth")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()
