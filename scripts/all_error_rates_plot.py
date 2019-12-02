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
# df = (df.groupby(['error_type'])['dataset']
#                      .value_counts(normalize=True)
#                      .rename('percentage')
#                      .mul(100)
#                      .reset_index()
#                      .sort_values('error_type'))
# indata = df.loc[df['k'] == 9]
print(df)

# y=sys.argv[3]
g = sns.catplot(x="dataset", y='percentage', col="type",
            data=df, hue="error_type", hue_order= ["subs", "ins", "del", "matches"],
            kind="bar",  aspect=1)

g.set(ylim=(0.1,100))
plt.yscale('log')
g.set_ylabels("Error rate (%)")
g.set_xlabels("Dataset")


plt.savefig(sys.argv[2])
plt.close()

# indata = df.loc[df['k'] == df['w']]

# y=sys.argv[3]
# g = sns.catplot(x="Depth", y=y,  col="k",col_wrap=1, col_order=[9],
#             data=indata, hue="type", hue_order= ["exact", "approx", "original"],
#             kind="box", aspect=1)

# g.set(ylim=(0,10))
# g.set_ylabels("Error rate (%)")
# g.set_xlabels("Fraction exon present")

# plt.savefig(sys.argv[2])
# plt.close()

# indata2 = df.loc[df['k'] + 2 == df['w'] ]
# g = sns.catplot(x="Depth", y=y,  col="k",col_wrap=3,
#             data=indata2, hue="type", hue_order= ["exact", "approx", "original"],
#             kind="box", aspect=1)

# g.set(ylim=(0,10))
# g.set_ylabels("Error rate (%)")
# g.set_xlabels("Fraction exon present")
# plt.savefig(sys.argv[2]+str('2.pdf'))
# plt.close()