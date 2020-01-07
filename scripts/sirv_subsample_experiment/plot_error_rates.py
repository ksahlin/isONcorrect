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
# indata = df.loc[df['k'] == 9]

y=sys.argv[3]
g = sns.catplot(x="Depth", y=y, col="w", row = "k",
            data=df, hue="type", hue_order= ["exact", "approx", "original"],
            kind="box", aspect=1)

g.set(ylim=(0,12))
g.set_ylabels("Error rate (%)")
g.set_xlabels("Read depth")
g.set_titles("{row_var} = {row_name}, {col_var} = {row_var} + {col_name}")
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