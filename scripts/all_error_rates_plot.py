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
g = sns.catplot(x="dataset", y='fraction', col="type", col_order=["original", "corrected"],
            data=df, hue="error_type", hue_order= ["subs", "ins", "del"],
            kind="bar",  aspect=1)

plt.yscale('log')
g.set_ylabels("Error rate")
g.set_xlabels("Dataset")

# # ax = g.ax
# for ax in g.axes:
#     print(type(ax))
#     for p in ax.patches:
#         ax.annotate('{:.2f}%'.format(100*p.get_height()/1000000.0), (p.get_x()+0.01, p.get_height()+1000), rotation=90, fontsize = 'x-small' )

g.set(ylim=(0.0001,0.1))
g.set(yticks= [0.0001, 0.0002, 0.0003, 0.0004, 0.0005,0.0006,0.0007,0.0008,0.0009, 0.001,0.002,0.003,0.004,0.005,0.006,0.007, 0.008,0.009, 0.01, 0.02,0.03,0.04,0.05] )
g.set(yticklabels= ["", '', '', '', '','', '','','', '0.001','','','','','','', '','', '0.01', '','','',''] )
# g.set(yticks = ["0", "0.001","0.002","0.003","0.004","0.005","0.006","0.007", "0.008","0.009", "0.01"])

g.set_titles("{col_name}")

# plt.setp(g.ax.get_yticklabels(["0.1","1.0", "10.0" ]))


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