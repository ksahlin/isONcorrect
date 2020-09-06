
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
df1 = df.loc[df['minor'] == 1] # select only reads that are of minor isoform
print(df1)
df2 = df1[['Depth','p', 'preserved']] # focus on three columns only as now all rows are about minor isoforms
print(df2)
# print(df2.groupby(['Depth','p']).count() / len(df2)) # get fraction of minor isoforms with preserved exon structure
print(df2.groupby(['Depth','p']).sum() / df2.groupby(['Depth','p']).count()) # get fraction of minor isoforms with preserved exon structure

# heatmap_data = df2.groupby(['Depth','p']).count() / len(df2)
heatmap_data = df2.groupby(['Depth','p']).sum() / df2.groupby(['Depth','p']).count()
heatmap_data = heatmap_data.reset_index()
print(heatmap_data)
heatmap_data = heatmap_data[['Depth','p', 'preserved']]
heatmap_data = heatmap_data.pivot("Depth", "p", "preserved")
heatmap_data = heatmap_data.sort_values('Depth', ascending=False)
heatmap_data = 100*heatmap_data

# heatmap_data = df1.groupby(['Depth', 'p']).sum()
# heatmap_data = df.groupby(['Depth', 'p']).sum().reset_index() 
# heatmap_data = heatmap_data[['Depth','p', 'preserved']]
# heatmap_data = heatmap_data.pivot("Depth", "p", "preserved")
# heatmap_data = heatmap_data.sort_values('Depth', ascending=False)
# heatmap_data = 10*heatmap_data
print(heatmap_data)

plt.clf()
# xticks= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,">15"]
plt.figure(figsize = (8,5))
ax = sns.heatmap(heatmap_data, cmap='coolwarm_r', annot=True, vmin=0, vmax=100,fmt='g', linewidths=.7, cbar_kws={'format': '%.0f%%'})
ax.set_ylabel("Read depth")
ax.set_xlabel("Minor isoform Frequency")

plt.savefig(sys.argv[2]+ '_heatmap.pdf')
# print("plotting", sys.argv[2]+ '_heatmap.pdf')
plt.close()


plt.savefig(sys.argv[2])
plt.close()
