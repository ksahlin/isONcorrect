
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

from matplotlib.colors import ListedColormap

# with sns.axes_style('white'):
#     sns.heatmap(pdf,
#                 cbar=False,
#                 square=False,
#                 annot=True,
#                 fmt='g',
#                 cmap=ListedColormap(['white']),
#                 linewidths=0.5)

def is_minor_after_corr(row):
    if row['preserved'] == 1 and row['minor'] == 1: # was minor and did not get overcorrected
        val = 1
    elif row['preserved'] == 0 and row['minor'] == 0: # was major but got overcorrected
        val = 1
    else:
        val = 0
    return val


data=sys.argv[1]

df = pd.read_csv(data)

# df1 = df.loc[df['minor'] == 1] # select only reads that are of minor isoform to study only the miscorrected minor isoforms. Decided to look at the total fraction instead because that is what is asked.
# print(df1)

df['is_minor_after_corr'] = df.apply(is_minor_after_corr, axis=1)

# df2 = df1[['Depth','p', 'preserved']] # focus on three relevant columns only
df2 = df[['Depth','p', 'is_minor_after_corr']] # focus on three relevant columns only
print(df2)

print(df2.groupby(['Depth','p']).sum() / df2.groupby(['Depth','p']).count()) # get fraction of minor isoforms after correction   # before: get fraction of minor isoforms with preserved exon structure

# heatmap_data = df2.groupby(['Depth','p']).count() / len(df2)
heatmap_data = df2.groupby(['Depth','p']).sum() / df2.groupby(['Depth','p']).count()
heatmap_data = heatmap_data.reset_index()
print(heatmap_data)
heatmap_data = heatmap_data[['Depth','p', 'is_minor_after_corr']]
heatmap_data = heatmap_data.pivot("Depth", "p", "is_minor_after_corr")
heatmap_data = heatmap_data.sort_values('Depth', ascending=False)
# heatmap_data = 100*heatmap_data
print(heatmap_data)

plt.clf()

#try 1 and 2
# plt.figure(figsize = (8,5))
# ax = sns.heatmap(heatmap_data, cmap='coolwarm_r', annot=True, vmin=0, vmax=100,fmt='g', linewidths=.7, cbar_kws={'format': '%.0f%%'}) # old working 
# ax = sns.heatmap(heatmap_data, cmap=ListedColormap(['white']), annot=True, vmin=0, vmax=100,fmt='g', square=True, linewidths=.7, cbar=False)
# ax.set_ylabel("Read depth")
# ax.set_xlabel("Minor isoform Frequency")

# try 3
# cm = ['coolwarm_r', 'coolwarm_r', 'coolwarm_r', 'coolwarm_r', 'coolwarm_r']
# f, axs = plt.subplots(1, heatmap_data.columns.size, figsize=(10, 3))
# for i, (s, a, c) in enumerate(zip(heatmap_data, axs, cm)):
# 	vmax = (i+1)*10
# 	sns.heatmap(np.array([heatmap_data[s].values]).T, yticklabels=heatmap_data.index, vmin=0, vmax=vmax, xticklabels=[s], annot=True, fmt='.2f', ax=a, cmap=c)
# 	if i>0:
# 		a.yaxis.set_ticks([])
# f.tight_layout()

# try 4
plot_data = df2.groupby(['Depth','p']).sum() / df2.groupby(['Depth','p']).count()
plot_data = plot_data.reset_index()
plot_data = plot_data[['Depth','p', 'is_minor_after_corr']]
plot_data = plot_data.pivot( "p","Depth", "is_minor_after_corr")
print(plot_data)
ax = plot_data.plot(xticks=[0.1,0.2,0.3,0.4,0.5], yticks=[0.1,0.2,0.3,0.4,0.5], xlim=[0.0,0.6], ylim=[0.0,0.6], xlabel= "Minor isoform fraction (before correction)", ylabel= "Minor isoform fraction (after correction)") 
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='k', ls='--')
plt.savefig(sys.argv[2]+ '_lineplot.pdf')
plt.close()



plt.savefig(sys.argv[2])
plt.close()
