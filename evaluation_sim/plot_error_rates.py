
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
# y=sys.argv[3]

df_corr = indata.loc[indata['type'] == 'corrected']
median_error =  df_corr['err_rate'].median()
# df_corr['subs_rate'] = df_corr['subs']/df_corr['read_length']
# df_corr['ins_rate'] = df_corr['ins']/df_corr['read_length']
# df_corr['del_rate'] = df_corr['del']/df_corr['read_length']

print("median error rate/subs/ins/del corrected:", median_error) #, 100*df_corr['subs_rate'].median(), 100*df_corr['ins_rate'].median(), 100*df_corr['del_rate'].median())
# sys.exit()
ab = 50
df_corr_ab = df_corr.loc[df_corr['transcript_cov'] == ab ]
median_error_ab =  df_corr_ab['err_rate'].median()
print("median error rate for ab:", ab, median_error_ab)
    

df_orig = indata.loc[indata['type'] == 'original']

median_error =  df_orig['err_rate'].median()
# df_orig['subs_rate'] = df_orig['subs']/df_orig['read_length']
# df_orig['ins_rate'] = df_orig['ins']/df_orig['read_length']
# df_orig['del_rate'] = df_orig['del']/df_orig['read_length']
# print(df_orig['del_rate'])

print("median error rate/subs/ins/del original:",median_error) #, 100*df_orig['subs_rate'].median(), 100*df_orig['ins_rate'].median(), 100*df_orig['del_rate'].median())


# g = sns.catplot(x="abundance_original", y="abundance_corrected", col="Depth", col_wrap=3,
#             data=indata, kind="point", aspect=1)
# ax = sns.scatterplot(x="abundance_original", y="abundance_corrected", #col="Depth", # col_wrap=3, #hue="time",
#                       data=indata)
# g = sns.lmplot(x="transcript_cov", y="err_rate",  hue="type", scatter= False, x_estimator = np.mean,  #x_jitter = 0.1, #col="Depth", 
#                 x_ci= 'sd', data=indata) #, col_wrap=2, height=3)

ax = sns.lineplot(x="transcript_cov", y="err_rate",  hue="type", 
                  ci = 'sd', estimator='median', data=indata)
ax.set_ylim(0,12)
ax.set_xscale('log')
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
ax.set_ylabel("Error rate (%)")
ax.set_xlabel("Reads per transcript")

# g.set(ylim=(0,100))
# ax.set_ylabel("Abundance after correction")
# ax.set_xlabel("Abundance before correction")

# g.set_ylabels("Error rate (%)")
# g.set_xlabels("Reads per transcript")

# ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
# ax.set_ylim(0,15)
# ax.set_ylabel("Error rate %")

plt.savefig(sys.argv[2])
plt.close()

