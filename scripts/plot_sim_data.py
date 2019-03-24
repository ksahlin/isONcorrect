
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
indata = pd.read_csv(data, sep="\t")

ax = sns.boxplot(x="p", y="rate",  data=indata)
plt.savefig(sys.argv[2])
plt.close()
