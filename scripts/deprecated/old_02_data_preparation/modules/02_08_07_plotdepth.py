import sys
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plotdir = os.path.dirname(sys.argv[1])
plotname = os.path.basename(sys.argv[2])
depthfile = pd.read_table(sys.argv[1], header=None)

mean = np.mean(depthfile[2])
maxval = int(mean) + 100
limit = str(maxval)

# Changes all values above a certain limit to the same value
# depthfile.loc[depthfile[2] >= (maxval + 1), 2] = (maxval + 1)

# Subsets dataframe to defined max value
depthfile = depthfile[depthfile[2] <= maxval]

plt.figure(figsize = (15, 8.5))
plt.style.use("ggplot")
# plt.hist(depthfile[2], bins = maxval + 1)
plt.hist(depthfile[2], bins = maxval)
plt.ylabel("Count")
plt.xlabel("Coverage")
# plt.xticks(np.arange(0, maxval + 1 , 5.0), rotation="vertical")
plt.xticks(np.arange(0, maxval, 5.0), rotation="vertical")
plt.axvline(mean, color = "k", linestyle = "dashed")
plt.savefig("{}/{}_filtered_depth<{}_histogram.png".format(plotdir, plotname, limit))