import sys
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plotdir = os.path.dirname(sys.argv[1])
plotname = os.path.basename(plotdir)
mapqfile = pd.read_table(sys.argv[1], header=None)

maxval = max(mapqfile[0])

plt.figure(figsize = (15, 8.5))
plt.style.use("ggplot")
plt.hist(mapqfile[0], bins = maxval)
plt.ylabel("Count")
plt.xlabel("MapQ")
plt.xticks(np.arange(0, maxval+1, 2.0))
plt.axvline(np.mean(mapqfile[0]), color = "k", linestyle = "dashed")
plt.savefig("{}/{}_mapq_histogram.png".format(plotdir, plotname))