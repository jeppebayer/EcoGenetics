import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("people/Jeppe_Bayer/steps/samtools/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.coverage", sep = "\t")

plt.scatter(data["#rname"], data["coverage"])
plt.ylabel("Coverage %"), plt.xlabel("Reference name")

plt.savefig("people/Jeppe_Bayer/plots/coverage.png")