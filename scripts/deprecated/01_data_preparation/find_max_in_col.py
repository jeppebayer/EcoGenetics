from numpy import maximum
import pandas as pd

df = pd.read_csv("people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.coverage",
                sep = "\t")

# print(df["meandepth"].max())

print(df.nlargest(n = 10, columns = ["meandepth"]))
