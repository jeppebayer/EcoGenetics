import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t', lineterminator='\n', header = None)

sumlist = []
for num in range(len(df.columns)):
    sumlist.append((df[num].sum()).astype(str))

with open("people/Jeppe_Bayer/scripts/tests/test4.txt", "w") as file:
    file.write("\t".join(sumlist) + "\n")

