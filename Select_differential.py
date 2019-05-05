import pandas as pd
file = pd.read_csv("asdfg.csv")
outfile = file[(abs(file['log2FoldChange'])>=1) & (file['pvalue']<0.05)]
outfile.to_csv('2hourVS0hour.csv', index=False)

