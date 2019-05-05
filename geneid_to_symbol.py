########################################################################################################################
# TAIRds are a number of serial numbers such as AT1G01010 etc.
# These are not easy to read and understand what their exactly function is.
# We need to transfer the ids to Symbol, which are the abbreviations of the gene name, contain the function
# DATE: 2019-4-10
########################################################################################################################

import pandas as pd

symbol = pd.read_csv("table.csv")


#print(symbol.head())

tairs = list(symbol["TAIR"])
symbs = list(symbol["SYMBOL"])
# # print(tair)
# # ditc = dict(zip(tair,symb))
# # print(ditc)
# setline = {}
# for key, value in pairs:
#     if key not in d:
#         d[key] = []
#     d[key].append(value)

for i in range(len(tairs)):
    if type(symbs[i]) == float:
        symbs[i] = tairs[i]
# print(type(symbs[17]))
# print(symbs)

file = pd.read_csv("choose_matrix12.csv")
# print(file.head())
genes = []
geneid = list(file["GeneID"])
# print(len(geneid))
for j in range(len(geneid)):
    for i in range(len(tairs)):
        if geneid[j] == tairs[i]:
            geneid[j] = symbs[i]
            continue
# print(geneid)
file.drop(['GeneID'],axis = 1, inplace= True)
file.insert(0,"symbol",geneid)
print(file.head(5))
file.to_csv("choose_matrix21.csv",index = None)
