########################################################################################################################
# This code will show gene related to response to water deprivation and genes related to ABA
# The analysis is based on DEGs analysis and GO annotation
# The functions have been annotated in GO term annotation in two files from TAIR and DAVID website
# This code was wrote to screen genes according to the databases
########################################################################################################################

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import os
from PIL import Image

file2 = open("ATH_GO_GOSLIM.txt")
informations = file2.readlines()
# print(informations)

list_b= []
list_c =[]
for information in informations:
    information = information.split("\t")
    if "response to water d"in information[4]:
        list_b.append(information[0])
        # output.write(information[0]+"\t" +information[4]+"\n")
    if "abscisic acid" in information[4] and "binding" not in information[4]:
        list_c.append(information[0])
list_b = list(set(list_b))
list_b.sort()
list_c.sort()


filenames = ["2hourVS0hour_upregulation.csv",
             "2hourVS0hour_downregulation.csv",
             "4hourVS0hour_upregulation.csv",
             "4hourVS0hour_downregulation.csv",
             "6hourVS0hour_upregulation.csv",
             "6hourVS0hour_downregulation.csv",
             "8hourVS0hour_upregulation.csv",
             "8hourVS0hour_downregulation.csv",
             "2hourVS0hour.csv",
            "4hourVS0hour.csv",
            "6hourVS0hour.csv",
            "8hourVS0hour.csv"]

for filename in filenames:
    file = pd.read_csv(filename)
    # print(file.head())
    genes = list(file["GeneID"])
    #print(list(genes))
    drought = []
    abscisic_acid = []
    for gene in genes:
        if gene in list_b:
            drought.append(gene)
        if gene in list_c:
            abscisic_acid.append(gene)
    venn2(subsets=[set(drought), set(abscisic_acid)], set_labels=("response to drought", "response to ABA"),
          set_colors=("r", "b"))
    # plt.show()
    if "regulation" in filename:
        name = filename.strip("regulation.csv")
    else:
        name = filename.strip(".csv")
    plt.axis('off')

    plt.savefig(name+".png")  # 保存
    plt.close()  # 关闭图 0



