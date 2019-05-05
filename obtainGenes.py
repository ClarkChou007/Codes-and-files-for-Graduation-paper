########################################################################################################################
# The genes are always with specific functions
# the functions have been annotated in GO term annotation in two files from TAIR and DAVID website
# this code was wrote to screen genes according to the databases
# date : 4th of May in 2019
########################################################################################################################
import matplotlib.pyplot as plt

data = [1, 3, 5, 6, 9, 1, 5, 8, 7, 2]

plt.plot(data, color='g')

plt.scatter(list(range(10)), data, color='r')

plt.savefig('000.png')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

file1 = pd.read_table("Ath_GO_annotation")
# print(file.head())
geneid = list(file1["Gene_id"])
goterm = list(file1["GO_term"])
list_a = []
output = open("drought_related_gene.txt",'w')
for i in range(len(geneid)):
    if "resonse to water d"in goterm[i]:
        output.write(geneid[i]+"\t"+ goterm[i]+"\n")
        list_a.append(geneid[i])
list_a = list(set(list_a))
list_a.sort()
print(list_a)
print(len(list_a))


file2 = open("ATH_GO_GOSLIM.txt")
informations = file2.readlines()
# print(informations)

list_b= []
list_c =[]
for information in informations:
    information = information.split("\t")
    if "response to water d"in information[4]:
        list_b.append(information[0])
        output.write(information[0]+"\t" +information[4]+"\n")
    if "abscisic acid" in information[4] and "binding" not in information[4]:
        list_c.append(information[0])
list_b = list(set(list_b))
list_b.sort()
list_c.sort()


print(list_b)
print(len(list_b))
# file2 = pd.read_csv("ATH_GO_GOSLIM.txt",sep = "\t")
# print(file2.head())
venn2(subsets = [set(list_b),set(list_c)],set_labels=("response to drought","response to ABA"), set_colors =("r","b"))
plt.show()
