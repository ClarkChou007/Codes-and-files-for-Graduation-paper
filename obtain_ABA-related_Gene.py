########################################################################################################################
# The genes are always with specific functions
# the functions have been annotated in GO term annotation in two files from TAIR and DAVID website
# this code was wrote to screen genes according to the databases
# date : 24th of April in 2019
########################################################################################################################

import pandas as pd
import numpy as np

file1 = pd.read_table("Ath_GO_annotation")
# print(file.head())
geneid = list(file1["Gene_id"])
goterm = list(file1["GO_term"])
list_a = []
output = open("ABA_related_gene.txt",'w')
for i in range(len(geneid)):
    if "abscisic acid" in goterm[i] and "pathway"not in goterm[i] and "binding" not in goterm[i]:
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
for information in informations:
    information = information.split("\t")
    if "abscisic acid" in information[4] and "pathway" not in information[4] and "binding" not in information[4]:
        list_b.append(information[0])
        output.write(information[0]+"\t" +information[4]+"\n")
list_b = list(set(list_b))
list_b.sort()
print(list_b)
print(len(list_b))
# file2 = pd.read_csv("ATH_GO_GOSLIM.txt",sep = "\t")
# print(file2.head())
