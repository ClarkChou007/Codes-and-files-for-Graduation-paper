# test code has been annotated
# #open the data from featureCounts
# file = open("drought_treat.txt")
# informations = file.readlines()
# file.close()
# length = []
# counts = []
# #print(informations[0],informations[1])
# information = informations[0].split(" ")
# zcv = informations[1].split("\t")
# print(information)
# print(zcv)
# # obtain the readcounts for one treatment of them
#
# for information in informations:
#     information = information.split("\t")
#     if "#" not in information[0] and "Geneid" not in information[0] :
#          length.append(information[5])
#          counts.append(information[6])
# # print(len(length))
# # print(len(counts))
# # rewrite the file
# # file1 = open("drought_treat.txt",'w')
# # file1.write(informations[0])
# # file1.write(informations[1]
# # calculate the EXPRESSION LEVEL
# count = [int(a) for a in counts]
# mappable = sum(count)
# print(type(count[1]))
# RPKMs = []
# for  i in range(len(counts)):
#     RPKM = count[i]/(int(length[i])/1000*mappable/1000000)
#     RPKMs.append(RPKM)
# print(RPKMs)

import pandas as pd

from pandas.core.frame import DataFrame
# open the file of expression matrix
file = pd.read_csv("drought_treat.csv")
# remove the gene and Chromosome information leave the length of gene and featureCounts
file1 = file.drop(columns = ['Chr','Start','End','Strand'])

# obtain the information row by row. Each row contains one sample

Geneid = list(file1["Geneid"])
Length = list(file1["Length"])
DR0_1 = [int(dor) for dor in list(file1["drought_0_1.bam"])]
DR0_2 = [int(dor) for dor in list(file1["drought_0_2.bam"])]
DR0_3 = [int(dor) for dor in list(file1["drought_0_3.bam"])]
DR2_1 = [int(dor) for dor in list(file1["drought_2_1.bam"])]
DR2_2 = [int(dor) for dor in list(file1["drought_2_2.bam"])]
DR2_3 = [int(dor) for dor in list(file1["drought_2_3.bam"])]
DR4_1 = [int(dor) for dor in list(file1["drought_4_1.bam"])]
DR4_2 = [int(dor) for dor in list(file1["drought_4_2.bam"])]
DR4_3 = [int(dor) for dor in list(file1["drought_4_3.bam"])]
DR6_1 = [int(dor) for dor in list(file1["drought_6_1.bam"])]
DR6_2 = [int(dor) for dor in list(file1["drought_6_2.bam"])]
DR6_3 = [int(dor) for dor in list(file1["drought_6_3.bam"])]
DR8_1 = [int(dor) for dor in list(file1["drought_8_1.bam"])]
DR8_2 = [int(dor) for dor in list(file1["drought_8_2.bam"])]
DR8_3 = [int(dor) for dor in list(file1["drought_8_3.bam"])]

# Calculate the relative expression level by RPKM
dr0_1 = []
mappable0_1 = sum(DR0_1)
for i in range(len(DR0_1)):
    RPKM = DR0_1[i] / (int(Length[i]) / 1000 * mappable0_1 / 1000000)
    dr0_1.append(RPKM)
dr0_2 = []
mappable0_2 = sum(DR0_2)
for i in range(len(DR0_2)):
    RPKM = DR0_2[i] / (int(Length[i]) / 1000 * mappable0_2 / 1000000)
    dr0_2.append(RPKM)
dr0_3= []
mappable0_3 = sum(DR0_3)
for i in range(len(DR0_3)):
    RPKM = DR0_3[i] / (int(Length[i]) / 1000 * mappable0_3 / 1000000)
    dr0_3.append(RPKM)
dr2_1 = []
mappable2_1 = sum(DR2_1)
for i in range(len(DR2_1)):
    RPKM = DR2_1[i] / (int(Length[i]) / 1000 * mappable2_1 / 1000000)
    dr2_1.append(RPKM)
dr2_2 = []
mappable2_2 = sum(DR2_2)
for i in range(len(DR2_2)):
    RPKM = DR2_2[i] / (int(Length[i]) / 1000 * mappable2_2 / 1000000)
    dr2_2.append(RPKM)
dr2_3 = []
mappable2_3 = sum(DR2_3)
for i in range(len(DR2_3)):
    RPKM = DR2_3[i] / (int(Length[i]) / 1000 * mappable2_3 / 1000000)
    dr2_3.append(RPKM)
dr4_1 = []
mappable4_1 = sum(DR4_1)
for i in range(len(DR4_1)):
    RPKM = DR4_1[i] / (int(Length[i]) / 1000 * mappable4_1 / 1000000)
    dr4_1.append(RPKM)
dr4_2 = []
mappable4_2 = sum(DR4_2)
for i in range(len(DR4_2)):
    RPKM = DR4_2[i] / (int(Length[i]) / 1000 * mappable4_2 / 1000000)
    dr4_2.append(RPKM)
dr4_3 = []
mappable4_3 = sum(DR4_3)
for i in range(len(DR4_3)):
    RPKM = DR4_3[i] / (int(Length[i]) / 1000 * mappable4_3 / 1000000)
    dr4_3.append(RPKM)
dr6_1 = []
mappable6_1 = sum(DR6_1)
for i in range(len(DR6_1)):
    RPKM = DR6_1[i] / (int(Length[i]) / 1000 * mappable6_1 / 1000000)
    dr6_1.append(RPKM)
dr6_2 = []
mappable6_2 = sum(DR6_2)
for i in range(len(DR6_2)):
    RPKM = DR6_2[i] / (int(Length[i]) / 1000 * mappable6_2 / 1000000)
    dr6_2.append(RPKM)
dr6_3 = []
mappable6_3 = sum(DR6_3)
for i in range(len(DR6_3)):
    RPKM = DR6_3[i] / (int(Length[i]) / 1000 * mappable6_3 / 1000000)
    dr6_3.append(RPKM)
dr8_1 = []
mappable8_1 = sum(DR8_1)
for i in range(len(DR8_1)):
    RPKM = DR8_1[i] / (int(Length[i]) / 1000 * mappable8_1 / 1000000)
    dr8_1.append(RPKM)

dr8_2 = []
mappable8_2 = sum(DR8_2)
for i in range(len(DR8_2)):
    RPKM = DR8_2[i] / (int(Length[i]) / 1000 * mappable8_2 / 1000000)
    dr8_2.append(RPKM)
dr8_3 = []
mappable8_3 = sum(DR8_3)
for i in range(len(DR8_3)):
    RPKM = DR8_3[i] / (int(Length[i]) / 1000 * mappable8_3 / 1000000)
    dr8_3.append(RPKM)
# get the information into a new .csv file
exprMatrix = {"Geneid":Geneid,
        "Length":Length,
        "DR0HS1":dr0_1,
        "DR0HS2":dr0_2,
        "DR0HS3":dr0_3,
        "DR2HS1":dr2_1,
        "DR2HS2":dr2_2,
        "DR2HS3":dr2_3,
        "DR4HS1":dr4_1,
        "DR4HS2":dr4_2,
        "DR4HS3":dr4_3,
        "DR6HS1":dr6_1,
        "DR6HS2":dr6_2,
        "DR6HS3":dr6_3,
        "DR8HS1":dr8_1,
        "DR8HS2":dr8_2,
        "DR8HS3":dr8_3,
}
data = DataFrame(exprMatrix)
data.to_csv("ExpressionMatrix.csv",sep=',', header=True, index=True)
