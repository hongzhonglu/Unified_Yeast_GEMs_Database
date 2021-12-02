import os
import re
import pandas as pd
import numpy as np
# set the directory
os.chdir('/Users/xluhon/Documents/GitHub/Unified_Yeast_GEMs_Database/code')
os.getcwd()
genesMatrix = pd.read_excel("../data/genesMatrix.xlsx")

#sample one
submatrix = genesMatrix.sample(1)
submatrix0 = submatrix.drop(['strain_name'], axis=1)
s1 = submatrix0.sum(axis=1)

#sample
def calculateCoreNum(sample, matrix=genesMatrix):
    sampelNum = sample
    submatrix5 = matrix.sample(sampelNum)
    submatrix50 =  submatrix5.drop(['strain_name'], axis=1)
    submatrix5t = submatrix50.T
    submatrix5t["sum"] = submatrix5t.sum(axis=1)
    ss = submatrix5t["sum"].tolist()
    coreGene = [x for x in ss if x == sampelNum]
    panGene =  [x for x in ss if x >= 1]
    coreNum = len(coreGene)
    panNum = len(panGene)
    return(coreNum)

def calculatePanNum(sample, matrix=genesMatrix):
    sampelNum = sample
    submatrix5 = matrix.sample(sampelNum)
    submatrix50 =  submatrix5.drop(['strain_name'], axis=1)
    submatrix5t = submatrix50.T
    submatrix5t["sum"] = submatrix5t.sum(axis=1)
    ss = submatrix5t["sum"].tolist()
    coreGene = [x for x in ss if x == sampelNum]
    panGene =  [x for x in ss if x >= 1]
    panNum = len(panGene)
    return(panNum)

def calculateAccessNum(sample, matrix=genesMatrix):
    ss = calculatePanNum(sample) - calculateCoreNum(sample)
    return(ss)



#loop to calculate the core and pangene number
coregene = list()
for i in range(1,1012):
    print(i)
    s = calculateCoreNum(i)
    print(s)
    coregene.append(s)


pangene = list()
for i in range(1,1012):
    print(i)
    s = calculatePanNum(i)
    print(s)
    pangene.append(s)

acessgene = list()
for i in range(1,1012):
    print(i)
    s = calculateAccessNum(i)
    print(s)
    acessgene.append(s)

strain_num = [i for i in range(1,1012)]
result=pd.DataFrame({'num':strain_num})
result['coregene'] = coregene
result['pangene'] = pangene
result['accessory_gene'] = acessgene
#save the data
writer = pd.ExcelWriter('../result/pangene and coregene.xlsx')
result.to_excel(writer,'Sheet1')
writer.save()