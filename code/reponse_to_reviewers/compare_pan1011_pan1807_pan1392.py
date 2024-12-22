# -*- coding: utf-8 -*-
# date : 2024/11/18 
# author : wangh
# file : compare_pan1011_pan1807_pan1392.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import json
from matplotlib import pyplot as plt


df_1011_ssGEMs=pd.read_csv(r'result/model_simulation/df_pan1011_ssGEMs_size.csv',index_col=0)
df_1392_ssGEMs=pd.read_csv(r'result/model_simulation/df_pan1392_ssGEMs_size.csv',index_col=0)
df_1807_ssGEMs=pd.read_csv(r'result/model_simulation/df_ssGEMs_size.csv',index_col=0)

growth_ratios=[
    len(df_1011_ssGEMs[df_1011_ssGEMs['aerobic_growth']>0.01])/len(df_1011_ssGEMs),
    len(df_1392_ssGEMs[df_1392_ssGEMs['aerobic_growth']>0.01])/len(df_1392_ssGEMs),
    len(df_1807_ssGEMs[df_1807_ssGEMs['aerobic_growth']>0.01])/len(df_1807_ssGEMs)
]

# plot the barplot
plt.bar(['Sce-pan1011','Sce-pan1392','Sce-pan1807'],growth_ratios)
plt.ylabel('Ratio of normal growth strains')
# set y 0-1
plt.ylim([0, 1])
plt.show()



# load gapfilling result
with open(r'code/5.build_ssGEMs/output/pan1011_gapfilling_solutions.json', 'r') as f:
    pan1011_gapfilling_solutions = json.load(f)

with open(r'code/5.build_ssGEMs/output/pan1011_gapfill_failed_strainList.json', 'r') as f:
    pan1011_gapfilling_failed_strainList = json.load(f)
count_freq=dict()
for strain in pan1011_gapfilling_solutions:
    for rxn in pan1011_gapfilling_solutions[strain]:
        if rxn in count_freq:
            count_freq[rxn]+=1
        else:
            count_freq[rxn]=1


df_count=pd.Series(count_freq)

# plot the distribution
df_count[df_count<5].plot(kind='hist')
plt.show()

# check if test strains in pan1011 failed strainList
test_strains=["GCA_019394525.1_ASM1939452v1_genomic.xml","GCA_019394085.1_ASM1939408v1_genomic.xml","GCA_019394815.1_ASM1939481v1_genomic.xml"]
for strain in test_strains:
    if strain in pan1011_gapfilling_failed_strainList:
        print(strain)
