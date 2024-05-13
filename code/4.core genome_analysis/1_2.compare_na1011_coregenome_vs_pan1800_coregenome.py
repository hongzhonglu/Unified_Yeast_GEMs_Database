# -*- coding: utf-8 -*-
# date : 2023/6/3 
# author : wangh
# file : 1_2.compare_na1011_coregenome_vs_pan1800_coregenome.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
from cobra.io import read_sbml_model

# find the core genome of na1011 pangenome
na1011_geneMatrix=pd.read_csv('data/geneMatrix/nature1011_geneMatrix.csv',index_col=0)
coregeneList=na1011_geneMatrix[na1011_geneMatrix.sum(axis=1)==1011].index.tolist()
df_na1011_core=pd.DataFrame(index=na1011_geneMatrix.index.tolist(),columns=["core100"])
df_na1011_core.loc[coregeneList,"core100"]=1
df_na1011_core.fillna(0,inplace=True)

# # save df_na1011_core
# df_na1011_core.to_excel("result/na1011_coregene.xlsx")
panYeast_geneList=[g.id for g in read_sbml_model('model/panYeast_v3.xml').genes]
na1011_core_metabolism_geneList=[g for g in coregeneList if g in panYeast_geneList]


# compare with pan1800 core genome
pan1800_coregenome=pd.read_excel("result/pan1800_coregene.xlsx",index_col=0)
pan1800_all_geneList=pan1800_coregenome.index.tolist()
pan1800_coregeneList=pan1800_coregenome[pan1800_coregenome["core98"]==1].index.tolist()
pan1800_core_metabolism_geneList=[g for g in pan1800_coregeneList if g in panYeast_geneList]

nacore100_nopan1800core=list(set(na1011_core_metabolism_geneList).difference(set(pan1800_core_metabolism_geneList)))
nacore100_nopan1800=list(set(na1011_core_metabolism_geneList).difference(set(pan1800_all_geneList)))







