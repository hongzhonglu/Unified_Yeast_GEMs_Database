# -*- coding: utf-8 -*-
# date : 2023/6/20 
# author : wangh
# file : convert2foldchange.py
# project : Unified_Yeast_GEMs_Database
'''convert absolute expression value to relateive expression value:
option1: s288c as reference
option2: mean value as reference
option3: median value as reference
'''
import matplotlib.pyplot as plt
import pandas as pd
from cobra.io import read_sbml_model
import math
import numpy as np


# load expression data
df_expression=pd.read_csv('code/7.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.csv',index_col=0)
# covert to log2
df_expression_log2=df_expression.applymap(lambda x:math.log2(x+1))

# option1: s288c as reference
df_expression_foldchange1=df_expression_log2.apply(lambda x: x/df_expression_log2['s288c'],axis=0)
# option2: mean value for rows as reference
# calculate mean value without 0 value in each row
df_expression_log2_mean=df_expression_log2.apply(lambda x: x[x!=0].mean(),axis=1)
df_expression_foldchange2=df_expression_log2.apply(lambda x: x/df_expression_log2_mean,axis=0)

# option3: median value as reference
df_expression_foldchange3=df_expression_log2.apply(lambda x: x/df_expression_log2.median(axis=1),axis=0)


# load panYeast
panYeast=read_sbml_model('model/panYeast_v4_5.xml')
geneList=[g.id for g in panYeast.genes]

df_expression_foldchange1=df_expression_foldchange1[df_expression_foldchange1.index.isin(geneList)]
df_expression_foldchange2=df_expression_foldchange2[df_expression_foldchange2.index.isin(geneList)]
df_expression_foldchange3=df_expression_foldchange3[df_expression_foldchange3.index.isin(geneList)]

# fill inf as 10, and nan as 0
df_expression_foldchange1.replace([np.inf, -np.inf], 10,inplace=True)
df_expression_foldchange2.replace([np.inf, -np.inf], 10,inplace=True)
df_expression_foldchange3.replace([np.inf, -np.inf], 10,inplace=True)
df_expression_foldchange1.fillna(0,inplace=True)
df_expression_foldchange2.fillna(0,inplace=True)
df_expression_foldchange3.fillna(0,inplace=True)

# reset index and name the original index as geneName
df_expression_foldchange1.reset_index(inplace=True)
df_expression_foldchange1.rename(columns={'re3':'geneName'},inplace=True)
df_expression_foldchange2.reset_index(inplace=True)
df_expression_foldchange2.rename(columns={'re3':'geneName'},inplace=True)
df_expression_foldchange3.reset_index(inplace=True)
df_expression_foldchange3.rename(columns={'re3':'geneName'},inplace=True)


# save
df_expression_foldchange1.to_excel('code/7.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_foldchange.xlsx',index=False)
df_expression_foldchange2.to_excel('code/7.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_foldchange_mean.xlsx',index=False)
df_expression_foldchange3.to_excel('code/7.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_foldchange_median.xlsx',index=False)

# check the value in df_expression
# check how many values are less than threshold


