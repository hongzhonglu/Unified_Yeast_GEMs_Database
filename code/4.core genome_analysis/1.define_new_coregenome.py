# -*- coding: utf-8 -*-
# date : 2023/5/11 
# author : wangh
# file : 1.define_new_coregenome.py
# project : Unified_Yeast_GEMs_Database
'''define core genome of pan1900_v4 genome
'''
from Bio import SeqIO
import pandas as pd
import os
import numpy as np

# load geneMatrix
geneMatrix = pd.read_csv('data/geneMatrix/pan1800_v2_blastp_50_70_cnvMatrix.csv',index_col=0)

# load all strain data
all_strain_info=pd.read_excel("data/1897_strains_info.xlsx", index_col=0)
kept_strainList=all_strain_info[all_strain_info["remove"]==False]["genome_id"].tolist()
kept_strainList=[s+".fa" for s in kept_strainList]
kept_strainList.append('s288c_R64.fa')


strainList=[s for s in geneMatrix.columns if s in kept_strainList]
geneMatrix2=geneMatrix.loc[:,strainList]
gene_ratio2=geneMatrix2.sum(axis=1)/geneMatrix2.shape[1]
len(gene_ratio2[gene_ratio2==1])
geneMatrix2.sum().describe()

core99_idlist=gene_ratio2[gene_ratio2>0.99].index.tolist()
core98_idlist=gene_ratio2[gene_ratio2>0.98].index.tolist()

df_core=pd.DataFrame(index=geneMatrix2.index.tolist(),columns=["core99","core98"])
df_core.loc[core99_idlist,"core99"]=1
df_core.loc[core98_idlist,"core98"]=1
df_core.fillna(0,inplace=True)
df_core.to_excel("result/pan1800_coregene.xlsx")





