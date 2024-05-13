# -*- coding: utf-8 -*-
# date : 2023/6/8 
# author : wangh
# file : 1.check_trancriptome_data.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import tqdm
from Bio import SeqIO

# load transcriptome data with .tab format
df_test=pd.read_csv('data/transcriptomics/test_data')

import chardet
# 检测文件文字编码格式
with open("data/transcriptomics/test_data", 'rb') as f:
    data = f.read()
    print(chardet.detect(data))

with open('data/transcriptomics/final_data_annotated_merged_04052022.tab', encoding='utf-8',errors='ignore') as f:
    # dataset=pd.read_csv(f,error_bad_lines=False)
    dataset=pd.read_csv(f)

# collect each gene family's TPM/count value in each strains
allstrainList=list(set(dataset['Strain']))
allgeneList=list(set(dataset['systematic_name'].dropna()))


# for each gene family, collect its TPM value in each strains
df_genetpmMatrix=pd.DataFrame(index=allgeneList,columns=allstrainList)
for gene in allgeneList:
    print(gene)
    df_gene_data=dataset[dataset['systematic_name']==gene]
    # set strain as index
    df_gene_data=df_gene_data.set_index('Strain')
    df_tpm=pd.Series(index=allstrainList)
    # fill the tpm value from df_gene_data to df_tpm according to the index
    # df_tpm.loc[df_gene_data.index]=df_gene_data['tpm']
    df_tpm.loc[df_gene_data.index]=df_gene_data['count']
    df_genetpmMatrix.loc[gene]=df_tpm

df_genetpmMatrix.fillna(0,inplace=True)

# rename the gene name in allgeneList
df_rename_genes=pd.DataFrame(index=allgeneList)
# round 1, remove X in the start of gene name
df_rename_genes['re1']=[g.lstrip('X') for g in allgeneList]
# round 2, if the gene name sart with number, replace . with - in the gene name
df_rename_genes['re2']=[g.replace('.','-') if g[0].isdigit() else g for g in df_rename_genes['re1']]

# check weather all gene id have been renamed through check if all gene exist in pan1011
pan1011_geneList=[g.id for g in SeqIO.parse("data/genome/na_pan1011.fasta","fasta")]
len([g for g in df_rename_genes['re2'] if g not in pan1011_geneList])   # 254 genes not include in pan1011 pan-genome

# round 3, repalce geneID with panID
df_panID_geneID=pd.read_excel("result/panID_pan1011.xlsx",index_col=0)
# if the gene in df_rename_genes['re2'] exist in df_panID_geneID['pan1011_geneID'], replace it with panID
panID_pan1011_geneList=df_panID_geneID['pan1011_geneID'].values.tolist()
re3_geneList=[]
for g in df_rename_genes['re2'].tolist():
    if g in panID_pan1011_geneList:
        panID=df_panID_geneID[df_panID_geneID['pan1011_geneID']==g].index[0]
        print(g,panID)
        re3_geneList.append(panID)
    else:
        re3_geneList.append(g)
df_rename_genes['re3']=re3_geneList
df_genetpmMatrix.index=df_rename_genes['re3']


# rename the strain name in allstrainList
all_strain_info=pd.read_excel("data/1897_strains_info.xlsx",index_col=0)
df_reanme_strain=pd.DataFrame(index=allstrainList)
# round 1 , remove SACE_ in the start of strain name
df_reanme_strain['re1']=[s.replace("SACE_","") for s in allstrainList]

# round 2  ,replace the strain name in re1 with genome_id in all_strain_info
all_strain_info['genome_id']=all_strain_info.index
# modify index by spliting the strain name by _
all_strain_info.index=[str(s).split('_')[0].replace('.re','') for s in all_strain_info.index]

df_reanme_strain['re2']=[all_strain_info.loc[s]['genome_id'] if s in all_strain_info.index else s for s in df_reanme_strain['re1']]
# count how many rows which re1 and re2 are different
df_genetpmMatrix.columns=df_reanme_strain['re2']


# save the df_genetpmMatrix
# df_genetpmMatrix.to_csv("code/7.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.csv")
df_genetpmMatrix.to_csv("code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix.csv")


