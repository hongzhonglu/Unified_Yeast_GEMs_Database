# -*- coding: utf-8 -*-
# date : 2024/11/2 
# author : wangh
# file : build_ssGEMs_by_orther_pan.py
# project : Unified_Yeast_GEMs_Database
import os
import sys
from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd
from cobra.flux_analysis import gapfill
import tqdm

sys.path.append("code")
from mainFunction import *
pan_model=read_sbml_model('model/panYeast.xml')

#produce the dataframe for the metabolites and the rxn based on panYeast
pan_gene = produceGeneList(pan_model)
# read geneMatrix
geneMatrix = pd.read_csv('data/geneMatrix/lg_pan1392_blastp_geneMatrix.csv',index_col=0)
output_dir="model/pan1392_ssGEMs/"



df_coregeneMatrix=pd.read_csv("code/5.build_ssGEMs/output/coregene_tblastn_geneMatrix.csv",index_col=0)
# change columns by repalace .fna to .fa
df_coregeneMatrix.columns=[s.replace(".fna",".fa") for s in df_coregeneMatrix.columns]


# get all strain list
all_strain_info=pd.read_excel("data/1897_strains_info.xlsx", index_col=0)
kept_strainList=all_strain_info[all_strain_info["remove"]==False].index.tolist()
kept_strainList=[s+".fa" for s in kept_strainList]
kept_strainList.append('s288c_R64.fa')
strainList=[s for s in geneMatrix.columns if s in kept_strainList]
geneMatrix2=geneMatrix.loc[:,strainList]


geneMatrix_core100=geneMatrix.copy()
for g in df_coregeneMatrix.index:
    if g in geneMatrix_core100.index.tolist():
        geneMatrix_core100.loc[g,:]=df_coregeneMatrix.loc[g,:]+geneMatrix_core100.loc[g,:]
    else:
        geneMatrix_core100.loc[g,:]=df_coregeneMatrix.loc[g,:]
        print(g)
geneMatrix_core100['s288c_R64.fa'].fillna(1,inplace=True)
# all more than 1 are 1
geneMatrix_core100[geneMatrix_core100>1]=1

# only keep kept_strainList strain columns
geneMatrix_core100=geneMatrix_core100.loc[:,strainList]


# if not exist create it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
# build strain-specific GEMs
geneMatrix2['geneID'] = geneMatrix2.index
geneMatrix_core100['geneID'] = geneMatrix_core100.index
count_normally_growable=0
count_half_growable=0
# # randomly sample 200 strains from the strainList
# import random
# random.seed(2021)
# strainList=random.sample(strainList,100)

add_rxnList=[pan_model.reactions.get_by_id(id) for id in ["r_1603","r_0226","r_0439","r_0438"]]

exist_ssGEMs=os.listdir(output_dir)
for strain in tqdm.tqdm(strainList):
    ssGEM_name=strain.rstrip(".fa")+'.xml'
    if ssGEM_name in exist_ssGEMs:
        continue
    # print("%s ssGEM is building"%strain)
    NEW = getStrainGEM(s0=strain, geneMatrix0=geneMatrix_core100, templateGEM=pan_model, templateGene= pan_gene)
    # NEW = getStrainGEM(s0=strain, geneMatrix0=geneMatrix2, templateGEM=pan_model, templateGene= pan_gene)
    NEW.id=strain
    NEW.description='S.cerevisiae_strain_%s_strain-specific_genome-scale_metabolic_model'%strain
    NEW.name=ssGEM_name
    # add rxns:
    NEW.add_reactions(add_rxnList)
    if NEW.slim_optimize()>0.04:
        count_normally_growable+=1
        print(NEW.slim_optimize())
    elif NEW.slim_optimize()>0:
        count_half_growable+=1
        print(NEW.slim_optimize())

    write_sbml_model(NEW, output_dir + ssGEM_name)
