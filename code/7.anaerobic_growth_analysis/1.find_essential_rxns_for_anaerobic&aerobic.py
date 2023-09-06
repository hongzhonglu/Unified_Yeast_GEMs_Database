# -*- coding: utf-8 -*-
# date : 2023/6/8 
# author : wangh
# file : 1.find_essential_rxns_for_anaerobic&aerobic.py
# project : Unified_Yeast_GEMs_Database
# use yeast8 & panYeast tofind essential rxns for anaerobic & aerobic condition respectively
from cobra.io import read_sbml_model
from cobra.flux_analysis import single_reaction_deletion
import pandas as pd
import os
import sys
sys.path.append('code')

# load model
yeast8=read_sbml_model("model/yeastGEM.xml")
panYeast=read_sbml_model("model/panYeast_v4_5.xml")

# anaerobic condition simulation
from model_modifications import anaerobic_simulation

anaerobic_yeast8=anaerobic_simulation(yeast8)
anaerobic_panYeast=anaerobic_simulation(panYeast)

# single reactions deletion in minimal media
df_anaerobic_single_rxn_deletion_panYeast=single_reaction_deletion(anaerobic_panYeast)

max_anaerobic_growth=df_anaerobic_single_rxn_deletion_panYeast['growth'].max()
len(df_anaerobic_single_rxn_deletion_panYeast[df_anaerobic_single_rxn_deletion_panYeast['growth']<max_anaerobic_growth*0.5])
anaerobic_essential_rxnList=df_anaerobic_single_rxn_deletion_panYeast[df_anaerobic_single_rxn_deletion_panYeast['growth']<max_anaerobic_growth*0.5]['ids'].astype(str).tolist()

# aerobic condition simulation
df_aerobic_single_rxn_deletion_panYeast=single_reaction_deletion(panYeast)
max_aerobic_growth=df_aerobic_single_rxn_deletion_panYeast['growth'].max()
len(df_aerobic_single_rxn_deletion_panYeast[df_aerobic_single_rxn_deletion_panYeast['growth']<max_aerobic_growth*0.5])
aerobic_essential_rxnList=df_aerobic_single_rxn_deletion_panYeast[df_aerobic_single_rxn_deletion_panYeast['growth']<max_aerobic_growth*0.5]['ids'].astype(str).tolist()

# remove '{\' and \'}' in rxn ids
import re
# extract rxn ids from list
id_pattern=re.compile(r'r_\d{4}')
aerobic_essential_rxnList=[id_pattern.findall(i)[0] for i in aerobic_essential_rxnList]
anaerobic_essential_rxnList=[id_pattern.findall(i)[0] for i in anaerobic_essential_rxnList]


# find aerobic specific essential rxns and anaerobic specific essential rxns
aerobic_specific_essential_rxnList=list(set(aerobic_essential_rxnList)-set(anaerobic_essential_rxnList))
anaerobic_specific_essential_rxnList=list(set(anaerobic_essential_rxnList)-set(aerobic_essential_rxnList))

allrxnList=[rxn.id for rxn in panYeast.reactions]
df_aerobic_anaerobic_essential_rxnList=pd.DataFrame(index=allrxnList,columns=['aerobic_essential','anaerobic_essential'])
df_aerobic_anaerobic_essential_rxnList.loc[aerobic_essential_rxnList,'aerobic_essential']=1
df_aerobic_anaerobic_essential_rxnList.loc[anaerobic_essential_rxnList,'anaerobic_essential']=1
df_aerobic_anaerobic_essential_rxnList.fillna(0,inplace=True)
df_aerobic_anaerobic_essential_rxnList.to_excel('result/ssGEM_simulation/panYeast_anaerobic&aerobic_essential_rxns.xlsx')


# find essential rxns related genes
aerobic_related_genes=[]
for rxnID in aerobic_essential_rxnList:
    rxn=panYeast.reactions.get_by_id(rxnID)
    for gene in rxn.genes:
        aerobic_related_genes.append(gene.id)
aerobic_related_genes=list(set(aerobic_related_genes))

anaerobic_related_genes=[]
for rxnID in anaerobic_essential_rxnList:
    rxn=panYeast.reactions.get_by_id(rxnID)
    for gene in rxn.genes:
        anaerobic_related_genes.append(gene.id)
anaerobic_related_genes=list(set(anaerobic_related_genes))

# check how many aerobic specifically related genes and anaerobic specifically related genes
aerobic_specific_related_genes=list(set(aerobic_related_genes)-set(anaerobic_related_genes))
anaerobic_specific_related_genes=list(set(anaerobic_related_genes)-set(aerobic_related_genes))

allgeneList=[gene.id for gene in panYeast.genes]
df_aerobic_anaerobic_essential_geneList=pd.DataFrame(index=allgeneList,columns=['aerobic_related','anaerobic_related'])
df_aerobic_anaerobic_essential_geneList.loc[aerobic_related_genes,'aerobic_related']=1
df_aerobic_anaerobic_essential_geneList.loc[anaerobic_related_genes,'anaerobic_related']=1
df_aerobic_anaerobic_essential_geneList.fillna(0,inplace=True)
df_aerobic_anaerobic_essential_geneList.to_excel('result/ssGEM_simulation/panYeast_anaerobic&aerobic_related_genes.xlsx')



