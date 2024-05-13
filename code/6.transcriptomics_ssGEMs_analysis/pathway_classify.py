# -*- coding: utf-8 -*-
# date : 2023/7/5 
# author : wangh
# file : pathway_classify.py
# project : Unified_Yeast_GEMs_Database
'''find all reactions genes involved in different metabolic pathways.
Divide glycolysis pathway into upstream and downstream.
Get the CNV and TPM information for each reaction according to the GPR rules. max for OR, min for AND'''
import pandas as pd
import os
from cobra.io import read_sbml_model,write_sbml_model


# load panModel
panYeast=read_sbml_model('model/panYeast_v4_5.xml')
pathwaylist=['Glycolysis','TCA cycle','Pentose phosphate pathway','Galactose','Oxidative phosphorylation']

# find pathway reactions according Yeast8 data by chengyu
yeast8_data=pd.read_excel('model/yeast-GEM_zcy.xlsx',sheet_name='reactions')
# set ID as index
yeast8_data=yeast8_data.set_index('ID')
# fill NA with 'unknown'
yeast8_data=yeast8_data.fillna('unknown')

pathway_rxndict=dict()
for pathway in pathwaylist:
    pathway_rxndict[pathway]=list(yeast8_data[yeast8_data['subsystem'].str.contains(pathway)].index)


# divide into upstream and downstream of glycolysis from glucose pathway

glycolysis_upstream=['r_4288','r_4285','r_4284','r_0322','r_0449','r_0450','r_0467','r_0886']
glycolysis_downstream=['r_0356','r_0486','r_0534','r_0884','r_0892','r_0893','r_0962','r_1054','r_0366']
anaerobic_fermentation=['r_0165','r_0173','r_0174','r_2115','r_2116','r_0163']

pathway_rxndict['Glycolysis_upstream']=glycolysis_upstream
pathway_rxndict['Glycolysis_downstream']=glycolysis_downstream
pathway_rxndict['Anaerobic_fermentation']=anaerobic_fermentation


# save the pathway_rxndict as dict
import json
with open('model/model_pathway_rxndict.json','w') as f:
    json.dump(pathway_rxndict,f)


# load the pathway_rxndict
with open('model/model_pathway_rxndict.json','r') as f:
    pathway_rxndict=json.load(f)

# add Biosynthesis of secondary metabolites and Biosynthesis of amino acids
add_pathways=['Biosynthesis of secondary metabolites','Biosynthesis of amino acids']
for pathway in add_pathways:
    rxnList=[rxn.id for rxn in panYeast.reactions if pathway in rxn.subsystem]
    pathway_rxndict[pathway]=rxnList

