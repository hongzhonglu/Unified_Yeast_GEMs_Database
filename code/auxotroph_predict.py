# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：auxotroph_predict.py
# 2022/4/27
'''fast predict auxotrophic phenotype for large amounts of strain-specific GEMs'''


import pandas as pd
from cobra.io import read_sbml_model
from cobra import flux_analysis
import os
import sys

sys.path.append('code')
from mainFunction import *


file_path='result/1011_ssGEMs/'
strains_list=os.listdir(file_path)

# 1、minimal medium growth simulation
mini_ungrowable_strainlist=[]
for strain in strains_list:
    model=read_sbml_model(file_path+strain)
    aerobic_value=model.slim_optimize()
    anaerobic_model=anaerobicsimulation(model)
    anaerobic_value=anaerobic_model.slim_optimize()
    if aerobic_value>0.001 and anaerobic_value>0.001:
        print('no auxotroph for %s'%strain)
    else:
        mini_ungrowable_strainlist.append(strain)


# 2、rich medium grow simulation
# get rich medium
nutrients_rxnlist_0=pd.read_excel('data/panYeast_exchangeRXNs.xlsx',sheet_name='auxotroph_analusis')
nutrients_list=nutrients_rxnlist_0['rxn_ID'].values.tolist()
rich_medium=model.medium
for nutrient in nutrients_list:
    rich_medium[nutrient]=1000


# rich medium simulation
rich_growable_strainlist=[]
problem_strainlist=[]       # the strains that might have other auxotroph instead of AA
for strain in mini_ungrowable_strainlist:
    model=read_sbml_model(file_path+strain)
    model.medium=rich_medium
    aerobic_value=model.slim_optimize()
    anaerobic_model=anaerobicsimulation(model)
    anaerobic_value=anaerobic_model.slim_optimize()
    if aerobic_value>0.001 and anaerobic_value>0.001:
        rich_growable_strainlist.append(strain)
    else:
        print('%s might has some specific auxotroph not included in 20 amido acids'%strain)
        problem_strainlist.append(strain)

# find the specific nutrient for auxotrophic phenotype
auxo_strainlist={}
unsolved_strains=[]
for strain in rich_growable_strainlist:
    model=read_sbml_model(file_path+strain)
    minimal_medium=model.medium
    auxo_nutrients=[]
    for nutrient in nutrients_list:
        selected_medium=minimal_medium
        selected_medium[nutrient]=1000
        with model:
            model.medium=selected_medium
            aerobic_value = model.slim_optimize()
            anaerobic_model = anaerobicsimulation(model)
            anaerobic_value = anaerobic_model.slim_optimize()
            if aerobic_value > 0.001 and anaerobic_value > 0.001:
                auxo_nutrients.append(nutrient)
    if len(auxo_nutrients):
        auxo_strainlist[strain]=auxo_nutrients
    else:
        print('%s has some problems')
        unsolved_strains.append(strain)


# save result
auxo_pred_result=pd.DataFrame(dict([(k, pd.Series(v)) for k, v in auxo_strainlist.items()]))
auxo_pred_result.to_csv('result/auxo_pred_result.csv')

with open('result/auxo_pred_unsolved_strains.txt','w') as f:
    f.write('the strainlist that can\'t grow on minimal medium and minimal medium+20 AAs.\n')
    for strain in problem_strainlist:
        f.write(strain+',')





