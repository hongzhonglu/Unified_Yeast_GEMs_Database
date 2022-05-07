# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：growthsimulation.py
# 2022/4/26

'''evaluate the growth simulation condition for ssGEM on glucose basic medium'''

import os
import sys
from cobra.io import read_sbml_model
import pandas as pd

sys.path.append('code/model_modification')
import anerobic_change

file_path='result/1011_ssGEMs/'
strains_list=os.listdir(file_path)
col=['gene_numb','rxn_numb','aerobic growth','aerobic_rxns','anaerobic_growth','anaerobic_rxns']
ssGEMs_eval= {}
for strain in strains_list:
    print('%s is caculating'%strain)
    model=read_sbml_model(file_path+strain)

    # 有氧生长表型
    aerobic_growth_value = model.slim_optimize()

    # 统计有氧激活的反应
    aerobic_active_rxn_count=0
    for x in model.reactions:
        if x.flux!=0:
            aerobic_active_rxn_count+=1

    # 统计无氧生长表型和激活反应数量
    try:
        anaerobic_model=anerobic_change.anaerobicsimulation(model)

        anaerobic_growth_value = anaerobic_model.slim_optimize()

        anaerobic_active_rxn_count = 0
        for x in anaerobic_model.reactions:
            if x.flux != 0:
                anaerobic_active_rxn_count += 1
        anaerobic_growth=[anaerobic_growth_value,anaerobic_active_rxn_count]
    except:
        print('%s has trouble in anerobicsimulation' %strain)
        anaerobic_growth=[0,0]


    ssGEMs_eval[strain] = str(len(model.genes)), str(len(model.reactions)),aerobic_growth_value,aerobic_active_rxn_count\
        ,anaerobic_growth[0],anaerobic_growth[1]

ssGEMs_matrix = pd.DataFrame(ssGEMs_eval, index=col).transpose()

ssGEMs_matrix.to_csv('result/growth_simulation_result.csv')

