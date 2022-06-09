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

sys.path.append('code')
from model_modifications import *
from mainFunction import *

file_path='model/new600_ssGEMs/'
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
        anaerobic_model=anaerobic_simulation(model)

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

ssGEMs_matrix.to_csv('result/new600_growth_simulation_result.csv')


# result analysis
na1011_growth_simulation_result=pd.read_csv('result/growth_simulation_result.csv')
na1011_growth_simulation_result['aerobic growth']=na1011_growth_simulation_result['aerobic growth'].round(6)
na1011_growth_simulation_result['anaerobic_growth']=na1011_growth_simulation_result['anaerobic_growth'].round(6)
na1011_growth_simulation_result['aerobic growth'].value_counts()
na1011_growth_simulation_result['anaerobic_growth'].value_counts()


lg_growth_simulation_result=pd.read_csv('result/lg1392_growth_simulation_result.csv')
lg_growth_simulation_result['aerobic growth']=lg_growth_simulation_result['aerobic growth'].round(6)
lg_growth_simulation_result['anaerobic_growth']=lg_growth_simulation_result['anaerobic_growth'].round(6)
lg_growth_simulation_result['aerobic growth'].value_counts()
lg_growth_simulation_result['anaerobic_growth'].value_counts()


lg_growth_simulation_result.rename(columns={'Unnamed: 0':'ssGEMs'},inplace=True)
lg_growth_simulation_result['rxn_numb'].describe()

half_growth_strains=lg_growth_simulation_result[(lg_growth_simulation_result['aerobic growth']>0.01)&
                         (lg_growth_simulation_result['aerobic growth']<0.02)]['ssGEMs'].tolist()



halfgrowth_strains_info=pd.DataFrame()
for strain in half_growth_strains:
    strain=strain.replace(".re","")
    strain=strain.split("_")[0]
    strain_info=check_strain_info(strain)
    halfgrowth_strains_info = pd.concat([halfgrowth_strains_info, strain_info], axis=0, join='outer', ignore_index=True)

halfgrowth_strains_info['type'].value_counts()


with open('result/auxo_pred_unsolved_strains_0505.txt') as f:
    unsolved_strainlist_0=f.readlines()[1].split(',')
unsolved_strainlist_0.remove('')
unsolved_strains_size=pd.DataFrame()
for strain in unsolved_strainlist_0:
    growth_result=na1011_growth_simulation_result[na1011_growth_simulation_result['ssGEMs']==strain]
    unsolved_strains_size = pd.concat([unsolved_strains_size, growth_result], axis=0, join='outer', ignore_index=True)
unsolved_strains_size['rxn_numb'].describe()
unsolved_strains_info=pd.DataFrame()
for strain in unsolved_strainlist_0:
    strain_info=check_strain_info(strain)
    unsolved_strains_info = pd.concat([unsolved_strains_info, strain_info], axis=0, join='outer', ignore_index=True)

# p-value calculate
import scipy.stats as stats
import matplotlib.pyplot as plt
t_stat, p_val = stats.ttest_ind(unsolved_strains_size['rxn_numb'], na1011_growth_simulation_result['rxn_numb'], equal_var=False)
p_val


all_strains_info=pd.read_excel('data/1897_strains_info.xlsx')
all_strains_info=all_strains_info[all_strains_info['source']=='1011_nature']
all_classif=all_strains_info['type'].value_counts()
unsolved_classf=unsolved_strains_info['type'].value_counts()
halfgrowth_classf=halfgrowth_strains_info['type'].value_counts()

all_classif.plot.pie(colors=['r','b','g','c','m'])
unsolved_classf.plot.pie(colors=['g','b','r','c','m'])
halfgrowth_classf.plot.pie(colors=['r','b','g','c','m'])
plt.show()


# find common lost rxn
# half growth strain
panYeast=read_sbml_model('model/panYeast_v3.xml')
unsolved_strains_lostrxns={}
pan_rxnList_0=panYeast.reactions
for strain in half_growth_strains:
    model=read_sbml_model('model/1011_ssGEMs/'+
                          strain)
    model_rxnList_0=model.reactions
    diff_rxnList=[]
    i=0
    for rxn in pan_rxnList_0:
        if rxn in model_rxnList_0:
            i+=1
        else:
            diff_rxnList.append(rxn)
    print(strain,':',i)
    diff_rxns_id=[i.id for i in diff_rxnList]
    unsolved_strains_lostrxns[strain]=diff_rxns_id

