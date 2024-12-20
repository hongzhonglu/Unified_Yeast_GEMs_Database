# -*- coding: utf-8 -*-
# date : 2024/12/6 
import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database/code')
from cobra.io import read_sbml_model,write_sbml_model
from model_modifications import set_SCmedium
import tqdm


def ssGEM_simulation(strain,file_path):
    sys.path.append('code')
    from model_modifications import anaerobic_simulation
    print('%s is caculating' % strain)
    model = read_sbml_model(file_path + strain)
    # set SCmedium
    model=set_SCmedium(model)
    # 有氧生长表型
    aerobic_growth_value = model.slim_optimize()
    # 统计有氧激活的反应
    aerobic_active_rxn_count = 0
    for x in model.reactions:
        if x.flux != 0:
            aerobic_active_rxn_count += 1
    # anaerobic growth phenotype and active reaction number
    try:
        anaerobic_model = anaerobic_simulation(model)
        anaerobic_growth_value = anaerobic_model.slim_optimize()

        anaerobic_active_rxn_count = 0
        for x in anaerobic_model.reactions:
            if x.flux != 0:
                anaerobic_active_rxn_count += 1
        anaerobic_growth = [anaerobic_growth_value, anaerobic_active_rxn_count]
    except:
        print('%s has trouble in anerobicsimulation' % strain)
        anaerobic_growth = [0, 0]
    result=[str(len(model.genes)), str(len(model.reactions)),aerobic_growth_value,aerobic_active_rxn_count\
            ,anaerobic_growth[0],anaerobic_growth[1]]
    return result



ssGEM_dir=r"E:\data\sce_paper\pan1392_ssGEMs/"
strains_list=os.listdir(ssGEM_dir)
col=['gene_numb','rxn_numb','aerobic_growth','aerobic_rxns','anaerobic_growth','anaerobic_rxns']
# pool=multiprocessing.Pool(processes=1)
ssGEMs_eval={}
for strain in tqdm.tqdm(strains_list):
    ssGEMs_eval[strain]=ssGEM_simulation(strain,ssGEM_dir)
ssGEM_simulation_result=pd.DataFrame(ssGEMs_eval,index=col).T
len(ssGEM_simulation_result[ssGEM_simulation_result['aerobic_growth']>0.1])
ssGEM_simulation_result.to_csv(r"result/model_simulation/df_pan1392_ssGEMs_SC_medium_size.csv")

