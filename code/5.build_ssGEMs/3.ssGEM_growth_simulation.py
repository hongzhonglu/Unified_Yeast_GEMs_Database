# -*- coding: utf-8 -*-
import os
import sys
from cobra.io import read_sbml_model
import pandas as pd
from tqdm import tqdm
# from mainFunction import *
import multiprocessing


def ssGEM_simulation(strain,file_path):
    sys.path.append('code')
    from model_modifications import anaerobic_simulation
    print('%s is caculating' % strain)
    model = read_sbml_model(file_path + strain)
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

# sleep 2 hours
# if __name__ == '__main__':
ssGEM_dir="model/pan1800_ssGEMs/"
strains_list=os.listdir(ssGEM_dir)
col=['gene_numb','rxn_numb','aerobic_growth','aerobic_rxns','anaerobic_growth','anaerobic_rxns']
# pool=multiprocessing.Pool(processes=1)
ssGEMs_eval={}
for strain in strains_list:
    ssGEMs_eval[strain]=ssGEM_simulation(strain,ssGEM_dir)
ssGEM_simulation_result=pd.DataFrame(ssGEMs_eval,index=col).T
len(ssGEM_simulation_result[ssGEM_simulation_result['aerobic_growth']>0.04])
ssGEM_simulation_result.to_csv(r"result/model_simulation/df_ssGEMs_size.csv")



