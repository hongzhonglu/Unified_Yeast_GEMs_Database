# -*- coding: utf-8 -*-
import pandas as pd
from cobra.flux_analysis import gapfill
from cobra.io import read_sbml_model
import os
import cobra
from tqdm import tqdm


cobra_config = cobra.Configuration()
cobra_config.solver="gurobi"
# cobra_config.solver="cplex"
def strains_gapfilling(template_model,ssGEM_dir,strain_list=list()):
    '''run gapfilling '''
    if not strain_list:
        print("have not provide a specific strain_list,get all strains in the ssGEM_dir")
        strain_list=os.listdir(ssGEM_dir)
    gapfill_solutions={}
    gapfilling_failed_strain=[]
    for strain in tqdm(strain_list):
        try:
            model=read_sbml_model(ssGEM_dir+strain)
        except:
            print("cant't find model for strain:%s"%strain)
            continue
        else:
            if model.slim_optimize()>0.05:
                print("%s could simulate normally"%strain)
            elif len(model.genes) <1000:
                print("%s is unnormal with a too small ssGEM size"%strain)
                gapfilling_failed_strain.append(strain)
            else:
                try:
                    solution = gapfill(model, template_model, exchange_reactions=False, demand_reactions=False, lower_bound=0.04)
                    solution0=solution[0]
                    rxnList=[rxn.id for rxn in solution0]
                except:
                    print('%s can\'t find a solution' %strain)
                    gapfilling_failed_strain.append(strain)
                else:
                    print('%s successfully get solution' %strain)
                    gapfill_solutions[strain] = rxnList
    return gapfilling_failed_strain,gapfill_solutions

# set 5 hours sleep time
# import time
# time.sleep(5*60*60)

panYeast=read_sbml_model("model/panYeast.xml")
ssGEM_dir="model/pan1800_ssGEMs/"
df_ssGEMs_simulate=pd.read_csv("result/model_simulation/df_ssGEMs_size.csv",index_col=0)
gapfill_strainList=df_ssGEMs_simulate[~(df_ssGEMs_simulate["aerobic_growth"]>0.04)].index.tolist()
# gapfill_strainList=os.listdir(ssGEM_dir)
# add .xml suffix
gapfill_strainList=[i+".xml" for i in gapfill_strainList]
gapf_failed_strainList,gapf_solutions=strains_gapfilling(template_model=panYeast,
                                                         ssGEM_dir=ssGEM_dir,
                                                         strain_list=gapfill_strainList)

# save result
import json

with open("code/5.build_ssGEMs/output/pan1800_gapfilling_solutions.json","w") as f:
    json.dump(gapf_solutions,f)

# save failed strain list
with open("code/5.build_ssGEMs/output/pan1800_gapfill_failed_strainList.json","w") as f:
    json.dump(gapf_failed_strainList,f)



