# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：auxofind_t1.py
# 2022/5/4


from cobra.io import read_sbml_model
import pandas as pd
from cobra.flux_analysis import gapfill

with open('result/auxo_pred_unsolved_strains_0505.txt') as f:
    unsolved_strainlist_0=f.readlines()[1].split(',')
unsolved_strainlist_0.remove('')

panYeast=read_sbml_model('model/panYeast_v3.xml')
yeastGEM=read_sbml_model('model/yeastGEM.xml')
ssGEM_dir='model/1011_ssGEMs/'
gapfill_solutions={}

# gapfill function
for strain in unsolved_strainlist_0:
    model=read_sbml_model('%s%s'%(ssGEM_dir,strain))
    with panYeast:
        try:
            solution=gapfill(model,panYeast,exchange_reactions=True,demand_reactions=False,lower_bound=0.04,iterations=3)
        except:
            print('%s can\'t find a solution'%strain)
        else:
            print('%s successfully get solution'%strain)
            gapfill_solutions[strain]=solution

gapfill_results={}
for strain in gapfill_solutions.keys():
    solution=[rxn.name for rxn in gapfill_solutions[strain][0]]
    gapfill_results[strain]=solution







