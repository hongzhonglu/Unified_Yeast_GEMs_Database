# -*- coding: utf-8 -*-
# date : 2024/12/4 
import pandas as pd
from cobra.io import read_sbml_model,write_sbml_model
import os
import random
import cobra
from cobra.flux_analysis import gapfill

panmodel=read_sbml_model("model/panYeast.xml")
# test_strains
pan1807_ssGEM_dir="model/ssGEMs/"
pan1011_ssGEM_dir=r"E:\data\sce_paper\pan1011_ssGEMs/"
pan1392_ssGEM_dir=r"E:\data\sce_paper\pan1392_ssGEMs/"
test_strains=["GCA_019394525.1_ASM1939452v1_genomic.xml","GCA_019394085.1_ASM1939408v1_genomic.xml","GCA_019394815.1_ASM1939481v1_genomic.xml"]

for strain in test_strains:
    model=read_sbml_model(pan1011_ssGEM_dir+strain)
    try:
        sol=gapfill(model,panmodel,lower_bound=0.01)
        print(sol)
        # model.add_reactions(sol[0])
    except:
        print('%s has trouble in gapfill' %strain)
    print(model.slim_optimize())

# according to gapfilling, r_0061,r_1838 need to be added to the 3 ssGEMs
add_rxnidList=['r_0061','r_1838']
add_rxnList=[panmodel.reactions.get_by_id(id) for id in add_rxnidList]

# 61种不同碳源
carbon_sources_0=pd.read_excel('data/panYeast_exchangeRXNs.xlsx',sheet_name='carbon_sources')
carbon_sources=carbon_sources_0['carbon sources'].values.tolist()
carbob_IDlist=carbon_sources_0['panID'].values.tolist()


growth_simulation=pd.DataFrame(index=carbon_sources)
# predict 3 ssGEMs
for strain in test_strains:
    model=read_sbml_model(pan1011_ssGEM_dir+strain)
    model.add_reactions(add_rxnList)
    growth_values=[]
    # model = anerobicsimulation(model)
    print(model.slim_optimize())
    for carbon in carbob_IDlist:
        with model:
            medium = model.medium
            medium['r_1714'] = 0.0  # 移除原培养基中的glucose
            medium[carbon]=1.0
            model.medium=medium
            growth_value=model.slim_optimize()
            growth_values.append(growth_value)
    growth_simulation[strain]=growth_values
strain_id_list=[i.strip(".xml") for i in test_strains]
