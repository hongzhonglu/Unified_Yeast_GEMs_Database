# -*- coding: utf-8 -*-
# date : 2023/7/12 
# author : wangh
# file : 1_anaerobic_vs_aerobic_simulation.py
# project : Unified_Yeast_GEMs_Database
'''comparion the flux distribution between wildtype strains which grow in aerobic condition and bioethanol strains which grow in anaerobic condition.
1. Simulate the wild type fluxes distribution under aerobic condition by pFBA.
2. Simulate the 3.Brazilian bioethanol strains fluxes distribution under anaerobic condition by pFBA.
'''
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
import pandas as pd
import os
from model_modifications import anaerobic_simulation
import tqdm


# load strain information data
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)

ssGEM_dir='model/pan1800_tblastn4_coregene_ssGEMs/'

# set medium: mimimal medium with glucose as only carbon source
gluc_uptake=-5
gluc_rxn='r_1714'

growth_id='r_2111'

# load panYeast
panYeast=read_sbml_model('model/panYeast_v4_5.xml')
rxnList=[i.id for i in panYeast.reactions]

# 1. Simulate the wild type fluxes distribution under aerobic condition by pFBA.
# extract wildtype strainList
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[df_strain_info['nature_clade'].isin(wildtypelist)].index.tolist()
wt_strainList=[i+'.xml' for i in wt_strainList]

df_wt_flux=pd.DataFrame(index=['growth']+rxnList)

# set wild type ethanol production rate
wt_ethanol=gluc_uptake*-0.2
ethanol_exchange='r_1761'

for strain in wt_strainList[:1]:
    model=read_sbml_model(ssGEM_dir+strain)
    model.reactions.get_by_id(gluc_rxn).lower_bound=gluc_uptake
    model.reactions.get_by_id(ethanol_exchange).lower_bound=wt_ethanol
    # block rxn r_4235
    model.reactions.get_by_id('r_4235').bounds=(0,0)
    model.objective=growth_id
    sol=pfba(model,fraction_of_optimum=0.9)
    gr=sol.fluxes[growth_id]
    df_wt_flux[strain]=sol.fluxes
    df_wt_flux.loc['growth',strain]=gr

# save result
df_wt_flux.to_csv('code/7.anaerobic_growth_analysis/output/wt_pfba_fluxes.csv')


# 2. Simulate the 3.Brazilian bioethanol strains fluxes distribution under anaerobic condition by pFBA.
bioethanol_strainList=df_strain_info[df_strain_info['nature_clade']=='3. Brazilian bioethanol '].index.tolist()
bioethanol_strainList=[i+'.xml' for i in bioethanol_strainList]

df_bioethanol_flux=pd.DataFrame(index=['growth']+rxnList)

# set wild type ethanol production rate
bioethanol_ethanol=gluc_uptake*-0.9
ethanol_exchange='r_1761'

for strain in tqdm.tqdm(bioethanol_strainList[:1]):
    model=read_sbml_model(ssGEM_dir+strain)
    model.reactions.get_by_id(gluc_rxn).lower_bound=gluc_uptake
    model.objective=growth_id
    model.reactions.get_by_id(ethanol_exchange).lower_bound=bioethanol_ethanol
    # model=anaerobic_simulation(model)
    sol=pfba(model,fraction_of_optimum=0.9)
    gr=sol.fluxes[growth_id]
    df_bioethanol_flux[strain]=sol.fluxes
    df_bioethanol_flux.loc['growth',strain]=gr

# check tca flux
# citrate：s_0524[m]
print(model.metabolites.get_by_id('s_0524[m]').summary())
# isocitrite:s_0941[m]
print(model.metabolites.get_by_id('s_0941[m]').summary())
# 2-oxoglutarate:s_0182[m]
print(model.metabolites.get_by_id('s_0182[m]').summary())
# succinate: s_1460[m]
print(model.metabolites.get_by_id('s_1460[m]').summary())
# fumarate:s_0727[m]
print(model.metabolites.get_by_id('s_0727[m]').summary())
# malate: s_0068[m]
print(model.metabolites.get_by_id('s_0068[m]').summary())

# save result
df_bioethanol_flux.to_csv('code/7.anaerobic_growth_analysis/output/bioethanol_pfba_fluxes_v2.csv')


