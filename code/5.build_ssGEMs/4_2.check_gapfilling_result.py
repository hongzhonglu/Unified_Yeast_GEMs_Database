# -*- coding: utf-8 -*-
# date : 2023/6/5 
# author : wangh
# file : 4_2.check_gapfilling_result.py
# project : Unified_Yeast_GEMs_Database
# check the gapfilling unsolved strains and the solved gapfilling solutions
import json
import pandas as pd
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model
import seaborn as sns

#load gaofilling result
with open(r'code/5.build_ssGEMs/output/pan1800_gapfilling_solutions.json','r') as f:
    pan1800_gapfilling_sol=json.load(f)
with open(r'code/5.build_ssGEMs/output/pan1800_gapfill_failed_strainList.json','r') as f:
    pan1800_gapfilling_failed_strainList=json.load(f)

gapfilling_unsolved=[x.strip('.xml') for x in pan1800_gapfilling_failed_strainList]

# load all strains information
all_strain_info=pd.read_excel(r'data\1897_strains_info.xlsx',index_col=0)
gapfilling_unsolved_info=all_strain_info[all_strain_info.index.isin(gapfilling_unsolved)]


# 1. check the unsolved strains
df_combine_type_count=all_strain_info['type'].value_counts()/len(all_strain_info)
df_combine_type_count=pd.DataFrame(df_combine_type_count).reset_index()
df_combine_type_count.columns=['type','relative ratio']
df_combine_type_count['hue']='all 1897 strains'
# add gapfilling unsolved info to the df
df_unsolved_type_count=gapfilling_unsolved_info['type'].value_counts()/len(gapfilling_unsolved_info)
df_unsolved_type_count=pd.DataFrame(df_unsolved_type_count).reset_index()
df_unsolved_type_count.columns=['type','relative ratio']
df_unsolved_type_count['hue']='gapfilling unsolved strains'
df_combine_type_count=df_combine_type_count.append(df_unsolved_type_count)

# check assemble level
df_combine_assemble_count=all_strain_info['assemble_level'].value_counts()/len(all_strain_info)
df_combine_assemble_count=pd.DataFrame(df_combine_assemble_count).reset_index()
df_combine_assemble_count.columns=['assemble_level','relative ratio']
df_combine_assemble_count['hue']='all 1897 strains'
# add gapfilling unsolved info to the df
df_unsolved_assemble_count=gapfilling_unsolved_info['assemble_level'].value_counts()/len(gapfilling_unsolved_info)
df_unsolved_assemble_count=pd.DataFrame(df_unsolved_assemble_count).reset_index()
df_unsolved_assemble_count.columns=['assemble_level','relative ratio']
df_unsolved_assemble_count['hue']='gapfilling unsolved strains'
df_combine_assemble_count=df_combine_assemble_count.append(df_unsolved_assemble_count)

# check source distribution
df_combine_source_count=all_strain_info['source'].value_counts()/len(all_strain_info)
df_combine_source_count=pd.DataFrame(df_combine_source_count).reset_index()
df_combine_source_count.columns=['source','relative ratio']
df_combine_source_count['hue']='all 1897 strains'
# add gapfilling unsolved info to the df
df_unsolved_source_count=gapfilling_unsolved_info['source'].value_counts()/len(gapfilling_unsolved_info)
df_unsolved_source_count=pd.DataFrame(df_unsolved_source_count).reset_index()
df_unsolved_source_count.columns=['source','relative ratio']
df_unsolved_source_count['hue']='gapfilling unsolved strains'
df_combine_source_count=df_combine_source_count.append(df_unsolved_source_count)


fig,ax=plt.subplots(1,3,figsize=(15,4))
sns.barplot(x='type',y='relative ratio',hue='hue',data=df_combine_type_count,ax=ax[0])
sns.barplot(x='assemble_level',y='relative ratio',hue='hue',data=df_combine_assemble_count,ax=ax[1])
sns.barplot(x='source',y='relative ratio',hue='hue',data=df_combine_source_count,ax=ax[2])
# rotate the xticklabels
for i in range(3):
    for tick in ax[i].get_xticklabels():
        tick.set_rotation(30)
plt.tight_layout()
plt.show()

# check number of rxns in gapfilling solutions
gapsol_rxnnumbs={}
for strain in pan1800_gapfilling_sol:
    gapsol_rxnnumbs[strain]=len(pan1800_gapfilling_sol[strain])
df_count=pd.Series(gapsol_rxnnumbs).value_counts()
# plot the barplot
plt.bar(df_count.index,df_count)
plt.show()
# save df_count


# check the gapfilling solutions
def count_rxns_in_sol(sol_dict):
    '''count the rxn frequency in the gapfilling solutions'''
    rxns=[]
    for i in sol_dict.values():
        rxns.extend(i)
    rxns=pd.Series(rxns)
    return rxns.value_counts()

df_gapfilling_rxn_count=count_rxns_in_sol(pan1800_gapfilling_sol)
df_gapfilling_rxn_count=pd.DataFrame(df_gapfilling_rxn_count)

panYeast=read_sbml_model('model/panYeast.xml')
rxnNamelist=list()
gprList=list()
for id in df_gapfilling_rxn_count.index.tolist():
    if 'DM' not in id:
        rxnNamelist.append(panYeast.reactions.get_by_id(id).name)
        gprList.append(panYeast.reactions.get_by_id(id).gene_reaction_rule)
    else:
        rxnNamelist.append(id)
        gprList.append('')

df_gapfilling_rxn_count['rxnName']=rxnNamelist
df_gapfilling_rxn_count['gpr']=gprList
df_gapfilling_rxn_count.rename(columns={0:'strain number'},inplace=True)

# save result
df_gapfilling_rxn_count.to_csv(r'code/5.build_ssGEMs/output/gapfilling_rxn_information.csv')


# check the strain that lost r_1838
r_1838_lost_strainList=[]
for strain,list in pan1800_gapfilling_sol.items():
    if 'r_1838' in list:
        r_1838_lost_strainList.append(strain.rstrip('.xml'))
df_r1838_lost_strain=all_strain_info[all_strain_info.index.isin(r_1838_lost_strainList)]