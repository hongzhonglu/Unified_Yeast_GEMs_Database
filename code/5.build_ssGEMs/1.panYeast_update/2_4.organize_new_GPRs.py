# -*- coding: utf-8 -*-
# date : 2023/5/31 
# author : wangh
# file : 2_4.organize_new_GPRs.py
# project : Unified_Yeast_GEMs_Database
from cobra.io import read_sbml_model
import pandas as pd
from Bio import SeqIO


panYeast=read_sbml_model("model/panYeast_v3.xml")
# 1. update GPR reactions
df_update_gpr=pd.read_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_rxns.xlsx",index_col=0)
# reset index
df_update_gpr=df_update_gpr.reset_index(drop=True)

# organize new GPR string
# 1.1 replace pan1011ID with panID into GPR: updated 31 reactions.
for i in df_update_gpr[df_update_gpr['panID'].isna()==False].index.tolist():
    panID=df_update_gpr.loc[i,'panID']
    pan1011ID=df_update_gpr.loc[i,'pan1011ID']
    old_gpr=df_update_gpr.loc[i,'original_gpr']
    new_gpr=old_gpr.replace(pan1011ID,panID)
    df_update_gpr.loc[i,'new_gpr']=new_gpr
    print(old_gpr+" --> "+new_gpr)

# 1.2 add new isoenzyme into GPR: updated 64 reactions.
for i in df_update_gpr[df_update_gpr['panID'].isna()].index.tolist():
    old_gpr=df_update_gpr.loc[i,'original_gpr']
    new_gpr=old_gpr
    if not 'and' in old_gpr:
        add_panIDlist=df_update_gpr.loc[i,'add_pangene'].split(';')
        for panID in add_panIDlist:
            new_gpr=new_gpr+' or '+panID
        df_update_gpr.loc[i,'new_gpr']=new_gpr
        print(old_gpr+" --> "+new_gpr)

# 1.3 add new isoenzyme into GPR from add new rxns filtering: 2 genes, 3 reactions.
rxnIDlist=['r_0362','r_0959','r_0960']
panIDlist=['scepan0396','scepan1240','scepan1240']
oldgprlist=list()
for i in range(len(rxnIDlist)):
    rxnID=rxnIDlist[i]
    panID=panIDlist[i]
    rxn=panYeast.reactions.get_by_id(rxnID)
    old_gpr=rxn.gene_reaction_rule
    oldgprlist.append(old_gpr)

df_update_gpr_from_new_rxns=pd.DataFrame({'rxnID':rxnIDlist,'panID':panIDlist,'original_gpr':oldgprlist})
# merge df_update_gpr_from_new_rxns into df_update_gpr
df_update_gpr=pd.concat([df_update_gpr,df_update_gpr_from_new_rxns],axis=0,ignore_index=True)
# new GPR organized in excel
df_update_gpr.to_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_rxns.xlsx")