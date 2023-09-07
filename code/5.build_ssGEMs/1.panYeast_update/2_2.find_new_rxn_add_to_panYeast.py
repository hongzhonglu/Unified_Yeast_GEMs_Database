# -*- coding: utf-8 -*-
# date : 2023/5/15 
# author : wangh
# file : add_rxnmet_to_panYeast.py
# project : Unified_Yeast_GEMs_Database
from cobra.io import read_sbml_model
import pandas as pd

panYeast=read_sbml_model("model/panYeast_v3.xml")
panrxnList=[rxn for rxn in panYeast.reactions]
panrxn_keggID=list()
panID_keggID_dict=dict()
for rxn in panrxnList:
    if 'kegg.reaction' in rxn.annotation.keys():
        panrxn_keggID.append(rxn.annotation['kegg.reaction'])
        panID_keggID_dict[rxn.id]=rxn.annotation['kegg.reaction']

df_add_rxn=pd.read_excel("code/5.build_ssGEMs/output/panYeast_add_gene&reaction.xlsx",index_col=0)

# extract GPR update reactions from add_rxn
df_update_gpr_rxns=pd.DataFrame(columns=['rxnID','add_pangene','keggID'])
rxnIDlist=list()
pangeneIDlist=list()
keggIDlist=list()
for panID,keggID in panID_keggID_dict.items():
    if keggID in df_add_rxn.index.tolist():
        df=df_add_rxn.loc[keggID,:]
        if isinstance(df,pd.DataFrame):
            gene=';'.join(df['panID'].tolist())
        elif isinstance(df,pd.Series):
            gene=df['panID']
        rxnIDlist.append(panID)
        pangeneIDlist.append(gene)
        keggIDlist.append(keggID)
df_update_gpr_rxns['rxnID']=rxnIDlist
df_update_gpr_rxns['add_pangene']=pangeneIDlist
df_update_gpr_rxns['keggID']=keggIDlist

rxnNameList=list()
rxnFormulaList=list()
gprList=list()
for rxnID in df_update_gpr_rxns['rxnID']:
    rxn=panYeast.reactions.get_by_id(rxnID)
    rxnName=rxn.name
    rxnFormula=rxn.reaction
    gpr=rxn.gene_reaction_rule
    rxnNameList.append(rxnName)
    rxnFormulaList.append(rxnFormula)
    gprList.append(gpr)
df_update_gpr_rxns['rxnName']=rxnNameList
df_update_gpr_rxns['rxnFormula']=rxnFormulaList
df_update_gpr_rxns['original_gpr']=gprList

df_update_gpr_rxns.to_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_rxns.xlsx")

#remove those update gpr reaction from add rxn
df_add_rxn=df_add_rxn[~df_add_rxn.index.isin(df_update_gpr_rxns['keggID'])]
df_add_rxn.to_excel("code/5.build_ssGEMs/output/panYeast_add_new_rxn.xlsx")


