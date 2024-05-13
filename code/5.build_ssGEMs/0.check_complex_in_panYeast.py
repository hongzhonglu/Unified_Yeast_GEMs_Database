# -*- coding: utf-8 -*-
# date : 2023/6/5 
# author : wangh
# file : 0.check_complex_in_panYeast.py
# project : Unified_Yeast_GEMs_Database
'''check the reactions involved complex enzyme which means need more than 1 genes simultaneously to activate the reaction
'''
from cobra.io import read_sbml_model

panYeast=read_sbml_model("model/panYeast_v4_5.xml")

complex_rxnList=[]
for rxn in panYeast.reactions:
    if 'and' in rxn.gene_reaction_rule:
        complex_rxnList.append(rxn)






