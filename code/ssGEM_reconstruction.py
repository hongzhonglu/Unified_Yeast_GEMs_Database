# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：ssGEM_reconstruction.py
# 2022/4/25

'''construction strain-specified GEM for S.cerevisiae strains'''

import os
import sys
from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd

sys.path.append("code")
from mainFunction import *
os.chdir('code')
# template model
pan_model=read_sbml_model('../model/panYeast_v3.xml')

#produce the dataframe for the metabolites and the rxn based on panYeast
pan_met_nov = produceMetaboliteList(pan_model)
pan_rxn_nov = produceRxnList(pan_model)
pan_gene = produceGeneList(pan_model)

# read geneMatrix
geneMatrix = pd.read_csv('../data/new600_geneMatrix.csv')
geneMatrix.rename(columns={"Unnamed: 0":"geneID"},inplace=True)
geneMatrix_1011 = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
strainList = list(geneMatrix.columns)
strainList0 = strainList[1:]


# build strain-specific GEMs
for strain in strainList0:
    strain_name=strain.replace(".fa.dmnd","")+'.xml'
    strain_list=os.listdir('../model/new600_ssGEMs')
    if strain_name in strain_list:
        print(strain_name+' already exists')
    else:
        print("%s ssGEM is building"%strain)
        NEW = getStrainGEM(s0=strain, geneMatrix0=geneMatrix, templateGEM=pan_model, templateGene= pan_gene)
        write_sbml_model(NEW, "../model/new600_ssGEMs/" + strain_name)