# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：ssGEM_reconstruction.py
# 2022/4/25

'''construction strain-specified GEM for S.cerevisiae strains'''

import os
import sys
from cobra.io import read_sbml_model, write_sbml_model


sys.path.append("code")
from mainFunction import *
os.chdir('code')
# template model
pan_model=read_sbml_model('../data/panYeast_v3.xml')

#produce the dataframe for the metabolites and the rxn based on panYeast
pan_met_nov = produceMetaboliteList(pan_model)
pan_rxn_nov = produceRxnList(pan_model)
pan_gene = produceGeneList(pan_model)

# read geneMatrix
geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
strainList = list(geneMatrix.columns)
strainList0 = strainList[2:]


# build strain-specific GEMs
for i in range(len(strainList0)):
    strain_name=strainList0[i]+'.xml'
    strain_list=os.listdir('../result/1011_ssGEMs')
    if strain_name in strain_list:
        print(strain_name+' already exist')
    else:
        print(i)
        NEW = getStrainGEM(s0=strainList0[i], geneMatrix0=geneMatrix, templateGEM=pan_model, templateGene= pan_gene)
        write_sbml_model(NEW, "../result/1011_ssGEMs/" + strainList0[i]+ ".xml")