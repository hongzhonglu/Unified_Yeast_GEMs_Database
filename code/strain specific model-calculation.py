# -*- coding: utf-8 -*-
'''this code is to carry out the strain specific model reconstruction
5th, Nov, 2018'''



# Import packages
import pandas as pd
import os
import sys
import pprint
from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import delete_model_genes, find_gene_knockout_reactions
from __future__ import absolute_import
import numpy as np
from functools import reduce


from cobra.core import Metabolite, Model, Reaction
from cobra.manipulation import *
from cobra.flux_analysis.variability import (
    find_blocked_reactions, find_essential_genes, find_essential_reactions,
    flux_variability_analysis)

os.chdir('/Users/luho/PycharmProjects/model/cobrapy/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")
pprint.pprint(sys.path)
# import self function
from mainFunction import *


# input the panYeast
#pan_nov = read_sbml_model('../data/panYeast_nov.xml')
pan_nov = read_sbml_model('../data/panYeast_v3.xml')
#pan_nov= correctSomeWrongFormat(pan_nov)

#produce the dataframe for the metabolites and the rxn based on panYeast
pan_met_nov = produceMetaboliteList(pan_nov)
pan_rxn_nov = produceRxnList(pan_nov)
pan_gene = produceGeneList(pan_nov)


#produce the strains specific model
#input the gene matrix
geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
strainList = list(geneMatrix.columns)
strainList0 = strainList[2:]
geneMatrix.loc[:,'geneID'] = geneMatrix.loc[:,'geneID'].str.replace('-', '__45__')
geneMatrix.loc[:,'geneID'] = geneMatrix.loc[:,'geneID'].str.replace('.', '__46__')


### output the result
# test print the rxn reduced for each strain specific models
for i in range(5):
    print(i)
    NEW = getStrainGEM(s0=strainList0[i], geneMatrix0=geneMatrix, templateGEM=pan_nov, templateGene= pan_gene)
    write_sbml_model(NEW, "result/" + strainList0[i]+ ".xml")




