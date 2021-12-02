# -*- coding: utf-8 -*-
'''this code is to carry out the strain specific model reconstruction
5th, Nov, 2018'''



# Import packages
import pandas as pd
import os    ##for directory
import sys
import pprint
from cobra.io import read_sbml_model
from cobra.manipulation.delete import delete_model_genes, find_gene_knockout_reactions
#from __future__ import absolute_import
import numpy as np
from functools import reduce

from cobra.core import Metabolite, Model, Reaction
#from cobra.manipulation import *


os.chdir('/Users/xluhon/Documents/GitHub/Unified_Yeast_GEMs_Database/code')
sys.path.append(r"/Users/xluhon/Documents/GitHub/Unified_Yeast_GEMs_Database/code")
pprint.pprint(sys.path)
# import self function
from mainFunction import *

# input the subsystem information
gem_dataframe = pd.read_excel('../data/yeastGEM_with subsystem.xlsx')

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
out9 = open('../result/strain specific model-reduced rxn number', 'w')
out10 = open('../result/strain specific model-rxn number', 'w')
out11 = open('../result/strain specific model-reduced gene list from strain specific model', 'w')
out12 = open('../result/common rxn from all strain specific model', 'w')
out13 = open('../result/strain specific model-gene number', 'w')
out14 = open('../result/common gene from all strain specific model', 'w')

rxn_list0 = []
gene_list0 = []
for i in range(len(strainList0)):
    print(i)
    rxn_list = getStrainGEMrxn(s0=strainList0[i], geneMatrix0=geneMatrix, templateGEM=pan_nov, templateGene= pan_gene)
    print(strainList0[i], list(set(pan_rxn_nov['rxnID'])-set(rxn_list)))
    rxn_list0.append(rxn_list)
    #save the reduced rxn
    rxn_info = list(set(pan_rxn_nov['rxnID'])-set(rxn_list))
    rxn_info.insert(0, strainList0[i])
    print(rxn_info)
    out9.write(';'.join(rxn_info) + '\n')
    #save the rxn number analysis
    rxn_num = [strainList0[i], len(rxn_list), 4040-len(rxn_list)]
    rxn_num1 = [str(x) for x in rxn_num]
    out10.write(';'.join(rxn_num1) + '\n')
    #save the removed gene list for each strain specific model
    reduced_gene_list =getRemoveGeneList(s0=strainList0[i], geneMatrix0=geneMatrix, templateGEM=pan_nov, templateGene= pan_gene)
    gene_list = list(set(pan_gene)-set(reduced_gene_list))
    gene_list0.append(gene_list)
    gene_num = [strainList0[i], len(gene_list), len(reduced_gene_list)]
    gene_num1 = [str(x) for x in gene_num]
    out13.write(';'.join(gene_num1) + '\n')
    reduced_gene_list.insert(0, strainList0[i])
    out11.write(';'.join(reduced_gene_list) + '\n')

# calculate the common rxn
common_rxn = list(reduce(np.intersect1d, rxn_list0))
common_gene = list(reduce(np.intersect1d, gene_list0))
out12.write(';'.join(common_rxn) + '\n')
out14.write(';'.join(common_gene) + '\n')

out9.close()
out10.close()
out11.close()
out12.close()
out13.close()
out14.close()
#example
# reduce(np.intersect1d, [[1, 3, 4, 3], [3, 1, 2, 1], [6, 3, 4, 2]])

#gene classification
gene_analysis = pd.DataFrame({'gene':pan_gene})
gene_analysis['common_gene'] = gene_analysis['gene'].isin(common_gene)
#standard the gene name
gene_analysis.loc[:,'gene'] = gene_analysis.loc[:,'gene'].str.replace('__45__','-')
gene_analysis.loc[:,'gene'] = gene_analysis.loc[:,'gene'].str.replace('__46__','.')

saveExcel(gene_analysis,'../result/gene_analysis_panYeast.xlsx')






#other code used only for panYeast of early version
#update the gene relation due to some wrong format
#this error is corrected in the version2
special_gene = pan_nov.genes.get_by_id("49__45__MEL1__45__R__46__14846__45____59__1083__45__augustus_masked__45__BPL_2__45__7584")
special_gene.reactions
rxn_need_check = []
for i in special_gene.reactions:
    print(i.id)
    rxn_need_check.append(i.id)
#change the gpr for the 11 reactions:
for i in rxn_need_check:
    print(i)
    pan_nov.reactions.get_by_id(i).gene_reaction_rule = "49__45__MEL1__45__R__46__14846__45__ or 1083__45__augustus_masked__45__BPL_2__45__7584"


####one example compared with COBRA3
s1 = ['geneID', 'BFC']
geneList = geneMatrix.loc[:, s1]
gene_exist = singleMapping(geneList.loc[:, 'BFC'].tolist(), geneList.loc[:, 'geneID'].tolist(), pan_gene,
                           dataframe=False)
gene_exist = [0 if v is None else v for v in gene_exist]
gene_remove = [x for x, y in zip(pan_gene, gene_exist) if y < 1]

newModel = pan_nov.copy()
remove_genes(newModel, gene_remove,  remove_reactions=True)
rxn = []
for x in newModel.reactions:
    rxn.append(x.id)

reducedRxnHong = list(set(pan_rxn_nov['rxnID'])-set(rxn))
reducedRxnFei = pd.read_excel('../data/removed_rxn_BFC.xlsx')
reducedRxnFei = reducedRxnFei['removed_rxn_BFC'].tolist()

# it can be found that 18 new reactions were removed
rxn_hong = list(set(reducedRxnHong) - set(reducedRxnFei))
rxn_fei = set(reducedRxnFei) -set(reducedRxnHong)
pan_rxn_nov['check_sign'] = pan_rxn_nov['rxnID'].isin(rxn_hong)

#obtain the removed gene for each reaction
#fistly establish the relation between gene and reaction
rxn_gene = getRXNgeneMapping(rxn0=pan_rxn_nov['rxnID'], gpr0=pan_rxn_nov['GPR'])
rxn_gene.loc[:,'gene'] = rxn_gene.loc[:,'gene'].str.replace('-', '__45__')
rxn_gene.loc[:,'gene'] = rxn_gene.loc[:,'gene'].str.replace('.', '__46__')
rxn_gene['removed_sign'] = rxn_gene['gene'].isin(gene_remove)
rxn_gene0 = rxn_gene[rxn_gene['removed_sign']==True]
pan_rxn_nov['remove_gene'] = multiMapping(rxn_gene0['gene'],rxn_gene0['rxnID'],pan_rxn_nov['rxnID'])
pan_rxn_nov_check = pan_rxn_nov[pan_rxn_nov['check_sign']==True]





#part3 add the panID for the collapsed gene in panYeast, which don't have the panID from pangenome
#it can be found that some genes from SGD were removed from pan. It is essential for us to add them into the panYeast
#firstly compare the yeastGEM and panYeast in gene
# input the coreYeast, panYeast and yeastGEM
geneAll = geneMatrix.loc[:,'geneID']
gene_collapse = set(pan_gene) - set(geneAll)
len(gene_collapse)
#save the results
gene_collapse0 = pd.DataFrame({'geneID': list(gene_collapse)})

#input the mapping analysis results: sgd-panID
sgd_panid = pd.read_excel('../data/sgd_panid.xlsx')
gene_collapse0['panID'] = singleMapping(sgd_panid['panID'], sgd_panid['Systematic_name'], gene_collapse0['geneID'])
gene_collapse0['Remark'] = singleMapping(sgd_panid['Remark'], sgd_panid['Systematic_name'], gene_collapse0['geneID'])

#save the result for Feiran
saveExcel(gene_collapse0, "../result/gene_collapse0_panID.xlsx")



#GR191W should be r_1201 YGR191W?????
special_gene1 = pan_nov.genes.get_by_id('GR191W')
pan_nov.reactions.get_by_id('r_4592').gene_reaction_rule = "YGR191W"

#so how to relace the genes?????
gene_collapse0['panID'] = gene_collapse0['panID'].str.replace('-', '__45__')
gene_collapse0['panID'] = gene_collapse0['panID'].str.replace('.', '__46__')
gene_collapse1 = gene_collapse0.dropna()

#obtain the rxn list contains the above genes
rxn_gene['update_sign'] = rxn_gene['gene'].isin(gene_collapse1['geneID'])
rxn_gene1=rxn_gene[rxn_gene['update_sign']==True]
rxn_update = rxn_gene1['rxnID'].tolist()
for i in rxn_update:
    print(i)
    old_gpr=pan_nov.reactions.get_by_id(i).gene_reaction_rule
    print(old_gpr)
    new_gpr=updateGPR(old_gpr, nameMapping=gene_collapse1)
    print(new_gpr)
    pan_nov.reactions.get_by_id(i).gene_reaction_rule = new_gpr




####one example compared with COBRA3 after we replace the new gene with
s1 = ['geneID', 'BFC']
geneList = geneMatrix.loc[:, s1]
gene_exist = singleMapping(geneList.loc[:, 'BFC'].tolist(), geneList.loc[:, 'geneID'].tolist(), pan_gene,dataframe=False)
gene_exist = [0 if v is None else v for v in gene_exist]
gene_remove = [x for x, y in zip(pan_gene, gene_exist) if y < 1]

#produce the rxn list for the new model
newModel = pan_nov.copy()
remove_genes(newModel, gene_remove,  remove_reactions=True)
rxn = []
for x in newModel.reactions:
    rxn.append(x.id)

#compare the result between cobrapy and cobrav3
reducedRxnHong = list(set(pan_rxn_nov['rxnID'])-set(rxn))
reducedRxnFei = pd.read_excel('../data/removed_rxn_BFC.xlsx')
reducedRxnFei = reducedRxnFei['removed_rxn_BFC'].tolist()
# it can be found that 18 new reactions were removed
rxn_hong = list(set(reducedRxnHong) - set(reducedRxnFei))
rxn_fei = set(reducedRxnFei) -set(reducedRxnHong)
pan_rxn_nov['check_sign'] = pan_rxn_nov['rxnID'].isin(rxn_hong)

#obtain the removed gene for each reaction
#fistly establish the relation between gene and reaction
rxn_gene = getRXNgeneMapping(rxn0=pan_rxn_nov['rxnID'], gpr0=pan_rxn_nov['GPR'])
rxn_gene.loc[:,'gene'] = rxn_gene.loc[:,'gene'].str.replace('-', '__45__')
rxn_gene.loc[:,'gene'] = rxn_gene.loc[:,'gene'].str.replace('.', '__46__')
rxn_gene['removed_sign'] = rxn_gene['gene'].isin(gene_remove)
rxn_gene0 = rxn_gene[rxn_gene['removed_sign']==True]
pan_rxn_nov['remove_gene'] = multiMapping(rxn_gene0['gene'],rxn_gene0['rxnID'],pan_rxn_nov['rxnID'])
pan_rxn_nov_check = pan_rxn_nov[pan_rxn_nov['check_sign']==True]
