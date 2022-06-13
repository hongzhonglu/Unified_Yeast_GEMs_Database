# -*- coding: utf-8 -*-
'''this code is to carry out the model comparative analysis among yeastGEM, panYeast and coreYeast
5th, Nov, 2018'''


# Import packages
import pandas as pd
import os    ##for directory
import sys
import pprint
from cobra.io import read_sbml_model


os.chdir('/Users/xluhon/Documents/GitHub/Unified_Yeast_GEMs_Database/code')
sys.path.append(r"/Users/xluhon/Documents/GitHub/Unified_Yeast_GEMs_Database/code")
pprint.pprint(sys.path)

# import self function
from mainFunction import *


# input the subsystem information
gem_dataframe = pd.read_excel('../data/yeastGEM_with subsystem.xlsx')


# input the coreYeast, panYeast and yeastGEM
pan_nov = read_sbml_model('../data/panYeast_v2.xml')

pan_nov = correctSomeWrongFormat(pan_nov)
pan_gene = produceGeneList(pan_nov)
pan_rxn = produceRxnList(pan_nov)

# input the core reaction and core genes from strains specific model reconstruction
with open('../result/common gene from all strain specific model', 'r') as common_gene0:
    for line in common_gene0:
        s0 = line.split(";")
        print(s0)
common_gene = [x.replace("\n","") for x in s0]
common_gene = [x.replace("__45__","-") for x in common_gene]

# obtain the accessory gene
accessory_gene = list(set(pan_gene)-set(common_gene))


#here we establish the relation between the gene and subsystem
#one gene could be located in different subsystems
rxn_gene = getRXNgeneMapping(gem_dataframe['rxnID'], gem_dataframe['GPR'])
rxn_gene['subsystem'] = singleMapping(gem_dataframe['subsytem_automatic'],gem_dataframe['rxnID'],rxn_gene.loc[:,'rxnID'])
#further establish the single mapping between gene and subsystem

rxn_gene['gene_id'] = [None]*len(rxn_gene['subsystem'])
for i in range(len(rxn_gene['gene_id'])):
    rxn_gene['gene_id'][i] = rxn_gene['gene'][i] + '@@' + str(i)

#split the subsystem
gene_subsystem = splitAndCombine(rxn_gene['subsystem'], rxn_gene['gene_id'], sep0= "//", moveDuplicate=True)
gene_subsystem['gene'] = [None]*len(gene_subsystem['V2'])
gene_subsystem['gene'] = gene_subsystem['V1'].str.split('@@', 1, expand=True)
gene_subsystem0 = gene_subsystem.loc[:,['gene','V2']]
gene_subsystem0['gene'] = gene_subsystem0['gene'].str.strip()
gene_subsystem0['V2'] = gene_subsystem0['V2'].str.strip()
gene_subsystem1 = gene_subsystem0.drop_duplicates()
gene_subsystem2 = gene_subsystem1[gene_subsystem1['V2'] != "None"]
gene_subsystem3 = gene_subsystem2[gene_subsystem2['gene'] != "NA"]
# standard the transport subsystem
subsystem_newRxn = gene_subsystem3['V2'].tolist()
subsystem_newRxn0 = []
for x in subsystem_newRxn:
    print(x)
    if x is None:
        subsystem_newRxn0.append(None)
    elif 'Transport' in x:
        subsystem_newRxn0.append('Transport')
    else:
        subsystem_newRxn0.append(x)
#calculate the frequency
all_gene_subsystem = calculateFrequency(list0=subsystem_newRxn0, item0="all_gene_number")


#for the accessory gene
gene_subsystem_accessory = gene_subsystem3[gene_subsystem3['gene'].isin(accessory_gene)]
subsystem_newRxn = gene_subsystem_accessory['V2'].tolist()
subsystem_newRxn0 = []
for x in subsystem_newRxn:
    print(x)
    if x is None:
        subsystem_newRxn0.append(None)
    elif 'Transport' in x:
        subsystem_newRxn0.append('Transport')
    else:
        subsystem_newRxn0.append(x)
#calculate the frequency
accessory_gene_subsystem = calculateFrequency(list0=subsystem_newRxn0, item0="all_gene_number")
all_gene_subsystem['accessory_number'] = singleMapping(accessory_gene_subsystem['number'],accessory_gene_subsystem['all_gene_number'],all_gene_subsystem['all_gene_number'])

#here we can mainly summarize the subsystem with more than 8 genes in total gene number
#save the results
saveExcel(all_gene_subsystem, "../result/subsystem summary for gene from panYeast.xlsx")


#analysis the rxn
#input the accessory reaction
#here for the reaction we don't have a very reliable classification, thus it may be not right based on such a analysis
with open('../result/common rxn from all strain specific model', 'r') as common_rxn0:
    for line in common_rxn0:
        s0 = line.split(";")
        print(s0)
common_rxn = [x.replace("\n","") for x in s0]

# obtain the subsystem
accessory_rxn = list(set(pan_rxn['rxnID'])-set(common_rxn))
gem_dataframe_accessory = gem_dataframe[gem_dataframe['rxnID'].isin(accessory_rxn)]














































