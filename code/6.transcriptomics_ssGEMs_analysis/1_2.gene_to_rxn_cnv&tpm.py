# -*- coding: utf-8 -*-
# date : 2023/7/6 
# author : wangh
# file : gene_to_rxn_cnv&tpm.py
# project : Unified_Yeast_GEMs_Database
'''Get the CNV and TPM information for each reaction according to the GPR rules.max for OR, min for AND'''

import pandas as pd
import os
from cobra.io import read_sbml_model

def gene_to_rxn_calculator(gene_dataMatrix,rxn):
    '''calculate the related copy number or related transcriptional level for each reaction according to the GPR rules and gene cnvMatrix data.
    calculating
    rules: max for OR, min for AND
    '''
    # get related gene cnvMatrix
    geneIDlist=[gene.id for gene in rxn.genes]
    rxn_gene_dataMatrix=gene_dataMatrix.loc[gene_dataMatrix.index.isin(geneIDlist),:]

    # calculate the related copy number for each reaction according to the GPR rules.
    gpr=rxn.gene_reaction_rule

    # if there is no and in gpr, then the copy number is the max of related genes
    if 'and' not in gpr:
        # rxn_data=rxn_gene_dataMatrix.max(axis=0) # max for OR
        rxn_data=rxn_gene_dataMatrix.sum(axis=0)   # sum for OR
    # if there is no or in gpr, then the copy number is the min of related genes
    elif 'or' not in gpr:
        rxn_data=rxn_gene_dataMatrix.min(axis=0)
    # if there is both and and or in gpr, then firstly, divede the gpr into several parts according to or,
    # then calculate the copy number for each part by min, finally, the copy number is the max of these parts.
    else:
        df_complexes_dataMatrix=pd.DataFrame()
        complexlist=gpr.split('or')
        for complex in complexlist:
            # remove the space in both sides
            complex=complex.strip()
            # remove the bracket in both sides
            complex=complex.strip('(')
            complex=complex.strip(')')
            geneIDlist=complex.split('and')
            # remove the space in both sides
            geneIDlist=[geneID.strip() for geneID in geneIDlist]
            complex_data=rxn_gene_dataMatrix.loc[rxn_gene_dataMatrix.index.isin(geneIDlist),:].min(axis=0)
            df_complexes_dataMatrix=pd.concat([df_complexes_dataMatrix,complex_data],axis=1)
        # rxn_data=df_complexes_dataMatrix.max(axis=1)
        rxn_data=df_complexes_dataMatrix.sum(axis=1)
    return rxn_data


# load panModel
panYeast=read_sbml_model('model/panYeast.xml')


# get rxn_cnvMatrix for each strain according to geneMatrix and GPR rules
# load copy number data
cnvMatrix=pd.read_csv('data/geneMatrix/pan1800_v2_blastp_50_70_cnvMatrix.csv',index_col=0)
# calculate rxn_cnvMatrix
df_rxn_cnvMatrix=pd.DataFrame(index=cnvMatrix.columns)
for rxn in panYeast.reactions:
    id=rxn.id
    rxn_data=gene_to_rxn_calculator(cnvMatrix,rxn)
    df_rxn_cnvMatrix[id]=rxn_data

df_rxn_cnvMatrix=df_rxn_cnvMatrix.T
# fillna with 0
df_rxn_cnvMatrix.fillna(0,inplace=True)

# save rxn_cnvMatrix
df_rxn_cnvMatrix.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce1800_rxn_cnvMatrix.csv')

# get rxn_tpmMatrix for each strain according to geneMatrix and GPR rules
# load tpm data
tpmMatrix=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.csv',index_col=0)
# calculate rxn_tpmMatrix
df_rxn_tpmMatrix=pd.DataFrame(index=tpmMatrix.columns)
for rxn in panYeast.reactions:
    id=rxn.id
    rxn_data=gene_to_rxn_calculator(tpmMatrix,rxn)
    df_rxn_tpmMatrix[id]=rxn_data

df_rxn_tpmMatrix=df_rxn_tpmMatrix.T
# fillna with 0
df_rxn_tpmMatrix.fillna(0,inplace=True)

# save rxn_tpmMatrix
df_rxn_tpmMatrix.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_rxn_tpmMatrix.csv')



