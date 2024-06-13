# -*- coding: utf-8 -*-
# date : 2024/2/16 
# author : wangh
import os
import tqdm
import pandas as pd
import cobra

# load panmodel
panmodel=cobra.io.read_sbml_model('model/panYeast.xml')

geneList=[gene.id for gene in panmodel.genes]
rxnList=[rxn.id for rxn in panmodel.reactions]

models_dir=r'code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_shrinked'
ssGEMsList=os.listdir(models_dir)

geneMatrix=pd.DataFrame(index=geneList,columns=ssGEMsList)
rxnMatrix=pd.DataFrame(index=rxnList,columns=ssGEMsList)

for ssGEM in tqdm.tqdm(ssGEMsList):
    model=cobra.io.read_sbml_model(os.path.join(models_dir,ssGEM))
    strain_geneList=[gene.id for gene in model.genes]
    strain_rxnList=[rxn.id for rxn in model.reactions]
    geneMatrix.loc[strain_geneList,ssGEM]=1
    rxnMatrix.loc[strain_rxnList,ssGEM]=1

geneMatrix.fillna(0,inplace=True)
rxnMatrix.fillna(0,inplace=True)

# remove rows with all 0
geneMatrix=geneMatrix.loc[geneMatrix.sum(axis=1)>0]
rxnMatrix=rxnMatrix.loc[rxnMatrix.sum(axis=1)>0]

# remove .xml in column names
geneMatrix.columns=[i.replace('.xml','') for i in geneMatrix.columns]
rxnMatrix.columns=[i.replace('.xml','') for i in rxnMatrix.columns]

# save result
geneMatrix.to_csv('code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_geneMatrix.csv')
rxnMatrix.to_csv('code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_rxnMatrix.csv')