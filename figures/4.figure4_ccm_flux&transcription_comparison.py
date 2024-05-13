# -*- coding: utf-8 -*-

import pandas as pd


downgly_rxnlist=['r_1054','r_0486','r_0892','r_0893','r_0366','r_0962','r_0961']
upgly_rxnlist=['r_0534','r_0467','r_0886','r_0450']
af_rxnlist=['r_0959','r_2115','r_0173']
tca_rxnlist=['r_0300','r_0280','r_0658','r_1022','r_1021','r_0451','r_0713']
op_rxnlist=['r_0773','r_0439','r_0438','r_0226']
df_flux_transcription_pathway=pd.DataFrame(index=[upgly_rxnlist+downgly_rxnlist+af_rxnlist+tca_rxnlist+op_rxnlist],columns=['wt_flux','bioethanol_flux','wt_transcription','bioethanol_transcription'])

# Get mean PTM transcription level for wt , bioethanol and beer strain.
df_bioethanol_vs_wt_tpmMatrix=pd.read_csv(r'code/6.anaerobic_growth_analysis/output/bioethanol_vs_wt_rxn_tpmMatrix.csv',index_col=0)
df_rxn_tpmMatrix=pd.read_csv(r'code/7.transcriptomics_ssGEMs_analysis/output/sce969_rxn_tpmMatrix.csv',index_col=0).T
df_rxn_tpmMatrix['type']=df_bioethanol_vs_wt_tpmMatrix['type']
df_rxn_tpmMatrix=df_rxn_tpmMatrix.loc[:,['type']+upgly_rxnlist+downgly_rxnlist+af_rxnlist+tca_rxnlist+op_rxnlist]

wt_transcription=df_rxn_tpmMatrix[df_rxn_tpmMatrix['type']=='wildtype'].mean(axis=0)
bioethanol_transcription=df_rxn_tpmMatrix[df_rxn_tpmMatrix['type']=='bioethanol'].mean(axis=0)

df_flux_transcription_pathway['wt_transcription']=wt_transcription.values
df_flux_transcription_pathway['bioethanol_transcription']=bioethanol_transcription.values


# get flux data
# load data
df_wt_flux=pd.read_csv(r'code/6.anaerobic_growth_analysis/output/wt_pfba_fluxes.csv',index_col=0)
df_bioethanol_flux=pd.read_csv(r'code/6.anaerobic_growth_analysis/output/bioethanol_pfba_fluxes_v2.csv',index_col=0)

# if the column has 0 value in growth index, then delete it
df_wt_flux=df_wt_flux.loc[:,df_wt_flux.loc['growth',:]>0]
df_bioethanol_flux=df_bioethanol_flux.loc[:,df_bioethanol_flux.loc['growth',:]>0]

# remove growth row
df_wt_flux=df_wt_flux.drop('growth',axis=0)
df_bioethanol_flux=df_bioethanol_flux.drop('growth',axis=0)

# only keep the flux of the rxn in the pathway
df_wt_flux=df_wt_flux.loc[upgly_rxnlist+downgly_rxnlist+af_rxnlist+tca_rxnlist+op_rxnlist,:]
df_bioethanol_flux=df_bioethanol_flux.loc[upgly_rxnlist+downgly_rxnlist+af_rxnlist+tca_rxnlist+op_rxnlist,:]

wt_fluxes=df_wt_flux.mean(axis=1)
bioethanol_fluxes=df_bioethanol_flux.mean(axis=1)

df_flux_transcription_pathway['wt_flux']=wt_fluxes.values
df_flux_transcription_pathway['bioethanol_flux']=bioethanol_fluxes.values

df_flux_transcription_pathway['flux_foldchange']=df_flux_transcription_pathway['bioethanol_flux']/df_flux_transcription_pathway['wt_flux']
df_flux_transcription_pathway['transcription_foldchange']=df_flux_transcription_pathway['bioethanol_transcription']/df_flux_transcription_pathway['wt_transcription']

# save result
df_flux_transcription_pathway.to_csv(r'figures/output/figure4_ccm_flux_transcription_comparison_v2.csv')