# -*- coding: utf-8 -*-
# date : 2023/7/12 
# author : wangh
# file : 2.different_flux_analysis.py
# project : Unified_Yeast_GEMs_Database
'''Different flux analysis of the model'''
import omicverse as ov
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import json



# load data
df_wt_flux=pd.read_csv(r'code/7.anaerobic_growth_analysis/output/wt_pfba_fluxes.csv',index_col=0)
df_bioethanol_flux=pd.read_csv(r'code/7.anaerobic_growth_analysis/output/bioethanol_pfba_fluxes.csv',index_col=0)

# if the column has 0 value in growth index, then delete it
df_wt_flux=df_wt_flux.loc[:,df_wt_flux.loc['growth',:]>0]
df_bioethanol_flux=df_bioethanol_flux.loc[:,df_bioethanol_flux.loc['growth',:]>0]

# remove growth row
df_wt_flux=df_wt_flux.drop('growth',axis=0)
df_bioethanol_flux=df_bioethanol_flux.drop('growth',axis=0)

# extract the strainList of different type
bioethanol_strainlist=df_bioethanol_flux.columns.tolist()
wt_strainlist=df_wt_flux.columns.tolist()

# combine the fluxes of different type according to index
df_combine_flux=pd.concat([df_wt_flux,df_bioethanol_flux],axis=1)
# fill nan with 0
df_combine_flux=df_combine_flux.fillna(0)

# calculate the difference of fluxes
dds_flux=dds=ov.bulk.pyDEG(df_combine_flux)

# normalize the fluxes
dds_flux.normalize()

# calculate the different flux reaction
result=dds_flux.deg_analysis(bioethanol_strainlist,wt_strainlist,method='ttest')

# This function automatically calculates the appropriate threshold based on the log2FC distribution
# -1 means automatically calculates
dds_flux.foldchange_set(fc_threshold=1.2,
                   pval_threshold=0.001,
                   logp_max=50)



# filter according to baseMean
filter_result=result.loc[result['BaseMean']>0.005,:]

# plot the volcano plot
dds_flux.result=filter_result
dds_flux.plot_volcano(title='Differential Flux Analysis',
                      titlefont={'weight':'bold','size':15},
                      figsize=(5,4),
                      plot_genes=['r_0893','r_0486','r_0501','r_0569','r_0438'],
                      plot_genes_fontsize=12,
                      legend_bbox=(1.05,1),)
plt.savefig('figures/output/additional5_diff_flux_analysis.svg',bbox_inches='tight',dpi=400,transparent=True)
plt.show()

# load YeastGEM reaction information
df_reaction_info=pd.read_excel(r'model/yeast-GEM8.7.xlsx',index_col=1)
filter_result=filter_result.loc[filter_result.index.isin(df_reaction_info.index),:]

# save result
filter_result.to_csv(r'code/6.anaerobic_growth_analysis/output/bioethanol_vs_wt_different_flux_analysis.csv')


# pathway enrichment analysis
# group df_reaction_info index by subsystem
df_reaction_info_group=df_reaction_info.groupby('SUBSYSTEM').groups
# convert to dict, key is subsystem, value is reaction list
pathway_dict={k:list(v) for k,v in df_reaction_info_group.items()}

# add subpathway information
with open(r'model/model_pathway_rxndict.json','r') as f:
    subpathway_dict=json.load(f)
# add Glycolysis_upstream,Glycolysis_downstream,'Anaerobic_fermentation' to pathway_dict
pathway_dict['Glycolysis_upstream']=subpathway_dict['Glycolysis_upstream']
pathway_dict['Glycolysis_downstream']=subpathway_dict['Glycolysis_downstream']
pathway_dict['Anaerobic_fermentation']=subpathway_dict['Anaerobic_fermentation']


# pathway enrichment analysis
sig_rxns=filter_result[filter_result['sig']!='normal'].index.tolist()
enr=ov.bulk.geneset_enrichment(gene_list=sig_rxns,
                                pathways_dict=pathway_dict,
                                pvalue_type='auto',
                                background=df_combine_flux.index.tolist(),
                               description='None')

# save result
enr.to_csv(r'code/6.anaerobic_growth_analysis/output/bioethanol_vs_wt_pathway_enrichment_analysis.csv')
