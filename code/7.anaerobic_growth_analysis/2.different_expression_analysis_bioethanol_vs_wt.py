# -*- coding: utf-8 -*-
# date : 2023/7/13 
# author : wangh
# file : 2.different_expression_analysis_bioethanol_vs_wt.py
# project : Unified_Yeast_GEMs_Database
import omicverse as ov
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import json


# load strain information data
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)

wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()


# load rxn_tpmMatrix
rxn_tpmMatrix=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_rxn_tpmMatrix.csv',index_col=0)

wt_strainList=[strain for strain in wt_strainList if strain in rxn_tpmMatrix.columns]
bioethanol_strainList=[strain for strain in bioethanol_strainList if strain in rxn_tpmMatrix.columns]


# do differential expression analysis
# fill na with 0
rxn_tpmMatrix=rxn_tpmMatrix.fillna(0)

# calculate the difference of fluxes
dds_exp=dds=ov.bulk.pyDEG(rxn_tpmMatrix)

# normalize the fluxes
dds_exp.normalize()

# calculate the different flux reaction
result=dds_exp.deg_analysis(bioethanol_strainList,wt_strainList,method='ttest')

dds_exp.foldchange_set(fc_threshold=0.6,
                   pval_threshold=0.05,
                   logp_max=50)

len(result[result['sig']!='normal'])
# plot the volcano plot
dds_exp.plot_volcano(title='Different Expression Analysis',
                      titlefont={'weight':'bold','size':15},
                      figsize=(5,4),
                      plot_genes_fontsize=12,
                    legend_bbox=(1.05,1),)
# save the figure
plt.savefig('figures/output/additional5_diff_expression_analysis.svg',bbox_inches='tight',dpi=400,transparent=True)
plt.show()

# save result
result.to_csv('code/7.anaerobic_growth_analysis/output/bioethanol_vs_wt_diff_exp.csv')

# pathway enrichment analysis
# group df_reaction_info index by subsystem
# load YeastGEM reaction information
df_reaction_info=pd.read_excel(r'model/yeast-GEM8.7.xlsx',index_col=1)
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
sig_rxns=result[result['sig']!='normal'].index.tolist()
enr=ov.bulk.geneset_enrichment(gene_list=sig_rxns,
                                pathways_dict=pathway_dict,
                                pvalue_type='auto',
                                background=rxn_tpmMatrix.index.tolist(),
                               description='None')

ov.bulk.geneset_plot(enr,figsize=(8,3),fig_title='Wiki Pathway enrichment',
                        cmap='Reds')
plt.show()

# save result
enr.to_csv(r'code/7.anaerobic_growth_analysis/output/bioethanol_vs_wt_diff_exp_pathway_enrichment_analysis.csv')
