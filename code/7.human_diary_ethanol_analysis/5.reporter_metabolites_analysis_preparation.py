'''Prepare the data for reporter metabolites analysis of diary ,human, bioethanol and wild type clade strains.
1.Relative gene expression of diary ,human, bioethanol vs wild type clade strains.
2.p value data
'''
import pandas as pd
import os
import numpy as np

# load RNAseq data
tpm_data=pd.read_excel('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.xlsx',index_col=0)
# remove rows with all 0
tpm_data=tpm_data.loc[(tpm_data!=0).any(axis=1)]

#load strainlist
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()
human_strainList=df_strain_info[(df_strain_info['nature_clade']=='10. French Guiana human ')&(df_strain_info['type']=='Human')].index.tolist()
diary_strainList=df_strain_info[(df_strain_info['nature_clade']=='5. French dairy ')&(df_strain_info['type']=='Fermentation')].index.tolist()

bioethanol_tpmMatrix=tpm_data[tpm_data.columns.intersection(bioethanol_strainList)]
human_tpmMatrix=tpm_data[tpm_data.columns.intersection(human_strainList)]
diary_tpmMatrix=tpm_data[tpm_data.columns.intersection(diary_strainList)]
wt_tpmMatrix=tpm_data[tpm_data.columns.intersection(wt_strainList)]

# calculate p value of each gene of each clade vs each gene of wild type strains
from scipy.stats import ttest_ind

df_clade_pvalue=pd.DataFrame(index=tpm_data.index)
# for each gene
for gene in tpm_data.index:
    # for each clade
    for clade in ['bioethanol','human','diary']:
        # calculate p value
        pvalue=ttest_ind(wt_tpmMatrix.loc[gene],eval(clade+'_tpmMatrix').loc[gene])[1]
        df_clade_pvalue.loc[gene,clade]=pvalue


# calculate the foldchange of each mean value to the mean value of wild type strains
df_clade_foldchange=pd.DataFrame(index=tpm_data.index)
df_clade_foldchange['bioethanol']=bioethanol_tpmMatrix.mean(axis=1)/wt_tpmMatrix.mean(axis=1)
df_clade_foldchange['human']=human_tpmMatrix.mean(axis=1)/wt_tpmMatrix.mean(axis=1)
df_clade_foldchange['diary']=diary_tpmMatrix.mean(axis=1)/wt_tpmMatrix.mean(axis=1)
# remove rows which include inf values
df_clade_foldchange=df_clade_foldchange.replace([np.inf, -np.inf], np.nan).dropna(how='all')
# remove rows with all nan values
df_clade_foldchange=df_clade_foldchange.dropna(how='all')


df_clade_pvalue=df_clade_pvalue.loc[df_clade_foldchange.index]

# save result
df_clade_pvalue.to_csv('code/human_diary_analysis/output/human_diary_bioethanol_pvalue.csv')
df_clade_foldchange.to_csv('code/human_diary_analysis/output/human_diary_bioethanol_genefoldchange.csv')