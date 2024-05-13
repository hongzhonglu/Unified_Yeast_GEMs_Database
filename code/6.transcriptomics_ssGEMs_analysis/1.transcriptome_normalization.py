'''Normalize the read count expression data by median of ratios method from DESeq2'''
from pydeseq2.preprocessing import deseq2_norm
import pandas as pd

# load expression data
df_exp=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix.csv',index_col=0)

# normalize the expression data
df_norm_exp,size_factors=deseq2_norm(df_exp.T)

df_norm_exp=df_norm_exp.T
size_factors=pd.Series(size_factors,index=df_norm_exp.columns)

# save the normalized expression data
df_norm_exp.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix_normalized.csv')
size_factors.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix_normalized_factor.csv')

