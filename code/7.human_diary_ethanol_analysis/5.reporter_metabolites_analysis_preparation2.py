import pandas as pd

# organize the differential expresion analysis result by DEseq2

df_bioethanol=pd.read_csv('code/human_diary_analysis/output/diff_expression_analysis_bioethanol.csv',index_col=0)
df_human=pd.read_csv('code/human_diary_analysis/output/diff_expression_analysis_human.csv',index_col=0)
df_diary=pd.read_csv('code/human_diary_analysis/output/diff_expression_analysis_diary.csv',index_col=0)


df_foldchange=pd.DataFrame(index=df_bioethanol.index)
df_foldchange['bioethanol']=df_bioethanol['log2FoldChange']
df_foldchange['human']=df_human['log2FoldChange']
df_foldchange['diary']=df_diary['log2FoldChange']

df_pvalue=pd.DataFrame(index=df_bioethanol.index)
df_pvalue['bioethanol']=df_bioethanol['padj']
df_pvalue['human']=df_human['padj']
df_pvalue['diary']=df_diary['padj']

# remove all nan rows in foldchange
df_foldchange=df_foldchange.dropna(how='all')
df_pvalue=df_pvalue.loc[df_foldchange.index]

# save result
df_foldchange.to_csv('code/human_diary_analysis/output/log2foldchange_DEseqs2.csv')
df_pvalue.to_csv('code/human_diary_analysis/output/pvalue_DEseqs2.csv')