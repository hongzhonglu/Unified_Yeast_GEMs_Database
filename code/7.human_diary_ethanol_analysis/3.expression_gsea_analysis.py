# -*- coding: utf-8 -*-
# date : 2024/3/5 
# author : wangh
import gseapy as gp
import json
import pandas as pd

# get database
kegg_gene_dict=json.load(open('data/sce_kegg_geneset.json'))

# load expression data
df_expression=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix_normalized.csv',index_col=0)

#load strainlist
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()
human_strainList=df_strain_info[(df_strain_info['nature_clade']=='10. French Guiana human ')&(df_strain_info['type']=='Human')].index.tolist()
diary_strainList=df_strain_info[(df_strain_info['nature_clade']=='5. French dairy ')&(df_strain_info['type']=='Fermentation')].index.tolist()

wt_expression=df_expression[df_expression.columns.intersection(wt_strainList)]
bioethanol_expression=df_expression[df_expression.columns.intersection(bioethanol_strainList)]
human_expression=df_expression[df_expression.columns.intersection(human_strainList)]
diary_expression=df_expression[df_expression.columns.intersection(diary_strainList)]

# run gsea for human, diary, bioethanol to wt

# human vs wt
from gseapy import GSEA

def run_gsea(df_expression_pos, df_expression_neg, gene_sets,pos_name,neg_name='WT'):

    class_vector = [pos_name] * df_expression_pos.shape[1] + [neg_name] * df_expression_neg.shape[1]
    df_expression=pd.concat([df_expression_pos,df_expression_neg],axis=1)

    gs = GSEA(data=df_expression,
              gene_sets=gene_sets,
              classes=class_vector,  # cls=class_vector
              # set permutation_type to phenotype if samples >=15
              permutation_type='phenotype',
              permutation_num=1000,  # reduce number to speed up test
              outdir=None,
              method='signal_to_noise',
              threads=4, seed=8)
    gs.pheno_pos = pos_name
    gs.pheno_neg = neg_name
    gs.run()
    result=gs.res2d
    return gs,result

human_gs,human_result=run_gsea(df_expression_pos=human_expression,df_expression_neg=wt_expression,gene_sets=kegg_gene_dict,
                     pos_name='human',neg_name='wt')

bioethanol_gs,bioethanol_result=run_gsea(df_expression_pos=bioethanol_expression,df_expression_neg=wt_expression,gene_sets=kegg_gene_dict,
                                         pos_name='bioethanol',neg_name='wt')

diary_gs,diary_result=run_gsea(df_expression_pos=diary_expression,df_expression_neg=wt_expression,gene_sets=kegg_gene_dict,
                                 pos_name='diary',neg_name='wt')


