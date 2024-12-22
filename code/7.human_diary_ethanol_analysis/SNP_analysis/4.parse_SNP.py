# -*- coding: utf-8 -*-
# date : 2024/11/18 
# author : wangh
# file : 4.parse_SNP.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd

df=pd.read_csv('code/SNP_analysis/output/human_bioethanol_dairy_wt_genotype.tsv',sep='\t')

# count nan for each rows
count_null=df.isnull().sum(axis=1)
len(count_null[count_null<1])