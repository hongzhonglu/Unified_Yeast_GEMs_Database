# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：predicted_allcds_eval.py
# 2022/5/27

'''collect information including genome size、assemble level、the number of genes，ssGEM size：the number of
rxns & genes'''

import pandas as pd
import numpy as np

# genome size statistic
new_strain_info=pd.DataFrame()
maker_genome_info=pd.read_csv('result/1900_allcds_maker_stats.txt',sep='\s+')
flo_genome_info=pd.read_csv("result/1900_allcds_flo_stats.txt",sep="\s+")
# na1011_strain_genome_size=pd.read_csv('data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
# strainlist=lg_strain_genome_size["file"].values[:10].tolist()
# for strain in strainlist:
#     strain=strain.split("_")[0]
#     gene_numb=len(na1011_strain_genome_size[na1011_strain_genome_size[strain]==1])
#     print(strain,":",gene_numb)
# genes_number
maker_genome_size=maker_genome_info['num_seqs']
maker_genome_size=maker_genome_size.str.replace(',','')
maker_genome_size=maker_genome_size.astype(int)

flo_genome_size=flo_genome_info['num_seqs']
flo_genome_size=flo_genome_size.str.replace(',','')
flo_genome_size=flo_genome_size.astype(int)
# flo_genome_size=list(flo_genome_size)

# genome_length
maker_genome_len=maker_genome_info['sum_len']
maker_genome_len=maker_genome_len.str.replace(',','')
maker_genome_len=maker_genome_len.astype(int)

flo_genome_len=flo_genome_info['sum_len']
flo_genome_len=flo_genome_len.str.replace(',','')
flo_genome_len=flo_genome_len.astype(int)


# figure
# libraries & dataset
import seaborn as sns
import matplotlib.pyplot as plt


# all_strain_genome_info=pd.read_csv('S.cerevisiae_ssGEMs_auto-construction/data/all_strain_genome_info.csv')
# new_strain_info=all_strain_genome_info[all_strain_genome_info['source']!='nature_1011']
# nature_1011_info=all_strain_genome_info[all_strain_genome_info['source']=='nature_1011']

# set a grey background (use sns.set_theme() if seaborn version 0.11.0 or above)
sns.set(style="darkgrid")
# plotting both distibutions on the same figure
# predicted_genes_number
fig_maker_gene_numb=sns.kdeplot(maker_genome_size,shade=True)
fig_flo_gene_numb=sns.kdeplot(flo_genome_size,shade=True)
plt.legend(loc='upper left',labels=['maker_annotation','liftOver_annotation'],title='result from different annotation processes')
plt.title('Numbers of genes in different strains')
plt.show()
plt.savefig('S.cerevisiae_ssGEMs_auto-construction/result/numb_of_genes_evaluation.png')
