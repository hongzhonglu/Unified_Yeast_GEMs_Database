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
new_strain_genome_size=pd.read_csv('result/new_allcds_info.txt',sep='\s+')
lg_strain_genome_size=pd.read_csv("result/lg1392_allcds_info.txt",sep="\s+")

# genes_number
new_genome_size=new_strain_genome_size['num_seqs']
new_genome_size=new_genome_size.str.replace(',','')
new_genome_size=new_genome_size.astype(int)

lg_genome_size=lg_strain_genome_size['num_seqs']
lg_genome_size=lg_genome_size.str.replace(',','')
lg_genome_size=lg_genome_size.astype(int)


# genome_length
new_genome_len=new_strain_genome_size['sum_len']
new_genome_len=new_genome_len.str.replace(',','')
new_genome_len=new_genome_len.astype(int)

lg_genome_len=lg_strain_genome_size['sum_len']
lg_genome_len=lg_genome_len.str.replace(',','')
lg_genome_len=lg_genome_len.astype(int)


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
# fig_all_genome_size =sns.kdeplot(all_strain_genome_info['genome_size'], shade=True, color="r")
# plt.title('gene_numb distribution of different strains')

# fig_all_gene_numb=sns.kdeplot(all_strain_genome_info['number_of_gene'],shade=True,color='g')

# predicted_genes_number
fig_new_gene_numb=sns.kdeplot(new_genome_size,shade=True)
fig_lg_gene_numb=sns.kdeplot(lg_genome_size,shade=True)
plt.legend(loc='upper left',labels=['new_strains','lg1392_strains'],title='strains from different annotation processes')
plt.title('Numbers of genes in different strains')
plt.show()
plt.savefig('S.cerevisiae_ssGEMs_auto-construction/result/numb_of_genes_evaluation.png')

t_stat, p_val = stats.ttest_ind(new_genome_size, lg_genome_size, equal_var=False)
p_val

# genome size
fig_new_genome_size =sns.kdeplot(new_genome_len, shade=True)
fig_lg_genome_size =sns.kdeplot(lg_genome_len, shade=True)
plt.legend(loc="upper left",labels=['new_strains_genome','lg1392_strains_genome'],title='strains from different sources')
plt.title('The genome size distribution of different strains')
# plt.savefig('S.cerevisiae_ssGEMs_auto-construction/result/genome_size_distribution.png')
# fig = sns.kdeplot(df['sepal_length'], shade=True, color="b")
plt.show()

import scipy.stats as stats
import matplotlib.pyplot as plt
t_stat, p_val = stats.ttest_ind(new_genome_len, lg_genome_len, equal_var=False)
p_val