# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：ssGEMs_analysis1.py
# 2022/5/31
'''analyse the growth simulation results of ssGEMs with strains classification data'''

import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

from mainFunction import *

# read strains classification info
allstrain_info=pd.read_excel('data/1897_strains_info.xlsx')

# read lg1392 ssGEMs growth simulation result
lg1392_ssGEM_result=pd.read_csv('result/lg1392_growth_simulation_result.csv')
lg1392_ssGEM_result.rename(columns={'Unnamed: 0':'ssGEMs'},inplace=True)

# read nature_1011 ssGEMs growth simulation result
na1011_ssGEM_result=pd.read_csv('result/na1011_growth_simulation_result.csv')
na1011_ssGEM_result.rename(columns={'Unnamed: 0':'ssGEMs'},inplace=True)


# new600_ssGEMs growth simulation result
new600_ssGEM_result=pd.read_csv("result/new600_growth_simulation_result.csv")
new600_ssGEM_result.rename(columns={'Unnamed: 0':'ssGEMs'},inplace=True)


# modify ssGEMs id
# nature_1011 strains: eg,AAA_6.re.xml -> AAA
ssGEM_idlist=lg1392_ssGEM_result['ssGEM_id'].tolist()
ssGEM_idlist_new_v1=[]
for ssGEM_id in ssGEM_idlist:
    ssGEM_id_new=re.sub("\.re\.xml","",ssGEM_id)
    ssGEM_idlist_new_v1.append(ssGEM_id_new)

# other strains: only use the genome assemble accession number as ssGEM id，
# eg. GCA_001738115.1_ASM173811v1_genomic.xml > GCA_001738115.1
ssGEM_idlist_new_v2=[]
for ssGEM_id in ssGEM_idlist_new_v1:
    ssGEM_id_new=re.match("GCA_\d{9}\.\d",ssGEM_id)
    if ssGEM_id_new:
        print(ssGEM_id_new.group())
        ssGEM_idlist_new_v2.append(ssGEM_id_new.group())
    else:
        ssGEM_idlist_new_v2.append(ssGEM_id)

lg1392_ssGEM_result['ssGEM_id']=ssGEM_idlist_new_v1

# update ssGEM result
# lg1392_ssGEM_result.to_csv('result/lg1392_growth_simulation_result.csv')
# new600_ssGEM_result.to_csv('result/new600_growth_simulation_result.csv')


# lg1392 ssGEM simulation result analysis
ungrowable_strains=lg1392_ssGEM_result[~(lg1392_ssGEM_result['aerobic growth']>0)]
anaerobic_strains=ungrowable_strains[ungrowable_strains['anaerobic_growth']>0]
anaerobic_strainslist=anaerobic_strains.ssGEM_id.tolist()
anaerobic_strains_info=pd.DataFrame()
for strain in anaerobic_strainslist:
    strain_info=allstrain_info[allstrain_info['ssGEM']==strain]
    anaerobic_strains_info = pd.concat([anaerobic_strains_info, strain_info], axis=0, join='outer', ignore_index=True)

lg1392_ssGEM_result['aerobic growth'].value_counts()

len(lg1392_ssGEM_result)

na1011_strainlist=allstrain_info[allstrain_info['source']=='1011_nature']['ssGEM'].tolist()
na1011_strainlist=[i.replace("SACE_","") for i in na1011_strainlist]
lg_1011_ssGEM_rxn_counts=[]
lg_1011_ssGEM_gene_count=[]
for strain in na1011_strainlist:
    ssGEM_gene_numb=lg1392_ssGEM_result[lg1392_ssGEM_result['ssGEM_id']==strain]['gene_numb']
    ssGEM_rxn_numb=lg1392_ssGEM_result[lg1392_ssGEM_result['ssGEM_id']==strain]['rxn_numb']
    if len(ssGEM_rxn_numb):
        lg_1011_ssGEM_rxn_counts.append(ssGEM_rxn_numb.values[0])
        lg_1011_ssGEM_gene_count.append(ssGEM_gene_numb.values[0])
    else:
        print("%s can't find"%strain)

# compare the ssGEM size among lg_originated ssGEMs and nature1011_originated ssGEMs and new600_ssGEMs
na1011_ssGEM_result['rxn_numb'].describe()
na1011_ssGEM_result['gene_numb'].describe()
lg_1011_ssGEM_rxn_counts=pd.DataFrame(lg_1011_ssGEM_rxn_counts)
lg_1011_ssGEM_gene_count=pd.DataFrame(lg_1011_ssGEM_gene_count)
lg_1011_ssGEM_rxn_counts.describe()
lg_1011_ssGEM_gene_count.describe()
lg1392_ssGEM_result['rxn_numb'].describe()
lg1392_ssGEM_result['gene_numb'].describe()
new600_ssGEM_result['rxn_numb'].describe()
new600_ssGEM_result['gene_numb'].describe()
# compare ssGEMs gene number
fig_new_gene_numb=sns.kdeplot(new600_ssGEM_result['gene_numb'],shade=True)
fig_lg_gene_numb=sns.kdeplot(lg1392_ssGEM_result['gene_numb'],shade=True)
fig_na1011_gene_numb=sns.kdeplot(na1011_ssGEM_result['gene_numb'],shade=True)
plt.legend(loc='upper left',labels=['new_strains','lg1392_strains','nature1011_strains'],title='strains from different annotation processes')
plt.title('Numbers of genes in ssGEMs from different genome annotation results')
plt.show()

# compare ssGEMs rxn number
fig_new_rxn_numb=sns.kdeplot(new600_ssGEM_result['rxn_numb'],shade=True)
fig_lg_rxn_numb=sns.kdeplot(lg1392_ssGEM_result['rxn_numb'],shade=True)
fig_na1011_rxn_numb=sns.kdeplot(na1011_ssGEM_result['rxn_numb'],shade=True)
plt.legend(loc='upper left',labels=['new_strains','lg1392_strains','nature1011_strains'],title='strains from different annotation processes')
plt.title('Numbers of rxn in ssGEMs from different genome annotation results ')
plt.show()

# evaluate difference in ssGEMs size between new600 and lg1392
rxn_numb_t_stat, rxn_numb_p_val = stats.ttest_ind(new600_ssGEM_result['rxn_numb'], lg1392_ssGEM_result['rxn_numb'], equal_var=False)
gene_numb_t_stat,gene_numb_p_val=stats.ttest_ind(new600_ssGEM_result['gene_numb'], lg1392_ssGEM_result['gene_numb'], equal_var=False)
rxn_numb_p_val
gene_numb_p_val

#
ungrowable_strains=lg1392_ssGEM_result[(~(lg1392_ssGEM_result['aerobic growth']>0))&(~(ungrowable_strains['anaerobic_growth']>0))]
ungrowable_strainslist=ungrowable_strains.ssGEM_id.tolist()
ungrowable_strains_info=pd.DataFrame()
for strain in ungrowable_strainslist:
    strain_info = allstrain_info[allstrain_info['ssGEM'] == strain]
    ungrowable_strains_info = pd.concat([ungrowable_strains_info, strain_info], axis=0, join='outer', ignore_index=True)


lg_allstrain_classify=allstrain_info[(allstrain_info['source']=='1011_nature')
                                     |(allstrain_info['source']=='lg_others')]['type'].value_counts()
lg_ungrowable_strain_classify=ungrowable_strains_info['type'].value_counts()

new600_strain_classify=allstrain_info[(allstrain_info['source']!='1011_nature')
                                      &(allstrain_info['source']!='lg_others')]['type'].value_counts()

lg_allstrain_classify.plot.pie(colors=['r','b','g','c','m'],)
lg_ungrowable_strain_classify.plot.pie()
new600_strain_classify.plot.pie(colors=['b','r','g','c','m'])
plt.show()


# new600_ssGEMs
new600_ssGEM_result['aerobic growth']=new600_ssGEM_result['aerobic growth'].round(6)
new600_ssGEM_result['anaerobic_growth']=new600_ssGEM_result['anaerobic_growth'].round(6)
new600_ssGEM_result['aerobic growth'].value_counts()
new600_ssGEM_result['anaerobic_growth'].value_counts()

new600_ungrowable_strainlist=list(new600_ssGEM_result[new600_ssGEM_result['anaerobic_growth']>0]['ssGEM_id'])
new600_ungrowable_strains_info=pd.DataFrame()
for strain in new600_ungrowable_strainlist:
    strain_info = allstrain_info[allstrain_info['ssGEM'] == strain]
    new600_ungrowable_strains_info = pd.concat([new600_ungrowable_strains_info, strain_info], axis=0, join='outer', ignore_index=True)
new600_ungrowable_strains_calassify=new600_ungrowable_strains_info['type'].value_counts()
new600_ungrowable_strains_calassify.plot.pie(colors=['b','r','g','c','m'])
plt.show()


# strain classification
