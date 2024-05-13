# -*- coding: utf-8 -*-
# date : 2023/5/23 
# author : wangh
# file : lgcontamination_check.py
# project : Unified_Yeast_GEMs_Database
'''parse the blastp result of pan1900_new vs protein database'''
import pandas as pd
import os
# set work dir
os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database\code\3.pan-genome_construction\contamination')

def parse_blastp_result(blastp_file,file_dir):
    #query,target,pident,bits,qcov,tcov,qlen,tlen,evalue,taxid,taxname,taxlineage
    col=['query','target','pident','bits','qcov','tcov','qlen','tlen','evalue','taxid','taxname','taxlineage']
    df=pd.read_csv(file_dir+blastp_file,sep='\t',names=col)
    return df

# parse all pan1900_new vs nr blastp result
blastp_result_dir='./output/'
blastp_file='lgremove900_vs_nr.txt'
df=parse_blastp_result(blastp_file,blastp_result_dir)
# df=df[(df['tcov']>0.5)&(df['pident']>0.6)&(df['qcov']>0.5)]
# group by query
df_group=df.groupby('query')
# get all group name list
group_name_List=df_group.groups.keys()
best_taxanameList=[]
best_taxaIDList=[]
best_taxalineageList=[]
best_taxa_pidList=[]
best_taxa_tcovList=[]
best_taxa_qcovList=[]
sce_indexList=[]



for name in group_name_List:
    group=df_group.get_group(name)
    best_taxaname = group['taxname'].value_counts().index[0]
    best_taxanameID = group['taxid'].value_counts().index[0]
    best_taxalineage = group['taxlineage'].value_counts().index[0]
    best_tax_group=group[group['taxname']==group['taxname']]
    # best_tax_group 根据bits排序，去bits最大的一行
    best_taxa_group_best=best_tax_group.sort_values(by='bits',ascending=False).iloc[0]
    best_taxa_group_hit_pid=best_taxa_group_best['pident']
    best_taxa_group_hit_qcov=best_taxa_group_best['qcov']
    best_taxa_group_hit_tcov=best_taxa_group_best['tcov']
    sce_index=len(group[group['taxname'].str.contains('Saccharomyces')])/len(group)
    sce_indexList.append(sce_index)
    best_taxanameList.append(best_taxaname)
    best_taxaIDList.append(best_taxanameID)
    best_taxa_pidList.append(best_taxa_group_hit_pid)
    best_taxa_tcovList.append(best_taxa_group_hit_tcov)
    best_taxa_qcovList.append(best_taxa_group_hit_qcov)
    best_taxalineageList.append(best_taxalineage)


# save result
df_result=pd.DataFrame({'query':group_name_List,'taxname':best_taxanameList,'taxid':best_taxaIDList,'pident':best_taxa_pidList,'qcov':best_taxa_qcovList,'tcov':best_taxa_tcovList,'sce_index':sce_indexList,'taxlineage':best_taxalineageList})


len(df_result[df_result['taxlineage'].str.contains('Bacteria')])