# -*- coding: utf-8 -*-
# date : 2023/5/18 
# author : wangh
# file : 2.parse_blastp_result.py
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
blastp_file='pan1900_new_v4_vs_nr.txt'
df=parse_blastp_result(blastp_file,blastp_result_dir)
# remove rows which taxname is root
df=df[~df['taxname'].str.contains('root')]
df=df[(df['pident']>0.5)&(df['qcov']>0.4)&(df['tcov']>0.4)]
# group by query
df_group=df.groupby('query')
# get all group name list
group_name_List=df_group.groups.keys()
besthit_idList=[]
besthit_bitscoreList=[]
besthit_taxnameList=[]
besthit_taxidList=[]
besthit_taxlineageList=[]
besthit_qcovList=[]
besthit_tcovList=[]
besthit_pidentList=[]
sce_indexList=[]


for name in group_name_List:
    group=df_group.get_group(name)
    # sort by bitscore
    group=group.sort_values(by='bits',ascending=False)
    # get the best hit
    best_hit=group.iloc[0]
    # get the sce index
    sce_index=len(group[group['taxname'].str.contains('Saccharomyces')])/len(group)
    sce_indexList.append(sce_index)
    # get the best hit info
    besthit_idList.append(best_hit['target'])
    besthit_bitscoreList.append(best_hit['bits'])
    besthit_taxnameList.append(best_hit['taxname'])
    besthit_taxidList.append(best_hit['taxid'])
    besthit_taxlineageList.append(best_hit['taxlineage'])
    besthit_qcovList.append(best_hit['qcov'])
    besthit_tcovList.append(best_hit['tcov'])
    besthit_pidentList.append(best_hit['pident'])


# create a dataframe to store the result
df_result=pd.DataFrame({'query':group_name_List,
                        'besthit_id':besthit_idList,
                        'besthit_bitscore':besthit_bitscoreList,
                        'besthit_taxname':besthit_taxnameList,
                        'besthit_taxid':besthit_taxidList,
                        'besthit_taxlineage':besthit_taxlineageList,
                        'qcov':besthit_qcovList,
                        'tcov':besthit_tcovList,
                        'pident':besthit_pidentList,
                        'sce_index':sce_indexList
                        })

# check Bacteria contamination
df_contamination=df_result[df_result['besthit_taxlineage'].str.contains('Bacteria')]
df_contamination[df_contamination['query'].str.contains('CDH.re')]
len(df_result[df_result['query'].str.contains('CDH.re')])

# check the completeness of pan1900_new_all genes
len(df_result[(df_result['pident']>0.9)&(df_result['tcov']<0.4)&(df_result['qcov']>0.9)])

# check difference with keep and remove gene result
df_remove_result=df_result[df_result['query'].str.contains('remove')]
df_keep_result=df_result[df_result['query'].str.contains('keep')]
len(df_keep_result[(df_keep_result['pident']>0.9)&(df_keep_result['tcov']<0.4)&(df_keep_result['qcov']>0.9)])
len(df_remove_result[(df_remove_result['pident']>0.9)&(df_remove_result['tcov']<0.4)&(df_remove_result['qcov']>0.9)])

# check the contamination
df_rmconta=df_result[df_result['besthit_taxname'].str.contains('Saccharomyces')]

# remove incomplete genes
df_rmconta_complete=df_rmconta[~((df_rmconta['pident']>0.9)&(df_rmconta['tcov']<0.4)&(df_rmconta['qcov']>0.9))]

# check dose it include all s288c gene
pan_rename_index=pd.read_csv(r'output/pan1900_new_all9323_index.csv',index_col=0)
s288cIDlist=pan_rename_index[pan_rename_index['geneID'].str.startswith('Y')]['renameID'].tolist()

len(set(s288cIDlist).intersection(set(df_rmconta_complete['query'].tolist())))
len(set(s288cIDlist).intersection(set(df_result['query'].tolist())))

df_s288c_result=df_result[df_result['query'].isin(s288cIDlist)]
len(df_s288c_result[df_s288c_result['taxname'].str.contains('Saccharomyces')])
len(df_s288c_result[df_s288c_result['sce_index']<0.005])



