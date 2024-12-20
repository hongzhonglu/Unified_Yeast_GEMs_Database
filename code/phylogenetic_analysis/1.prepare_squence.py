# -*- coding: utf-8 -*-
# date : 2024/10/11 
# author : wangh
# file : 1.prepare_squence.py
# project : Unified_Yeast_GEMs_Database
'''Prepare sequence data for building phylogenetic tree.
1. pan-genome presence/absence data.
2. Core genes protein sequences'''
import os
import pandas as pd
from Bio import SeqIO
import tqdm

def get_sequence_by_id(fasta_file, sequence_id):
    # 读取FASTA文件并遍历序列
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == sequence_id:
            return str(record.seq)  # 返回序列内容
    return None  # 如果没有找到匹配的ID，则返回None

# 1. pan-genome presence/absence data as fasta format
cnvMatrix=pd.read_csv(r'data/geneMatrix/pan1800_v2_blastp_50_70_cnvMatrix.csv',index_col=0)
# fill nan as 0
cnvMatrix.fillna(0,inplace=True)
# set as int
cnvMatrix=cnvMatrix.astype(int)


# ignore genes with max gene copy number more than 9
n=0
for i in cnvMatrix.index:
    if cnvMatrix.loc[i].max()>9:
        n+=1
        # print(cnvMatrix.loc[i].describe())
        # remove i row
        cnvMatrix.drop(i,inplace=True)

# ignore single-copy core genes and unique genes -- 90% is 1 or 0
n=0
threshold=0.95
for i in cnvMatrix.index:
    if cnvMatrix.loc[i].value_counts().max()>cnvMatrix.shape[1]*threshold:
        n+=1
        # print(cnvMatrix.loc[i].describe())
        # remove i row
        cnvMatrix.drop(i,inplace=True)

with open(r'code/phylogenetic_analysis/output/pangenome_existence.fa','w') as f:
    for i in cnvMatrix.columns:
        strainName=i.rstrip('.fa')
        # convert cnvMatrix[i] to string
        # seq=''.join(cnvMatrix[i].map(mapping).tolist())
        seq=''.join(cnvMatrix[i].astype(str).tolist())
        if len(seq)!=len(cnvMatrix):
            print(f'{strainName}:{len(seq)}')
            print(cnvMatrix[i].value_counts())
        f.write('>'+strainName+'\n')
        f.write(seq+'\n')


# 2.extract gene sequence for each ORFs
proteome_dir=r"data/genome/sce/pep/"
output_dir=r'code/phylogenetic_analysis/output/core_ORFs/'
# define the range of genes
threshold=0.5
# load genenameMatrix
geneMatrix=pd.read_csv(r'code/3.pan-genome_construction/2.build_geneMatrix/output/pan1800_new_50_70_blastp_genenameMatrix.csv',index_col=0)
panID_geneID=pd.read_excel(r'result/panID_geneID.xlsx')
panID_geneID_dict=dict(zip(panID_geneID['geneID'],panID_geneID['panID']))
geneMatrix.index=geneMatrix.index.map(panID_geneID_dict)

strainList=os.listdir(proteome_dir)
# single_copy_core_geneList=cnvMatrix[(cnvMatrix.max(axis=1)<=1)&(cnvMatrix.isin([0]).sum(axis=1)<=cnvMatrix.shape[1]*(1-core_threshold))].index.tolist()
core_geneList=cnvMatrix[(cnvMatrix.isin([0]).sum(axis=1)<=cnvMatrix.shape[1]*(1-threshold))].index.tolist()
# for orf in tqdm.tqdm(core_geneList):
#     orf_genes=geneMatrix.loc[orf]
#     # remove nan
#     orf_genes=orf_genes[~orf_genes.isna()]
#     # check if the file orf.fa exists in the dir
#     if os.path.exists(output_dir+orf+'.fa'):
#         continue
#     with open(output_dir+orf+'.fa','w') as f:
#         for g in orf_genes:
#             strainName=g.split('|')[0]
#             fileName=strainName+'.fa'
#             if fileName in strainList:
#                 geneName=g.split('|')[1]
#                 seq=get_sequence_by_id(proteome_dir+fileName,g)
#                 if seq:
#                     f.write('>'+strainName+'\n')
#                     f.write(seq+'\n')
#                 else:
#                     print(f"Can not find seq {g}")

from concurrent.futures import ThreadPoolExecutor
def process_orf(orf):
    orf_genes = geneMatrix.loc[orf].dropna()
    output_path = f"{output_dir}{orf}.fa"

    if os.path.exists(output_path):
        return  # 文件已存在，跳过

    with open(output_path, 'w') as f:
        for g in orf_genes:
            strainName = g.split('|')[0]
            fileName = f"{strainName}.fa"
            if fileName in strainList:
                geneName = g.split('|')[1]
                seq = get_sequence_by_id(proteome_dir + fileName, g)
                if seq:
                    f.write(f'>{strainName}\n{seq}\n')
                else:
                    print(f"Can not find seq {g}")

num_threads = 10
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    list(tqdm.tqdm(executor.map(process_orf, core_geneList), total=len(core_geneList)))






