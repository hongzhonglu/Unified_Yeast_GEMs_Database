# -*- coding: utf-8 -*-
# date : 2023/5/8 
# author : wangh
# file : 6.build_blastp_geneMatrix.py
# project : Unified_Yeast_GEMs_Database
import os
import pandas as pd
import sys
import numpy as np
from Bio import SeqIO
from tqdm import tqdm
sys.path.append(r"D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database\code")
from mainFunction import get_gene_lens
import re


def parse_blastp_result(blastp_file,blastp_dir,query,query_dir,subject,subject_dir):
    """
    parse the blastp result file and calculate COV for query and subject
    :param blastp_file: blastp result file
    :param query: query file name
    :param query_dir: query file directory
    :param subject: subject file name
    :param subject_dir: subject file directory
    :return: blastp_file with COV for query and subject
    """
    # load blastp result file
    df_blastp_file = pd.read_csv(blastp_dir+blastp_file, sep="\t", header=None, index_col=0)
    columns = ["subject", "identity", "alignment length", "mismatches", "gap opens", "q_start", "q_end", "s_start",
               "s_end", "evalue", "bit score"]
    df_blastp_file.columns = columns

    # get gene length for query and subject
    query_lens=get_gene_lens(query,in_folder=query_dir)
    subject_lens=get_gene_lens(subject,in_folder=subject_dir)
    query_lens.set_index("gene",inplace=True)
    subject_lens.set_index("gene",inplace=True)

    # map query lens to blastp_file and name the column as "query_lens"
    df_blastp_file=df_blastp_file.join(query_lens,how="left")
    df_blastp_file.rename(columns={"gene_length":"query_lens"},inplace=True)
    # map subject lens to blastp_file and name the column as "subject_lens"
    df_blastp_file=df_blastp_file.join(subject_lens,how="left",on="subject")
    df_blastp_file.rename(columns={"gene_length":"subject_lens"},inplace=True)

    # calculate COV for query and subject
    df_blastp_file["query_cov"]=(df_blastp_file["q_end"]-df_blastp_file["q_start"]+1)/df_blastp_file["query_lens"]
    df_blastp_file["subject_cov"]=(df_blastp_file["s_end"]-df_blastp_file["s_start"]+1)/df_blastp_file["subject_lens"]

    return df_blastp_file


def build_geneMatrix(strainList,geneList,pangenome,pangenome_dir,proteome_dir,bbh_file_dir,cov_cutoff=0.6,pid_cutoff=0.8):
    """
    build gene matrix for each strain
    :param strainList: strain list
    :param geneList: gene list
    :param bbh_file_dir: bbh file directory
    :return: gene matrix for each strain
    """
    cnvMatrix=pd.DataFrame(index=geneList,columns=strainList)
    # pan_re=r'_new_cov.*\.fasta$'
    for strain in tqdm(strainList):
        try:
            blastp_file=parse_blastp_result(blastp_file=strain.rstrip(".fa")+"_vs_"+"pan1800_v2"+".txt",
                                  blastp_dir=bbh_file_dir,
                                  query=strain,
                                  query_dir=proteome_dir,
                                  subject=pangenome,
                                  subject_dir=pangenome_dir)
        except:
            continue
        blastp_file=blastp_file[(blastp_file["query_cov"]>=cov_cutoff)&(blastp_file["subject_cov"]>=cov_cutoff)&(blastp_file["identity"]>=pid_cutoff)]
        # blastp_file中的subject根据index进行groupby，取identity最大的那个,如果最大值有多行，只取一行
        blastp_file["query"]=blastp_file.index.str.split("/").map(lambda x:x[0])
        df_blastp_group=blastp_file.groupby(by="query")
        blastp_file2=df_blastp_group.apply(lambda x:x.sort_values(by="identity",ascending=False).head(1))
        # blastp_file2.set_index("query",inplace=True)
        df_gene_counts = blastp_file2['subject'].value_counts()
        # geneList=list(set(blastp_file2["subject"].tolist()))
        cnvMatrix.loc[geneList,strain]=df_gene_counts
    cnvMatrix.fillna(0,inplace=True)
    cnvMatrix=cnvMatrix.astype(int)
    geneMatrix=cnvMatrix.copy()
    geneMatrix[geneMatrix>0]=1
    return geneMatrix,cnvMatrix


# get geneMatrix
import os
from Bio.SeqIO import parse
strainList=os.listdir(r"data/genome/predicted_proteomes/combined_proteomes_old")
pangenome_dir=r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/"
pangenome="pan1900_new_v4_50_70_v2_filter.fasta"
blastp_file_dir=r"code/3.pan-genome_construction/2.build_geneMatrix/output/blastp_vs_pan1900_new_v4_50_70_v2_filter/"
geneList=[g.id for g in parse(pangenome_dir+pangenome,"fasta")]
proteome_dir=r"data/genome/predicted_proteomes/combined_proteomes_old/"


# # set 3 hours sleep time
# import time
# time.sleep(3*60*60)

# rename all files name in blastp_file_dir by changing "pan1900_new_v4_50_70_v2_filter" to "pan1800_v2"
for file in os.listdir(blastp_file_dir):
    os.rename(blastp_file_dir+file,blastp_file_dir+file.replace("pan1900_new_v4_50_70_v2_filter","pan1800_v2"))

# get geneMatrix and cnvMatrix
geneMatrix,cnvMatrix=build_geneMatrix(strainList=strainList,
                            geneList=geneList,
                            pangenome=pangenome,
                            pangenome_dir=pangenome_dir,
                            proteome_dir=proteome_dir,
                            bbh_file_dir=blastp_file_dir,
                            cov_cutoff=0.5,
                            pid_cutoff=0.7)

# rebuild some lost strain which the sum of geneMatrix is 0
strainList_lost=geneMatrix.columns[geneMatrix.sum()==0].tolist()
geneMatrix_lost,cnvMatrix_lost=build_geneMatrix(strainList=strainList_lost,
                            geneList=geneList,
                            pangenome=pangenome,
                            pangenome_dir=pangenome_dir,
                            proteome_dir=proteome_dir,
                            bbh_file_dir=blastp_file_dir,
                            cov_cutoff=0.5,
                            pid_cutoff=0.7)

# replace columns in geneMatrix and cnvMatrix by columns in geneMatrix_lost and cnvMatrix_lost
geneMatrix[geneMatrix.columns[geneMatrix.sum()==0]]=geneMatrix_lost
cnvMatrix[cnvMatrix.columns[cnvMatrix.sum()==0]]=cnvMatrix_lost

# replace geneID with panID
df_panID_geneID_dict=pd.read_excel("result/panID_geneID.xlsx")
df_panID_geneID_dict.set_index("geneID",inplace=True)
geneMatrix.index=geneMatrix.index.map(lambda x:df_panID_geneID_dict.loc[x,"panID"])
cnvMatrix.index=cnvMatrix.index.map(lambda x:df_panID_geneID_dict.loc[x,"panID"])

geneMatrix.to_csv("data/geneMatrix/pan1900_new_50_70_v2_blastp_geneMatrix.csv")
cnvMatrix.to_csv("data/geneMatrix/pan1900_new_50_70_v2_blastp_cnvMatrix.csv")

# check which column has all 0: AKB_2.re.fa
geneMatrix.loc[:,geneMatrix.sum()==0]

# check gene_ratio
gene_ratio=geneMatrix.sum(axis=1)/geneMatrix.shape[1]
len(gene_ratio[gene_ratio==1])

# load all strain data
all_strain_info=pd.read_excel("data/1897_strains_info.xlsx", index_col=0)
kept_strainList=all_strain_info[all_strain_info["remove"]==False]["genome_id"].tolist()
kept_strainList=[s+".fa" for s in kept_strainList]
kept_strainList.append('s288c_R64.fa')


strainList2=[s for s in strainList if s in kept_strainList]
geneMatrix2=geneMatrix.loc[:,strainList2]
cnvMatrix2=cnvMatrix.loc[:,strainList2]
# remove column AKB_2.re.fa
# geneMatrix2.drop(columns="AKB_2.re.fa",inplace=True)
# cnvMatrix2.drop(columns="AKB_2.re.fa",inplace=True)
gene_ratio2=geneMatrix2.sum(axis=1)/geneMatrix2.shape[1]
len(gene_ratio2[gene_ratio2==1])
geneMatrix2.sum().describe()
for i in np.linspace(0.98,1,21):
    print(i,len(gene_ratio2[gene_ratio2>=i]))



