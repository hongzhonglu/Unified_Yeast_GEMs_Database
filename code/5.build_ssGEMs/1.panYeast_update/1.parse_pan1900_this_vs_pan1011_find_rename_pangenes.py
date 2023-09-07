# -*- coding: utf-8 -*-
# date : 2023/5/15 
# author : wangh
# file : pan1900_this_vs_pan1011.py
# project : Unified_Yeast_GEMs_Database
import sys
sys.path.append(r"D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database\code")
import pandas as pd
from Bio import SeqIO
from mainFunction import get_gene_lens
from cobra.io import read_sbml_model
import time


def parse_bbh_result(bbh_file,bbh_dir,query,query_dir,subject,subject_dir):
    """
    parse the bbh result file and calculate COV for query and subject
    :param bbh_file: bbh result file
    :param query: query file name
    :param query_dir: query file directory
    :param subject: subject file name
    :param subject_dir: subject file directory
    :return: bbh_file with COV for query and subject
    """
    # load bbh result file
    df_bbh_file = pd.read_csv(bbh_dir+bbh_file, sep="\t", header=None, index_col=0)
    columns = ["subject", "identity", "alignment length", "mismatches", "gap opens", "q_start", "q_end", "s_start",
               "s_end", "evalue", "bit score"]
    df_bbh_file.columns = columns

    # get gene length for query and subject
    query_lens=get_gene_lens(query,in_folder=query_dir)
    subject_lens=get_gene_lens(subject,in_folder=subject_dir)
    query_lens.set_index("gene",inplace=True)
    subject_lens.set_index("gene",inplace=True)

    # map query lens to bbh_file and name the column as "query_lens"
    df_bbh_file=df_bbh_file.join(query_lens,how="left")
    df_bbh_file.rename(columns={"gene_length":"query_lens"},inplace=True)
    # map subject lens to bbh_file and name the column as "subject_lens"
    df_bbh_file=df_bbh_file.join(subject_lens,how="left",on="subject")
    df_bbh_file.rename(columns={"gene_length":"subject_lens"},inplace=True)

    # calculate COV for query and subject
    df_bbh_file["query_cov"]=(df_bbh_file["q_end"]-df_bbh_file["q_start"]+1)/df_bbh_file["query_lens"]
    df_bbh_file["subject_cov"]=(df_bbh_file["s_end"]-df_bbh_file["s_start"]+1)/df_bbh_file["subject_lens"]

    return df_bbh_file


df_bbh=parse_bbh_result(bbh_file="pan1800_v2_vs_pan1011_bbh.txt",
                        bbh_dir=r"code/3.pan-genome_construction/3.pan-genome_comparison/output/",
                        query="pan1800_50_70_v2.fasta",
                        query_dir=r"data/genome/",
                        subject="pan1011_v1.fasta",
                        subject_dir=r"data/genome/")


cov_cutoff=0.5
pid_cutoff=0.7
df_bbh=df_bbh[(df_bbh["query_cov"]>=cov_cutoff)&(df_bbh["subject_cov"]>=cov_cutoff)&(df_bbh["identity"]>=pid_cutoff)]

# remove S288c genes from pan1900 index
s288c_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/S288c_R64.fasta","fasta")]
df_bbh=df_bbh[~df_bbh.index.isin(s288c_geneIDlist)]
df_bbh=df_bbh[~df_bbh["subject"].isin(s288c_geneIDlist)]
df_panID_pan1011=df_bbh["subject"]
df_panID_pan1011.to_excel(r"result/panID_pan1011.xlsx")

# only keep genes include in panYeast
panYeast=read_sbml_model("model/panYeast_v3.xml")
pangeneIDlist=[i.id for i in panYeast.genes]
df_rename_bbh=df_bbh[df_bbh['subject'].isin(pangeneIDlist)]  # 7 gene need to be renamed
df_panID=pd.read_excel(r"result/panID_geneID.xlsx")
# set geneID as index
df_panID.set_index("geneID",inplace=True)

df_rename_gene=pd.DataFrame(columns=["pan1011","panID","geneID"])
df_rename_gene["pan1011"]=df_rename_bbh["subject"].tolist()
df_rename_gene["geneID"]=df_rename_bbh.index.tolist()
df_rename_gene["panID"]=df_rename_gene["geneID"].map(df_panID["panID"].to_dict())
df_rename_gene["geneName"]=df_rename_gene["pan1011"].apply(lambda x:panYeast.genes.get_by_id(x).name)
df_rename_gene.to_excel(r"code/5.build_ssGEMs/output/panYeast_rename_panID.xlsx",index=False)




