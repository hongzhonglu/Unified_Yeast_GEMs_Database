# -*- coding: utf-8 -*-
# date : 2023/6/1 
# author : wangh
# file : rename_pangenome_v2.py
# project : Unified_Yeast_GEMs_Database
from Bio import SeqIO
import pandas as pd

df_panID_geneID_dict=pd.read_excel("result/panID_geneID.xlsx")
df_panID_geneID_dict.set_index("geneID",inplace=True)

pan1800_v2=[g for g in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2_filter.fasta","fasta")]
with open("data/genome/pan1800_50_70_v2.fasta","w") as f:
    for g in pan1800_v2:
        geneID=g.id
        panID=df_panID_geneID_dict.loc[geneID,"panID"]
        f.write(">"+panID+"\n"+str(g.seq)+"\n")



