# -*- coding: utf-8 -*-
# date : 2023/6/1 
# author : wangh
# file : extract_seq_pan1800_v2_vs_v1.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
from Bio import SeqIO

geneMatrix=pd.read_csv(r"code/3.pan-genome_construction/2.build_geneMatrix/output/pan1800_new_v4_50_70_blastp_geneMatrix.csv",index_col=0)
pan1800_v1_geneIDlist=geneMatrix.index.tolist()

pan1800_v2=[g for g in SeqIO.parse(r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1800_new_v4_50_70_v2_filter.fasta","fasta")]

pan1800_v2_new=[g for g in pan1800_v2 if g.id not in pan1800_v1_geneIDlist]

with open("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1800_v2_new.fasta","w+") as f:
    for g in pan1800_v2_new:
        f.write(">"+g.id+"\n"+str(g.seq)+"\n")




