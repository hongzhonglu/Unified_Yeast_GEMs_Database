# -*- coding: utf-8 -*-
# date : 2023/5/29 
# author : wangh
# file : 6.rename_panID.py
# project : Unified_Yeast_GEMs_Database
'''rename nonref geneID to panID
'''
import pandas as pd
from Bio import SeqIO

# load panID_geneID_dict
df_panID_geneID_dict=pd.read_excel(r"result/panID_geneID.xlsx")
#load all s288c geneID
s288c_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/S288c_R64.fasta","fasta")]

#load pangenome
pangenome=[g for g in SeqIO.parse(r"data/genome/pan1900_new_v4_50_70_filter.fasta","fasta")]
pangenome_v2=[g for g in SeqIO.parse(r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2_filter.fasta","fasta")]
pangenome_v2_newgeneIDlist=[g.id for g in pangenome_v2 if g.id not in df_panID_geneID_dict["geneID"].tolist()]

# add pangenome_v2_new into df_panID_geneID_dict
df_new=pd.DataFrame({"panID":pangenome_v2_newgeneIDlist,"geneID":pangenome_v2_newgeneIDlist})
df_panID_geneID_dict=df_panID_geneID_dict.append(df_new)
df_panID_geneID_dict.to_excel(r"result/panID_geneID.xlsx",index=False)

pangenome_rename=[]
panID_dict={}
i=0
#rename panID as sce0001,sce0002,sce0003...
for g in pangenome:
    if g.id in s288c_geneIDlist:
        a=1
        # pangenome_rename.append(g)
        # panID_dict[g.id]=g.id
    else:
        i+=1
        id = g.id
        numb='0'*(4-len(str(i)))+str(i)
        panID="scepan"+numb
        g.id=panID
        pangenome_rename.append(g)
        panID_dict[id]=panID

with open(r"data/genome/pan1900_new_v4_50_70_filter_rename.fasta","w+") as f:
    for g in pangenome_rename:
        f.write(">"+g.id+"\n"+str(g.seq)+"\n")

df_panID=pd.DataFrame.from_dict(panID_dict,orient="index",columns=["panID"])
df_panID["geneID"]=df_panID.index
df_panID.to_excel(r"result/panID_geneID.xlsx",index=False)

# rename geneMatrix & cnvMatrix
df_geneMatrix=pd.read_csv(r"code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_geneMatrix2.csv",index_col=0)
df_cnvMatrix=pd.read_csv(r"code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_cnvMatrix2.csv",index_col=0)

# rename geneID in index to panID
df_geneMatrix.index=[panID_dict[i] for i in df_geneMatrix.index]
df_cnvMatrix.index=[panID_dict[i] for i in df_cnvMatrix.index]

# save result
df_geneMatrix.to_csv(r"code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_geneMatrix2.csv")
df_cnvMatrix.to_csv(r"code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_cnvMatrix2.csv")