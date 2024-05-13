# -*- coding: utf-8 -*-
# date : 2023/5/7 
# author : wangh
# file : 2_3.filter_nonref_genes.py
# project : Unified_Yeast_GEMs_Database
#filter the nonref genes by remove bad strains and bad gene in 3 steps:
# 1. remove strains with mapped s288c gene number < 0.85 * total mapped s288c gene number or nonref gene number > 0.15 * total mapped s288c gene number
# 2. remove strain according genome completeness ecaluated by BUSCO
# 3. remove incomplted genes which have 'X' in the sequence

import pandas as pd
from Bio import SeqIO

#load
all_nonref=[g for g in SeqIO.parse(r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900_s288c_nonref_v4.fasta","fasta")]

# 1. remove strains with mapped s288c gene number < 0.85 * median mapped s288c gene number
df_sce_nonref = pd.read_csv(r"code\3.pan-genome_construction\1.new_pan-genome_pipeline\output\1900_s288c_nonref_gene_count_v4.csv",index_col=0)
step1_toremove_strainList=df_sce_nonref[(df_sce_nonref['s288c_gene_numb']<6716*0.85)|(df_sce_nonref['nonref_gene_numb']>6716*0.15)].index.tolist()
step1_toremove_strainList=[g.strip(".fa") for g in step1_toremove_strainList]
step1_toremove_geneList=[]
for gene in all_nonref:
    if gene.id.split("|")[0] in step1_toremove_strainList:
        step1_toremove_geneList.append(gene.id)


# 2. remove strain according genome completeness ecaluated by BUSCO
all_strain_info=pd.read_excel("data/1897_strains_info.xlsx", index_col=0)
kept_strainList=all_strain_info[all_strain_info["remove"]==False]["genome_id"].tolist()
# add s288c
kept_strainList.append("s288c")
step2_toremove_geneList=[]
for gene in all_nonref:
    if gene.id.split("|")[0] not in kept_strainList:
        step2_toremove_geneList.append(gene.id)


# 3. remove incomplted genes which have 'X' in the sequence
step3_toremove_geneList=[]
for gene in all_nonref:
    length=len(gene.seq)
    if gene.seq.count("X")/length>0:
        step3_toremove_geneList.append(gene.id)


# 4. remove genes from strain CDH.re which got amounts of contaminated genes
step4_toremove_geneList=[]
for gene in all_nonref:
    if gene.id.startswith("CDH.re"):
        step4_toremove_geneList.append(gene.id)



# write the filtered_nonref_genes .fasta file
all_toremove_geneList=list(set(step1_toremove_geneList+step2_toremove_geneList+step3_toremove_geneList+step4_toremove_geneList))
with open("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900_s288c_nonref_v4_filtered.fasta","w+") as f:
    for gene in all_nonref:
        if gene.id not in all_toremove_geneList:
            f.write(">"+gene.id+"\n")
            f.write(str(gene.seq)+"\n")



