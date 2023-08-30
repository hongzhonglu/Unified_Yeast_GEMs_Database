# -*- coding: utf-8 -*-
# date : 2023/5/7 
# author : wangh
# file : 4_2.chek_clustering_result.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import os
from Bio import SeqIO

cluster_result_path="code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900sce_s288c+nonref_v4_filtered_cov50_pid70_cluster.tsv"

# load all representative sequences
all_pan_genome=[g for g in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900sce_s288c+nonref_v4_filtered_cov50_pid70_rep_seq.fasta","fasta")]

df_mmseqs_clu= pd.read_csv(cluster_result_path, sep="\t", header=None)
df_mmseqs_clu.columns = ["representive_id", "gene_id"]
strainlist=list(set([g.split("|")[0] for g in df_mmseqs_clu["gene_id"].tolist()]))

# build 1811 strains geneMatrix
df_mmseqs_clu_group=df_mmseqs_clu.groupby("representive_id")
df_mmseqs_clu_group=df_mmseqs_clu_group.apply(lambda x: x["gene_id"].tolist())
df_mmseqs_clu_group=df_mmseqs_clu_group.apply(lambda x: [g.split("|")[0] for g in x])
# df_nonref_includ_s288c=df_mmseqs_clu_group[(~df_mmseqs_clu_group.index.str.startswith("s288c"))&(df_mmseqs_clu_group.apply(lambda x: "s288c" in x))]
# df_mmseqs_clu_group_2=df_mmseqs_clu_group[(df_mmseqs_clu_group.apply(lambda x: len(set(x)))==2)&(~df_mmseqs_clu_group.index.str.startswith("s288c"))]

clustr_strain_count= {}
for i in df_mmseqs_clu_group.index.tolist():
    clustr_strain_count[i]=len(set(df_mmseqs_clu_group[i]))
df_clustr_strain_count=pd.Series(clustr_strain_count)
unique_clustr_list=df_clustr_strain_count[df_clustr_strain_count==1].index.tolist()
s288c_all=[clustr for clustr in df_mmseqs_clu_group if "s288c" in clustr]
s288c_unique=[g for g in unique_clustr_list if g.startswith("s288c")]
non_s288c_unique=[g for g in unique_clustr_list if not g.startswith("s288c")]
print("total pangenome size: ", len(df_clustr_strain_count[(df_clustr_strain_count>1)|(df_clustr_strain_count.index.str.startswith("s288c"))]))
filter_rep_idList=[g for g in df_clustr_strain_count[(df_clustr_strain_count>1)|(df_clustr_strain_count.index.str.startswith("s288c"))].index.tolist()]
with open("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_50_70_ori.fasta", "w") as f:
    for g in all_pan_genome:
        if g.id in filter_rep_idList:
            f.write(">"+g.id+"\n"+str(g.seq)+"\n")

# build geneMatrix
df_geneMatrix=pd.DataFrame(index=df_clustr_strain_count.index.tolist(), columns=strainlist)
for i in df_clustr_strain_count.index.tolist():
    df_geneMatrix.loc[i, list(set(df_mmseqs_clu_group[i]))]=1
df_geneMatrix=df_geneMatrix.fillna(0)
gene_ratios=df_geneMatrix.sum(axis=1)
len(df_geneMatrix[df_geneMatrix.index.str.startswith("s288c")])
len(gene_ratios[gene_ratios>1800])

# write pan-genome sequence
from Bio import SeqIO
pan_geneIDlist=df_clustr_strain_count[(df_clustr_strain_count>1)|(df_clustr_strain_count.index.str.startswith("s288c"))].index.tolist()
all_records=[g for g in SeqIO.parse("code/3.pan-genome_construction/new_pan-genome_pipeline/output/1900sce_s288c+nonref_v2_filtered_cov50_pid60_mode0_cluster.tsv", "fasta")]
with open("code/3.pan-genome_construction/new_pan-genome_pipeline/output/pan1900_new_cov60_pid60.fasta", "w") as f:
    for g in all_records:
        if g.id in pan_geneIDlist:
            f.write(">"+g.id+"\n"+str(g.seq)+"\n")


