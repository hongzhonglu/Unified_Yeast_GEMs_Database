# -*- coding: utf-8 -*-
# date : 2023/5/23 
# author : wangh
# file : 5_2.parse_allblastp_result.py
# project : Unified_Yeast_GEMs_Database
'''parse all vs all blastp result, and collapse those similar genes and incomplet genes:
1.similar genes: if the 2 representive sequence matched with pid>70% and qcov>50% and tcov>50%, the 2 clusters will be merged
2.incomplet genes: if the representive sequence matched with pid>95% and qcov>90% and tcov<50% & the cluster only exist in less than 5 genes, the query cluster will be
regardes as the fragemented gene of target cluster,and will be merged into the target cluster
'''
import pandas as pd
import os
from Bio import SeqIO


# query,target,pident,bits,qcov,tcov,qlen,tlen,evalue
col=['query','target','pident','bits','qcov','tcov','qlen','tlen','evalue']

# df_blastp=pd.read_csv('code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_all9323_vs_all9323.txt',sep='\t',names=col)
df_blastp=pd.read_csv('code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_v4_50_70_v2_all_vs_all.txt',sep='\t',names=col)

# 1. check incomplet gene: 13 genes
df_blastp=df_blastp[df_blastp['query']!=df_blastp['target']]
#group by query
df_blastp_group=df_blastp.groupby('query')
# choose the best hit for each query
df_blastp_besthit=df_blastp_group.apply(lambda x:x.sort_values('bits',ascending=False).head(1))
df_incomplet_gene=df_blastp_besthit[(df_blastp_besthit['pident']>95)&(df_blastp_besthit['tcov']<0.50)&(df_blastp_besthit['qcov']>0.9)]
# ignore s288c gene in query gene list
df_incomplet_gene=df_incomplet_gene[~df_incomplet_gene['query'].str.startswith('Y')]

# chek the existence of genes in all stains
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
# count how many element in each group
df_strain_count=df_mmseqs_clu_group.apply(lambda x: len(set(x)))
df_strain_count_less_than5=df_strain_count[(df_strain_count<5)&(~df_strain_count.index.str.startswith('s288c'))]
# check how many incoplete genes in less than 5 strains
df_incomplet_gene=df_incomplet_gene[df_incomplet_gene['query'].isin(df_strain_count_less_than5.index.tolist())]
# reset query as index
df_incomplet_gene=df_incomplet_gene.set_index('query')


# 2.check similar genes
df_similar_repseqs=df_blastp[(df_blastp['pident']>70)&(df_blastp['tcov']>0.50)&(df_blastp['qcov']>0.5)]
len(set(df_similar_repseqs['target'].tolist()).intersection(set(df_similar_repseqs['query'].tolist())))
# re clustering those similar genes and choose the longest one as the representive sequence
similar_reclu_dict=dict()
checked_genes = list()
repr_info=dict()
for index,row in df_similar_repseqs.iterrows():
    query=row['query']
    target=row['target']
    i=0
    for old_rep,group in similar_reclu_dict.items():
        if query in group:
            i=1
            break
        elif target in group:
            i=1
            break
    if row['qlen']>row['tlen']:
        representive=query
        rep_len=row['qlen']
    else:
        representive=target
        rep_len=row['tlen']
    if i==0:
        similar_reclu_dict[representive]=[query,target]
        repr_info[representive]=rep_len
    elif i==1:
        if rep_len>repr_info[old_rep]:
            similar_reclu_dict[representive]=similar_reclu_dict.pop(old_rep)
            similar_reclu_dict[representive]=list(set(group+[query,target]))
            repr_info[representive]=rep_len
            repr_info.pop(old_rep)
        else:
            similar_reclu_dict[old_rep]=list(set(group+[query,target]))

remove_similar_seqs=[g for group in similar_reclu_dict.values() for g in group if g not in similar_reclu_dict.keys()]
# discard yeast8 genes from remove_similar_seqs
from cobra.io import read_sbml_model
panYeast_geneIDlist=[g.id for g in read_sbml_model("model/panYeast_v3.xml").genes]
remove_similar_seqs=[g for g in remove_similar_seqs if g not in panYeast_geneIDlist]

remove_incomplete_seqs=df_incomplet_gene.index.tolist()
remove_all_seqs=list(set(remove_similar_seqs+remove_incomplete_seqs))

# update the pan-genome representative sequence
pangenome_filter=[]
for g in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2.fasta","fasta"):
    if g.id not in remove_all_seqs:
        pangenome_filter.append(g)
with open("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2_filter.fasta","w") as f:
    for g in pangenome_filter:
        f.write(">"+g.id+"\n"+str(g.seq)+"\n")


# compare pan1900_new_v4_50_70_v2_filter vs pan1900_new_v4_50_70_filter
pan1800_v1=[g.id for g in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_filter.fasta","fasta")]
pan1800_v2=[g.id for g in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2_filter.fasta","fasta")]

# check difference
len(set(pan1800_v1).difference(set(pan1800_v2)))
len(set(pan1800_v2).difference(set(pan1800_v1)))



# update geneMatrix & cnvMatrix
geneMatrix=pd.read_csv("code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_geneMatrix.csv",index_col=0)
cnvMatrix=pd.read_csv("code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_cnvMatrix.csv",index_col=0)

# collapse incomplete genes
geneMatrix2=geneMatrix.copy()
cnvMatrix2=cnvMatrix.copy()
for g in remove_incomplete_seqs:
    target=df_incomplet_gene.loc[g]['target']
    #update target row by merging incomplete row and target row
    geneMatrix2.loc[target]=geneMatrix2.loc[target]+geneMatrix2.loc[g]
    cnvMatrix2.loc[target]=cnvMatrix2.loc[target]+cnvMatrix2.loc[g]

# collapse similar genes
for g in remove_similar_seqs:
    target=[k for k,v in similar_reclu_dict.items() if g in v][0]
    #update target row by merging incomplete row and target row
    geneMatrix2.loc[target]=geneMatrix2.loc[target]+geneMatrix2.loc[g]
    cnvMatrix2.loc[target]=cnvMatrix2.loc[target]+cnvMatrix2.loc[g]

# remove similar row
geneMatrix2=geneMatrix2.drop(remove_all_seqs)
cnvMatrix2=cnvMatrix2.drop(remove_all_seqs)
# set all value >1 as 1 in geneMatrix2
geneMatrix2[geneMatrix2>1.0]=1.0

geneMatrix.sum().describe()
geneMatrix2.sum().describe()

cnvMatrix.sum().describe()
cnvMatrix2.sum().describe()

# save geneMatrix2 & cnvMatrix2
geneMatrix2.to_csv("code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_geneMatrix2.csv")
cnvMatrix2.to_csv("code/3.pan-genome_construction/2.build_geneMatrix/output/pan1900_new_v4_50_70_blastp_cnvMatrix2.csv")


