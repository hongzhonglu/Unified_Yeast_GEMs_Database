# -*- coding: utf-8 -*-
# date : 2023/5/30 
# author : wangh
# file : update_gpr_for_s288c_isoenzyme.py
# project : Unified_Yeast_GEMs_Database
'''some gene in YeastGEM has been collapsed into one gene clusters according to sequence similarity.
So, it should be found those redundant isoenyzme genes to update GPR with panID'''
from Bio import SeqIO
from cobra.io import read_sbml_model
import pandas as pd
import pickle

panYeast=read_sbml_model("model/panYeast_v3.xml")


s288c_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/S288c_R64.fasta","fasta")]
pan1011_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/pan1011_v1.fasta","fasta")]
pan1900_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/pan1800_50_70_v2.fasta","fasta")]
pan1011_nons288c=list(set(pan1011_geneIDlist)-set(s288c_geneIDlist))
pan1011_non1900=list(set(pan1011_geneIDlist)-set(pan1900_geneIDlist))

remove_geneList=[gene.id for gene in panYeast.genes if (gene.id in pan1011_non1900)&~(gene.id in pan1011_nons288c)]

# parse mmseqs clustering result
df_clu=pd.read_csv(r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900sce_s288c+nonref_v4_filtered_cov50_pid70_cluster.tsv",sep='\t',header=None)
# remove 's288c|' in all df_clu
df_clu[0]=df_clu[0].str.replace('s288c\|','')
df_clu[1]=df_clu[1].str.replace('s288c\|','')
df_clu_s288c_add=df_clu[(df_clu[1].isin(remove_geneList))&(~df_clu[0].isin(pan1900_geneIDlist))]
df_clu_s288c_add.columns=['panID','s288cID']
# load genome cluster info
with open("data/genome/pan1800_50_70_v2_cluster.pkl","rb") as f:
    pan1800_v2_cluster=pickle.load(f)

panIDlist=list()
for geneID in df_clu_s288c_add['s288cID'].tolist():
    for repID in pan1800_v2_cluster.keys():
        if geneID in pan1800_v2_cluster[repID]:
            panIDlist.append(repID)
            break
df_clu_s288c_add['panID']=panIDlist
# remove rows with same panID and s288cID
df_clu_s288c_add=df_clu_s288c_add.drop_duplicates(subset=['panID','s288cID'])
df_clu_s288c_add.to_excel(r"code/5.build_ssGEMs/output/df_clu_s288c_add.xlsx",index=False)

# find all rxn need to be updated GPR
rxnIDlist=list()
panIDlist=list()
s288cgeneIDlist=list()
oldgprlist=list()
rxnNameList=list()

for geneID in df_clu_s288c_add['s288cID'].tolist():
    rxns=panYeast.genes.get_by_id(geneID).reactions
    panID=df_clu_s288c_add[df_clu_s288c_add['s288cID']==geneID]['panID'].values[0]
    for rxn in rxns:
        rxnIDlist.append(rxn.id)
        panIDlist.append(panID)
        s288cgeneIDlist.append(geneID)
        oldgprlist.append(rxn.gene_reaction_rule)
        rxnNameList.append(rxn.name)

panYeast_update_gpr_s288c=pd.DataFrame({'rxnID':rxnIDlist,'panID':panIDlist,'s288cID':s288cgeneIDlist,'oldGPR':oldgprlist,'rxnName':rxnNameList})
panYeast_update_gpr_s288c.to_excel(r"code/5.build_ssGEMs/output/panYeast_update_gpr_s288c.xlsx",index=False)

# recheck
panYeast_update_gpr_s288c=pd.read_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_s288c.xlsx")
panYeast_update_gpr_s288c=panYeast_update_gpr_s288c[panYeast_update_gpr_s288c['s288cID'].isin(remove_geneList)]
panYeast_geneidList=[gene.id for gene in panYeast.genes]
panYeast_update_gpr_s288c_check=panYeast_update_gpr_s288c[panYeast_update_gpr_s288c['panID'].isin(panYeast_geneidList)]
updateList=list()
# check does the replace gene and panID have the same GPR. If not same, the GPR should be updated
for i in panYeast_update_gpr_s288c_check.index:
    rm_geneID=panYeast_update_gpr_s288c.loc[i,'s288cID']
    panID=panYeast_update_gpr_s288c.loc[i,'panID']
    rm_gene_gpr=panYeast.genes.get_by_id(rm_geneID).reactions
    panID_gpr=panYeast.genes.get_by_id(panID).reactions
    if rm_gene_gpr!=panID_gpr:
        # print('not same')
        print(panID,":",rm_geneID)
        print(str(panID_gpr),":")
        print(str(rm_gene_gpr))
        updateList.append(i)
    else:
        # print('same')
        pass

# remove those s288cID in ignoreList
panYeast_update_gpr_s288c=panYeast_update_gpr_s288c.loc[updateList,:]
panYeast_update_gpr_s288c.to_excel(r"code/5.build_ssGEMs/output/panYeast_update_gpr_s288c.xlsx",index=False)
