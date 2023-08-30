# -*- coding: utf-8 -*-
# date : 2023/5/10 
# author : wangh
# file : 4_3.extract_representive_sequences.py
# project : Unified_Yeast_GEMs_Database
# extract representative sequences from clustering result:
# filter: if only 1 non-s288c strain exist in a cluster, then remove it.
# representive sequence selection rules:
# 1. In the cluster, if only 1 sequnce belong to the s288c, the s288c sequence is the representive sequence
# 2. In the cluster, if more than 1 sequnce belong to the s288c, firstly,the gene already existed in Yeast8 will be choosen, and then the longest s288c sequence is the representive sequence
# 3. In the cluster, if no s288c sequence, the original representive sequence from mmseqs2 result is selected
from Bio import SeqIO
import pandas as pd
from cobra.io import read_sbml_model

all_seq=[record for record in SeqIO.parse("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900sce_s288c+nonref_v4_filtered_cov50_pid70_all_seqs.fasta","fasta")]
s288c_geneIDs=[g.id for g in SeqIO.parse("data/genome/S288c_R64.fasta","fasta")]
panYeast_s288c_geneIDs=[g.id for g in read_sbml_model("model/yeastGEM.xml").genes if g.id in s288c_geneIDs]
# check YHR208W in panYeast_s288c_geneIDs


# get all cluster dict: original representive sequence ID as key, all sequences in the cluster as value formed as list
all_seq_dict={}
for g in all_seq:
    if len(g.seq)==0:
        ori_repID=g.id
        all_seq_dict[g.id]=[]
        continue
    else:
        all_seq_dict[ori_repID].append(g)

# new_rep_seq_dict: new representive sequence ID as key, all sequences in the cluster as value formed as list
count_panYeast1=0
count_panYeast_more1=0
panYeast_gene_collapse={}
new_all_seq_dict={}
new_rep_seq_dict={}
for ori_repID,seqList in all_seq_dict.items():
    strainList=[seq.id.split("|")[0] for seq in seqList]
    seqIDlist=[seq.id for seq in seqList]

    if len(seqIDlist)==1:
        # if only one sequence in the cluster, it is the representive sequence
        new_rep_seq_dict[ori_repID]=seqList[0]
        new_all_seq_dict[ori_repID]=seqList
        continue

    elif len(seqIDlist)>1:
        # count the number of "s288c" in the strainList
        s288c_geneList=[g for g in seqList if "s288c" in g.id]
        if len(s288c_geneList)==1:
            # if only 1 s288c gene exist, choose it as representive sequence
            new_rep_gene=s288c_geneList[0]
            new_rep_id=new_rep_gene.id
            new_rep_seq_dict[new_rep_id]=new_rep_gene
            new_all_seq_dict[new_rep_id]=seqList
            continue
        elif len(s288c_geneList)>1:
            # count the number of s288c gene already existed in Yeast8
            panYeast_s288c_geneList=[g for g in s288c_geneList if g.id.strip("s288c|") in panYeast_s288c_geneIDs]
            if len(panYeast_s288c_geneList)==0:
                # if non Yeast8 gene exist, choose the longest s288c sequence as representive sequence
                s288c_geneList.sort(key=lambda x:len(str(x.seq)),reverse=True)
                new_rep_gene=s288c_geneList[0]
                new_rep_id=new_rep_gene.id
                new_rep_seq_dict[new_rep_id]=new_rep_gene
                new_all_seq_dict[new_rep_id]=seqList
                continue
            elif len(panYeast_s288c_geneList)==1:
                count_panYeast1+=1
                # if only 1 Yeast8 gene exist, choose it as representive sequence
                new_rep_gene=panYeast_s288c_geneList[0]
                new_rep_id=new_rep_gene.id
                new_rep_seq_dict[new_rep_id]=new_rep_gene
                new_all_seq_dict[new_rep_id]=seqList
                continue
            else:
                count_panYeast_more1+=1
                # if more than 1 Yeast8 gene exist, choose the longest one as representive sequence
                panYeast_s288c_geneList.sort(key=lambda x:len(str(x.seq)),reverse=True)
                new_rep_gene=panYeast_s288c_geneList[0]
                new_rep_id=new_rep_gene.id
                new_rep_seq_dict[new_rep_id]=new_rep_gene
                new_all_seq_dict[new_rep_id]=seqList
                panYeast_gene_collapse[new_rep_id]=panYeast_s288c_geneList
                continue
        else:
            # if no s288c sequence, the original representive sequence from mmseqs2 result is selected
            new_rep_id=ori_repID
            new_rep_gene=[g for g in seqList if g.id==new_rep_id][0]
            new_rep_seq_dict[new_rep_id]=new_rep_gene
            new_all_seq_dict[new_rep_id]=seqList
            continue


#according to GPR, genes in the same cluster but have different GPR will be add again to the pan-genome
model=read_sbml_model("model/yeastGEM.xml")
add_panYeast_geneIDList=list()
for rep_gene,geneList in panYeast_gene_collapse.items():
    # print(rep_gene+"----------------------------")
    rep_geneID=rep_gene.strip("s288c|")
    rep_gpr=model.genes.get_by_id(rep_geneID).reactions
    # print(str(rep_gpr))
    for g in geneList:
        id=g.id.strip("s288c|")
        gene=model.genes.get_by_id(id)
        gpr=gene.reactions
        if gpr==rep_gpr:
            continue
        else:
            print(rep_gene + "----------------------------")
            # print(str(rep_gpr))
            print(id+" : "+str(gene.reactions)+str(len(str(g.seq))))
            add_panYeast_geneIDList.append(id)
# remove YOL086C,YOL156W
add_panYeast_geneIDList.remove("YOL086C")
add_panYeast_geneIDList.remove("YOL156W")
add_geneList=[g for g in all_seq if (g.id.strip("s288c|") in add_panYeast_geneIDList)&(len(g.seq)>0)]
for g in add_geneList:
    new_rep_seq_dict[g.id]=g
    new_all_seq_dict[g.id]=[g]

new_all_seqID_dict=dict()
for id,seqList in new_all_seq_dict.items():
    rep_id=id.strip("s288c|")
    seqIDlist=[seq.id.strip("s288c|") for seq in seqList]
    new_all_seqID_dict[rep_id]=seqIDlist

# save new_all_seqID_dict
import pickle
with open("data/genome/pan1800_50_70_v2_cluster.pkl","wb") as f:
    pickle.dump(new_all_seqID_dict,f)


# write the pan-genome representive sequence
with open("code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2.fasta","w") as f:
    for id,seq in new_rep_seq_dict.items():
        id=id.lstrip('s288c|')
        f.write(">"+id+"\n"+str(seq.seq)+"\n")



