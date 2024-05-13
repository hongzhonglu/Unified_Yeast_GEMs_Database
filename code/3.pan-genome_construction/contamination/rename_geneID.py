# -*- coding: utf-8 -*-
# date : 2023/5/17 
# author : wangh
# file : id_to_tax.py
# project : Unified_Yeast_GEMs_Database
# rename geneID to shorten it
import pandas as pd
import os
from Bio import SeqIO

# rename pan1900_removed1600
records=[record for record in SeqIO.parse('data/genome/pan1900_new_50_70_removed1600.fasta','fasta')]
removeIDlist=[record.id for record in records]
count=0
geneID_dict={}
with open("data/genome/pan1900_new_50_70_removed1600_rename.fasta",'w+') as f:
    for record in records:
        count+=1
        reID='remove%s'%count
        seq=str(record.seq)
        geneID_dict[record.id]=reID
        f.write('>'+reID+'\n'+seq+'\n')

# save the geneID_dict
df=pd.DataFrame.from_dict(geneID_dict,orient='index')
df.columns=['renameID']
df['geneID']=df.index
df.to_csv('code/3.pan-genome_construction/contamination/output/pan1900_new_remove1900_index.csv')


# rename pan1900_new_50_70_all9323
records=[record for record in SeqIO.parse('data/genome/pan1900_new_50_70_all9323.fasta','fasta')]
count=0
geneID_dict={}
with open("data/genome/pan1900_new_50_70_all9323_rename.fasta",'w+') as f:
    for record in records:
        count+=1
        if record.id in removeIDlist:
            reID='remove%s'%count
        else:
            reID='keep%s'%count
        seq=str(record.seq)
        geneID_dict[record.id]=reID
        f.write('>'+reID+'\n'+seq+'\n')

# save the geneID_dict
df=pd.DataFrame.from_dict(geneID_dict,orient='index')
df.columns=['renameID']
df['geneID']=df.index
df.to_csv('code/3.pan-genome_construction/contamination/output/pan1900_new_all9323_index.csv')

