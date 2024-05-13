# -*- coding: utf-8 -*-
# date : 2023/5/6 
# author : wangh
# file : rename_proteomes_geneID.py
# project : Unified_Yeast_GEMs_Database_from_13pro
from Bio.SeqIO import parse
import os
import tqdm

proteomes_dir="data/genome/predicted_allcds/combined_proteomes/"
rename_dir="Unified_Yeast_GEMs_Database/data/genome/predicted_proteomes/combine_proteomes/"
strainList=os.listdir(proteomes_dir)
for strain in tqdm.tqdm(strainList):
    proteome=[g for g in parse(proteomes_dir+strain,"fasta")]
    with open(rename_dir+strain,'w+') as f:
        count=1
        for gene in proteome:
            id=str(count)+"_"+gene.id
            f.write(">"+id+"\n"+str(gene.seq)+"\n")
            count+=1