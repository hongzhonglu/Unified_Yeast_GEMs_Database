# -*- coding: utf-8 -*-
# date : 2023/5/7 
# author : wangh
# file : 3.combine_s288c&nonref.py
# project : Unified_Yeast_GEMs_Database
# combine two fasta file into one fasta file

s288c_path="data/genome/S288c_R64.fasta"
nonref_path=r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900_s288c_nonref_v4_filtered.fasta"
output_path=r"code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/1900sce_s288c+nonref_v4_filtered.fasta"
with open(output_path,"w") as f:
    with open(s288c_path,"r") as f1:
        for line in f1:
            if line.startswith(">"):
                line=line.replace(">",">s288c|")
            f.write(line)
    with open(nonref_path,"r") as f2:
        for line in f2:
            f.write(line)

