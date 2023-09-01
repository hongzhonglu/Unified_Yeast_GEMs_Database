# -*- coding: utf-8 -*-
# date : 2023/5/6 
# author : wangh
# file : 3.combine_nonref_s288c.py
# project : Unified_Yeast_GEMs_Database
# combine two fasta file into one fasta file

s288c_path="data/genome/S288c_R64.fasta"
nonref_path=r"code/4.pan-genome_analysis/compare_with_closed_yeasts/output/7yeasts_s288c_nonref.fasta"
output_path=r"code/4.pan-genome_analysis/compare_with_closed_yeasts/output/s288c+7yeasts_nonref.fasta"
with open(output_path,"w") as f:
    with open(s288c_path,"r") as f1:
        for line in f1:
            f.write(line)
    with open(nonref_path,"r") as f2:
        for line in f2:
            f.write(line)
