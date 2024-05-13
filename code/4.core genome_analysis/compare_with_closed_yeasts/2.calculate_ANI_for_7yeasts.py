# -*- coding: utf-8 -*-
# date : 2023/6/12 
# author : wangh
# file : 2.calculate_ANI_for_7yeasts.py
# project : Unified_Yeast_GEMs_Database
import os

strainList=[i.replace(".fasta",".fna") for i in os.listdir("data/genome/343_yeast_genomes/") if i.endswith(".fasta")][1:]
strainList.append("saccharomyces_cerevisiae.fna")

#data\343_yeast\0_332yeast_genomes\332_genome_assemblies\332_genome_assemblies
wsl_assembled_dir=r"/mnt/e/data/343_yeast/0_332yeast_genomes/332_genome_assemblies/332_genome_assemblies/"

with open("code/4.pan-genome_analysis/compare_with_closed_yeasts/7yeasts_assemblies_path.txt","w") as f:
    for i in strainList:
        f.write(wsl_assembled_dir+i+"\n")


# run fastani in WSL
# fastANI --ql code/4.pan-genome_analysis/compare_with_closed_yeasts/7yeasts_assemblies_path.txt --rl code/4.pan-genome_analysis/compare_with_closed_yeasts/7yeasts_assemblies_path.txt -o code/4.pan-genome_analysis/compare_with_closed_yeasts/output/7yeasts_ANI --matrix --threads 4


