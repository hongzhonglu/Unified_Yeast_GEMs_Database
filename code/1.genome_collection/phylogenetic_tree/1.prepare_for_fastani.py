# -*- coding: utf-8 -*-
# date : 2023/6/10 
# author : wangh
# file : 4_1.prepare_for_fastani.py
# project : Unified_Yeast_GEMs_Database
# prepare assembled genome paths list
import os

assemble_genome_dir=r"/dssg/home/acct-clslhz/clslhz/why_ssGEM/data/assembled_genome/"
strainList=os.listdir(r"E:\data\sce_ssGEMs\genomes\assembled_sequence\1900_assembled_genome")

with open('code/1.genome_collection/outputs/1900sce_assebled_genome_paths.txt','w') as f:
    for strain in strainList:
        path=assemble_genome_dir+strain
        f.write(path+'\n')