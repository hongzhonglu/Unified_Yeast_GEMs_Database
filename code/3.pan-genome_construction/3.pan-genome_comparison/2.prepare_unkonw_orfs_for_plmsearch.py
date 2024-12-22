# -*- coding: utf-8 -*-
# date : 2024/10/13 
# author : wangh
# file : 2.prepare_unkonw_orfs_for_plmsearch.py
# project : Unified_Yeast_GEMs_Database
'''Prepare representative seqs of functionally unkonw ORFs for functional prediction by PLMsearch.
'''
import pandas as pd
from Bio import SeqIO

pan_file=r'data/genome/pan1800.fasta'
annot_file=r'data/genome/pan1800_functional_annotations.xlsx'

# new ORFs
new_orfList=[i for i in SeqIO.parse(pan_file,'fasta') if 'scepan' in i.id]

# write in fasta format
with open(r'code/3.pan-genome_construction/3.pan-genome_comparison/output/new_sceorfs.fa','w') as f:
    for i in new_orfList:
        f.write('>'+i.id+'\n')
        f.write(str(i.seq)+'\n')