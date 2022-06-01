# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：build_geneMatrix.py
# 2022/5/10


import os  ##for directory
import sys
import pandas as pd
import numpy as np
from glob import glob
from Bio import SeqIO,SeqRecord


sys.path.append('code')
from mainFunction import *



# get all strains' all_cds file
strain_list=os.listdir('test/build_geneMatrix/')
strains_dir='test/build_geneMatrix/'
# 将基因组全部mRNA sequence翻译为protein sequence
for strain in strain_list:
    trans_genome(strain)


#以1011_pan-genome作为reference,同样将其翻译成蛋白序列
records=SeqIO.parse('nature_1011_pangenome.fasta','fasta')
with open('extra_strain_genome_trans/ref_pan1011.fa','w') as output:
    for record in records:
        seq = record.seq.translate()
        seqid0 = record.id
        seqid=seqid0.split('-',1)
        output.write('>%s\n%s\n' % (seqid[1], seq))


# 分别对new strain以及reference genome建库
ref_id='ref_pan1011'
make_blast_db(ref_id, folder='extra_strain_genome_trans', db_type='prot')
for strain in extra_strainlist:
    make_blast_db(strain, folder='extra_strain_genome_trans', db_type='prot')


# execute bi-directional blast for each target strain against ref_strain,and save result to 'bi_blastp_result'folder
for strain in extra_strainlist:
    get_bbh(strain,ref_id,in_folder='bi_blastp_result')


# Parse the BLAST Results into one Homology Matrix of the Reconstruction Genes
# test blast result
blast_files=glob('bi_blastp_result/*_parsed.csv')
for blast in blast_files:
    bbh=pd.read_csv(blast)
    print(blast,bbh.shape)

# 以1011geneMatrix中的genelist作为reference
geneMatrix0=pd.read_csv('../geneMatrix0 of 1011 yeast strains.txt',sep='\t')
listGeneIDs=list(geneMatrix0['geneID'])
# 以panGEM中的基因作为genelist
# listGeneIDs=[]
# import cobra
# model = cobra.io.read_sbml_model('../panYeast_v3.xml')
# # ref_geneome=SeqIO.parse('1011_pan_genome.fasta.txt','fasta')
# for gene in model.genes:
#     listGeneIDs.append(gene.id)

ortho_matrix=pd.DataFrame(index=listGeneIDs,columns=extra_strainlist)
geneIDs_matrix=pd.DataFrame(index=listGeneIDs,columns=extra_strainlist)
for blast in blast_files:
    bbh=pd.read_csv(blast)
    listIDs = []
    listPID = []
    for r in ortho_matrix.iterrows():
        print(r[0])
        try:
            currentOrtholog = bbh[bbh['subject'] == r[0]].reset_index()
            listIDs.append(currentOrtholog.iloc[0]['subject'])
            listPID.append(currentOrtholog.iloc[0]['PID'])
        except:
            listIDs.append('None')
            listPID.append(0)
    for col in ortho_matrix.columns:
        if col in blast:
            ortho_matrix[col] = listPID
            geneIDs_matrix[col] = listIDs

for column in ortho_matrix:
    ortho_matrix.loc[ortho_matrix[column]<=80.0,column]=0
    ortho_matrix.loc[ortho_matrix[column]>80.0,column]=1

ortho_matrix.to_csv('extra_strain_geneMatrix_v0.csv')
geneIDs_matrix.to_csv('extra_strains_geneMatrix_v0.csv')