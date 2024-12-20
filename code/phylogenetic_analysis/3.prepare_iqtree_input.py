# -*- coding: utf-8 -*-
# date : 2024/10/13 
# author : wangh
# file : 3.prepare_iqtree_input.py
# project : Unified_Yeast_GEMs_Database
'''Prepare ORFs used for building concatenated phylogenetic tree.
rules:
1. >50% strain occupancy
2. alignment length > 200 aa
3. locolize in cytosol
'''
import os
import pandas as pd
import numpy as np
from Bio import SeqIO

# set workdir
# os.chdir(r'code/phylogenetic_analysis/')
pan_geneList=[i for i in SeqIO.parse(r'data/genome/pan1800.fasta','fasta')]
occupancy_threshold=0.5
align_length_threshold=60

# 1. strain occupancy filter
geneMatrix=pd.read_csv(r'data/geneMatrix/pan1800_v2_blastp_50_70_geneMatrix.csv',index_col=0)
df_count=geneMatrix.sum(axis=1)/geneMatrix.shape[1]
to_remove=df_count[df_count<occupancy_threshold].index.tolist()
keep_geneList=[i for i in pan_geneList if i.id not in to_remove]
print(f'{len(keep_geneList)} genes left after {occupancy_threshold} occupancy filter')

# 2. alignment length filter
from Bio import AlignIO
align_dir=r'code/phylogenetic_analysis/output/core_ORFs_msa_trimmed/'
fileList=os.listdir(align_dir)
length_dict={}
for file in fileList:
    name=file.replace('.fa','')
    aln=AlignIO.read(align_dir+file,'fasta')
    aln_length=aln.get_alignment_length()
    length_dict[name]=aln_length

to_remove2=[i for i in length_dict.keys() if length_dict[i]<align_length_threshold]
keep_geneList=[i for i in keep_geneList if i.id not in to_remove2]
print(f'{len(keep_geneList)} genes left after {align_length_threshold} length filter')

# 3. cytosol filter
uniprot_annot=pd.read_excel(r'data/s288c_localization.xlsx',index_col=0)
# remove Gene Ontology (cellular component) is nan
# uniprot_annot=uniprot_annot[~uniprot_annot['Gene Ontology (cellular component)'].isna()]
uniprot_annot=uniprot_annot[~uniprot_annot['Subcellular location [CC]'].isna()]
# set Gene Names (ordered locus) as index
uniprot_annot=uniprot_annot.set_index('Gene Names (ordered locus)')
# keep gene contain Cytoplasm in Subcellular location [CC]
cyctosal_geneIDlist=uniprot_annot[uniprot_annot['Subcellular location [CC]'].str.contains('Cytoplasm')].index.tolist()
keep_geneList=[i for i in keep_geneList if i.id in cyctosal_geneIDlist]

# predict localization by deeploc2(https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) for unknown genes
# to_predict_geneList=[i for i in keep_geneList if i.id not in uniprot_annot.index]
# # write to fasta
# with open(r'code/phylogenetic_analysis/output/to_predict_localization.fasta','w') as f:
#     for i in to_predict_geneList:
#         f.write('>'+i.id+'\n')
#         f.write(str(i.seq)+'\n')

# parse deeploc2 output
deeploc_threshold=0.7
deeploc2_result=pd.read_csv(r'code/phylogenetic_analysis/output/predict_localization_by_deeploc2.csv',index_col=0)
predict_cyctosal_geneIDlist=deeploc2_result[(deeploc2_result['Localizations']=='Cytoplasm')&(deeploc2_result['Cytoplasm']>deeploc_threshold)].index.tolist()

cyctosal_geneIDlist=cyctosal_geneIDlist+predict_cyctosal_geneIDlist

keep_geneList=[i for i in keep_geneList if i.id in cyctosal_geneIDlist]
keep_geneIDlist=[i.id for i in keep_geneList]

print(f'{len(keep_geneList)} genes left after cytosol filter')
# extract all msa_trimmed ORFs into iqtree input directory
import shutil
iqtree_input_dir=r'code/phylogenetic_analysis/iqtree_input/'
all_alignments_dir=r'code/phylogenetic_analysis/output/core_ORFs_msa_trimmed/'

# if input dir not exist, create it
if not os.path.exists(iqtree_input_dir):
    os.mkdir(iqtree_input_dir)

for id in keep_geneIDlist:
    fileName=id+'.fa'
    if not os.path.exists(all_alignments_dir+fileName):
        print(f'Can not find {fileName}')
        continue
    # copy to iqtree input dir
    if not os.path.exists(iqtree_input_dir+fileName):
        shutil.copy(all_alignments_dir+fileName,iqtree_input_dir)
