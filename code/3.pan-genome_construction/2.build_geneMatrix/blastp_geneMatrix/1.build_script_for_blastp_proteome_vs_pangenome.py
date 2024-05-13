# -*- coding: utf-8 -*-
# date : 2023/5/8 
# author : wangh
# file : 5.build_script_for_proteome_vs_pangenome.py
# project : Unified_Yeast_GEMs_Database
import os

# tblastn script parameters:
ref_genome="pan1800.fasta"
ref_name=ref_genome.rstrip(".fasta")
output_script="code/3.pan-genome_construction/2.build_geneMatrix/blastp_geneMatrix/1.mmseqs2_blastp_vs_%s.sh"%ref_name
target_dir = "/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/"
blastp_output_dir="/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/2.build_geneMatrix/output/blastp_vs_%s/"%ref_name
proteomes_dir_wsl="/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/predicted_proteomes/combined_proteomes_old/"
blastp_template_cmd='mmseqs easy-search QUERY TARGET OUTPUT tmp --min-seq-id 0.6 -e 1e-5'

# building code parameters:
genomes_dir_win=r"data/genome/predicted_proteomes/combined_proteomes"
genomeList=os.listdir(genomes_dir_win)

# build scripts for 100 genomes to do tblastn process by mmeseqs2
# write cmd to set the workdir
with open(output_script,'w') as f:
    f.write('mkdir '+blastp_output_dir+"\n")
    f.write('cd '+blastp_output_dir+"\n")


# write cmd to do tblastn process for each genomes
for genome in genomeList:
    query=proteomes_dir_wsl+genome
    subject=target_dir+ref_genome
    output_file="%s_vs_%s_blastp.txt"%(genome.rstrip(".fa"),ref_name)
    blastp_output=blastp_output_dir+output_file
    blastp_cmd=blastp_template_cmd.replace("QUERY",query).replace("TARGET",subject).replace("OUTPUT",blastp_output)
    with open(output_script,'a') as f:
        f.write(blastp_cmd+"\n")
        f.write("\n")


