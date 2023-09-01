# -*- coding: utf-8 -*-
# date : 2023/5/10 
# author : wangh
# file : 5.build_script_for_blastp_proteome_vs_pangenome.py
# project : Unified_Yeast_GEMs_Database
import os

# tblastn script parameters:
ref_genome="s288c_7yeastsnonref_cov50_pid70_rep_seq.fasta"
ref_name=ref_genome.rstrip(".fasta")
output_script="code/4.pan-genome_analysis/compare_with_closed_yeasts/new_pipeline_pan-genome/5.mmseqs2_blastp_vs_%s.sh"%ref_name
target_dir = "/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/4.pan-genome_analysis/compare_with_closed_yeasts/output/"
blastp_output_dir="/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/4.pan-genome_analysis/compare_with_closed_yeasts/output/blastp_vs_%s/"%ref_name
proteomes_dir_wsl="/mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/343_yeast_genomes/"
blastp_template_cmd='mmseqs easy-search QUERY TARGET OUTPUT tmp --min-seq-id 0.6 -e 1e-5'

# building code parameters:
genomes_dir_win=r"data/genome/343_yeast_genomes/"
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
    output_file="%s_vs_%s_blastp.txt"%(genome.rstrip(".fasta"),ref_name)
    blastp_output=blastp_output_dir+output_file
    blastp_cmd=blastp_template_cmd.replace("QUERY",query).replace("TARGET",subject).replace("OUTPUT",blastp_output)
    with open(output_script,'a') as f:
        f.write(blastp_cmd+"\n")
        f.write("\n")

