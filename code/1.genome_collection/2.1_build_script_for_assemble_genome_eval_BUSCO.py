# -*- coding: utf-8 -*-
# date : 2022/11/23 
# author : wangh
# file : build_script_for_assemble_genome_eval.py
# project : Unified_Yeast_GEMs_Database_from_13pro
'''evaluate assemble genome quality:
用法：busco -m genome -i INPUT.nucleotides -o OUTPUT -l LINEAGE'''
import os

output_file="code_standardize/1.genome_sequence_collection/1900_assemblies_busco_completeness_eval.slurm"
predicted_allcds_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome"
busco_output_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_genome_completeness"

head='''#!/bin/bash

#SBATCH --job-name=1900_busco_check_genome_completeness
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# module purge
# module load miniconda3
# source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/busco
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_assembly_completeness'''

with open(output_file,"w") as f:
    f.write(head+"\n\n\n")

k=0
strainlist=os.listdir(r'data/genome/1900_assembled_genome')
for strain in strainlist:
    k+=1
    # strain=strain.replace(".fna",".fa")
    busco_cmd="(busco -i %s/%s -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o %s.out -m genome -e 1e-3"%(predicted_allcds_dir,strain,strain)
    mv_cmd="mv %s.out/short_summary.specific.saccharomycetes_odb10.%s.out.txt 1900_assemblies_BUSCO_results/"%(strain,strain)
    clean_cmd="rm -r %s.out)&"%strain
    with open(output_file,"a+") as f:
        f.write(busco_cmd+"\n")
        f.write(mv_cmd+"\n")
        f.write(clean_cmd+"\n")
        if k%10==0:
            f.write("wait\n\n")



