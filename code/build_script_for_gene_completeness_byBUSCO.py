# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：build_script_for_gene_completeness_byBUSCO.py
# 2022/6/5
'''evaluate gene completeness by BUSCO'''
import os

output_file="code/new_strain_analysis/new600_busco_part2.slurm"
assembled_genomes_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/new600_allproteomes_maker"
busco_output_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_genome_completeness"

head='''#!/bin/bash

#SBATCH --job-name=busco_check_genome_completeness
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
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_genome_completeness'''

with open(output_file,"w") as f:
    f.write(head+"\n\n\n")

k=0
strainlist=os.listdir('data/genome/new_strains_assembled_genome')
for strain in strainlist[30:]:
    k+=1
    strain=strain.replace(".fna",".fa")
    busco_cmd="(busco -i %s/%s -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o %s.out -m prot -e 1e-3"%(assembled_genomes_dir,strain,strain)
    mv_cmd="mv %s.out/short_summary.specific.saccharomycetes_odb10.%s.out.txt new600_busco_results/"%(strain,strain)
    clean_cmd="rm -r %s.out)&"%strain
    with open(output_file,"a+") as f:
        f.write(busco_cmd+"\n")
        f.write(mv_cmd+"\n")
        f.write(clean_cmd+"\n")
        if k%20==0:
            f.write("wait\n\n")