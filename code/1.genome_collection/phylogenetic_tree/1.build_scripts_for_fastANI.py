# -*- coding: utf-8 -*-
# date : 2023/6/11 
# author : wangh
# file : build_scripts_for_fastANI.py
# project : Unified_Yeast_GEMs_Database
'''用法：fastANI --ql <QUERY_LIST> --rl <REFERENCE_LIST> -o <OUTPUT_FILE> --matrix --threads 3
'''
import os

script_file="code/1.genome_collection/4_2.fastANI_1900sce.slurm"
query_list_path="/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/1.genome_collection/outputs/1900sce_assebled_genome_paths.txt"
ref_list_dir="/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/1.genome_collection/outputs/split_1900sce_genome_path/"
output_dir="/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/1.genome_collection/outputs/fastANI_split_result/"

head='''#!/bin/bash

#SBATCH --job-name=fastani_1900
#SBATCH --partition=64c512g
#SBATCH -N 2
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /dssg/home/acct-clslhz/clslhz/miniconda3/envs/fastani
'''

with open(script_file,'w') as f:
    f.write(head+"\n\n\n")

k=0
ref_list=os.listdir("code/1.genome_collection/outputs/split_1900sce_genome_path")
for ref in ref_list:
    ref_list_path=ref_list_dir+ref
    output=output_dir+"1900sce_fastANI_split%s"%k
    fastani_cmd="fastANI --ql %s --rl %s -o %s --matrix --threads 6 &"%(query_list_path,ref_list_path,output)
    k+=1
    with open(script_file,'a+') as f:
        f.write(fastani_cmd+"\n")
        if k%20==0:
            f.write("wait\n\n")
