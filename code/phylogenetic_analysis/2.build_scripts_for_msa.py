# -*- coding: utf-8 -*-
# date : 2024/10/12 
# author : wangh
# file : 2.build_scripts_for_msa.py
# project : Unified_Yeast_GEMs_Database
'''Build scripts for core ORFs MSA by maffat and trimming by trimAl
** commands:
# MSA
mafft --thread 4 --auto --maxiterate 1000 input.fa > output_msa.fa
# trimming
trimal -in test -out aligned_trimed.fa -gappyout

'''
import os

script_path=r'code/phylogenetic_analysis/2.core_ORFs_msa.slurm'
parallel_num=20

wkdir=r'/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database'
input_dir=r'code/phylogenetic_analysis/output/core_ORFs'
mas_output_dir=r'code/phylogenetic_analysis/output/core_ORFs_msa'
trim_output_dir=r'code/phylogenetic_analysis/output/core_ORFs_msa_trimmed'

msa_command_template='mafft --thread 3 --auto --maxiterate 1000 INPUT > OUTPUT'
trimming_command_template='trimal -in INPUT -out OUTPUT -automated1'


head=f'''#!/bin/bash

#SBATCH --job-name=coreORFs_msa
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /dssg/home/acct-clslhz/clslhz/miniconda3/envs/hyena-d/iqtree

cd {wkdir}
# mkdir
if [ ! -d {mas_output_dir} ]; then
    mkdir {mas_output_dir}
fi
if [ ! -d {trim_output_dir} ]; then
    mkdir {trim_output_dir}
fi
'''

with open(script_path,'w+') as f:
    f.write(head)

with open(script_path,'a') as f:
    i=0
    for file in os.listdir(input_dir):
        if file.endswith('.fa'):
            i+=1
            input_file=input_dir+'/'+file
            output_msa_file=mas_output_dir+'/'+file
            output_trim_file=trim_output_dir+'/'+file
            msa_command='('+msa_command_template.replace('INPUT',input_file).replace('OUTPUT',output_msa_file)
            trimming_command=trimming_command_template.replace('INPUT',output_msa_file).replace('OUTPUT',output_trim_file)+')&'
            f.write(msa_command+'\n')
            f.write(trimming_command+'\n')
            f.write('\n')

            if i%parallel_num==0:
                f.write('wait\n')
                f.write('\n')
