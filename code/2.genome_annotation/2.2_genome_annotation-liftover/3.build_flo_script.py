# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：build_flo_script.py
# 2022/6/27
# build scripts to annotate assembled genome sequence by liftover

import os
output_script="flo_1900.slurm"

head='''#!/bin/bash

#SBATCH --job-name=flo_new600
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/flo
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/flo_sce

'''
word_dir="code/genome_annotation-liftover"
strainList="../../data/genome/1900_assembled_genome/"
input_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/"
cmds = open('code/genome_annotation-liftover/liftover_template').read()
flo_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo"

os.chdir(word_dir)

existed_list_0=os.listdir("../../data/genome/predicted_allcds/combined_proteomes")
existed_list=[i.replace(".fa",".fna") for i in existed_list_0]

fhand = open(output_script, 'a')
fhand.write(head)
for name in os.listdir(strainList):
    if not name.endswith('fna'): continue
    if name in existed_list:continue
    new_cmds = cmds.replace('INPUT',os.path.join(input_dir,name))
    new_cmds = new_cmds.replace('OUT',name.strip(".fna"))
    new_cmds=new_cmds.replace("FLO_PATH",flo_dir)
    fhand.write(new_cmds)
    fhand.write('\n\n')
fhand.close()

# rename assembled sequence
# filelist=os.listdir("data/genome/1900_assembled_genome/")
# for file in filelist:
#     file_1=file.replace(".fa",".fna")
#     if file_1!=file:
#         os.rename("data/genome/1900_assembled_genome/"+file,"data/genome/1900_assembled_genome/"+file_1)


