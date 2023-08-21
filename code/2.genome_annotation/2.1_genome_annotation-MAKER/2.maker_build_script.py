# -*- coding: utf-8 -*-
# braker_AAB_6_allcds.fa
# why
# file：maker_build_script.py
# 2022/5/18
'''build genome annotation pipleline by MAKER'''

import os

output_file='code/genome_annotation-MAKER/maker_failstrain_20220921.slurm'
allstrain_dir="data/genome/1900_assembled_genome"

head='''#!/bin/bash

#SBATCH --job-name=maker_fail
#SBATCH --partition=cpu
#SBATCH -n 80
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/maker
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker'''

with open(output_file,'w+') as f:
    f.write(head+'\n')


masked_genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/masked_genome"
proteomes_result_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker_proteomes"
gff_result_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker_gff"
# get all strains list
# strainlist=os.listdir(allstrain_dir)
# new600_strainlist=os.listdir("data/genome/new600_assembled_genome")

# get strain list which failed to maker process last time
maker_list=os.listdir("data/genome/predicted_allcds/maker_proteomes")
flo_list=os.listdir("data/genome/predicted_allcds/flo_proteomes")
maker_list=[i.replace(".fa","") for i in maker_list]
flo_list=[i.replace(".fasta","") for i in flo_list]
maker_fail_list=[]
for strain in flo_list:
    if strain not in maker_list:maker_fail_list.append(strain)

i=0
for strain in maker_fail_list:
    i+=1
    file=str(i%20+1)
    # strain_name=strain.strip('.fna')
    maker1_cmd ="( cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker/%s" % file
    maker2_cmd="sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/%s\.fna\.masked/g' maker_opts.ctl"%strain
    maker3_cmd="maker -cpus 4 -fix_nucleotides "
    maker4_cmd="fasta_merge -d %s.fna.maker.output/%s.fna_master_datastore_index.log"%(strain,strain)
    maker5_cmd="cp %s.fna.all.maker.proteins.fasta %s/%s.fa"%(strain,proteomes_result_dir,strain)
    maker6_cmd="gff3_merge -d %s.fna.maker.output/%s.fna_master_datastore_index.log"%(strain,strain)
    maker7_cmd="cp %s.fna.all.gff %s/%s.gff"%(strain,gff_result_dir,strain)
    maker8_cmd="rm -r %s.fna.maker.output"%strain
    maker9_cmd="cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker )&"
    with open(output_file,'a+') as f:
        f.write(maker1_cmd + '\n')
        f.write(maker2_cmd + '\n')
        f.write(maker3_cmd + '\n')
        f.write(maker4_cmd + '\n')
        f.write(maker5_cmd + '\n')
        f.write(maker6_cmd + '\n')
        f.write(maker7_cmd + '\n')
        f.write(maker8_cmd + '\n')
        f.write(maker9_cmd + '\n')
        f.write('\n')
        if i%20==0:
            f.write('wait'+'\n')

