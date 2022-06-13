# -*- coding: utf-8 -*-
# braker_AAB_6_allcds.fa
# why
# file：maker_build_script.py
# 2022/5/18
'''build genome annotation pipleline by MAKER'''

import os

output_file='code/maker_script/maker_new600.slurm'

head='''#!/bin/bash

#SBATCH --job-name=maker_new600
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


rm_output_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/repeatmasker"     #repeatmasker output
genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome"
masked_genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/masked_genome"
allcds_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker"
# get all strains list
strainlist=os.listdir('data/genome/new_strains_assembled_genome')
i=0
for strain in strainlist[50:]:
    i+=1
    file=str(i%20+1)
    strain_name=strain.strip('.fna')
    maker1_cmd ="(cd %s" % file
    maker2_cmd="sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/%s\.fna\.masked/g' maker_opts.ctl"%strain_name
    maker3_cmd="maker -cpus 4 -fix_nucleotides "
    maker4_cmd="fasta_merge -d %s.fna.maker.output/%s.fna_master_datastore_index.log"%(strain_name,strain_name)
    maker5_cmd="cp %s.fna.all.maker.proteins.fasta %s/%s.fa"%(strain_name,allcds_dir,strain_name)
    maker6_cmd="rm %s.fna.all*"%strain_name
    maker7_cmd="cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &"
    with open(output_file,'a+') as f:
        f.write(maker1_cmd + '\n')
        f.write(maker2_cmd + '\n')
        f.write(maker3_cmd + '\n')
        f.write(maker4_cmd + '\n')
        f.write(maker5_cmd + '\n')
        f.write(maker6_cmd + '\n')
        f.write(maker7_cmd + '\n')
        f.write('\n')
        if i%20==0:
            f.write('wait'+'\n')

