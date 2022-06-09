# -*- coding: utf-8 -*-
# why
# file：build_bbh_script.py
# 2022/5/21
'''build bash script to build blustDB and blastp by diamond for lg_1392 annotated genomes'''
import os

output_file="code/bbh_new.slurm"
allcds_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes"
blastdb_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db"
bbh_result_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result"
ref_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs"

head='''#!/bin/bash

#SBATCH --job-name=lg_bbh
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/maker
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh'''


ref_id="pan1011_cds"
# build blastdb for ref_pan1011
ref_blastdb_cmd="diamond makedb --in %s/%s.fa -d %s/%s.dmnd"%(ref_dir,ref_id,blastdb_dir,ref_id)

with open(output_file,"w") as f:
    f.write(head+"\n\n\n")
    f.write("# build blastdb for ref_pan1011\n")
    f.write(ref_blastdb_cmd+"\n\n")

strainlist=os.listdir('D:\code\Github\prac_auto_ssGEMs\data\strain genome sequences\pan_genome_gang\Combined_CDS(genomes)')
i=0
for strain in strainlist:
    strain_name=strain.replace(".fa","")
    i+=1
    blastdb_cmd="( diamond makedb --in %s/%s -d %s/%s.dmnd"%(allcds_dir,strain,blastdb_dir,strain)
    cmd1_out="%s_vs_%s.txt"%(strain_name,ref_id)
    cmd2_out="%s_vs_%s.txt"%(ref_id,strain_name)
    bbh_cmd1="diamond blastp -d %s/%s.dmnd -q %s/%s -o %s/%s --more-sensitive --evalue 0.001 " \
             "--max-target-seqs 20"%(blastdb_dir,ref_id,allcds_dir,strain,bbh_result_dir,cmd1_out)
    bbh_cmd2="diamond blastp -d %s/%s.dmnd -q %s/%s.fa -o %s/%s --more-sensitive --evalue 0.001 " \
             "--max-target-seqs 20 ) &"%(blastdb_dir,strain,ref_dir,ref_id,bbh_result_dir,cmd2_out)
    with open(output_file,"a+") as f:
        f.write(blastdb_cmd+"\n")
        f.write(bbh_cmd1+"\n")
        f.write(bbh_cmd2+"\n\n")
        if i%10==0:
            f.write("wait"+"\n")

