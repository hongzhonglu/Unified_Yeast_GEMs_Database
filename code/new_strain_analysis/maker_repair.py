# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：maker_repair.py
# 2022/5/23
import os


output_file="code/new_strain_maker_repaire.sh"
rm_output_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/repeatmasker"     #repeatmasker output
genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome"
masked_genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/masked_genome"
allcds_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker"

with open(output_file,"w+") as f:
    f.write("cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker\n")
# get all strains list
strainlist=os.listdir('data/genome/new_strains_assembled_genome')
for i in range(1,11):
    for strain in strainlist:
        strain_name=strain.strip('.fna')
        maker_rp_cmd1="fasta_merge -d %s/%s.fna.maker.output/%s.fna_master_datastore_index.log"%(i,strain_name,strain_name)
        maker_rp_cmd2="cp %s.fna.all.maker.proteins.fasta %s/%s.fa"%(strain_name,allcds_dir,strain_name)
        maker_rp_cmd3="rm %s.fna.all* "%strain_name
        with open(output_file,'a+') as f:
            f.write(maker_rp_cmd1 + '\n')
            f.write(maker_rp_cmd2 + '\n')
            f.write(maker_rp_cmd3 + '\n')
            f.write('\n')