# -*- coding: utf-8 -*-
# braker_AAB_6_allcds.fa
# why
# file：build_seqtranslate_script.py
# 2022/5/23

import os
output_file="code/translate.sh"
genomes_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_1392_allcds"
proteomes_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes"

strainlist=os.listdir('D:\code\Github\prac_auto_ssGEMs\data\strain genome sequences\pan_genome_gang\Combined_CDS(genomes)')
for strain in strainlist:
    trans_cmd="seqkit translate %s/%s -o %s/%s"%(genomes_dir,strain,proteomes_dir,strain)
    with open(output_file,"a+") as f:
        f.write(trans_cmd+"\n")
