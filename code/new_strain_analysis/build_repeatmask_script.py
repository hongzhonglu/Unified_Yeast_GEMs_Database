# -*- coding: utf-8 -*-
# braker_AAB_6_allcds.fa
# why
# file：build_repeatmask_script.py
# 2022/5/18
import os

output_file='code/repeatmask_scripts/repeatmask_other.sh'
rm_output_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/repeatmasker"     #repeatmasker output
genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome"
masked_genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/masked_genome"
cd_cmd="cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker"
with open(output_file,'w+') as f:
    f.write(cd_cmd+'\n')
i=0
strainlist=os.listdir('data/genome/assembled_genome')
for strain in strainlist[100:]:
    i+=1
    strain_name=strain.strip('.fna')
    rm1_cmd="(RepeatMasker -species \"saccharomyces cerevisiae\" -xsmall -e ncbi -dir %s/%s %s/%s.fna"\
           %(rm_output_dir,strain_name,genome_dir,strain_name)
    rm2_cmd='cp %s/%s/%s.fna.masked %s/ ) &'%(rm_output_dir,strain_name,strain_name,masked_genome_dir)
    with open(output_file,'a') as f:
        f.write(rm1_cmd+'\n')
        f.write(rm2_cmd+'\n')
        f.write('\n')
        if i%20==0:
            f.write('wait'+'\n')