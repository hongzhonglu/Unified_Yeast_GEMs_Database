# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：getBUSCO_scores.py
# 2022/6/7
'''get BUSCO_scores of all genomes
ref: lg_code :3rd_genomeCompleteness/2nd_getBUSCO_Scores.py'''

import os
import re

head='''#Complete BUSCOs (C)
#Complete and single-copy BUSCOs (S)
#Complete and duplicated BUSCOs (D)
#Fragmented BUSCOs (F)
#Missing BUSCOs (M)
#Total BUSCO groups searched (n)
'''

input_dir = r'/lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_genome_completeness/1900_busco_results/'
os.chdir(input_dir)

busco_results=open('../1900_busco_scores.txt','w')
busco_results.write('#GenomeName,C(%),S(%),D(%),F(%),M(%),n\n')



for name in os.listdir(input_dir):
    if not name.endswith('.txt'):continue
    fhand=open(input_dir+name,'r')
    for line in fhand.readlines():
        line=line.strip()
        if 'C:' not in line:continue
        C=re.findall('C:(.*?)%',line)[0]
        S=re.findall('S:(.*?)%',line)[0]
        D=re.findall('D:(.*?)%',line)[0]
        F=re.findall('F:(.*?)%',line)[0]
        M=re.findall('M:(.*?)%',line)[0]
        n=re.findall('n:(.*)',line)[0]

        GenomeName=name.replace('short_summary.specific.saccharomycetes_odb10.','').replace('.fa.out.txt','')
        busco_results.write('{},{},{},{},{},{},{}\n'.format(GenomeName,C,S,D,F,M,n))
busco_results.close()