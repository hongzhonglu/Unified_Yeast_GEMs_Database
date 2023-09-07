# -*- coding: utf-8 -*-
# date : 2023/5/15 
# author : wangh
# file : extract_busco_score.py
# project : Unified_Yeast_GEMs_Database
import os
import re
import os
os.chdir(r"D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database")

head='''#Complete BUSCOs (C)
#Complete and single-copy BUSCOs (S)
#Complete and duplicated BUSCOs (D)
#Fragmented BUSCOs (F)
#Missing BUSCOs (M)
#Total BUSCO groups searched (n)
'''

input_dir = r'code/3.pan-genome_construction/3.pan-genome_comparison/output/BUSCO_summary/'

busco_results=open(r'code/3.pan-genome_construction/3.pan-genome_comparison/output/pangenomes_busco_scores.txt','w+')
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

        GenomeName=name.replace(".txt","")
        busco_results.write('{},{},{},{},{},{},{}\n'.format(GenomeName,C,S,D,F,M,n))
busco_results.close()