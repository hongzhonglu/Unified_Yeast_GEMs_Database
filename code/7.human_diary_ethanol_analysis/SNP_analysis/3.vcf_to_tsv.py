# -*- coding: utf-8 -*-
# date : 2024/11/16 
# author : wangh
# file : 3.vcf_to_tsv.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import vcfpy


# extract all SNPs and indels respectively and output as tsv
reader = vcfpy.Reader.from_path('code/SNP_analysis/output/ssnpeff.ann.vcf')


cols = ['chromosome', 'position','QUAL','type', 'ref', 'alt','gene','effect','impact_level','LOF'] + reader.header.samples.names
df = pd.DataFrame(columns=cols)
for record in reader:
    chrosome = record.CHROM
    position = record.POS
    QUAL = record.QUAL
    type = record.ALT[0].type
    ref = record.REF
    alt = record.ALT[0].value
    try:
        ann=record.INFO['ANN'][0].split('|')
        geneID=ann[3]
        effect=ann[1]
        putative_impact=ann[2]
    except:
        geneID=None
        effect=None
        putative_impact=None
    try:
        record.INFO['LOF']
        lof=1
    except:
        lof=0
    line = [chrosome, position,QUAL,type, ref, alt,geneID,effect,putative_impact,lof]
    for sample in record.calls:
        gt=sample.data['GT']
        if gt =='0/0':
            gt=0
        elif gt=='1/1':
            gt=1
        else:
            gt=None
        line.append(gt)
    df.loc[len(df)] = line
    print('Parse %s:%s'%(chrosome,position))

# save result
df.to_csv('code/SNP_analysis/output/human_bioethanol_dairy_wt_genotype.tsv', sep='\t', index=False)