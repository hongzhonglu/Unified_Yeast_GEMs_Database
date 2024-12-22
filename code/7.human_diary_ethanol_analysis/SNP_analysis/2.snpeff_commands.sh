# build database
cp data/genome/s288c_R64-5-1/* /home/why/miniconda3/envs/data/share/snpeff-5.2-1/./data/s288c_R64-5-1
echo -e "\ns288c_R64-5-1.genome : s288c_R64-5-1" >> /home/why/miniconda3/envs/data/share/snpeff-5.2-1/snpEff.config
snpEff build -gff3 -v s288c_R64-5-1

## check detabase
#snpEff databases|grep 'cerevisiae'
#snpEff download -v R64-1-1.75

# annotate
snpEff ann -v s288c_R64-5-1 -ud 500 -csvStats code/SNP_analysis/output/snpeff_stats.csv -s code/SNP_analysis/output/snpeff_variant.html -lof data/human_bioethanol_dairy_wt_SNP_indel.gvcf > code/SNP_analysis/output/ssnpeff.ann.vcf



