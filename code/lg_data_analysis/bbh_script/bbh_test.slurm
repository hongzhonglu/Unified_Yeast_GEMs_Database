#!/bin/bash

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
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh


# build blastdb for ref_pan1011
diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAA_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAA_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAA_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAA_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAA_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAA_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAB_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAB_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAB_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAB_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAB_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAB_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAC_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAC_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAC_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAC_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAC_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAC_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAD_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAD_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAD_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAD_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAD_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAD_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAE_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAE_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAE_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAE_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAE_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAE_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAG_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAG_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAG_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAG_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAG_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAG_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAH_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAH_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAH_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAH_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAH_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAH_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAI_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAI_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAI_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAI_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAI_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAI_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAK_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAK_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAK_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAK_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAK_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAK_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAL_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAL_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAL_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAL_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAL_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAL_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

wait
( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAM_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAM_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAM_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAM_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAM_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAM_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAN_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAN_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAN_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAN_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAN_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAN_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAP_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAP_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAP_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAP_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAP_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAP_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAQ_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAQ_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAQ_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAQ_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAQ_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAQ_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAR_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAR_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAR_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAR_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAR_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAR_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAS_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAS_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAS_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAS_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAS_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAS_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAT_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAT_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAT_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAT_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAT_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAT_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAV_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAV_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/AAV_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/AAV_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/AAV_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_AAV_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABA_3.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABA_3.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABA_3.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABA_3.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABA_3.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABA_3.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABB_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABB_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABB_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABB_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABB_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABB_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

wait
( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABC_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABC_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABC_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABC_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABC_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABC_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABD_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABD_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABD_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABD_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABD_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABD_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABE_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABE_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABE_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABE_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABE_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABE_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABF_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABF_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABF_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABF_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABF_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABF_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABG_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABG_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABG_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABG_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABG_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABG_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABH_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABH_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABH_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABH_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABH_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABH_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABI_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABI_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABI_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABI_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABI_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABI_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABK_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABK_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABK_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABK_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABK_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABK_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABL_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABL_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABL_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABL_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABL_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABL_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABM_5.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABM_5.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABM_5.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABM_5.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABM_5.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABM_5.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

wait
( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABP_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABP_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABP_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABP_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABP_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABP_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABQ_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABQ_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABQ_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABQ_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABQ_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABQ_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABR_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABR_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABR_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABR_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABR_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABR_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABS_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABS_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABS_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABS_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABS_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABS_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABT_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABT_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABT_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABT_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABT_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABT_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABV_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABV_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ABV_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ABV_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ABV_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ABV_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACA_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACA_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACA_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACA_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACA_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACA_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACB_6.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACB_6.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACB_6.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACB_6.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACB_6.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACB_6.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACC_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACC_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACC_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACC_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACC_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACC_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACD_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACD_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACD_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACD_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACD_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACD_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

wait
( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACF_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACF_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACF_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACF_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACF_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACF_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACG_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACG_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACG_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACG_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACG_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACG_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACH_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACH_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACH_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACH_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACH_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACH_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACI_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACI_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACI_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACI_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACI_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACI_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACK_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACK_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACK_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACK_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACK_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACK_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACL_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACL_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACL_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACL_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACL_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACL_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACM_7.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACM_7.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACM_7.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACM_7.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACM_7.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACM_7.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACN_8.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACN_8.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACN_8.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACN_8.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACN_8.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACN_8.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACP_8.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACP_8.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACP_8.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACP_8.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACP_8.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACP_8.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

( diamond makedb --in /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACQ_8.re.fa -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACQ_8.re.fa.dmnd
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/pan1011_cds.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/lg_allproteomes/ACQ_8.re.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/ACQ_8.re_vs_pan1011_cds.txt --more-sensitive --evalue 0.001 --max-target-seqs 20
diamond blastp -d /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/blast_db/ACQ_8.re.fa.dmnd -q /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa -o /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/bbh_result/pan1011_cds_vs_ACQ_8.re.txt --more-sensitive --evalue 0.001 --max-target-seqs 20 ) &

wait
