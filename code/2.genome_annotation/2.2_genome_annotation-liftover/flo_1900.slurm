#!/bin/bash

#SBATCH --job-name=flo_new600
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/flo
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/flo_sce

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAD_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAD_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAD_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAG_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAG_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAG_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAH_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAH_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAH_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAI_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAI_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAI_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAK_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAK_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAK_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAN_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAN_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAN_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AAV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AAV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABA_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABA_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABA_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABB_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABB_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABB_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABC_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABC_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABC_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABE_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABE_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABE_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABG_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABG_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABG_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABK_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABK_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABK_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABL_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABL_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABL_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABM_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABM_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABM_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABP_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABP_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABP_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABQ_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABQ_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABQ_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABR_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABR_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABR_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABS_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABS_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABS_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABT_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABT_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABT_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABV_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ABV_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ABV_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACC_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACC_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACC_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACD_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACD_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACD_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACF_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACF_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACF_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACG_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACG_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACG_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACH_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACH_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACH_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACI_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACI_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACI_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACK_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACK_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACK_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACL_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACL_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACL_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACM_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACM_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACM_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACN_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACN_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACN_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACP_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACP_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACP_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACQ_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACQ_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACQ_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACS_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACS_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACS_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ACV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ACV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADK_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADK_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADK_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ADV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ADV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEA_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEA_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEA_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEB_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEB_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEB_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEC_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEC_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEC_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEE_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEE_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEE_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEF_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEF_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEF_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEG_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEG_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEG_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEI_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEI_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEI_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEK_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEK_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEK_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEN_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEN_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEN_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AER_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AER_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AER_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AES_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AES_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AES_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AET_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AET_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AET_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AEV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AEV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFK_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFK_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFK_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFL_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFL_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFL_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFM_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFM_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFM_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFR_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFR_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFR_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFV_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AFV_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AFV_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGB_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGB_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGB_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGG_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGG_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGG_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGI_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGI_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGI_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGK_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGK_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGK_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGP_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGP_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGP_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGS_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGS_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGS_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AGV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AGV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHD_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHD_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHD_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHE_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHE_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHE_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHF_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHF_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHF_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHI_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHI_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHI_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHK_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHK_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHK_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHN_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHN_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHN_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHV_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AHV_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AHV_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIA_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIA_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIA_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIB_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIB_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIB_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AID_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AID_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AID_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AII_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AII_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AII_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIN_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIN_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIN_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIP_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIP_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIP_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AIT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AIT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKE_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKE_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKE_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKN_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKN_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKN_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKQ_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKQ_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKQ_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AKV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AKV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALB_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALB_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALB_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALG_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALG_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALG_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ALV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ALV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMB_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMB_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMB_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMC_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMC_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMC_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AML_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AML_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AML_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMQ_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMQ_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMQ_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AMT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AMT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AND_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AND_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AND_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANP_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANP_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANP_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANQ_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANQ_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANQ_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANS_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANS_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANS_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ANV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ANV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/API_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/API_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/API_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APT_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APT_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APT_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/APV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/APV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQB_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQB_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQB_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQE_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQE_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQE_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQH_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQH_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQH_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQN_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQN_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQN_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AQT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AQT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARC_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARC_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARC_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARE_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARE_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARE_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARG_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARG_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARG_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARK_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARK_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARK_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARM_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARM_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARM_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARN_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARN_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARN_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARP_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARP_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARP_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARQ_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARQ_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARQ_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARR_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARR_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARR_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ART_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ART_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ART_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARV_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ARV_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ARV_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASL_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASL_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASL_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASM_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASM_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASM_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASN_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASN_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASN_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASP_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASP_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASP_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AST_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AST_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AST_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASV_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ASV_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ASV_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATA_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATA_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATA_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATB_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATB_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATB_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATD_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATD_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATD_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATQ_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATQ_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATQ_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/ATV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/ATV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVB_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVB_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVB_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVD_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVD_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVD_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVK_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVK_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVK_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVL_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVL_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVL_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVM_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVM_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVM_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVQ_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVQ_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVQ_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVR_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVR_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVR_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVT_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVT_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVT_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVV_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/AVV_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/AVV_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAB_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAB_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAB_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAC_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAC_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAC_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAD_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAD_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAD_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAG_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAG_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAG_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAH_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAH_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAH_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAI_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAI_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAI_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAK_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAK_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAK_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAL_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAL_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAL_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAM_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAM_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAM_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAN_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAN_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAN_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAP_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAP_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAP_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAQ_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAQ_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAQ_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAR_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAR_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAR_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAS_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAS_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAS_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAT_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAT_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAT_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAV_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BAV_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BAV_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBA_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBA_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBA_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBC_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBC_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBC_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBD_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBD_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBD_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBG_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBG_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBG_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBH_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBH_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBH_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBI_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBI_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBI_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBK_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBK_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBK_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBL_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBL_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBL_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBM_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBM_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBM_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BBV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BBV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCA_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCA_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCA_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCB_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCB_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCB_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCD_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCD_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCD_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCE_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCE_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCE_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCF_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCF_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCF_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCG_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCG_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCG_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCI_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCI_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCI_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCK_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCK_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCK_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCL_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCL_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCL_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCM_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCM_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCM_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCN_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCN_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCN_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCP_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCP_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCP_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCQ_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCQ_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCQ_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCR_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCR_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCR_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCS_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCS_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCS_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCT_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCT_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCT_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCV_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BCV_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BCV_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDA_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDA_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDA_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDB_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDB_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDB_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDC_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDC_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDC_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDD_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDD_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDD_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDE_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDE_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDE_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDF_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDF_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDF_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDG_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDG_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDG_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDI_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDI_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDI_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDK_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDK_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDK_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDL_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDL_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDL_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDM_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDM_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDM_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDQ_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDQ_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDQ_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDR_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDR_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDR_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDS_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDS_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDS_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDT_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDT_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDT_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDV_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BDV_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BDV_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BED_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BED_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BED_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEF_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEF_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEF_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEG_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEG_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEG_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEH_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEH_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEH_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEI_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEI_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEI_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEK_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEK_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEK_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEM_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEM_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEM_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEN_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEN_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEN_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BER_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BER_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BER_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BES_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BES_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BES_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BET_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BET_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BET_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BEV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BEV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFA_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFA_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFA_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFD_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFD_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFD_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFQ_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFQ_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFQ_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFR_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFR_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFR_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFS_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFS_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFS_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFT_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFT_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFT_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFV_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BFV_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BFV_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGA_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGA_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGA_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGB_7.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGB_7.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGB_7.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGC_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGC_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGC_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGD_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGD_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGD_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGE_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGE_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGE_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGF_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGF_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGF_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGG_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGG_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGG_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGI_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGI_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGI_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGK_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGK_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGK_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGL_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGL_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGL_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGN_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGN_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGN_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BGV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BGV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHF_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHF_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHF_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHH_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHH_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHH_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHT_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHT_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHT_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHV_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BHV_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BHV_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIA_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIA_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIA_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIB_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIB_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIB_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BID_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BID_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BID_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIG_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIG_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIG_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BII_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BII_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BII_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIK_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIK_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIK_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIN_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIN_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIN_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIS_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIS_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIS_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BIT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BIT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKD_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKD_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKD_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKF_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKF_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKF_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKI_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKI_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKI_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKL_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKL_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKL_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKQ_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKQ_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKQ_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BKT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BKT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLD_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLD_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLD_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLF_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLF_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLF_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLG_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLG_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLG_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLI_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLI_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLI_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLL_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLL_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLL_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BLV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BLV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMB_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMB_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMB_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMC_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMC_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMC_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BME_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BME_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BME_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMK_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMK_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMK_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BML_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BML_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BML_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMQ_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMQ_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMQ_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BMV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BMV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BND_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BND_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BND_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BNV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BNV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BPV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BPV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BQT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BQT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BRV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BRV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSI_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSI_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSI_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BSS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BSS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BST_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BST_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BST_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTB_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTB_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTB_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTD_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTD_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTD_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTE_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTE_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTE_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTH_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTH_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTH_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BTV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BTV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVF_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVF_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVF_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVM_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVM_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVM_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/BVV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/BVV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAF_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAF_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAF_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAI_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAI_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAI_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAT_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAT_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAT_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CAV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CAV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBA_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBA_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBA_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBD_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBD_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBD_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBE_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBE_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBE_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBF_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBF_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBF_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBG_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBG_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBG_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBH_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBH_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBH_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBI_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBI_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBI_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBR_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CBR_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CBR_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCL_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCL_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCL_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCN_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCN_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCN_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCS_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCS_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCS_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCV_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CCV_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CCV_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDH.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDH.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDH.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDI_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDI_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDI_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDL_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDL_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDL_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDP_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDP_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDP_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDQ_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDQ_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDQ_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDT_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDT_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDT_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CDV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CDV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEB_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEB_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEB_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CED_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CED_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CED_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEI_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEI_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEI_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CER_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CER_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CER_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CES_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CES_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CES_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CET_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CET_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CET_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEV_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CEV_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CEV_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFH.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFH.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFH.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CFV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CFV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGB_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGB_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGB_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGE_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGE_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGE_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGH_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGH_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGH_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGN_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGN_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGN_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGP_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGP_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGP_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGR_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGR_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGR_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CGV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CGV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHD_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHD_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHD_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CHV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CHV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CID_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CID_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CID_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIF_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIF_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIF_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIH_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIH_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIH_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CII_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CII_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CII_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIL_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIL_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIL_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIM_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIM_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIM_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIN_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIN_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIN_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIP_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIP_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIP_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIQ_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIQ_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIQ_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIR_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIR_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIR_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIS_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIS_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIS_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIT_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIT_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIT_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIV_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CIV_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CIV_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKA_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKA_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKA_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKB_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKB_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKB_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKC_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKC_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKC_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKD_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKD_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKD_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKE_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKE_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKE_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKF_3.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKF_3.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKF_3.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKS_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKS_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKS_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CKV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CKV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLA_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLA_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLA_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLB_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLB_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLB_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLD_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLD_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLD_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLT_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLT_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLT_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLV_5.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CLV_5.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CLV_5.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMD_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMD_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMD_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CME_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CME_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CME_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMG_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMG_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMG_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMH_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMH_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMH_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMI_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMI_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMI_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMK_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMK_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMK_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CML_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CML_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CML_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMM_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMM_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMM_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMN_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMN_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMN_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMP_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMP_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMP_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMQ_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMQ_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMQ_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMS_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMS_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMS_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMT_8.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CMT_8.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CMT_8.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNE_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNE_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNE_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNF_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNF_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNF_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNG_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNG_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNG_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNH_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNH_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNH_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNM_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNM_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNM_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNN_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNN_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNN_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNP_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNP_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNP_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNQ_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNQ_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNQ_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNR_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNR_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNR_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CNV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CNV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPA_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPA_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPA_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPB_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPB_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPB_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPC_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPC_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPC_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPD_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPD_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPD_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPE_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPE_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPE_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPF_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPF_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPF_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPG_1.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPG_1.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPG_1.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPI_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPI_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPI_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CPT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CPT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQA_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQA_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQA_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQB_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQB_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQB_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQC_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQC_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQC_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQD_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQD_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQD_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQE_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQE_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQE_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQF_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQF_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQF_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQG_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQG_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQG_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQH_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQH_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQH_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQI_6.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQI_6.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQI_6.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQK_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQK_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQK_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQL_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQL_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQL_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQM_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQM_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQM_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQN_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQN_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQN_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQP_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQP_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQP_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQQ_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQQ_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQQ_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQR_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQR_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQR_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQS_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQS_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQS_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQT_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQT_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQT_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQV_4.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CQV_4.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CQV_4.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRA_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRA_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRA_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRB_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRB_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRB_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRC_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRC_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRC_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRD_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRD_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRD_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRE_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRE_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRE_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRG_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRG_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRG_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRH_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRH_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRH_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRI_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRI_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRI_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRK_2.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRK_2.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRK_2.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRL_.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/CRL_.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/CRL_.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GAL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GAL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAP.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GAP.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GAP.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAS.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GAS.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GAS.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAT.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GAT.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GAT.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAV.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GAV.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GAV.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000151485.1_ASM15148v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000151485.1_ASM15148v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000151485.1_ASM15148v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000181435.1_ASM18143v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000181435.1_ASM18143v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000181435.1_ASM18143v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182075.1_M22_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000182075.1_M22_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000182075.1_M22_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182175.1_ASM18217v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000182175.1_ASM18217v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000182175.1_ASM18217v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182315.2_ASM18231v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000182315.2_ASM18231v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000182315.2_ASM18231v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190195.1_ScAWRI796_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000190195.1_ScAWRI796_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000190195.1_ScAWRI796_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190215.1_ScVIN13_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000190215.1_ScVIN13_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000190215.1_ScVIN13_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190235.1_ScVL3_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000190235.1_ScVL3_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000190235.1_ScVL3_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190255.1_ScFostersB_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000190255.1_ScFostersB_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000190255.1_ScFostersB_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000234495.1_EC9-8_1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000234495.1_EC9-8_1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000234495.1_EC9-8_1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000260735.1_K7_1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000260735.1_K7_1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000260735.1_K7_1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000275665.1_YJSH1_1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000275665.1_YJSH1_1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000275665.1_YJSH1_1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000308935.1_ZTW1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000308935.1_ZTW1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000308935.1_ZTW1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000325965.1_ScQA23_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000325965.1_ScQA23_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000325965.1_ScQA23_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000326005.1_ScFostersO_v01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000326005.1_ScFostersO_v01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000326005.1_ScFostersO_v01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000416405.1_S._cerevisiae_NY1308_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000416405.1_S._cerevisiae_NY1308_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000508805.2_ScNAM34-4C_assembly01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000508825.2_ScIR-2_assembly01_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000508825.2_ScIR-2_assembly01_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000508825.2_ScIR-2_assembly01_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568005.1_R008_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000568005.1_R008_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000568005.1_R008_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568055.1_p301_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000568055.1_p301_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000568055.1_p301_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568295.1_P283_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000568295.1_P283_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000568295.1_P283_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568365.1_R103_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000568365.1_R103_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000568365.1_R103_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000723645.1_ZhangWP_SECN85_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000723645.1_ZhangWP_SECN85_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000723645.1_ZhangWP_SECN85_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000733235.3_ASM73323v3_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000733235.3_ASM73323v3_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000733235.3_ASM73323v3_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000756235.1_ASM75623v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000756235.1_ASM75623v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000756235.1_ASM75623v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000756245.1_ASM75624v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000756245.1_ASM75624v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000756245.1_ASM75624v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766165.2_ASM76616v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766165.2_ASM76616v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766165.2_ASM76616v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766175.2_ASM76617v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766175.2_ASM76617v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766175.2_ASM76617v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766185.2_ASM76618v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766185.2_ASM76618v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766185.2_ASM76618v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766245.2_ASM76624v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766245.2_ASM76624v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766245.2_ASM76624v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766265.2_ASM76626v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766265.2_ASM76626v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766265.2_ASM76626v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766305.2_ASM76630v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766305.2_ASM76630v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766305.2_ASM76630v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766315.2_ASM76631v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766315.2_ASM76631v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766315.2_ASM76631v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766375.2_ASM76637v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766375.2_ASM76637v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766375.2_ASM76637v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766395.2_ASM76639v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766395.2_ASM76639v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766395.2_ASM76639v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766415.2_ASM76641v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766415.2_ASM76641v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766415.2_ASM76641v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766455.2_ASM76645v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766455.2_ASM76645v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766455.2_ASM76645v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766495.2_ASM76649v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766495.2_ASM76649v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766495.2_ASM76649v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766515.2_ASM76651v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766515.2_ASM76651v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766515.2_ASM76651v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766535.2_ASM76653v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766535.2_ASM76653v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766535.2_ASM76653v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766555.2_ASM76655v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766555.2_ASM76655v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766555.2_ASM76655v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766575.2_ASM76657v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766575.2_ASM76657v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766575.2_ASM76657v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766595.2_ASM76659v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766595.2_ASM76659v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766595.2_ASM76659v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766615.2_ASM76661v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766615.2_ASM76661v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766615.2_ASM76661v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766635.2_ASM76663v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000766635.2_ASM76663v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000766635.2_ASM76663v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000767965.1_ASM76796v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000767965.1_ASM76796v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000767965.1_ASM76796v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000768095.1_ASM76809v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000768095.1_ASM76809v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000768095.1_ASM76809v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000773925.1_ASM77392v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000773925.1_ASM77392v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000773925.1_ASM77392v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000775555.1_ASM77555v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000775555.1_ASM77555v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000775555.1_ASM77555v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976785.2_Sc_YJM1202_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976785.2_Sc_YJM1202_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976785.2_Sc_YJM1202_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976815.3_Sc_YJM1208_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976815.3_Sc_YJM1208_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976815.3_Sc_YJM1208_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976845.3_Sc_YJM1242_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976845.3_Sc_YJM1242_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976845.3_Sc_YJM1242_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976875.3_Sc_YJM1244_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976875.3_Sc_YJM1244_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976875.3_Sc_YJM1244_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976905.2_Sc_YJM1248_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976905.2_Sc_YJM1248_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976905.2_Sc_YJM1248_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976935.2_Sc_YJM1250_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976935.2_Sc_YJM1250_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976935.2_Sc_YJM1250_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976965.3_Sc_YJM1252_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976965.3_Sc_YJM1252_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976965.3_Sc_YJM1252_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976995.2_Sc_YJM1273_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000976995.2_Sc_YJM1273_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000976995.2_Sc_YJM1273_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977025.2_Sc_YJM1304_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000977025.2_Sc_YJM1304_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000977025.2_Sc_YJM1304_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977055.3_Sc_YJM1307_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_000977055.3_Sc_YJM1307_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_000977055.3_Sc_YJM1307_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001006245.1_ScSp60_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001006245.1_ScSp60_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001006245.1_ScSp60_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001029075.1_ASM102907v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001029075.1_ASM102907v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001029075.1_ASM102907v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001029995.1_ASM102999v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001029995.1_ASM102999v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001029995.1_ASM102999v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001051215.1_ASM105121v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001051215.1_ASM105121v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001051215.1_ASM105121v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001078085.1_ASM107808v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001078085.1_ASM107808v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001078085.1_ASM107808v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001078105.1_ASM107810v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001078105.1_ASM107810v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001078105.1_ASM107810v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001282415.1_Y-12632v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001282415.1_Y-12632v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001282415.1_Y-12632v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001584535.1_ASM158453v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001584535.1_ASM158453v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001584535.1_ASM158453v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001592655.1_ASM159265v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001592655.1_ASM159265v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001592655.1_ASM159265v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001611835.1_ASM161183v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001611835.1_ASM161183v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001611835.1_ASM161183v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001634645.1_GLBRCY22-3_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001634645.1_GLBRCY22-3_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001634645.1_GLBRCY22-3_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001735825.1_ASM173582v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001735825.1_ASM173582v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001735825.1_ASM173582v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737145.2_ASM173714v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737145.2_ASM173714v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737145.2_ASM173714v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737155.1_ASM173715v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737155.1_ASM173715v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737155.1_ASM173715v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737165.1_ASM173716v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737165.1_ASM173716v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737165.1_ASM173716v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737175.1_ASM173717v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737175.1_ASM173717v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737175.1_ASM173717v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737265.1_ASM173726v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737265.1_ASM173726v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737265.1_ASM173726v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737275.1_ASM173727v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737275.1_ASM173727v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737275.1_ASM173727v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737285.1_ASM173728v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737285.1_ASM173728v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737285.1_ASM173728v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737295.2_ASM173729v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737295.2_ASM173729v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737295.2_ASM173729v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737385.2_ASM173738v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737385.2_ASM173738v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737385.2_ASM173738v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737395.2_ASM173739v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737395.2_ASM173739v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737395.2_ASM173739v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737415.2_ASM173741v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737415.2_ASM173741v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737415.2_ASM173741v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737445.2_ASM173744v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737445.2_ASM173744v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737445.2_ASM173744v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737505.2_ASM173750v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737505.2_ASM173750v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737505.2_ASM173750v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737515.2_ASM173751v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737515.2_ASM173751v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737515.2_ASM173751v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737535.2_ASM173753v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737535.2_ASM173753v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737535.2_ASM173753v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737545.2_ASM173754v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737545.2_ASM173754v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737545.2_ASM173754v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737625.2_ASM173762v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737625.2_ASM173762v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737625.2_ASM173762v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737635.2_ASM173763v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737635.2_ASM173763v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737635.2_ASM173763v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737645.2_ASM173764v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737645.2_ASM173764v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737645.2_ASM173764v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737655.2_ASM173765v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737655.2_ASM173765v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737655.2_ASM173765v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737745.2_ASM173774v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737745.2_ASM173774v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737745.2_ASM173774v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737755.2_ASM173775v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737755.2_ASM173775v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737755.2_ASM173775v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737765.2_ASM173776v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737765.2_ASM173776v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737765.2_ASM173776v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737785.2_ASM173778v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737785.2_ASM173778v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737785.2_ASM173778v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737865.2_ASM173786v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737865.2_ASM173786v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737865.2_ASM173786v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737875.1_ASM173787v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737875.1_ASM173787v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737875.1_ASM173787v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737885.1_ASM173788v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737885.1_ASM173788v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737885.1_ASM173788v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737895.1_ASM173789v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737895.1_ASM173789v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737895.1_ASM173789v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737985.1_ASM173798v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737985.1_ASM173798v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737985.1_ASM173798v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737995.1_ASM173799v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001737995.1_ASM173799v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001737995.1_ASM173799v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738005.1_ASM173800v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738005.1_ASM173800v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738005.1_ASM173800v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738015.1_ASM173801v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738015.1_ASM173801v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738015.1_ASM173801v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738105.1_ASM173810v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738105.1_ASM173810v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738105.1_ASM173810v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738115.1_ASM173811v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738115.1_ASM173811v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738115.1_ASM173811v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738125.1_ASM173812v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738125.1_ASM173812v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738125.1_ASM173812v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738135.2_ASM173813v2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738135.2_ASM173813v2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738135.2_ASM173813v2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738225.1_ASM173822v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738225.1_ASM173822v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738225.1_ASM173822v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738235.1_ASM173823v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738235.1_ASM173823v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738235.1_ASM173823v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738255.1_ASM173825v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738255.1_ASM173825v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738255.1_ASM173825v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738265.1_ASM173826v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738265.1_ASM173826v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738265.1_ASM173826v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738345.1_ASM173834v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738345.1_ASM173834v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738345.1_ASM173834v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738355.1_ASM173835v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738355.1_ASM173835v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738355.1_ASM173835v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738375.1_ASM173837v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738375.1_ASM173837v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738375.1_ASM173837v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738395.1_ASM173839v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738395.1_ASM173839v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738395.1_ASM173839v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738465.1_ASM173846v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738465.1_ASM173846v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738465.1_ASM173846v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738485.1_ASM173848v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738485.1_ASM173848v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738485.1_ASM173848v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738495.1_ASM173849v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738495.1_ASM173849v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738495.1_ASM173849v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738515.1_ASM173851v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738515.1_ASM173851v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738515.1_ASM173851v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738585.1_ASM173858v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738585.1_ASM173858v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738585.1_ASM173858v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738595.1_ASM173859v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738595.1_ASM173859v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738595.1_ASM173859v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738605.1_ASM173860v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738605.1_ASM173860v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738605.1_ASM173860v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738615.1_ASM173861v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738615.1_ASM173861v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738615.1_ASM173861v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738705.1_ASM173870v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738705.1_ASM173870v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738705.1_ASM173870v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738715.1_ASM173871v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738715.1_ASM173871v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738715.1_ASM173871v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738725.1_ASM173872v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738725.1_ASM173872v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738725.1_ASM173872v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738755.1_ASM173875v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738755.1_ASM173875v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738755.1_ASM173875v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738825.1_ASM173882v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738825.1_ASM173882v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738825.1_ASM173882v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738835.1_ASM173883v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738835.1_ASM173883v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738835.1_ASM173883v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738845.1_ASM173884v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738845.1_ASM173884v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738845.1_ASM173884v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738855.1_ASM173885v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738855.1_ASM173885v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738855.1_ASM173885v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738945.1_ASM173894v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738945.1_ASM173894v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738945.1_ASM173894v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738955.1_ASM173895v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738955.1_ASM173895v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738955.1_ASM173895v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738965.1_ASM173896v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738965.1_ASM173896v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738965.1_ASM173896v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738985.1_ASM173898v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001738985.1_ASM173898v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001738985.1_ASM173898v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739065.1_ASM173906v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739065.1_ASM173906v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739065.1_ASM173906v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739075.1_ASM173907v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739075.1_ASM173907v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739075.1_ASM173907v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739095.1_ASM173909v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739095.1_ASM173909v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739095.1_ASM173909v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739105.1_ASM173910v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739105.1_ASM173910v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739105.1_ASM173910v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739185.1_ASM173918v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739185.1_ASM173918v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739185.1_ASM173918v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739215.1_ASM173921v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739215.1_ASM173921v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739215.1_ASM173921v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739225.1_ASM173922v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739225.1_ASM173922v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739225.1_ASM173922v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739235.1_ASM173923v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739235.1_ASM173923v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739235.1_ASM173923v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739305.1_ASM173930v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739305.1_ASM173930v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739305.1_ASM173930v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739335.1_ASM173933v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739335.1_ASM173933v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739335.1_ASM173933v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739345.1_ASM173934v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739345.1_ASM173934v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739345.1_ASM173934v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739355.1_ASM173935v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739355.1_ASM173935v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739355.1_ASM173935v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739425.1_ASM173942v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739425.1_ASM173942v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739425.1_ASM173942v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739455.1_ASM173945v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739455.1_ASM173945v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739455.1_ASM173945v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739475.1_ASM173947v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739475.1_ASM173947v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739475.1_ASM173947v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739495.1_ASM173949v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739495.1_ASM173949v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739495.1_ASM173949v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739545.1_ASM173954v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739545.1_ASM173954v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739545.1_ASM173954v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739575.1_ASM173957v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739575.1_ASM173957v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739575.1_ASM173957v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739595.1_ASM173959v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739595.1_ASM173959v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739595.1_ASM173959v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739625.1_ASM173962v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739625.1_ASM173962v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739625.1_ASM173962v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739665.1_ASM173966v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739665.1_ASM173966v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739665.1_ASM173966v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739695.1_ASM173969v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739695.1_ASM173969v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739695.1_ASM173969v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739705.1_ASM173970v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739705.1_ASM173970v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739705.1_ASM173970v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739735.1_ASM173973v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739735.1_ASM173973v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739735.1_ASM173973v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739775.1_ASM173977v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739775.1_ASM173977v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739775.1_ASM173977v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739815.1_ASM173981v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739815.1_ASM173981v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739815.1_ASM173981v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739825.1_ASM173982v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739825.1_ASM173982v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739825.1_ASM173982v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739845.1_ASM173984v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739845.1_ASM173984v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739845.1_ASM173984v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739865.1_ASM173986v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739865.1_ASM173986v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739865.1_ASM173986v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739935.1_ASM173993v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739935.1_ASM173993v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739935.1_ASM173993v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739945.1_ASM173994v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739945.1_ASM173994v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739945.1_ASM173994v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739965.1_ASM173996v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739965.1_ASM173996v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739965.1_ASM173996v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739995.1_ASM173999v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001739995.1_ASM173999v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001739995.1_ASM173999v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740055.1_ASM174005v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740055.1_ASM174005v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740055.1_ASM174005v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740075.1_ASM174007v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740075.1_ASM174007v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740075.1_ASM174007v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740095.1_ASM174009v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740095.1_ASM174009v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740095.1_ASM174009v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740115.1_ASM174011v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740115.1_ASM174011v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740115.1_ASM174011v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740175.1_ASM174017v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740175.1_ASM174017v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740175.1_ASM174017v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740185.1_ASM174018v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740185.1_ASM174018v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740185.1_ASM174018v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740205.1_ASM174020v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740205.1_ASM174020v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740205.1_ASM174020v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740225.1_ASM174022v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740225.1_ASM174022v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740225.1_ASM174022v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740295.1_ASM174029v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740295.1_ASM174029v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740295.1_ASM174029v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740305.1_ASM174030v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740305.1_ASM174030v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740305.1_ASM174030v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740325.1_ASM174032v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740325.1_ASM174032v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740325.1_ASM174032v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740345.1_ASM174034v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740345.1_ASM174034v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740345.1_ASM174034v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740415.1_ASM174041v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740415.1_ASM174041v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740415.1_ASM174041v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740425.1_ASM174042v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740425.1_ASM174042v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740425.1_ASM174042v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740445.1_ASM174044v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740445.1_ASM174044v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740445.1_ASM174044v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740455.1_ASM174045v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740455.1_ASM174045v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740455.1_ASM174045v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740535.1_ASM174053v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740535.1_ASM174053v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740535.1_ASM174053v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740545.1_ASM174054v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740545.1_ASM174054v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740545.1_ASM174054v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740565.1_ASM174056v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740565.1_ASM174056v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740565.1_ASM174056v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740575.1_ASM174057v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740575.1_ASM174057v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740575.1_ASM174057v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740655.1_ASM174065v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740655.1_ASM174065v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740655.1_ASM174065v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740665.1_ASM174066v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740665.1_ASM174066v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740665.1_ASM174066v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740675.1_ASM174067v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740675.1_ASM174067v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740675.1_ASM174067v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740705.1_ASM174070v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740705.1_ASM174070v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740705.1_ASM174070v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740775.1_ASM174077v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740775.1_ASM174077v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740775.1_ASM174077v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740785.1_ASM174078v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740785.1_ASM174078v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740785.1_ASM174078v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740795.1_ASM174079v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740795.1_ASM174079v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740795.1_ASM174079v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740805.1_ASM174080v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740805.1_ASM174080v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740805.1_ASM174080v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740895.1_ASM174089v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740895.1_ASM174089v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740895.1_ASM174089v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740905.1_ASM174090v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740905.1_ASM174090v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740905.1_ASM174090v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740915.1_ASM174091v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740915.1_ASM174091v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740915.1_ASM174091v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740925.1_ASM174092v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001740925.1_ASM174092v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001740925.1_ASM174092v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741015.1_ASM174101v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741015.1_ASM174101v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741015.1_ASM174101v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741025.1_ASM174102v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741025.1_ASM174102v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741025.1_ASM174102v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741035.1_ASM174103v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741035.1_ASM174103v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741035.1_ASM174103v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741045.1_ASM174104v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741045.1_ASM174104v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741045.1_ASM174104v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741135.1_ASM174113v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741135.1_ASM174113v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741135.1_ASM174113v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741145.1_ASM174114v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741145.1_ASM174114v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741145.1_ASM174114v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741165.1_ASM174116v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741165.1_ASM174116v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741165.1_ASM174116v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741195.1_ASM174119v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741195.1_ASM174119v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741195.1_ASM174119v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741255.1_ASM174125v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741255.1_ASM174125v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741255.1_ASM174125v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741285.1_ASM174128v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741285.1_ASM174128v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741285.1_ASM174128v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741295.1_ASM174129v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741295.1_ASM174129v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741295.1_ASM174129v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741305.1_ASM174130v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741305.1_ASM174130v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741305.1_ASM174130v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741375.1_ASM174137v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741375.1_ASM174137v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741375.1_ASM174137v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741405.1_ASM174140v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741405.1_ASM174140v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741405.1_ASM174140v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741415.1_ASM174141v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741415.1_ASM174141v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741415.1_ASM174141v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741425.1_ASM174142v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741425.1_ASM174142v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741425.1_ASM174142v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741495.1_ASM174149v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741495.1_ASM174149v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741495.1_ASM174149v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741525.1_ASM174152v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741525.1_ASM174152v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741525.1_ASM174152v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741535.1_ASM174153v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741535.1_ASM174153v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741535.1_ASM174153v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741555.1_ASM174155v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741555.1_ASM174155v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741555.1_ASM174155v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741605.1_ASM174160v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741605.1_ASM174160v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741605.1_ASM174160v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741685.1_ASM174168v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741685.1_ASM174168v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741685.1_ASM174168v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741695.1_ASM174169v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741695.1_ASM174169v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741695.1_ASM174169v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741705.1_ASM174170v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741705.1_ASM174170v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741705.1_ASM174170v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741725.1_ASM174172v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741725.1_ASM174172v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741725.1_ASM174172v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741805.1_ASM174180v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741805.1_ASM174180v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741805.1_ASM174180v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741815.1_ASM174181v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001741815.1_ASM174181v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001741815.1_ASM174181v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001932575.1_CTBEBG1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001932575.1_CTBEBG1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001932575.1_CTBEBG1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001937245.1_ASM193724v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001937245.1_ASM193724v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001937245.1_ASM193724v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001983315.1_ASM198331v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_001983315.1_ASM198331v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_001983315.1_ASM198331v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057805.1_ASM205780v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002057805.1_ASM205780v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002057805.1_ASM205780v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057885.1_ASM205788v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002057885.1_ASM205788v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002057885.1_ASM205788v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057995.1_ASM205799v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002057995.1_ASM205799v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002057995.1_ASM205799v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002058095.1_ASM205809v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002058095.1_ASM205809v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002058095.1_ASM205809v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002163515.1_ASM216351v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002163515.1_ASM216351v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002163515.1_ASM216351v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002214855.1_ASM221485v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002214855.1_ASM221485v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002214855.1_ASM221485v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002804325.1_ASM280432v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002804325.1_ASM280432v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002804325.1_ASM280432v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002933895.1_ASM293389v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_002933895.1_ASM293389v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_002933895.1_ASM293389v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003020065.1_ASM302006v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_003020065.1_ASM302006v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_003020065.1_ASM302006v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046695.1_ASM304669v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_003046695.1_ASM304669v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_003046695.1_ASM304669v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046715.1_ASM304671v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_003046715.1_ASM304671v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_003046715.1_ASM304671v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046745.1_ASM304674v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_003046745.1_ASM304674v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_003046745.1_ASM304674v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003086655.1_ASM308665v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_003086655.1_ASM308665v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_003086655.1_ASM308665v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394885.1_ASM1939488v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_019394885.1_ASM1939488v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_019394885.1_ASM1939488v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394895.1_ASM1939489v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_019394895.1_ASM1939489v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_019394895.1_ASM1939489v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394915.1_ASM1939491v1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_019394915.1_ASM1939491v1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_019394915.1_ASM1939491v1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099165.1_FS2D_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099165.1_FS2D_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099165.1_FS2D_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099195.1_20B2_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099195.1_20B2_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099195.1_20B2_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099255.1_TS12-A7_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099255.1_TS12-A7_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099255.1_TS12-A7_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099265.1_F25_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099265.1_F25_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099265.1_F25_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099275.1_22A4_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099275.1_22A4_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099275.1_22A4_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099285.1_K1-280-2B_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099285.1_K1-280-2B_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099285.1_K1-280-2B_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099305.1_P3-D5_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099305.1_P3-D5_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099305.1_P3-D5_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099345.1_D47-60-3B_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099345.1_D47-60-3B_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099345.1_D47-60-3B_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099355.1_1014-F5_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099355.1_1014-F5_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099355.1_1014-F5_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099375.1_GE7-4A_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099375.1_GE7-4A_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099375.1_GE7-4A_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099395.1_F12-3B_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099395.1_F12-3B_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099395.1_F12-3B_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099425.1_GUF54-A1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099425.1_GUF54-A1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099425.1_GUF54-A1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099445.1_6320-A7_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900099445.1_6320-A7_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900099445.1_6320-A7_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177835.1_L560_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177835.1_L560_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177835.1_L560_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177845.1_L175_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177845.1_L175_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177845.1_L175_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177855.1_GMY321_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177855.1_GMY321_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177855.1_GMY321_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177865.1_L761_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177865.1_L761_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177865.1_L761_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177875.1_L328_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177875.1_L328_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177875.1_L328_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177885.1_L709_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177885.1_L709_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177885.1_L709_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177905.1_L564_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177905.1_L564_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177905.1_L564_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177915.1_L461_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177915.1_L461_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177915.1_L461_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177925.1_L716_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177925.1_L716_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177925.1_L716_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177935.1_L710_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177935.1_L710_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177935.1_L710_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177945.1_L174_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177945.1_L174_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177945.1_L174_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177975.1_L292_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177975.1_L292_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177975.1_L292_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177985.1_L712_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177985.1_L712_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177985.1_L712_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177995.1_L723_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900177995.1_L723_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900177995.1_L723_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178015.1_L725_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178015.1_L725_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178015.1_L725_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178025.1_L561_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178025.1_L561_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178025.1_L561_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178035.1_L528_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178035.1_L528_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178035.1_L528_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178045.1_SD085_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178045.1_SD085_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178045.1_SD085_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178065.1_L711_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178065.1_L711_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178065.1_L711_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178075.1_L181_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178075.1_L181_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178075.1_L181_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178105.1_BA_Souche789_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178105.1_BA_Souche789_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178105.1_BA_Souche789_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178125.1_BH_Souche796_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178125.1_BH_Souche796_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178125.1_BH_Souche796_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178155.1_BD_Souche2333_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178155.1_BD_Souche2333_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178155.1_BD_Souche2333_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178205.1_BF_Souche2251_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178205.1_BF_Souche2251_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178205.1_BF_Souche2251_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178225.1_BG_Souche2331_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178225.1_BG_Souche2331_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178225.1_BG_Souche2331_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178275.1_BE_Souche2334_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178275.1_BE_Souche2334_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178275.1_BE_Souche2334_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178285.1_DG_SoucheMTF889_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178505.1_BC_Souche1152_V1_genomic.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/GCA_900178505.1_BC_Souche1152_V1_genomic.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/GCA_900178505.1_BC_Souche1152_V1_genomic.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAA.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/MAA.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/MAA.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAB.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/MAB.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/MAB.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/MAL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/MAL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAA.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAA.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAA.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAB.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAB.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAB.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAC.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAC.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAC.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAD.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAD.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAD.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAE.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAE.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAE.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAF.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAF.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAF.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAG.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAG.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAG.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAK.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAK.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAK.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAM.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAM.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAM.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAN.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAN.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAN.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAP.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAP.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAP.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAQ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAQ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAQ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAR.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAR.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAR.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAT.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAT.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAT.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAU.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAU.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAU.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAV.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAV.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAV.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAW.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAW.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAW.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAX.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAX.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAX.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAY.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAY.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAY.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAZ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YAZ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YAZ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBA.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBA.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBA.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBB.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBB.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBB.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBE.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBE.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBE.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBF.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBF.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBF.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBG.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBG.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBG.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBH.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBH.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBH.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBJ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBJ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBJ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBK.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBK.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBK.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBM.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBM.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBM.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBN.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBN.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBN.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBO.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBO.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBO.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBP.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBP.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBP.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBQ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBQ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBQ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBR.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBR.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBR.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBS.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBS.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBS.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBT.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBT.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBT.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBU.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBU.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBU.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBV.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBV.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBV.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBW.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBW.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBW.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBX.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBX.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBX.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBY.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBY.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBY.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBZ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YBZ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YBZ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCA.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCA.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCA.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCB.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCB.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCB.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCC.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCC.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCC.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCD.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCD.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCD.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCE.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCE.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCE.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCF.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCF.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCF.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCG.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCG.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCG.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCH.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCH.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCH.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCI.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCI.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCI.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCJ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCJ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCJ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCK.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCK.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCK.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCM.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCM.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCM.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCN.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCN.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCN.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCO.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCO.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCO.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCP.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCP.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCP.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCQ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCQ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCQ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCR.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCR.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCR.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCS.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCS.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCS.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCT.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCT.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCT.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCU.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCU.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCU.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCV.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCV.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCV.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCW.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCW.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCW.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCX.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCX.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCX.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCY.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCY.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCY.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCZ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YCZ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YCZ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDA.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDA.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDA.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDB.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDB.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDB.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDC.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDC.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDC.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDD.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDD.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDD.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDE.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDE.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDE.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDF.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDF.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDF.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDG.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDG.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDG.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDH.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDH.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDH.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDI.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDI.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDI.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDJ.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDJ.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDJ.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDK.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDK.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDK.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDL.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDL.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDL.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDM.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDM.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDM.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDN.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDN.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDN.re.gff
sh clean.sh

cat /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDO.re.fna|cut -f1 -d ' ' >input.fa
rake -f /lustre/home/acct-clslhz/clslhz/why_ssGEMs/biosoft/flo/Rakefile --trace
gffread run/R64_cleaned_parent/lifted_cleaned.gff -g input.fa -J -y protein.fasta -S
cp protein.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_proteomes/YDO.re.fasta
cp run/R64_cleaned_parent/lifted_cleaned.gff /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/flo_gff/YDO.re.gff
sh clean.sh

