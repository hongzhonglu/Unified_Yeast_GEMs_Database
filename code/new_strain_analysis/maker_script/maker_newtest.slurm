#!/bin/bash

#SBATCH --job-name=maker_test
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
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000662435.2_Sc_YJM993_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000662435.2_Sc_YJM993_v1_genomic.fna.maker.output/GCA_000662435.2_Sc_YJM993_v1_genomic.fna_master_datastore_index.log
cp GCA_000662435.2_Sc_YJM993_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000662435.2_Sc_YJM993_v1_genomic.fa
rm GCA_000662435.2_Sc_YJM993_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975585.2_Sc_YJM195_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975585.2_Sc_YJM195_v1_genomic.fna.maker.output/GCA_000975585.2_Sc_YJM195_v1_genomic.fna_master_datastore_index.log
cp GCA_000975585.2_Sc_YJM195_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975585.2_Sc_YJM195_v1_genomic.fa
rm GCA_000975585.2_Sc_YJM195_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975615.2_Sc_YJM244_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975615.2_Sc_YJM244_v1_genomic.fna.maker.output/GCA_000975615.2_Sc_YJM244_v1_genomic.fna_master_datastore_index.log
cp GCA_000975615.2_Sc_YJM244_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975615.2_Sc_YJM244_v1_genomic.fa
rm GCA_000975615.2_Sc_YJM244_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975645.3_Sc_YJM1078_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.maker.output/GCA_000975645.3_Sc_YJM1078_v1_genomic.fna_master_datastore_index.log
cp GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975645.3_Sc_YJM1078_v1_genomic.fa
rm GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975675.2_Sc_YJM1083_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.maker.output/GCA_000975675.2_Sc_YJM1083_v1_genomic.fna_master_datastore_index.log
cp GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975675.2_Sc_YJM1083_v1_genomic.fa
rm GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975705.2_Sc_YJM1129_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.maker.output/GCA_000975705.2_Sc_YJM1129_v1_genomic.fna_master_datastore_index.log
cp GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975705.2_Sc_YJM1129_v1_genomic.fa
rm GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975735.4_Sc_YJM189_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975735.4_Sc_YJM189_v1_genomic.fna.maker.output/GCA_000975735.4_Sc_YJM189_v1_genomic.fna_master_datastore_index.log
cp GCA_000975735.4_Sc_YJM189_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975735.4_Sc_YJM189_v1_genomic.fa
rm GCA_000975735.4_Sc_YJM189_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975765.2_Sc_YJM193_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975765.2_Sc_YJM193_v1_genomic.fna.maker.output/GCA_000975765.2_Sc_YJM193_v1_genomic.fna_master_datastore_index.log
cp GCA_000975765.2_Sc_YJM193_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975765.2_Sc_YJM193_v1_genomic.fa
rm GCA_000975765.2_Sc_YJM193_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975795.2_Sc_YJM248_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975795.2_Sc_YJM248_v1_genomic.fna.maker.output/GCA_000975795.2_Sc_YJM248_v1_genomic.fna_master_datastore_index.log
cp GCA_000975795.2_Sc_YJM248_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975795.2_Sc_YJM248_v1_genomic.fa
rm GCA_000975795.2_Sc_YJM248_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975825.2_Sc_YJM270_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975825.2_Sc_YJM270_v1_genomic.fna.maker.output/GCA_000975825.2_Sc_YJM270_v1_genomic.fna_master_datastore_index.log
cp GCA_000975825.2_Sc_YJM270_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975825.2_Sc_YJM270_v1_genomic.fa
rm GCA_000975825.2_Sc_YJM270_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975855.3_Sc_YJM271_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975855.3_Sc_YJM271_v1_genomic.fna.maker.output/GCA_000975855.3_Sc_YJM271_v1_genomic.fna_master_datastore_index.log
cp GCA_000975855.3_Sc_YJM271_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975855.3_Sc_YJM271_v1_genomic.fa
rm GCA_000975855.3_Sc_YJM271_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975885.2_Sc_YJM320_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975885.2_Sc_YJM320_v1_genomic.fna.maker.output/GCA_000975885.2_Sc_YJM320_v1_genomic.fna_master_datastore_index.log
cp GCA_000975885.2_Sc_YJM320_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975885.2_Sc_YJM320_v1_genomic.fa
rm GCA_000975885.2_Sc_YJM320_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975915.2_Sc_YJM326_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975915.2_Sc_YJM326_v1_genomic.fna.maker.output/GCA_000975915.2_Sc_YJM326_v1_genomic.fna_master_datastore_index.log
cp GCA_000975915.2_Sc_YJM326_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975915.2_Sc_YJM326_v1_genomic.fa
rm GCA_000975915.2_Sc_YJM326_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975945.2_Sc_YJM428_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975945.2_Sc_YJM428_v1_genomic.fna.maker.output/GCA_000975945.2_Sc_YJM428_v1_genomic.fna_master_datastore_index.log
cp GCA_000975945.2_Sc_YJM428_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975945.2_Sc_YJM428_v1_genomic.fa
rm GCA_000975945.2_Sc_YJM428_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000975975.3_Sc_YJM450_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000975975.3_Sc_YJM450_v1_genomic.fna.maker.output/GCA_000975975.3_Sc_YJM450_v1_genomic.fna_master_datastore_index.log
cp GCA_000975975.3_Sc_YJM450_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000975975.3_Sc_YJM450_v1_genomic.fa
rm GCA_000975975.3_Sc_YJM450_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976005.2_Sc_YJM451_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976005.2_Sc_YJM451_v1_genomic.fna.maker.output/GCA_000976005.2_Sc_YJM451_v1_genomic.fna_master_datastore_index.log
cp GCA_000976005.2_Sc_YJM451_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976005.2_Sc_YJM451_v1_genomic.fa
rm GCA_000976005.2_Sc_YJM451_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976035.2_Sc_YJM453_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976035.2_Sc_YJM453_v1_genomic.fna.maker.output/GCA_000976035.2_Sc_YJM453_v1_genomic.fna_master_datastore_index.log
cp GCA_000976035.2_Sc_YJM453_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976035.2_Sc_YJM453_v1_genomic.fa
rm GCA_000976035.2_Sc_YJM453_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976065.2_Sc_YJM456_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976065.2_Sc_YJM456_v1_genomic.fna.maker.output/GCA_000976065.2_Sc_YJM456_v1_genomic.fna_master_datastore_index.log
cp GCA_000976065.2_Sc_YJM456_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976065.2_Sc_YJM456_v1_genomic.fa
rm GCA_000976065.2_Sc_YJM456_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976095.2_Sc_YJM470_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976095.2_Sc_YJM470_v1_genomic.fna.maker.output/GCA_000976095.2_Sc_YJM470_v1_genomic.fna_master_datastore_index.log
cp GCA_000976095.2_Sc_YJM470_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976095.2_Sc_YJM470_v1_genomic.fa
rm GCA_000976095.2_Sc_YJM470_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976125.2_Sc_YJM541_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976125.2_Sc_YJM541_v1_genomic.fna.maker.output/GCA_000976125.2_Sc_YJM541_v1_genomic.fna_master_datastore_index.log
cp GCA_000976125.2_Sc_YJM541_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976125.2_Sc_YJM541_v1_genomic.fa
rm GCA_000976125.2_Sc_YJM541_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976155.4_Sc_YJM554_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976155.4_Sc_YJM554_v1_genomic.fna.maker.output/GCA_000976155.4_Sc_YJM554_v1_genomic.fna_master_datastore_index.log
cp GCA_000976155.4_Sc_YJM554_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976155.4_Sc_YJM554_v1_genomic.fa
rm GCA_000976155.4_Sc_YJM554_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976185.2_Sc_YJM555_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976185.2_Sc_YJM555_v1_genomic.fna.maker.output/GCA_000976185.2_Sc_YJM555_v1_genomic.fna_master_datastore_index.log
cp GCA_000976185.2_Sc_YJM555_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976185.2_Sc_YJM555_v1_genomic.fa
rm GCA_000976185.2_Sc_YJM555_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976215.2_Sc_YJM627_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976215.2_Sc_YJM627_v1_genomic.fna.maker.output/GCA_000976215.2_Sc_YJM627_v1_genomic.fna_master_datastore_index.log
cp GCA_000976215.2_Sc_YJM627_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976215.2_Sc_YJM627_v1_genomic.fa
rm GCA_000976215.2_Sc_YJM627_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976245.2_Sc_YJM681_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976245.2_Sc_YJM681_v1_genomic.fna.maker.output/GCA_000976245.2_Sc_YJM681_v1_genomic.fna_master_datastore_index.log
cp GCA_000976245.2_Sc_YJM681_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976245.2_Sc_YJM681_v1_genomic.fa
rm GCA_000976245.2_Sc_YJM681_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976275.4_Sc_YJM682_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976275.4_Sc_YJM682_v1_genomic.fna.maker.output/GCA_000976275.4_Sc_YJM682_v1_genomic.fna_master_datastore_index.log
cp GCA_000976275.4_Sc_YJM682_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976275.4_Sc_YJM682_v1_genomic.fa
rm GCA_000976275.4_Sc_YJM682_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976305.2_Sc_YJM683_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976305.2_Sc_YJM683_v1_genomic.fna.maker.output/GCA_000976305.2_Sc_YJM683_v1_genomic.fna_master_datastore_index.log
cp GCA_000976305.2_Sc_YJM683_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976305.2_Sc_YJM683_v1_genomic.fa
rm GCA_000976305.2_Sc_YJM683_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976335.2_Sc_YJM689_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976335.2_Sc_YJM689_v1_genomic.fna.maker.output/GCA_000976335.2_Sc_YJM689_v1_genomic.fna_master_datastore_index.log
cp GCA_000976335.2_Sc_YJM689_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976335.2_Sc_YJM689_v1_genomic.fa
rm GCA_000976335.2_Sc_YJM689_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976365.2_Sc_YJM693_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976365.2_Sc_YJM693_v1_genomic.fna.maker.output/GCA_000976365.2_Sc_YJM693_v1_genomic.fna_master_datastore_index.log
cp GCA_000976365.2_Sc_YJM693_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976365.2_Sc_YJM693_v1_genomic.fa
rm GCA_000976365.2_Sc_YJM693_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976395.2_Sc_YJM969_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976395.2_Sc_YJM969_v1_genomic.fna.maker.output/GCA_000976395.2_Sc_YJM969_v1_genomic.fna_master_datastore_index.log
cp GCA_000976395.2_Sc_YJM969_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976395.2_Sc_YJM969_v1_genomic.fa
rm GCA_000976395.2_Sc_YJM969_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976425.2_Sc_YJM972_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976425.2_Sc_YJM972_v1_genomic.fna.maker.output/GCA_000976425.2_Sc_YJM972_v1_genomic.fna_master_datastore_index.log
cp GCA_000976425.2_Sc_YJM972_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976425.2_Sc_YJM972_v1_genomic.fa
rm GCA_000976425.2_Sc_YJM972_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976455.3_Sc_YJM975_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976455.3_Sc_YJM975_v1_genomic.fna.maker.output/GCA_000976455.3_Sc_YJM975_v1_genomic.fna_master_datastore_index.log
cp GCA_000976455.3_Sc_YJM975_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976455.3_Sc_YJM975_v1_genomic.fa
rm GCA_000976455.3_Sc_YJM975_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976485.2_Sc_YJM978_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976485.2_Sc_YJM978_v1_genomic.fna.maker.output/GCA_000976485.2_Sc_YJM978_v1_genomic.fna_master_datastore_index.log
cp GCA_000976485.2_Sc_YJM978_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976485.2_Sc_YJM978_v1_genomic.fa
rm GCA_000976485.2_Sc_YJM978_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976515.2_Sc_YJM981_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976515.2_Sc_YJM981_v1_genomic.fna.maker.output/GCA_000976515.2_Sc_YJM981_v1_genomic.fna_master_datastore_index.log
cp GCA_000976515.2_Sc_YJM981_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976515.2_Sc_YJM981_v1_genomic.fa
rm GCA_000976515.2_Sc_YJM981_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976545.2_Sc_YJM984_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976545.2_Sc_YJM984_v1_genomic.fna.maker.output/GCA_000976545.2_Sc_YJM984_v1_genomic.fna_master_datastore_index.log
cp GCA_000976545.2_Sc_YJM984_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976545.2_Sc_YJM984_v1_genomic.fa
rm GCA_000976545.2_Sc_YJM984_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976575.4_Sc_YJM987_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976575.4_Sc_YJM987_v1_genomic.fna.maker.output/GCA_000976575.4_Sc_YJM987_v1_genomic.fna_master_datastore_index.log
cp GCA_000976575.4_Sc_YJM987_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976575.4_Sc_YJM987_v1_genomic.fa
rm GCA_000976575.4_Sc_YJM987_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976605.3_Sc_YJM990_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976605.3_Sc_YJM990_v1_genomic.fna.maker.output/GCA_000976605.3_Sc_YJM990_v1_genomic.fna_master_datastore_index.log
cp GCA_000976605.3_Sc_YJM990_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976605.3_Sc_YJM990_v1_genomic.fa
rm GCA_000976605.3_Sc_YJM990_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976665.2_Sc_YJM996_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976665.2_Sc_YJM996_v1_genomic.fna.maker.output/GCA_000976665.2_Sc_YJM996_v1_genomic.fna_master_datastore_index.log
cp GCA_000976665.2_Sc_YJM996_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976665.2_Sc_YJM996_v1_genomic.fa
rm GCA_000976665.2_Sc_YJM996_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976695.3_Sc_YJM1133_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.maker.output/GCA_000976695.3_Sc_YJM1133_v1_genomic.fna_master_datastore_index.log
cp GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976695.3_Sc_YJM1133_v1_genomic.fa
rm GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976725.3_Sc_YJM1190_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.maker.output/GCA_000976725.3_Sc_YJM1190_v1_genomic.fna_master_datastore_index.log
cp GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976725.3_Sc_YJM1190_v1_genomic.fa
rm GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976755.2_Sc_YJM1199_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.maker.output/GCA_000976755.2_Sc_YJM1199_v1_genomic.fna_master_datastore_index.log
cp GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976755.2_Sc_YJM1199_v1_genomic.fa
rm GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976785.2_Sc_YJM1202_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.maker.output/GCA_000976785.2_Sc_YJM1202_v1_genomic.fna_master_datastore_index.log
cp GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976785.2_Sc_YJM1202_v1_genomic.fa
rm GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976815.3_Sc_YJM1208_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.maker.output/GCA_000976815.3_Sc_YJM1208_v1_genomic.fna_master_datastore_index.log
cp GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976815.3_Sc_YJM1208_v1_genomic.fa
rm GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976845.3_Sc_YJM1242_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.maker.output/GCA_000976845.3_Sc_YJM1242_v1_genomic.fna_master_datastore_index.log
cp GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976845.3_Sc_YJM1242_v1_genomic.fa
rm GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976875.3_Sc_YJM1244_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.maker.output/GCA_000976875.3_Sc_YJM1244_v1_genomic.fna_master_datastore_index.log
cp GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976875.3_Sc_YJM1244_v1_genomic.fa
rm GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976905.2_Sc_YJM1248_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.maker.output/GCA_000976905.2_Sc_YJM1248_v1_genomic.fna_master_datastore_index.log
cp GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976905.2_Sc_YJM1248_v1_genomic.fa
rm GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976935.2_Sc_YJM1250_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.maker.output/GCA_000976935.2_Sc_YJM1250_v1_genomic.fna_master_datastore_index.log
cp GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976935.2_Sc_YJM1250_v1_genomic.fa
rm GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976965.3_Sc_YJM1252_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.maker.output/GCA_000976965.3_Sc_YJM1252_v1_genomic.fna_master_datastore_index.log
cp GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976965.3_Sc_YJM1252_v1_genomic.fa
rm GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000976995.2_Sc_YJM1273_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.maker.output/GCA_000976995.2_Sc_YJM1273_v1_genomic.fna_master_datastore_index.log
cp GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000976995.2_Sc_YJM1273_v1_genomic.fa
rm GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977025.2_Sc_YJM1304_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.maker.output/GCA_000977025.2_Sc_YJM1304_v1_genomic.fna_master_datastore_index.log
cp GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977025.2_Sc_YJM1304_v1_genomic.fa
rm GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977055.3_Sc_YJM1307_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.maker.output/GCA_000977055.3_Sc_YJM1307_v1_genomic.fna_master_datastore_index.log
cp GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977055.3_Sc_YJM1307_v1_genomic.fa
rm GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
