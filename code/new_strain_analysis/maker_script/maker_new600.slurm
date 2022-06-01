#!/bin/bash

#SBATCH --job-name=maker_new600
#SBATCH --partition=cpu
#SBATCH -n 80
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
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977085.2_Sc_YJM1311_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.maker.output/GCA_000977085.2_Sc_YJM1311_v1_genomic.fna_master_datastore_index.log
cp GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977085.2_Sc_YJM1311_v1_genomic.fa
rm GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977115.3_Sc_YJM1326_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.maker.output/GCA_000977115.3_Sc_YJM1326_v1_genomic.fna_master_datastore_index.log
cp GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977115.3_Sc_YJM1326_v1_genomic.fa
rm GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977145.2_Sc_YJM1332_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.maker.output/GCA_000977145.2_Sc_YJM1332_v1_genomic.fna_master_datastore_index.log
cp GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977145.2_Sc_YJM1332_v1_genomic.fa
rm GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977175.3_Sc_YJM1336_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.maker.output/GCA_000977175.3_Sc_YJM1336_v1_genomic.fna_master_datastore_index.log
cp GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977175.3_Sc_YJM1336_v1_genomic.fa
rm GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977205.2_Sc_YJM1338_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.maker.output/GCA_000977205.2_Sc_YJM1338_v1_genomic.fna_master_datastore_index.log
cp GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977205.2_Sc_YJM1338_v1_genomic.fa
rm GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977235.2_Sc_YJM1341_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.maker.output/GCA_000977235.2_Sc_YJM1341_v1_genomic.fna_master_datastore_index.log
cp GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977235.2_Sc_YJM1341_v1_genomic.fa
rm GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977265.3_Sc_YJM1342_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.maker.output/GCA_000977265.3_Sc_YJM1342_v1_genomic.fna_master_datastore_index.log
cp GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977265.3_Sc_YJM1342_v1_genomic.fa
rm GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977295.2_Sc_YJM1355_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.maker.output/GCA_000977295.2_Sc_YJM1355_v1_genomic.fna_master_datastore_index.log
cp GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977295.2_Sc_YJM1355_v1_genomic.fa
rm GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977325.2_Sc_YJM1356_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.maker.output/GCA_000977325.2_Sc_YJM1356_v1_genomic.fna_master_datastore_index.log
cp GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977325.2_Sc_YJM1356_v1_genomic.fa
rm GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977355.3_Sc_YJM1381_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.maker.output/GCA_000977355.3_Sc_YJM1381_v1_genomic.fna_master_datastore_index.log
cp GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977355.3_Sc_YJM1381_v1_genomic.fa
rm GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977385.2_Sc_YJM1383_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.maker.output/GCA_000977385.2_Sc_YJM1383_v1_genomic.fna_master_datastore_index.log
cp GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977385.2_Sc_YJM1383_v1_genomic.fa
rm GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977415.2_Sc_YJM1385_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.maker.output/GCA_000977415.2_Sc_YJM1385_v1_genomic.fna_master_datastore_index.log
cp GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977415.2_Sc_YJM1385_v1_genomic.fa
rm GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977445.2_Sc_YJM1386_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.maker.output/GCA_000977445.2_Sc_YJM1386_v1_genomic.fna_master_datastore_index.log
cp GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977445.2_Sc_YJM1386_v1_genomic.fa
rm GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977475.3_Sc_YJM1387_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.maker.output/GCA_000977475.3_Sc_YJM1387_v1_genomic.fna_master_datastore_index.log
cp GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977475.3_Sc_YJM1387_v1_genomic.fa
rm GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977505.2_Sc_YJM1388_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.maker.output/GCA_000977505.2_Sc_YJM1388_v1_genomic.fna_master_datastore_index.log
cp GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977505.2_Sc_YJM1388_v1_genomic.fa
rm GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977535.2_Sc_YJM1389_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.maker.output/GCA_000977535.2_Sc_YJM1389_v1_genomic.fna_master_datastore_index.log
cp GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977535.2_Sc_YJM1389_v1_genomic.fa
rm GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977565.2_Sc_YJM1399_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.maker.output/GCA_000977565.2_Sc_YJM1399_v1_genomic.fna_master_datastore_index.log
cp GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977565.2_Sc_YJM1399_v1_genomic.fa
rm GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977595.2_Sc_YJM1400_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.maker.output/GCA_000977595.2_Sc_YJM1400_v1_genomic.fna_master_datastore_index.log
cp GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977595.2_Sc_YJM1400_v1_genomic.fa
rm GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977625.2_Sc_YJM1401_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.maker.output/GCA_000977625.2_Sc_YJM1401_v1_genomic.fna_master_datastore_index.log
cp GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977625.2_Sc_YJM1401_v1_genomic.fa
rm GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977655.2_Sc_YJM1402_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.maker.output/GCA_000977655.2_Sc_YJM1402_v1_genomic.fna_master_datastore_index.log
cp GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977655.2_Sc_YJM1402_v1_genomic.fa
rm GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977685.2_Sc_YJM1415_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.maker.output/GCA_000977685.2_Sc_YJM1415_v1_genomic.fna_master_datastore_index.log
cp GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977685.2_Sc_YJM1415_v1_genomic.fa
rm GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977715.4_Sc_YJM1417_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.maker.output/GCA_000977715.4_Sc_YJM1417_v1_genomic.fna_master_datastore_index.log
cp GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977715.4_Sc_YJM1417_v1_genomic.fa
rm GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977745.2_Sc_YJM1418_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.maker.output/GCA_000977745.2_Sc_YJM1418_v1_genomic.fna_master_datastore_index.log
cp GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977745.2_Sc_YJM1418_v1_genomic.fa
rm GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977775.3_Sc_YJM1419_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.maker.output/GCA_000977775.3_Sc_YJM1419_v1_genomic.fna_master_datastore_index.log
cp GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977775.3_Sc_YJM1419_v1_genomic.fa
rm GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977805.2_Sc_YJM1433_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.maker.output/GCA_000977805.2_Sc_YJM1433_v1_genomic.fna_master_datastore_index.log
cp GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977805.2_Sc_YJM1433_v1_genomic.fa
rm GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977835.2_Sc_YJM1434_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.maker.output/GCA_000977835.2_Sc_YJM1434_v1_genomic.fna_master_datastore_index.log
cp GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977835.2_Sc_YJM1434_v1_genomic.fa
rm GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977865.2_Sc_YJM1439_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.maker.output/GCA_000977865.2_Sc_YJM1439_v1_genomic.fna_master_datastore_index.log
cp GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977865.2_Sc_YJM1439_v1_genomic.fa
rm GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977895.2_Sc_YJM1443_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.maker.output/GCA_000977895.2_Sc_YJM1443_v1_genomic.fna_master_datastore_index.log
cp GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977895.2_Sc_YJM1443_v1_genomic.fa
rm GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977925.2_Sc_YJM1444_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.maker.output/GCA_000977925.2_Sc_YJM1444_v1_genomic.fna_master_datastore_index.log
cp GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977925.2_Sc_YJM1444_v1_genomic.fa
rm GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977955.2_Sc_YJM1447_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.maker.output/GCA_000977955.2_Sc_YJM1447_v1_genomic.fna_master_datastore_index.log
cp GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977955.2_Sc_YJM1447_v1_genomic.fa
rm GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000977985.2_Sc_YJM1450_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.maker.output/GCA_000977985.2_Sc_YJM1450_v1_genomic.fna_master_datastore_index.log
cp GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000977985.2_Sc_YJM1450_v1_genomic.fa
rm GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978015.2_Sc_YJM1460_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.maker.output/GCA_000978015.2_Sc_YJM1460_v1_genomic.fna_master_datastore_index.log
cp GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978015.2_Sc_YJM1460_v1_genomic.fa
rm GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978045.2_Sc_YJM1463_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.maker.output/GCA_000978045.2_Sc_YJM1463_v1_genomic.fna_master_datastore_index.log
cp GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978045.2_Sc_YJM1463_v1_genomic.fa
rm GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978075.3_Sc_YJM1477_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.maker.output/GCA_000978075.3_Sc_YJM1477_v1_genomic.fna_master_datastore_index.log
cp GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978075.3_Sc_YJM1477_v1_genomic.fa
rm GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978105.2_Sc_YJM1478_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.maker.output/GCA_000978105.2_Sc_YJM1478_v1_genomic.fna_master_datastore_index.log
cp GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978105.2_Sc_YJM1478_v1_genomic.fa
rm GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978135.2_Sc_YJM1479_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.maker.output/GCA_000978135.2_Sc_YJM1479_v1_genomic.fna_master_datastore_index.log
cp GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978135.2_Sc_YJM1479_v1_genomic.fa
rm GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978165.2_Sc_YJM1526_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.maker.output/GCA_000978165.2_Sc_YJM1526_v1_genomic.fna_master_datastore_index.log
cp GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978165.2_Sc_YJM1526_v1_genomic.fa
rm GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978195.2_Sc_YJM1527_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.maker.output/GCA_000978195.2_Sc_YJM1527_v1_genomic.fna_master_datastore_index.log
cp GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978195.2_Sc_YJM1527_v1_genomic.fa
rm GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978225.2_Sc_YJM1549_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.maker.output/GCA_000978225.2_Sc_YJM1549_v1_genomic.fna_master_datastore_index.log
cp GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978225.2_Sc_YJM1549_v1_genomic.fa
rm GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978255.2_Sc_YJM1573_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.maker.output/GCA_000978255.2_Sc_YJM1573_v1_genomic.fna_master_datastore_index.log
cp GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978255.2_Sc_YJM1573_v1_genomic.fa
rm GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978285.2_Sc_YJM1574_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.maker.output/GCA_000978285.2_Sc_YJM1574_v1_genomic.fna_master_datastore_index.log
cp GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978285.2_Sc_YJM1574_v1_genomic.fa
rm GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978315.2_Sc_YJM1592_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.maker.output/GCA_000978315.2_Sc_YJM1592_v1_genomic.fna_master_datastore_index.log
cp GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978315.2_Sc_YJM1592_v1_genomic.fa
rm GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_000978345.2_Sc_YJM1615_v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.maker.output/GCA_000978345.2_Sc_YJM1615_v1_genomic.fna_master_datastore_index.log
cp GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_000978345.2_Sc_YJM1615_v1_genomic.fa
rm GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001669875.1_ASM166987v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001669875.1_ASM166987v1_genomic.fna.maker.output/GCA_001669875.1_ASM166987v1_genomic.fna_master_datastore_index.log
cp GCA_001669875.1_ASM166987v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001669875.1_ASM166987v1_genomic.fa
rm GCA_001669875.1_ASM166987v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001669905.1_ASM166990v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001669905.1_ASM166990v1_genomic.fna.maker.output/GCA_001669905.1_ASM166990v1_genomic.fna_master_datastore_index.log
cp GCA_001669905.1_ASM166990v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001669905.1_ASM166990v1_genomic.fa
rm GCA_001669905.1_ASM166990v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001669935.1_ASM166993v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001669935.1_ASM166993v1_genomic.fna.maker.output/GCA_001669935.1_ASM166993v1_genomic.fna_master_datastore_index.log
cp GCA_001669935.1_ASM166993v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001669935.1_ASM166993v1_genomic.fa
rm GCA_001669935.1_ASM166993v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001669965.1_ASM166996v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001669965.1_ASM166996v1_genomic.fna.maker.output/GCA_001669965.1_ASM166996v1_genomic.fna_master_datastore_index.log
cp GCA_001669965.1_ASM166996v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001669965.1_ASM166996v1_genomic.fa
rm GCA_001669965.1_ASM166996v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001669995.1_ASM166999v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001669995.1_ASM166999v1_genomic.fna.maker.output/GCA_001669995.1_ASM166999v1_genomic.fna_master_datastore_index.log
cp GCA_001669995.1_ASM166999v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001669995.1_ASM166999v1_genomic.fa
rm GCA_001669995.1_ASM166999v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670025.1_ASM167002v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670025.1_ASM167002v1_genomic.fna.maker.output/GCA_001670025.1_ASM167002v1_genomic.fna_master_datastore_index.log
cp GCA_001670025.1_ASM167002v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670025.1_ASM167002v1_genomic.fa
rm GCA_001670025.1_ASM167002v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670055.1_ASM167005v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670055.1_ASM167005v1_genomic.fna.maker.output/GCA_001670055.1_ASM167005v1_genomic.fna_master_datastore_index.log
cp GCA_001670055.1_ASM167005v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670055.1_ASM167005v1_genomic.fa
rm GCA_001670055.1_ASM167005v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670085.1_ASM167008v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670085.1_ASM167008v1_genomic.fna.maker.output/GCA_001670085.1_ASM167008v1_genomic.fna_master_datastore_index.log
cp GCA_001670085.1_ASM167008v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670085.1_ASM167008v1_genomic.fa
rm GCA_001670085.1_ASM167008v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670115.1_ASM167011v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670115.1_ASM167011v1_genomic.fna.maker.output/GCA_001670115.1_ASM167011v1_genomic.fna_master_datastore_index.log
cp GCA_001670115.1_ASM167011v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670115.1_ASM167011v1_genomic.fa
rm GCA_001670115.1_ASM167011v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670145.1_ASM167014v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670145.1_ASM167014v1_genomic.fna.maker.output/GCA_001670145.1_ASM167014v1_genomic.fna_master_datastore_index.log
cp GCA_001670145.1_ASM167014v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670145.1_ASM167014v1_genomic.fa
rm GCA_001670145.1_ASM167014v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670175.1_ASM167017v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670175.1_ASM167017v1_genomic.fna.maker.output/GCA_001670175.1_ASM167017v1_genomic.fna_master_datastore_index.log
cp GCA_001670175.1_ASM167017v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670175.1_ASM167017v1_genomic.fa
rm GCA_001670175.1_ASM167017v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670205.1_ASM167020v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670205.1_ASM167020v1_genomic.fna.maker.output/GCA_001670205.1_ASM167020v1_genomic.fna_master_datastore_index.log
cp GCA_001670205.1_ASM167020v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670205.1_ASM167020v1_genomic.fa
rm GCA_001670205.1_ASM167020v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670235.1_ASM167023v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670235.1_ASM167023v1_genomic.fna.maker.output/GCA_001670235.1_ASM167023v1_genomic.fna_master_datastore_index.log
cp GCA_001670235.1_ASM167023v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670235.1_ASM167023v1_genomic.fa
rm GCA_001670235.1_ASM167023v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670265.1_ASM167026v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670265.1_ASM167026v1_genomic.fna.maker.output/GCA_001670265.1_ASM167026v1_genomic.fna_master_datastore_index.log
cp GCA_001670265.1_ASM167026v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670265.1_ASM167026v1_genomic.fa
rm GCA_001670265.1_ASM167026v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670295.1_ASM167029v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670295.1_ASM167029v1_genomic.fna.maker.output/GCA_001670295.1_ASM167029v1_genomic.fna_master_datastore_index.log
cp GCA_001670295.1_ASM167029v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670295.1_ASM167029v1_genomic.fa
rm GCA_001670295.1_ASM167029v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670325.1_ASM167032v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670325.1_ASM167032v1_genomic.fna.maker.output/GCA_001670325.1_ASM167032v1_genomic.fna_master_datastore_index.log
cp GCA_001670325.1_ASM167032v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670325.1_ASM167032v1_genomic.fa
rm GCA_001670325.1_ASM167032v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670355.1_ASM167035v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670355.1_ASM167035v1_genomic.fna.maker.output/GCA_001670355.1_ASM167035v1_genomic.fna_master_datastore_index.log
cp GCA_001670355.1_ASM167035v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670355.1_ASM167035v1_genomic.fa
rm GCA_001670355.1_ASM167035v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670385.1_ASM167038v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670385.1_ASM167038v1_genomic.fna.maker.output/GCA_001670385.1_ASM167038v1_genomic.fna_master_datastore_index.log
cp GCA_001670385.1_ASM167038v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670385.1_ASM167038v1_genomic.fa
rm GCA_001670385.1_ASM167038v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670415.1_ASM167041v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670415.1_ASM167041v1_genomic.fna.maker.output/GCA_001670415.1_ASM167041v1_genomic.fna_master_datastore_index.log
cp GCA_001670415.1_ASM167041v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670415.1_ASM167041v1_genomic.fa
rm GCA_001670415.1_ASM167041v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670845.1_ASM167084v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670845.1_ASM167084v1_genomic.fna.maker.output/GCA_001670845.1_ASM167084v1_genomic.fna_master_datastore_index.log
cp GCA_001670845.1_ASM167084v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670845.1_ASM167084v1_genomic.fa
rm GCA_001670845.1_ASM167084v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670935.1_ASM167093v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670935.1_ASM167093v1_genomic.fna.maker.output/GCA_001670935.1_ASM167093v1_genomic.fna_master_datastore_index.log
cp GCA_001670935.1_ASM167093v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670935.1_ASM167093v1_genomic.fa
rm GCA_001670935.1_ASM167093v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670965.1_ASM167096v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670965.1_ASM167096v1_genomic.fna.maker.output/GCA_001670965.1_ASM167096v1_genomic.fna_master_datastore_index.log
cp GCA_001670965.1_ASM167096v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670965.1_ASM167096v1_genomic.fa
rm GCA_001670965.1_ASM167096v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001670995.1_ASM167099v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001670995.1_ASM167099v1_genomic.fna.maker.output/GCA_001670995.1_ASM167099v1_genomic.fna_master_datastore_index.log
cp GCA_001670995.1_ASM167099v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001670995.1_ASM167099v1_genomic.fa
rm GCA_001670995.1_ASM167099v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671025.1_ASM167102v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671025.1_ASM167102v1_genomic.fna.maker.output/GCA_001671025.1_ASM167102v1_genomic.fna_master_datastore_index.log
cp GCA_001671025.1_ASM167102v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671025.1_ASM167102v1_genomic.fa
rm GCA_001671025.1_ASM167102v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671055.1_ASM167105v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671055.1_ASM167105v1_genomic.fna.maker.output/GCA_001671055.1_ASM167105v1_genomic.fna_master_datastore_index.log
cp GCA_001671055.1_ASM167105v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671055.1_ASM167105v1_genomic.fa
rm GCA_001671055.1_ASM167105v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671085.1_ASM167108v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671085.1_ASM167108v1_genomic.fna.maker.output/GCA_001671085.1_ASM167108v1_genomic.fna_master_datastore_index.log
cp GCA_001671085.1_ASM167108v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671085.1_ASM167108v1_genomic.fa
rm GCA_001671085.1_ASM167108v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671115.1_ASM167111v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671115.1_ASM167111v1_genomic.fna.maker.output/GCA_001671115.1_ASM167111v1_genomic.fna_master_datastore_index.log
cp GCA_001671115.1_ASM167111v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671115.1_ASM167111v1_genomic.fa
rm GCA_001671115.1_ASM167111v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671165.1_ASM167116v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671165.1_ASM167116v1_genomic.fna.maker.output/GCA_001671165.1_ASM167116v1_genomic.fna_master_datastore_index.log
cp GCA_001671165.1_ASM167116v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671165.1_ASM167116v1_genomic.fa
rm GCA_001671165.1_ASM167116v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671195.1_ASM167119v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671195.1_ASM167119v1_genomic.fna.maker.output/GCA_001671195.1_ASM167119v1_genomic.fna_master_datastore_index.log
cp GCA_001671195.1_ASM167119v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671195.1_ASM167119v1_genomic.fa
rm GCA_001671195.1_ASM167119v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671225.1_ASM167122v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671225.1_ASM167122v1_genomic.fna.maker.output/GCA_001671225.1_ASM167122v1_genomic.fna_master_datastore_index.log
cp GCA_001671225.1_ASM167122v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671225.1_ASM167122v1_genomic.fa
rm GCA_001671225.1_ASM167122v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671255.1_ASM167125v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671255.1_ASM167125v1_genomic.fna.maker.output/GCA_001671255.1_ASM167125v1_genomic.fna_master_datastore_index.log
cp GCA_001671255.1_ASM167125v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671255.1_ASM167125v1_genomic.fa
rm GCA_001671255.1_ASM167125v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671285.1_ASM167128v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671285.1_ASM167128v1_genomic.fna.maker.output/GCA_001671285.1_ASM167128v1_genomic.fna_master_datastore_index.log
cp GCA_001671285.1_ASM167128v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671285.1_ASM167128v1_genomic.fa
rm GCA_001671285.1_ASM167128v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671315.1_ASM167131v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671315.1_ASM167131v1_genomic.fna.maker.output/GCA_001671315.1_ASM167131v1_genomic.fna_master_datastore_index.log
cp GCA_001671315.1_ASM167131v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671315.1_ASM167131v1_genomic.fa
rm GCA_001671315.1_ASM167131v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671345.1_ASM167134v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671345.1_ASM167134v1_genomic.fna.maker.output/GCA_001671345.1_ASM167134v1_genomic.fna_master_datastore_index.log
cp GCA_001671345.1_ASM167134v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671345.1_ASM167134v1_genomic.fa
rm GCA_001671345.1_ASM167134v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671375.1_ASM167137v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671375.1_ASM167137v1_genomic.fna.maker.output/GCA_001671375.1_ASM167137v1_genomic.fna_master_datastore_index.log
cp GCA_001671375.1_ASM167137v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671375.1_ASM167137v1_genomic.fa
rm GCA_001671375.1_ASM167137v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671405.1_ASM167140v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671405.1_ASM167140v1_genomic.fna.maker.output/GCA_001671405.1_ASM167140v1_genomic.fna_master_datastore_index.log
cp GCA_001671405.1_ASM167140v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671405.1_ASM167140v1_genomic.fa
rm GCA_001671405.1_ASM167140v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671435.1_ASM167143v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671435.1_ASM167143v1_genomic.fna.maker.output/GCA_001671435.1_ASM167143v1_genomic.fna_master_datastore_index.log
cp GCA_001671435.1_ASM167143v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671435.1_ASM167143v1_genomic.fa
rm GCA_001671435.1_ASM167143v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671565.1_ASM167156v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671565.1_ASM167156v1_genomic.fna.maker.output/GCA_001671565.1_ASM167156v1_genomic.fna_master_datastore_index.log
cp GCA_001671565.1_ASM167156v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671565.1_ASM167156v1_genomic.fa
rm GCA_001671565.1_ASM167156v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671595.1_ASM167159v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671595.1_ASM167159v1_genomic.fna.maker.output/GCA_001671595.1_ASM167159v1_genomic.fna_master_datastore_index.log
cp GCA_001671595.1_ASM167159v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671595.1_ASM167159v1_genomic.fa
rm GCA_001671595.1_ASM167159v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671625.1_ASM167162v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671625.1_ASM167162v1_genomic.fna.maker.output/GCA_001671625.1_ASM167162v1_genomic.fna_master_datastore_index.log
cp GCA_001671625.1_ASM167162v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671625.1_ASM167162v1_genomic.fa
rm GCA_001671625.1_ASM167162v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671655.1_ASM167165v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671655.1_ASM167165v1_genomic.fna.maker.output/GCA_001671655.1_ASM167165v1_genomic.fna_master_datastore_index.log
cp GCA_001671655.1_ASM167165v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671655.1_ASM167165v1_genomic.fa
rm GCA_001671655.1_ASM167165v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671685.1_ASM167168v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671685.1_ASM167168v1_genomic.fna.maker.output/GCA_001671685.1_ASM167168v1_genomic.fna_master_datastore_index.log
cp GCA_001671685.1_ASM167168v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671685.1_ASM167168v1_genomic.fa
rm GCA_001671685.1_ASM167168v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671715.1_ASM167171v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671715.1_ASM167171v1_genomic.fna.maker.output/GCA_001671715.1_ASM167171v1_genomic.fna_master_datastore_index.log
cp GCA_001671715.1_ASM167171v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671715.1_ASM167171v1_genomic.fa
rm GCA_001671715.1_ASM167171v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671745.1_ASM167174v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671745.1_ASM167174v1_genomic.fna.maker.output/GCA_001671745.1_ASM167174v1_genomic.fna_master_datastore_index.log
cp GCA_001671745.1_ASM167174v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671745.1_ASM167174v1_genomic.fa
rm GCA_001671745.1_ASM167174v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671775.1_ASM167177v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671775.1_ASM167177v1_genomic.fna.maker.output/GCA_001671775.1_ASM167177v1_genomic.fna_master_datastore_index.log
cp GCA_001671775.1_ASM167177v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671775.1_ASM167177v1_genomic.fa
rm GCA_001671775.1_ASM167177v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671805.1_ASM167180v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671805.1_ASM167180v1_genomic.fna.maker.output/GCA_001671805.1_ASM167180v1_genomic.fna_master_datastore_index.log
cp GCA_001671805.1_ASM167180v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671805.1_ASM167180v1_genomic.fa
rm GCA_001671805.1_ASM167180v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671835.1_ASM167183v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671835.1_ASM167183v1_genomic.fna.maker.output/GCA_001671835.1_ASM167183v1_genomic.fna_master_datastore_index.log
cp GCA_001671835.1_ASM167183v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671835.1_ASM167183v1_genomic.fa
rm GCA_001671835.1_ASM167183v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671865.1_ASM167186v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671865.1_ASM167186v1_genomic.fna.maker.output/GCA_001671865.1_ASM167186v1_genomic.fna_master_datastore_index.log
cp GCA_001671865.1_ASM167186v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671865.1_ASM167186v1_genomic.fa
rm GCA_001671865.1_ASM167186v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671895.1_ASM167189v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671895.1_ASM167189v1_genomic.fna.maker.output/GCA_001671895.1_ASM167189v1_genomic.fna_master_datastore_index.log
cp GCA_001671895.1_ASM167189v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671895.1_ASM167189v1_genomic.fa
rm GCA_001671895.1_ASM167189v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671925.1_ASM167192v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671925.1_ASM167192v1_genomic.fna.maker.output/GCA_001671925.1_ASM167192v1_genomic.fna_master_datastore_index.log
cp GCA_001671925.1_ASM167192v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671925.1_ASM167192v1_genomic.fa
rm GCA_001671925.1_ASM167192v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671955.1_ASM167195v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671955.1_ASM167195v1_genomic.fna.maker.output/GCA_001671955.1_ASM167195v1_genomic.fna_master_datastore_index.log
cp GCA_001671955.1_ASM167195v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671955.1_ASM167195v1_genomic.fa
rm GCA_001671955.1_ASM167195v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_001671985.1_ASM167198v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_001671985.1_ASM167198v1_genomic.fna.maker.output/GCA_001671985.1_ASM167198v1_genomic.fna_master_datastore_index.log
cp GCA_001671985.1_ASM167198v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_001671985.1_ASM167198v1_genomic.fa
rm GCA_001671985.1_ASM167198v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270095.1_ASM327009v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270095.1_ASM327009v1_genomic.fna.maker.output/GCA_003270095.1_ASM327009v1_genomic.fna_master_datastore_index.log
cp GCA_003270095.1_ASM327009v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270095.1_ASM327009v1_genomic.fa
rm GCA_003270095.1_ASM327009v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270135.1_ASM327013v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270135.1_ASM327013v1_genomic.fna.maker.output/GCA_003270135.1_ASM327013v1_genomic.fna_master_datastore_index.log
cp GCA_003270135.1_ASM327013v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270135.1_ASM327013v1_genomic.fa
rm GCA_003270135.1_ASM327013v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270145.1_ASM327014v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270145.1_ASM327014v1_genomic.fna.maker.output/GCA_003270145.1_ASM327014v1_genomic.fna_master_datastore_index.log
cp GCA_003270145.1_ASM327014v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270145.1_ASM327014v1_genomic.fa
rm GCA_003270145.1_ASM327014v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270175.1_ASM327017v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270175.1_ASM327017v1_genomic.fna.maker.output/GCA_003270175.1_ASM327017v1_genomic.fna_master_datastore_index.log
cp GCA_003270175.1_ASM327017v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270175.1_ASM327017v1_genomic.fa
rm GCA_003270175.1_ASM327017v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270185.1_ASM327018v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270185.1_ASM327018v1_genomic.fna.maker.output/GCA_003270185.1_ASM327018v1_genomic.fna_master_datastore_index.log
cp GCA_003270185.1_ASM327018v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270185.1_ASM327018v1_genomic.fa
rm GCA_003270185.1_ASM327018v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270205.1_ASM327020v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270205.1_ASM327020v1_genomic.fna.maker.output/GCA_003270205.1_ASM327020v1_genomic.fna_master_datastore_index.log
cp GCA_003270205.1_ASM327020v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270205.1_ASM327020v1_genomic.fa
rm GCA_003270205.1_ASM327020v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270245.1_ASM327024v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270245.1_ASM327024v1_genomic.fna.maker.output/GCA_003270245.1_ASM327024v1_genomic.fna_master_datastore_index.log
cp GCA_003270245.1_ASM327024v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270245.1_ASM327024v1_genomic.fa
rm GCA_003270245.1_ASM327024v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270305.1_ASM327030v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270305.1_ASM327030v1_genomic.fna.maker.output/GCA_003270305.1_ASM327030v1_genomic.fna_master_datastore_index.log
cp GCA_003270305.1_ASM327030v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270305.1_ASM327030v1_genomic.fa
rm GCA_003270305.1_ASM327030v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270315.1_ASM327031v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270315.1_ASM327031v1_genomic.fna.maker.output/GCA_003270315.1_ASM327031v1_genomic.fna_master_datastore_index.log
cp GCA_003270315.1_ASM327031v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270315.1_ASM327031v1_genomic.fa
rm GCA_003270315.1_ASM327031v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270325.1_ASM327032v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270325.1_ASM327032v1_genomic.fna.maker.output/GCA_003270325.1_ASM327032v1_genomic.fna_master_datastore_index.log
cp GCA_003270325.1_ASM327032v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270325.1_ASM327032v1_genomic.fa
rm GCA_003270325.1_ASM327032v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270335.1_ASM327033v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270335.1_ASM327033v1_genomic.fna.maker.output/GCA_003270335.1_ASM327033v1_genomic.fna_master_datastore_index.log
cp GCA_003270335.1_ASM327033v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270335.1_ASM327033v1_genomic.fa
rm GCA_003270335.1_ASM327033v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270425.1_ASM327042v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270425.1_ASM327042v1_genomic.fna.maker.output/GCA_003270425.1_ASM327042v1_genomic.fna_master_datastore_index.log
cp GCA_003270425.1_ASM327042v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270425.1_ASM327042v1_genomic.fa
rm GCA_003270425.1_ASM327042v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270435.1_ASM327043v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270435.1_ASM327043v1_genomic.fna.maker.output/GCA_003270435.1_ASM327043v1_genomic.fna_master_datastore_index.log
cp GCA_003270435.1_ASM327043v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270435.1_ASM327043v1_genomic.fa
rm GCA_003270435.1_ASM327043v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270455.1_ASM327045v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270455.1_ASM327045v1_genomic.fna.maker.output/GCA_003270455.1_ASM327045v1_genomic.fna_master_datastore_index.log
cp GCA_003270455.1_ASM327045v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270455.1_ASM327045v1_genomic.fa
rm GCA_003270455.1_ASM327045v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270465.1_ASM327046v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270465.1_ASM327046v1_genomic.fna.maker.output/GCA_003270465.1_ASM327046v1_genomic.fna_master_datastore_index.log
cp GCA_003270465.1_ASM327046v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270465.1_ASM327046v1_genomic.fa
rm GCA_003270465.1_ASM327046v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270535.1_ASM327053v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270535.1_ASM327053v1_genomic.fna.maker.output/GCA_003270535.1_ASM327053v1_genomic.fna_master_datastore_index.log
cp GCA_003270535.1_ASM327053v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270535.1_ASM327053v1_genomic.fa
rm GCA_003270535.1_ASM327053v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270565.1_ASM327056v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270565.1_ASM327056v1_genomic.fna.maker.output/GCA_003270565.1_ASM327056v1_genomic.fna_master_datastore_index.log
cp GCA_003270565.1_ASM327056v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270565.1_ASM327056v1_genomic.fa
rm GCA_003270565.1_ASM327056v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270575.1_ASM327057v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270575.1_ASM327057v1_genomic.fna.maker.output/GCA_003270575.1_ASM327057v1_genomic.fna_master_datastore_index.log
cp GCA_003270575.1_ASM327057v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270575.1_ASM327057v1_genomic.fa
rm GCA_003270575.1_ASM327057v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270595.1_ASM327059v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270595.1_ASM327059v1_genomic.fna.maker.output/GCA_003270595.1_ASM327059v1_genomic.fna_master_datastore_index.log
cp GCA_003270595.1_ASM327059v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270595.1_ASM327059v1_genomic.fa
rm GCA_003270595.1_ASM327059v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270665.1_ASM327066v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270665.1_ASM327066v1_genomic.fna.maker.output/GCA_003270665.1_ASM327066v1_genomic.fna_master_datastore_index.log
cp GCA_003270665.1_ASM327066v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270665.1_ASM327066v1_genomic.fa
rm GCA_003270665.1_ASM327066v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270695.1_ASM327069v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270695.1_ASM327069v1_genomic.fna.maker.output/GCA_003270695.1_ASM327069v1_genomic.fna_master_datastore_index.log
cp GCA_003270695.1_ASM327069v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270695.1_ASM327069v1_genomic.fa
rm GCA_003270695.1_ASM327069v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270705.1_ASM327070v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270705.1_ASM327070v1_genomic.fna.maker.output/GCA_003270705.1_ASM327070v1_genomic.fna_master_datastore_index.log
cp GCA_003270705.1_ASM327070v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270705.1_ASM327070v1_genomic.fa
rm GCA_003270705.1_ASM327070v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270715.1_ASM327071v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270715.1_ASM327071v1_genomic.fna.maker.output/GCA_003270715.1_ASM327071v1_genomic.fna_master_datastore_index.log
cp GCA_003270715.1_ASM327071v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270715.1_ASM327071v1_genomic.fa
rm GCA_003270715.1_ASM327071v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270765.1_ASM327076v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270765.1_ASM327076v1_genomic.fna.maker.output/GCA_003270765.1_ASM327076v1_genomic.fna_master_datastore_index.log
cp GCA_003270765.1_ASM327076v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270765.1_ASM327076v1_genomic.fa
rm GCA_003270765.1_ASM327076v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270815.1_ASM327081v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270815.1_ASM327081v1_genomic.fna.maker.output/GCA_003270815.1_ASM327081v1_genomic.fna_master_datastore_index.log
cp GCA_003270815.1_ASM327081v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270815.1_ASM327081v1_genomic.fa
rm GCA_003270815.1_ASM327081v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270825.1_ASM327082v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270825.1_ASM327082v1_genomic.fna.maker.output/GCA_003270825.1_ASM327082v1_genomic.fna_master_datastore_index.log
cp GCA_003270825.1_ASM327082v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270825.1_ASM327082v1_genomic.fa
rm GCA_003270825.1_ASM327082v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270835.1_ASM327083v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270835.1_ASM327083v1_genomic.fna.maker.output/GCA_003270835.1_ASM327083v1_genomic.fna_master_datastore_index.log
cp GCA_003270835.1_ASM327083v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270835.1_ASM327083v1_genomic.fa
rm GCA_003270835.1_ASM327083v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270855.1_ASM327085v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270855.1_ASM327085v1_genomic.fna.maker.output/GCA_003270855.1_ASM327085v1_genomic.fna_master_datastore_index.log
cp GCA_003270855.1_ASM327085v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270855.1_ASM327085v1_genomic.fa
rm GCA_003270855.1_ASM327085v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270935.1_ASM327093v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270935.1_ASM327093v1_genomic.fna.maker.output/GCA_003270935.1_ASM327093v1_genomic.fna_master_datastore_index.log
cp GCA_003270935.1_ASM327093v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270935.1_ASM327093v1_genomic.fa
rm GCA_003270935.1_ASM327093v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270945.1_ASM327094v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270945.1_ASM327094v1_genomic.fna.maker.output/GCA_003270945.1_ASM327094v1_genomic.fna_master_datastore_index.log
cp GCA_003270945.1_ASM327094v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270945.1_ASM327094v1_genomic.fa
rm GCA_003270945.1_ASM327094v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270955.1_ASM327095v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270955.1_ASM327095v1_genomic.fna.maker.output/GCA_003270955.1_ASM327095v1_genomic.fna_master_datastore_index.log
cp GCA_003270955.1_ASM327095v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270955.1_ASM327095v1_genomic.fa
rm GCA_003270955.1_ASM327095v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003270985.1_ASM327098v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003270985.1_ASM327098v1_genomic.fna.maker.output/GCA_003270985.1_ASM327098v1_genomic.fna_master_datastore_index.log
cp GCA_003270985.1_ASM327098v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003270985.1_ASM327098v1_genomic.fa
rm GCA_003270985.1_ASM327098v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271045.1_ASM327104v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271045.1_ASM327104v1_genomic.fna.maker.output/GCA_003271045.1_ASM327104v1_genomic.fna_master_datastore_index.log
cp GCA_003271045.1_ASM327104v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271045.1_ASM327104v1_genomic.fa
rm GCA_003271045.1_ASM327104v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271085.1_ASM327108v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271085.1_ASM327108v1_genomic.fna.maker.output/GCA_003271085.1_ASM327108v1_genomic.fna_master_datastore_index.log
cp GCA_003271085.1_ASM327108v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271085.1_ASM327108v1_genomic.fa
rm GCA_003271085.1_ASM327108v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271115.1_ASM327111v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271115.1_ASM327111v1_genomic.fna.maker.output/GCA_003271115.1_ASM327111v1_genomic.fna_master_datastore_index.log
cp GCA_003271115.1_ASM327111v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271115.1_ASM327111v1_genomic.fa
rm GCA_003271115.1_ASM327111v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271125.1_ASM327112v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271125.1_ASM327112v1_genomic.fna.maker.output/GCA_003271125.1_ASM327112v1_genomic.fna_master_datastore_index.log
cp GCA_003271125.1_ASM327112v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271125.1_ASM327112v1_genomic.fa
rm GCA_003271125.1_ASM327112v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271135.1_ASM327113v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271135.1_ASM327113v1_genomic.fna.maker.output/GCA_003271135.1_ASM327113v1_genomic.fna_master_datastore_index.log
cp GCA_003271135.1_ASM327113v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271135.1_ASM327113v1_genomic.fa
rm GCA_003271135.1_ASM327113v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271185.1_ASM327118v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271185.1_ASM327118v1_genomic.fna.maker.output/GCA_003271185.1_ASM327118v1_genomic.fna_master_datastore_index.log
cp GCA_003271185.1_ASM327118v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271185.1_ASM327118v1_genomic.fa
rm GCA_003271185.1_ASM327118v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271225.1_ASM327122v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271225.1_ASM327122v1_genomic.fna.maker.output/GCA_003271225.1_ASM327122v1_genomic.fna_master_datastore_index.log
cp GCA_003271225.1_ASM327122v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271225.1_ASM327122v1_genomic.fa
rm GCA_003271225.1_ASM327122v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271245.1_ASM327124v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271245.1_ASM327124v1_genomic.fna.maker.output/GCA_003271245.1_ASM327124v1_genomic.fna_master_datastore_index.log
cp GCA_003271245.1_ASM327124v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271245.1_ASM327124v1_genomic.fa
rm GCA_003271245.1_ASM327124v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271255.1_ASM327125v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271255.1_ASM327125v1_genomic.fna.maker.output/GCA_003271255.1_ASM327125v1_genomic.fna_master_datastore_index.log
cp GCA_003271255.1_ASM327125v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271255.1_ASM327125v1_genomic.fa
rm GCA_003271255.1_ASM327125v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271315.1_ASM327131v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271315.1_ASM327131v1_genomic.fna.maker.output/GCA_003271315.1_ASM327131v1_genomic.fna_master_datastore_index.log
cp GCA_003271315.1_ASM327131v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271315.1_ASM327131v1_genomic.fa
rm GCA_003271315.1_ASM327131v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271355.1_ASM327135v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271355.1_ASM327135v1_genomic.fna.maker.output/GCA_003271355.1_ASM327135v1_genomic.fna_master_datastore_index.log
cp GCA_003271355.1_ASM327135v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271355.1_ASM327135v1_genomic.fa
rm GCA_003271355.1_ASM327135v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271365.1_ASM327136v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271365.1_ASM327136v1_genomic.fna.maker.output/GCA_003271365.1_ASM327136v1_genomic.fna_master_datastore_index.log
cp GCA_003271365.1_ASM327136v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271365.1_ASM327136v1_genomic.fa
rm GCA_003271365.1_ASM327136v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271395.1_ASM327139v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271395.1_ASM327139v1_genomic.fna.maker.output/GCA_003271395.1_ASM327139v1_genomic.fna_master_datastore_index.log
cp GCA_003271395.1_ASM327139v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271395.1_ASM327139v1_genomic.fa
rm GCA_003271395.1_ASM327139v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271425.1_ASM327142v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271425.1_ASM327142v1_genomic.fna.maker.output/GCA_003271425.1_ASM327142v1_genomic.fna_master_datastore_index.log
cp GCA_003271425.1_ASM327142v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271425.1_ASM327142v1_genomic.fa
rm GCA_003271425.1_ASM327142v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271475.1_ASM327147v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271475.1_ASM327147v1_genomic.fna.maker.output/GCA_003271475.1_ASM327147v1_genomic.fna_master_datastore_index.log
cp GCA_003271475.1_ASM327147v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271475.1_ASM327147v1_genomic.fa
rm GCA_003271475.1_ASM327147v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271505.1_ASM327150v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271505.1_ASM327150v1_genomic.fna.maker.output/GCA_003271505.1_ASM327150v1_genomic.fna_master_datastore_index.log
cp GCA_003271505.1_ASM327150v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271505.1_ASM327150v1_genomic.fa
rm GCA_003271505.1_ASM327150v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271535.1_ASM327153v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271535.1_ASM327153v1_genomic.fna.maker.output/GCA_003271535.1_ASM327153v1_genomic.fna_master_datastore_index.log
cp GCA_003271535.1_ASM327153v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271535.1_ASM327153v1_genomic.fa
rm GCA_003271535.1_ASM327153v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271545.1_ASM327154v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271545.1_ASM327154v1_genomic.fna.maker.output/GCA_003271545.1_ASM327154v1_genomic.fna_master_datastore_index.log
cp GCA_003271545.1_ASM327154v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271545.1_ASM327154v1_genomic.fa
rm GCA_003271545.1_ASM327154v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271555.1_ASM327155v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271555.1_ASM327155v1_genomic.fna.maker.output/GCA_003271555.1_ASM327155v1_genomic.fna_master_datastore_index.log
cp GCA_003271555.1_ASM327155v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271555.1_ASM327155v1_genomic.fa
rm GCA_003271555.1_ASM327155v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271565.1_ASM327156v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271565.1_ASM327156v1_genomic.fna.maker.output/GCA_003271565.1_ASM327156v1_genomic.fna_master_datastore_index.log
cp GCA_003271565.1_ASM327156v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271565.1_ASM327156v1_genomic.fa
rm GCA_003271565.1_ASM327156v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271645.1_ASM327164v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271645.1_ASM327164v1_genomic.fna.maker.output/GCA_003271645.1_ASM327164v1_genomic.fna_master_datastore_index.log
cp GCA_003271645.1_ASM327164v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271645.1_ASM327164v1_genomic.fa
rm GCA_003271645.1_ASM327164v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271675.1_ASM327167v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271675.1_ASM327167v1_genomic.fna.maker.output/GCA_003271675.1_ASM327167v1_genomic.fna_master_datastore_index.log
cp GCA_003271675.1_ASM327167v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271675.1_ASM327167v1_genomic.fa
rm GCA_003271675.1_ASM327167v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271695.1_ASM327169v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271695.1_ASM327169v1_genomic.fna.maker.output/GCA_003271695.1_ASM327169v1_genomic.fna_master_datastore_index.log
cp GCA_003271695.1_ASM327169v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271695.1_ASM327169v1_genomic.fa
rm GCA_003271695.1_ASM327169v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271705.1_ASM327170v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271705.1_ASM327170v1_genomic.fna.maker.output/GCA_003271705.1_ASM327170v1_genomic.fna_master_datastore_index.log
cp GCA_003271705.1_ASM327170v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271705.1_ASM327170v1_genomic.fa
rm GCA_003271705.1_ASM327170v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271775.1_ASM327177v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271775.1_ASM327177v1_genomic.fna.maker.output/GCA_003271775.1_ASM327177v1_genomic.fna_master_datastore_index.log
cp GCA_003271775.1_ASM327177v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271775.1_ASM327177v1_genomic.fa
rm GCA_003271775.1_ASM327177v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271805.1_ASM327180v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271805.1_ASM327180v1_genomic.fna.maker.output/GCA_003271805.1_ASM327180v1_genomic.fna_master_datastore_index.log
cp GCA_003271805.1_ASM327180v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271805.1_ASM327180v1_genomic.fa
rm GCA_003271805.1_ASM327180v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271815.1_ASM327181v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271815.1_ASM327181v1_genomic.fna.maker.output/GCA_003271815.1_ASM327181v1_genomic.fna_master_datastore_index.log
cp GCA_003271815.1_ASM327181v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271815.1_ASM327181v1_genomic.fa
rm GCA_003271815.1_ASM327181v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271825.1_ASM327182v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271825.1_ASM327182v1_genomic.fna.maker.output/GCA_003271825.1_ASM327182v1_genomic.fna_master_datastore_index.log
cp GCA_003271825.1_ASM327182v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271825.1_ASM327182v1_genomic.fa
rm GCA_003271825.1_ASM327182v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271835.1_ASM327183v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271835.1_ASM327183v1_genomic.fna.maker.output/GCA_003271835.1_ASM327183v1_genomic.fna_master_datastore_index.log
cp GCA_003271835.1_ASM327183v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271835.1_ASM327183v1_genomic.fa
rm GCA_003271835.1_ASM327183v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271925.1_ASM327192v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271925.1_ASM327192v1_genomic.fna.maker.output/GCA_003271925.1_ASM327192v1_genomic.fna_master_datastore_index.log
cp GCA_003271925.1_ASM327192v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271925.1_ASM327192v1_genomic.fa
rm GCA_003271925.1_ASM327192v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271935.1_ASM327193v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271935.1_ASM327193v1_genomic.fna.maker.output/GCA_003271935.1_ASM327193v1_genomic.fna_master_datastore_index.log
cp GCA_003271935.1_ASM327193v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271935.1_ASM327193v1_genomic.fa
rm GCA_003271935.1_ASM327193v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271945.1_ASM327194v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271945.1_ASM327194v1_genomic.fna.maker.output/GCA_003271945.1_ASM327194v1_genomic.fna_master_datastore_index.log
cp GCA_003271945.1_ASM327194v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271945.1_ASM327194v1_genomic.fa
rm GCA_003271945.1_ASM327194v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003271955.1_ASM327195v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003271955.1_ASM327195v1_genomic.fna.maker.output/GCA_003271955.1_ASM327195v1_genomic.fna_master_datastore_index.log
cp GCA_003271955.1_ASM327195v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003271955.1_ASM327195v1_genomic.fa
rm GCA_003271955.1_ASM327195v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272045.1_ASM327204v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272045.1_ASM327204v1_genomic.fna.maker.output/GCA_003272045.1_ASM327204v1_genomic.fna_master_datastore_index.log
cp GCA_003272045.1_ASM327204v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272045.1_ASM327204v1_genomic.fa
rm GCA_003272045.1_ASM327204v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272055.1_ASM327205v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272055.1_ASM327205v1_genomic.fna.maker.output/GCA_003272055.1_ASM327205v1_genomic.fna_master_datastore_index.log
cp GCA_003272055.1_ASM327205v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272055.1_ASM327205v1_genomic.fa
rm GCA_003272055.1_ASM327205v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272065.1_ASM327206v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272065.1_ASM327206v1_genomic.fna.maker.output/GCA_003272065.1_ASM327206v1_genomic.fna_master_datastore_index.log
cp GCA_003272065.1_ASM327206v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272065.1_ASM327206v1_genomic.fa
rm GCA_003272065.1_ASM327206v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272075.1_ASM327207v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272075.1_ASM327207v1_genomic.fna.maker.output/GCA_003272075.1_ASM327207v1_genomic.fna_master_datastore_index.log
cp GCA_003272075.1_ASM327207v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272075.1_ASM327207v1_genomic.fa
rm GCA_003272075.1_ASM327207v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272165.1_ASM327216v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272165.1_ASM327216v1_genomic.fna.maker.output/GCA_003272165.1_ASM327216v1_genomic.fna_master_datastore_index.log
cp GCA_003272165.1_ASM327216v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272165.1_ASM327216v1_genomic.fa
rm GCA_003272165.1_ASM327216v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272175.1_ASM327217v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272175.1_ASM327217v1_genomic.fna.maker.output/GCA_003272175.1_ASM327217v1_genomic.fna_master_datastore_index.log
cp GCA_003272175.1_ASM327217v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272175.1_ASM327217v1_genomic.fa
rm GCA_003272175.1_ASM327217v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272185.1_ASM327218v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272185.1_ASM327218v1_genomic.fna.maker.output/GCA_003272185.1_ASM327218v1_genomic.fna_master_datastore_index.log
cp GCA_003272185.1_ASM327218v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272185.1_ASM327218v1_genomic.fa
rm GCA_003272185.1_ASM327218v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272195.1_ASM327219v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272195.1_ASM327219v1_genomic.fna.maker.output/GCA_003272195.1_ASM327219v1_genomic.fna_master_datastore_index.log
cp GCA_003272195.1_ASM327219v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272195.1_ASM327219v1_genomic.fa
rm GCA_003272195.1_ASM327219v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272285.1_ASM327228v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272285.1_ASM327228v1_genomic.fna.maker.output/GCA_003272285.1_ASM327228v1_genomic.fna_master_datastore_index.log
cp GCA_003272285.1_ASM327228v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272285.1_ASM327228v1_genomic.fa
rm GCA_003272285.1_ASM327228v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272295.1_ASM327229v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272295.1_ASM327229v1_genomic.fna.maker.output/GCA_003272295.1_ASM327229v1_genomic.fna_master_datastore_index.log
cp GCA_003272295.1_ASM327229v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272295.1_ASM327229v1_genomic.fa
rm GCA_003272295.1_ASM327229v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272305.1_ASM327230v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272305.1_ASM327230v1_genomic.fna.maker.output/GCA_003272305.1_ASM327230v1_genomic.fna_master_datastore_index.log
cp GCA_003272305.1_ASM327230v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272305.1_ASM327230v1_genomic.fa
rm GCA_003272305.1_ASM327230v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272315.1_ASM327231v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272315.1_ASM327231v1_genomic.fna.maker.output/GCA_003272315.1_ASM327231v1_genomic.fna_master_datastore_index.log
cp GCA_003272315.1_ASM327231v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272315.1_ASM327231v1_genomic.fa
rm GCA_003272315.1_ASM327231v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272405.1_ASM327240v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272405.1_ASM327240v1_genomic.fna.maker.output/GCA_003272405.1_ASM327240v1_genomic.fna_master_datastore_index.log
cp GCA_003272405.1_ASM327240v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272405.1_ASM327240v1_genomic.fa
rm GCA_003272405.1_ASM327240v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272415.1_ASM327241v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272415.1_ASM327241v1_genomic.fna.maker.output/GCA_003272415.1_ASM327241v1_genomic.fna_master_datastore_index.log
cp GCA_003272415.1_ASM327241v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272415.1_ASM327241v1_genomic.fa
rm GCA_003272415.1_ASM327241v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272425.1_ASM327242v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272425.1_ASM327242v1_genomic.fna.maker.output/GCA_003272425.1_ASM327242v1_genomic.fna_master_datastore_index.log
cp GCA_003272425.1_ASM327242v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272425.1_ASM327242v1_genomic.fa
rm GCA_003272425.1_ASM327242v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272445.1_ASM327244v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272445.1_ASM327244v1_genomic.fna.maker.output/GCA_003272445.1_ASM327244v1_genomic.fna_master_datastore_index.log
cp GCA_003272445.1_ASM327244v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272445.1_ASM327244v1_genomic.fa
rm GCA_003272445.1_ASM327244v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272525.1_ASM327252v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272525.1_ASM327252v1_genomic.fna.maker.output/GCA_003272525.1_ASM327252v1_genomic.fna_master_datastore_index.log
cp GCA_003272525.1_ASM327252v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272525.1_ASM327252v1_genomic.fa
rm GCA_003272525.1_ASM327252v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272535.1_ASM327253v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272535.1_ASM327253v1_genomic.fna.maker.output/GCA_003272535.1_ASM327253v1_genomic.fna_master_datastore_index.log
cp GCA_003272535.1_ASM327253v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272535.1_ASM327253v1_genomic.fa
rm GCA_003272535.1_ASM327253v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272575.1_ASM327257v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272575.1_ASM327257v1_genomic.fna.maker.output/GCA_003272575.1_ASM327257v1_genomic.fna_master_datastore_index.log
cp GCA_003272575.1_ASM327257v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272575.1_ASM327257v1_genomic.fa
rm GCA_003272575.1_ASM327257v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272595.1_ASM327259v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272595.1_ASM327259v1_genomic.fna.maker.output/GCA_003272595.1_ASM327259v1_genomic.fna_master_datastore_index.log
cp GCA_003272595.1_ASM327259v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272595.1_ASM327259v1_genomic.fa
rm GCA_003272595.1_ASM327259v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272625.1_ASM327262v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272625.1_ASM327262v1_genomic.fna.maker.output/GCA_003272625.1_ASM327262v1_genomic.fna_master_datastore_index.log
cp GCA_003272625.1_ASM327262v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272625.1_ASM327262v1_genomic.fa
rm GCA_003272625.1_ASM327262v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272665.1_ASM327266v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272665.1_ASM327266v1_genomic.fna.maker.output/GCA_003272665.1_ASM327266v1_genomic.fna_master_datastore_index.log
cp GCA_003272665.1_ASM327266v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272665.1_ASM327266v1_genomic.fa
rm GCA_003272665.1_ASM327266v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272705.1_ASM327270v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272705.1_ASM327270v1_genomic.fna.maker.output/GCA_003272705.1_ASM327270v1_genomic.fna_master_datastore_index.log
cp GCA_003272705.1_ASM327270v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272705.1_ASM327270v1_genomic.fa
rm GCA_003272705.1_ASM327270v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272715.1_ASM327271v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272715.1_ASM327271v1_genomic.fna.maker.output/GCA_003272715.1_ASM327271v1_genomic.fna_master_datastore_index.log
cp GCA_003272715.1_ASM327271v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272715.1_ASM327271v1_genomic.fa
rm GCA_003272715.1_ASM327271v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272725.1_ASM327272v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272725.1_ASM327272v1_genomic.fna.maker.output/GCA_003272725.1_ASM327272v1_genomic.fna_master_datastore_index.log
cp GCA_003272725.1_ASM327272v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272725.1_ASM327272v1_genomic.fa
rm GCA_003272725.1_ASM327272v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272755.1_ASM327275v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272755.1_ASM327275v1_genomic.fna.maker.output/GCA_003272755.1_ASM327275v1_genomic.fna_master_datastore_index.log
cp GCA_003272755.1_ASM327275v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272755.1_ASM327275v1_genomic.fa
rm GCA_003272755.1_ASM327275v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272815.1_ASM327281v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272815.1_ASM327281v1_genomic.fna.maker.output/GCA_003272815.1_ASM327281v1_genomic.fna_master_datastore_index.log
cp GCA_003272815.1_ASM327281v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272815.1_ASM327281v1_genomic.fa
rm GCA_003272815.1_ASM327281v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272835.1_ASM327283v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272835.1_ASM327283v1_genomic.fna.maker.output/GCA_003272835.1_ASM327283v1_genomic.fna_master_datastore_index.log
cp GCA_003272835.1_ASM327283v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272835.1_ASM327283v1_genomic.fa
rm GCA_003272835.1_ASM327283v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272855.1_ASM327285v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272855.1_ASM327285v1_genomic.fna.maker.output/GCA_003272855.1_ASM327285v1_genomic.fna_master_datastore_index.log
cp GCA_003272855.1_ASM327285v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272855.1_ASM327285v1_genomic.fa
rm GCA_003272855.1_ASM327285v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272905.1_ASM327290v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272905.1_ASM327290v1_genomic.fna.maker.output/GCA_003272905.1_ASM327290v1_genomic.fna_master_datastore_index.log
cp GCA_003272905.1_ASM327290v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272905.1_ASM327290v1_genomic.fa
rm GCA_003272905.1_ASM327290v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272925.1_ASM327292v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272925.1_ASM327292v1_genomic.fna.maker.output/GCA_003272925.1_ASM327292v1_genomic.fna_master_datastore_index.log
cp GCA_003272925.1_ASM327292v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272925.1_ASM327292v1_genomic.fa
rm GCA_003272925.1_ASM327292v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272945.1_ASM327294v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272945.1_ASM327294v1_genomic.fna.maker.output/GCA_003272945.1_ASM327294v1_genomic.fna_master_datastore_index.log
cp GCA_003272945.1_ASM327294v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272945.1_ASM327294v1_genomic.fa
rm GCA_003272945.1_ASM327294v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003272975.1_ASM327297v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003272975.1_ASM327297v1_genomic.fna.maker.output/GCA_003272975.1_ASM327297v1_genomic.fna_master_datastore_index.log
cp GCA_003272975.1_ASM327297v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003272975.1_ASM327297v1_genomic.fa
rm GCA_003272975.1_ASM327297v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273025.1_ASM327302v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273025.1_ASM327302v1_genomic.fna.maker.output/GCA_003273025.1_ASM327302v1_genomic.fna_master_datastore_index.log
cp GCA_003273025.1_ASM327302v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273025.1_ASM327302v1_genomic.fa
rm GCA_003273025.1_ASM327302v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273035.1_ASM327303v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273035.1_ASM327303v1_genomic.fna.maker.output/GCA_003273035.1_ASM327303v1_genomic.fna_master_datastore_index.log
cp GCA_003273035.1_ASM327303v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273035.1_ASM327303v1_genomic.fa
rm GCA_003273035.1_ASM327303v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273045.1_ASM327304v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273045.1_ASM327304v1_genomic.fna.maker.output/GCA_003273045.1_ASM327304v1_genomic.fna_master_datastore_index.log
cp GCA_003273045.1_ASM327304v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273045.1_ASM327304v1_genomic.fa
rm GCA_003273045.1_ASM327304v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273075.1_ASM327307v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273075.1_ASM327307v1_genomic.fna.maker.output/GCA_003273075.1_ASM327307v1_genomic.fna_master_datastore_index.log
cp GCA_003273075.1_ASM327307v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273075.1_ASM327307v1_genomic.fa
rm GCA_003273075.1_ASM327307v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273145.1_ASM327314v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273145.1_ASM327314v1_genomic.fna.maker.output/GCA_003273145.1_ASM327314v1_genomic.fna_master_datastore_index.log
cp GCA_003273145.1_ASM327314v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273145.1_ASM327314v1_genomic.fa
rm GCA_003273145.1_ASM327314v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273155.1_ASM327315v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273155.1_ASM327315v1_genomic.fna.maker.output/GCA_003273155.1_ASM327315v1_genomic.fna_master_datastore_index.log
cp GCA_003273155.1_ASM327315v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273155.1_ASM327315v1_genomic.fa
rm GCA_003273155.1_ASM327315v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273165.1_ASM327316v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273165.1_ASM327316v1_genomic.fna.maker.output/GCA_003273165.1_ASM327316v1_genomic.fna_master_datastore_index.log
cp GCA_003273165.1_ASM327316v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273165.1_ASM327316v1_genomic.fa
rm GCA_003273165.1_ASM327316v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273175.1_ASM327317v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273175.1_ASM327317v1_genomic.fna.maker.output/GCA_003273175.1_ASM327317v1_genomic.fna_master_datastore_index.log
cp GCA_003273175.1_ASM327317v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273175.1_ASM327317v1_genomic.fa
rm GCA_003273175.1_ASM327317v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273265.1_ASM327326v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273265.1_ASM327326v1_genomic.fna.maker.output/GCA_003273265.1_ASM327326v1_genomic.fna_master_datastore_index.log
cp GCA_003273265.1_ASM327326v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273265.1_ASM327326v1_genomic.fa
rm GCA_003273265.1_ASM327326v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273275.1_ASM327327v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273275.1_ASM327327v1_genomic.fna.maker.output/GCA_003273275.1_ASM327327v1_genomic.fna_master_datastore_index.log
cp GCA_003273275.1_ASM327327v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273275.1_ASM327327v1_genomic.fa
rm GCA_003273275.1_ASM327327v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273285.1_ASM327328v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273285.1_ASM327328v1_genomic.fna.maker.output/GCA_003273285.1_ASM327328v1_genomic.fna_master_datastore_index.log
cp GCA_003273285.1_ASM327328v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273285.1_ASM327328v1_genomic.fa
rm GCA_003273285.1_ASM327328v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273295.1_ASM327329v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273295.1_ASM327329v1_genomic.fna.maker.output/GCA_003273295.1_ASM327329v1_genomic.fna_master_datastore_index.log
cp GCA_003273295.1_ASM327329v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273295.1_ASM327329v1_genomic.fa
rm GCA_003273295.1_ASM327329v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273385.1_ASM327338v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273385.1_ASM327338v1_genomic.fna.maker.output/GCA_003273385.1_ASM327338v1_genomic.fna_master_datastore_index.log
cp GCA_003273385.1_ASM327338v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273385.1_ASM327338v1_genomic.fa
rm GCA_003273385.1_ASM327338v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273395.1_ASM327339v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273395.1_ASM327339v1_genomic.fna.maker.output/GCA_003273395.1_ASM327339v1_genomic.fna_master_datastore_index.log
cp GCA_003273395.1_ASM327339v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273395.1_ASM327339v1_genomic.fa
rm GCA_003273395.1_ASM327339v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273405.1_ASM327340v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273405.1_ASM327340v1_genomic.fna.maker.output/GCA_003273405.1_ASM327340v1_genomic.fna_master_datastore_index.log
cp GCA_003273405.1_ASM327340v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273405.1_ASM327340v1_genomic.fa
rm GCA_003273405.1_ASM327340v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273415.1_ASM327341v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273415.1_ASM327341v1_genomic.fna.maker.output/GCA_003273415.1_ASM327341v1_genomic.fna_master_datastore_index.log
cp GCA_003273415.1_ASM327341v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273415.1_ASM327341v1_genomic.fa
rm GCA_003273415.1_ASM327341v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273495.1_ASM327349v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273495.1_ASM327349v1_genomic.fna.maker.output/GCA_003273495.1_ASM327349v1_genomic.fna_master_datastore_index.log
cp GCA_003273495.1_ASM327349v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273495.1_ASM327349v1_genomic.fa
rm GCA_003273495.1_ASM327349v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273535.1_ASM327353v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273535.1_ASM327353v1_genomic.fna.maker.output/GCA_003273535.1_ASM327353v1_genomic.fna_master_datastore_index.log
cp GCA_003273535.1_ASM327353v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273535.1_ASM327353v1_genomic.fa
rm GCA_003273535.1_ASM327353v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273555.1_ASM327355v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273555.1_ASM327355v1_genomic.fna.maker.output/GCA_003273555.1_ASM327355v1_genomic.fna_master_datastore_index.log
cp GCA_003273555.1_ASM327355v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273555.1_ASM327355v1_genomic.fa
rm GCA_003273555.1_ASM327355v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273575.1_ASM327357v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273575.1_ASM327357v1_genomic.fna.maker.output/GCA_003273575.1_ASM327357v1_genomic.fna_master_datastore_index.log
cp GCA_003273575.1_ASM327357v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273575.1_ASM327357v1_genomic.fa
rm GCA_003273575.1_ASM327357v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273585.1_ASM327358v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273585.1_ASM327358v1_genomic.fna.maker.output/GCA_003273585.1_ASM327358v1_genomic.fna_master_datastore_index.log
cp GCA_003273585.1_ASM327358v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273585.1_ASM327358v1_genomic.fa
rm GCA_003273585.1_ASM327358v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273625.1_ASM327362v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273625.1_ASM327362v1_genomic.fna.maker.output/GCA_003273625.1_ASM327362v1_genomic.fna_master_datastore_index.log
cp GCA_003273625.1_ASM327362v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273625.1_ASM327362v1_genomic.fa
rm GCA_003273625.1_ASM327362v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273685.1_ASM327368v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273685.1_ASM327368v1_genomic.fna.maker.output/GCA_003273685.1_ASM327368v1_genomic.fna_master_datastore_index.log
cp GCA_003273685.1_ASM327368v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273685.1_ASM327368v1_genomic.fa
rm GCA_003273685.1_ASM327368v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273695.1_ASM327369v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273695.1_ASM327369v1_genomic.fna.maker.output/GCA_003273695.1_ASM327369v1_genomic.fna_master_datastore_index.log
cp GCA_003273695.1_ASM327369v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273695.1_ASM327369v1_genomic.fa
rm GCA_003273695.1_ASM327369v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273705.1_ASM327370v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273705.1_ASM327370v1_genomic.fna.maker.output/GCA_003273705.1_ASM327370v1_genomic.fna_master_datastore_index.log
cp GCA_003273705.1_ASM327370v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273705.1_ASM327370v1_genomic.fa
rm GCA_003273705.1_ASM327370v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273725.1_ASM327372v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273725.1_ASM327372v1_genomic.fna.maker.output/GCA_003273725.1_ASM327372v1_genomic.fna_master_datastore_index.log
cp GCA_003273725.1_ASM327372v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273725.1_ASM327372v1_genomic.fa
rm GCA_003273725.1_ASM327372v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273805.1_ASM327380v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273805.1_ASM327380v1_genomic.fna.maker.output/GCA_003273805.1_ASM327380v1_genomic.fna_master_datastore_index.log
cp GCA_003273805.1_ASM327380v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273805.1_ASM327380v1_genomic.fa
rm GCA_003273805.1_ASM327380v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273815.1_ASM327381v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273815.1_ASM327381v1_genomic.fna.maker.output/GCA_003273815.1_ASM327381v1_genomic.fna_master_datastore_index.log
cp GCA_003273815.1_ASM327381v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273815.1_ASM327381v1_genomic.fa
rm GCA_003273815.1_ASM327381v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273825.1_ASM327382v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273825.1_ASM327382v1_genomic.fna.maker.output/GCA_003273825.1_ASM327382v1_genomic.fna_master_datastore_index.log
cp GCA_003273825.1_ASM327382v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273825.1_ASM327382v1_genomic.fa
rm GCA_003273825.1_ASM327382v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273835.1_ASM327383v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273835.1_ASM327383v1_genomic.fna.maker.output/GCA_003273835.1_ASM327383v1_genomic.fna_master_datastore_index.log
cp GCA_003273835.1_ASM327383v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273835.1_ASM327383v1_genomic.fa
rm GCA_003273835.1_ASM327383v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273925.1_ASM327392v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273925.1_ASM327392v1_genomic.fna.maker.output/GCA_003273925.1_ASM327392v1_genomic.fna_master_datastore_index.log
cp GCA_003273925.1_ASM327392v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273925.1_ASM327392v1_genomic.fa
rm GCA_003273925.1_ASM327392v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273935.1_ASM327393v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273935.1_ASM327393v1_genomic.fna.maker.output/GCA_003273935.1_ASM327393v1_genomic.fna_master_datastore_index.log
cp GCA_003273935.1_ASM327393v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273935.1_ASM327393v1_genomic.fa
rm GCA_003273935.1_ASM327393v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273955.1_ASM327395v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273955.1_ASM327395v1_genomic.fna.maker.output/GCA_003273955.1_ASM327395v1_genomic.fna_master_datastore_index.log
cp GCA_003273955.1_ASM327395v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273955.1_ASM327395v1_genomic.fa
rm GCA_003273955.1_ASM327395v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003273985.1_ASM327398v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003273985.1_ASM327398v1_genomic.fna.maker.output/GCA_003273985.1_ASM327398v1_genomic.fna_master_datastore_index.log
cp GCA_003273985.1_ASM327398v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003273985.1_ASM327398v1_genomic.fa
rm GCA_003273985.1_ASM327398v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274025.1_ASM327402v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274025.1_ASM327402v1_genomic.fna.maker.output/GCA_003274025.1_ASM327402v1_genomic.fna_master_datastore_index.log
cp GCA_003274025.1_ASM327402v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274025.1_ASM327402v1_genomic.fa
rm GCA_003274025.1_ASM327402v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274055.1_ASM327405v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274055.1_ASM327405v1_genomic.fna.maker.output/GCA_003274055.1_ASM327405v1_genomic.fna_master_datastore_index.log
cp GCA_003274055.1_ASM327405v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274055.1_ASM327405v1_genomic.fa
rm GCA_003274055.1_ASM327405v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274085.1_ASM327408v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274085.1_ASM327408v1_genomic.fna.maker.output/GCA_003274085.1_ASM327408v1_genomic.fna_master_datastore_index.log
cp GCA_003274085.1_ASM327408v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274085.1_ASM327408v1_genomic.fa
rm GCA_003274085.1_ASM327408v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274095.1_ASM327409v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274095.1_ASM327409v1_genomic.fna.maker.output/GCA_003274095.1_ASM327409v1_genomic.fna_master_datastore_index.log
cp GCA_003274095.1_ASM327409v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274095.1_ASM327409v1_genomic.fa
rm GCA_003274095.1_ASM327409v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274125.1_ASM327412v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274125.1_ASM327412v1_genomic.fna.maker.output/GCA_003274125.1_ASM327412v1_genomic.fna_master_datastore_index.log
cp GCA_003274125.1_ASM327412v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274125.1_ASM327412v1_genomic.fa
rm GCA_003274125.1_ASM327412v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274195.1_ASM327419v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274195.1_ASM327419v1_genomic.fna.maker.output/GCA_003274195.1_ASM327419v1_genomic.fna_master_datastore_index.log
cp GCA_003274195.1_ASM327419v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274195.1_ASM327419v1_genomic.fa
rm GCA_003274195.1_ASM327419v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274205.1_ASM327420v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274205.1_ASM327420v1_genomic.fna.maker.output/GCA_003274205.1_ASM327420v1_genomic.fna_master_datastore_index.log
cp GCA_003274205.1_ASM327420v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274205.1_ASM327420v1_genomic.fa
rm GCA_003274205.1_ASM327420v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274215.1_ASM327421v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274215.1_ASM327421v1_genomic.fna.maker.output/GCA_003274215.1_ASM327421v1_genomic.fna_master_datastore_index.log
cp GCA_003274215.1_ASM327421v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274215.1_ASM327421v1_genomic.fa
rm GCA_003274215.1_ASM327421v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274225.1_ASM327422v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274225.1_ASM327422v1_genomic.fna.maker.output/GCA_003274225.1_ASM327422v1_genomic.fna_master_datastore_index.log
cp GCA_003274225.1_ASM327422v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274225.1_ASM327422v1_genomic.fa
rm GCA_003274225.1_ASM327422v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274315.1_ASM327431v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274315.1_ASM327431v1_genomic.fna.maker.output/GCA_003274315.1_ASM327431v1_genomic.fna_master_datastore_index.log
cp GCA_003274315.1_ASM327431v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274315.1_ASM327431v1_genomic.fa
rm GCA_003274315.1_ASM327431v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274325.1_ASM327432v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274325.1_ASM327432v1_genomic.fna.maker.output/GCA_003274325.1_ASM327432v1_genomic.fna_master_datastore_index.log
cp GCA_003274325.1_ASM327432v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274325.1_ASM327432v1_genomic.fa
rm GCA_003274325.1_ASM327432v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274345.1_ASM327434v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274345.1_ASM327434v1_genomic.fna.maker.output/GCA_003274345.1_ASM327434v1_genomic.fna_master_datastore_index.log
cp GCA_003274345.1_ASM327434v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274345.1_ASM327434v1_genomic.fa
rm GCA_003274345.1_ASM327434v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274365.1_ASM327436v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274365.1_ASM327436v1_genomic.fna.maker.output/GCA_003274365.1_ASM327436v1_genomic.fna_master_datastore_index.log
cp GCA_003274365.1_ASM327436v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274365.1_ASM327436v1_genomic.fa
rm GCA_003274365.1_ASM327436v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274425.1_ASM327442v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274425.1_ASM327442v1_genomic.fna.maker.output/GCA_003274425.1_ASM327442v1_genomic.fna_master_datastore_index.log
cp GCA_003274425.1_ASM327442v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274425.1_ASM327442v1_genomic.fa
rm GCA_003274425.1_ASM327442v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274465.1_ASM327446v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274465.1_ASM327446v1_genomic.fna.maker.output/GCA_003274465.1_ASM327446v1_genomic.fna_master_datastore_index.log
cp GCA_003274465.1_ASM327446v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274465.1_ASM327446v1_genomic.fa
rm GCA_003274465.1_ASM327446v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274475.1_ASM327447v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274475.1_ASM327447v1_genomic.fna.maker.output/GCA_003274475.1_ASM327447v1_genomic.fna_master_datastore_index.log
cp GCA_003274475.1_ASM327447v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274475.1_ASM327447v1_genomic.fa
rm GCA_003274475.1_ASM327447v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274485.1_ASM327448v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274485.1_ASM327448v1_genomic.fna.maker.output/GCA_003274485.1_ASM327448v1_genomic.fna_master_datastore_index.log
cp GCA_003274485.1_ASM327448v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274485.1_ASM327448v1_genomic.fa
rm GCA_003274485.1_ASM327448v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274495.1_ASM327449v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274495.1_ASM327449v1_genomic.fna.maker.output/GCA_003274495.1_ASM327449v1_genomic.fna_master_datastore_index.log
cp GCA_003274495.1_ASM327449v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274495.1_ASM327449v1_genomic.fa
rm GCA_003274495.1_ASM327449v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274585.1_ASM327458v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274585.1_ASM327458v1_genomic.fna.maker.output/GCA_003274585.1_ASM327458v1_genomic.fna_master_datastore_index.log
cp GCA_003274585.1_ASM327458v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274585.1_ASM327458v1_genomic.fa
rm GCA_003274585.1_ASM327458v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274595.1_ASM327459v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274595.1_ASM327459v1_genomic.fna.maker.output/GCA_003274595.1_ASM327459v1_genomic.fna_master_datastore_index.log
cp GCA_003274595.1_ASM327459v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274595.1_ASM327459v1_genomic.fa
rm GCA_003274595.1_ASM327459v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274605.1_ASM327460v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274605.1_ASM327460v1_genomic.fna.maker.output/GCA_003274605.1_ASM327460v1_genomic.fna_master_datastore_index.log
cp GCA_003274605.1_ASM327460v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274605.1_ASM327460v1_genomic.fa
rm GCA_003274605.1_ASM327460v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274615.1_ASM327461v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274615.1_ASM327461v1_genomic.fna.maker.output/GCA_003274615.1_ASM327461v1_genomic.fna_master_datastore_index.log
cp GCA_003274615.1_ASM327461v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274615.1_ASM327461v1_genomic.fa
rm GCA_003274615.1_ASM327461v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274705.1_ASM327470v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274705.1_ASM327470v1_genomic.fna.maker.output/GCA_003274705.1_ASM327470v1_genomic.fna_master_datastore_index.log
cp GCA_003274705.1_ASM327470v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274705.1_ASM327470v1_genomic.fa
rm GCA_003274705.1_ASM327470v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274715.1_ASM327471v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274715.1_ASM327471v1_genomic.fna.maker.output/GCA_003274715.1_ASM327471v1_genomic.fna_master_datastore_index.log
cp GCA_003274715.1_ASM327471v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274715.1_ASM327471v1_genomic.fa
rm GCA_003274715.1_ASM327471v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274725.1_ASM327472v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274725.1_ASM327472v1_genomic.fna.maker.output/GCA_003274725.1_ASM327472v1_genomic.fna_master_datastore_index.log
cp GCA_003274725.1_ASM327472v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274725.1_ASM327472v1_genomic.fa
rm GCA_003274725.1_ASM327472v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274735.1_ASM327473v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274735.1_ASM327473v1_genomic.fna.maker.output/GCA_003274735.1_ASM327473v1_genomic.fna_master_datastore_index.log
cp GCA_003274735.1_ASM327473v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274735.1_ASM327473v1_genomic.fa
rm GCA_003274735.1_ASM327473v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274825.1_ASM327482v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274825.1_ASM327482v1_genomic.fna.maker.output/GCA_003274825.1_ASM327482v1_genomic.fna_master_datastore_index.log
cp GCA_003274825.1_ASM327482v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274825.1_ASM327482v1_genomic.fa
rm GCA_003274825.1_ASM327482v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274845.1_ASM327484v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274845.1_ASM327484v1_genomic.fna.maker.output/GCA_003274845.1_ASM327484v1_genomic.fna_master_datastore_index.log
cp GCA_003274845.1_ASM327484v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274845.1_ASM327484v1_genomic.fa
rm GCA_003274845.1_ASM327484v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274855.1_ASM327485v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274855.1_ASM327485v1_genomic.fna.maker.output/GCA_003274855.1_ASM327485v1_genomic.fna_master_datastore_index.log
cp GCA_003274855.1_ASM327485v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274855.1_ASM327485v1_genomic.fa
rm GCA_003274855.1_ASM327485v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274865.1_ASM327486v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274865.1_ASM327486v1_genomic.fna.maker.output/GCA_003274865.1_ASM327486v1_genomic.fna_master_datastore_index.log
cp GCA_003274865.1_ASM327486v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274865.1_ASM327486v1_genomic.fa
rm GCA_003274865.1_ASM327486v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274905.1_ASM327490v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274905.1_ASM327490v1_genomic.fna.maker.output/GCA_003274905.1_ASM327490v1_genomic.fna_master_datastore_index.log
cp GCA_003274905.1_ASM327490v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274905.1_ASM327490v1_genomic.fa
rm GCA_003274905.1_ASM327490v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274975.1_ASM327497v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274975.1_ASM327497v1_genomic.fna.maker.output/GCA_003274975.1_ASM327497v1_genomic.fna_master_datastore_index.log
cp GCA_003274975.1_ASM327497v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274975.1_ASM327497v1_genomic.fa
rm GCA_003274975.1_ASM327497v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274985.1_ASM327498v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274985.1_ASM327498v1_genomic.fna.maker.output/GCA_003274985.1_ASM327498v1_genomic.fna_master_datastore_index.log
cp GCA_003274985.1_ASM327498v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274985.1_ASM327498v1_genomic.fa
rm GCA_003274985.1_ASM327498v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003274995.1_ASM327499v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003274995.1_ASM327499v1_genomic.fna.maker.output/GCA_003274995.1_ASM327499v1_genomic.fna_master_datastore_index.log
cp GCA_003274995.1_ASM327499v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003274995.1_ASM327499v1_genomic.fa
rm GCA_003274995.1_ASM327499v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275005.1_ASM327500v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275005.1_ASM327500v1_genomic.fna.maker.output/GCA_003275005.1_ASM327500v1_genomic.fna_master_datastore_index.log
cp GCA_003275005.1_ASM327500v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275005.1_ASM327500v1_genomic.fa
rm GCA_003275005.1_ASM327500v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275095.1_ASM327509v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275095.1_ASM327509v1_genomic.fna.maker.output/GCA_003275095.1_ASM327509v1_genomic.fna_master_datastore_index.log
cp GCA_003275095.1_ASM327509v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275095.1_ASM327509v1_genomic.fa
rm GCA_003275095.1_ASM327509v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275105.1_ASM327510v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275105.1_ASM327510v1_genomic.fna.maker.output/GCA_003275105.1_ASM327510v1_genomic.fna_master_datastore_index.log
cp GCA_003275105.1_ASM327510v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275105.1_ASM327510v1_genomic.fa
rm GCA_003275105.1_ASM327510v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275115.1_ASM327511v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275115.1_ASM327511v1_genomic.fna.maker.output/GCA_003275115.1_ASM327511v1_genomic.fna_master_datastore_index.log
cp GCA_003275115.1_ASM327511v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275115.1_ASM327511v1_genomic.fa
rm GCA_003275115.1_ASM327511v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275125.1_ASM327512v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275125.1_ASM327512v1_genomic.fna.maker.output/GCA_003275125.1_ASM327512v1_genomic.fna_master_datastore_index.log
cp GCA_003275125.1_ASM327512v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275125.1_ASM327512v1_genomic.fa
rm GCA_003275125.1_ASM327512v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275215.1_ASM327521v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275215.1_ASM327521v1_genomic.fna.maker.output/GCA_003275215.1_ASM327521v1_genomic.fna_master_datastore_index.log
cp GCA_003275215.1_ASM327521v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275215.1_ASM327521v1_genomic.fa
rm GCA_003275215.1_ASM327521v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275225.1_ASM327522v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275225.1_ASM327522v1_genomic.fna.maker.output/GCA_003275225.1_ASM327522v1_genomic.fna_master_datastore_index.log
cp GCA_003275225.1_ASM327522v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275225.1_ASM327522v1_genomic.fa
rm GCA_003275225.1_ASM327522v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275235.1_ASM327523v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275235.1_ASM327523v1_genomic.fna.maker.output/GCA_003275235.1_ASM327523v1_genomic.fna_master_datastore_index.log
cp GCA_003275235.1_ASM327523v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275235.1_ASM327523v1_genomic.fa
rm GCA_003275235.1_ASM327523v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275245.1_ASM327524v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275245.1_ASM327524v1_genomic.fna.maker.output/GCA_003275245.1_ASM327524v1_genomic.fna_master_datastore_index.log
cp GCA_003275245.1_ASM327524v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275245.1_ASM327524v1_genomic.fa
rm GCA_003275245.1_ASM327524v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275335.1_ASM327533v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275335.1_ASM327533v1_genomic.fna.maker.output/GCA_003275335.1_ASM327533v1_genomic.fna_master_datastore_index.log
cp GCA_003275335.1_ASM327533v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275335.1_ASM327533v1_genomic.fa
rm GCA_003275335.1_ASM327533v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275345.1_ASM327534v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275345.1_ASM327534v1_genomic.fna.maker.output/GCA_003275345.1_ASM327534v1_genomic.fna_master_datastore_index.log
cp GCA_003275345.1_ASM327534v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275345.1_ASM327534v1_genomic.fa
rm GCA_003275345.1_ASM327534v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275365.1_ASM327536v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275365.1_ASM327536v1_genomic.fna.maker.output/GCA_003275365.1_ASM327536v1_genomic.fna_master_datastore_index.log
cp GCA_003275365.1_ASM327536v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275365.1_ASM327536v1_genomic.fa
rm GCA_003275365.1_ASM327536v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275385.1_ASM327538v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275385.1_ASM327538v1_genomic.fna.maker.output/GCA_003275385.1_ASM327538v1_genomic.fna_master_datastore_index.log
cp GCA_003275385.1_ASM327538v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275385.1_ASM327538v1_genomic.fa
rm GCA_003275385.1_ASM327538v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275455.1_ASM327545v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275455.1_ASM327545v1_genomic.fna.maker.output/GCA_003275455.1_ASM327545v1_genomic.fna_master_datastore_index.log
cp GCA_003275455.1_ASM327545v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275455.1_ASM327545v1_genomic.fa
rm GCA_003275455.1_ASM327545v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275465.1_ASM327546v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275465.1_ASM327546v1_genomic.fna.maker.output/GCA_003275465.1_ASM327546v1_genomic.fna_master_datastore_index.log
cp GCA_003275465.1_ASM327546v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275465.1_ASM327546v1_genomic.fa
rm GCA_003275465.1_ASM327546v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275475.1_ASM327547v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275475.1_ASM327547v1_genomic.fna.maker.output/GCA_003275475.1_ASM327547v1_genomic.fna_master_datastore_index.log
cp GCA_003275475.1_ASM327547v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275475.1_ASM327547v1_genomic.fa
rm GCA_003275475.1_ASM327547v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275505.1_ASM327550v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275505.1_ASM327550v1_genomic.fna.maker.output/GCA_003275505.1_ASM327550v1_genomic.fna_master_datastore_index.log
cp GCA_003275505.1_ASM327550v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275505.1_ASM327550v1_genomic.fa
rm GCA_003275505.1_ASM327550v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275575.1_ASM327557v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275575.1_ASM327557v1_genomic.fna.maker.output/GCA_003275575.1_ASM327557v1_genomic.fna_master_datastore_index.log
cp GCA_003275575.1_ASM327557v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275575.1_ASM327557v1_genomic.fa
rm GCA_003275575.1_ASM327557v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275585.1_ASM327558v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275585.1_ASM327558v1_genomic.fna.maker.output/GCA_003275585.1_ASM327558v1_genomic.fna_master_datastore_index.log
cp GCA_003275585.1_ASM327558v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275585.1_ASM327558v1_genomic.fa
rm GCA_003275585.1_ASM327558v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275605.1_ASM327560v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275605.1_ASM327560v1_genomic.fna.maker.output/GCA_003275605.1_ASM327560v1_genomic.fna_master_datastore_index.log
cp GCA_003275605.1_ASM327560v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275605.1_ASM327560v1_genomic.fa
rm GCA_003275605.1_ASM327560v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275615.1_ASM327561v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275615.1_ASM327561v1_genomic.fna.maker.output/GCA_003275615.1_ASM327561v1_genomic.fna_master_datastore_index.log
cp GCA_003275615.1_ASM327561v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275615.1_ASM327561v1_genomic.fa
rm GCA_003275615.1_ASM327561v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275685.1_ASM327568v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275685.1_ASM327568v1_genomic.fna.maker.output/GCA_003275685.1_ASM327568v1_genomic.fna_master_datastore_index.log
cp GCA_003275685.1_ASM327568v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275685.1_ASM327568v1_genomic.fa
rm GCA_003275685.1_ASM327568v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275715.1_ASM327571v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275715.1_ASM327571v1_genomic.fna.maker.output/GCA_003275715.1_ASM327571v1_genomic.fna_master_datastore_index.log
cp GCA_003275715.1_ASM327571v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275715.1_ASM327571v1_genomic.fa
rm GCA_003275715.1_ASM327571v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275745.1_ASM327574v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275745.1_ASM327574v1_genomic.fna.maker.output/GCA_003275745.1_ASM327574v1_genomic.fna_master_datastore_index.log
cp GCA_003275745.1_ASM327574v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275745.1_ASM327574v1_genomic.fa
rm GCA_003275745.1_ASM327574v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275755.1_ASM327575v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275755.1_ASM327575v1_genomic.fna.maker.output/GCA_003275755.1_ASM327575v1_genomic.fna_master_datastore_index.log
cp GCA_003275755.1_ASM327575v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275755.1_ASM327575v1_genomic.fa
rm GCA_003275755.1_ASM327575v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275785.1_ASM327578v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275785.1_ASM327578v1_genomic.fna.maker.output/GCA_003275785.1_ASM327578v1_genomic.fna_master_datastore_index.log
cp GCA_003275785.1_ASM327578v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275785.1_ASM327578v1_genomic.fa
rm GCA_003275785.1_ASM327578v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275825.1_ASM327582v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275825.1_ASM327582v1_genomic.fna.maker.output/GCA_003275825.1_ASM327582v1_genomic.fna_master_datastore_index.log
cp GCA_003275825.1_ASM327582v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275825.1_ASM327582v1_genomic.fa
rm GCA_003275825.1_ASM327582v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275835.1_ASM327583v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275835.1_ASM327583v1_genomic.fna.maker.output/GCA_003275835.1_ASM327583v1_genomic.fna_master_datastore_index.log
cp GCA_003275835.1_ASM327583v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275835.1_ASM327583v1_genomic.fa
rm GCA_003275835.1_ASM327583v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275875.1_ASM327587v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275875.1_ASM327587v1_genomic.fna.maker.output/GCA_003275875.1_ASM327587v1_genomic.fna_master_datastore_index.log
cp GCA_003275875.1_ASM327587v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275875.1_ASM327587v1_genomic.fa
rm GCA_003275875.1_ASM327587v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275885.1_ASM327588v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275885.1_ASM327588v1_genomic.fna.maker.output/GCA_003275885.1_ASM327588v1_genomic.fna_master_datastore_index.log
cp GCA_003275885.1_ASM327588v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275885.1_ASM327588v1_genomic.fa
rm GCA_003275885.1_ASM327588v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275955.1_ASM327595v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275955.1_ASM327595v1_genomic.fna.maker.output/GCA_003275955.1_ASM327595v1_genomic.fna_master_datastore_index.log
cp GCA_003275955.1_ASM327595v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275955.1_ASM327595v1_genomic.fa
rm GCA_003275955.1_ASM327595v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003275965.1_ASM327596v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003275965.1_ASM327596v1_genomic.fna.maker.output/GCA_003275965.1_ASM327596v1_genomic.fna_master_datastore_index.log
cp GCA_003275965.1_ASM327596v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003275965.1_ASM327596v1_genomic.fa
rm GCA_003275965.1_ASM327596v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276005.1_ASM327600v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276005.1_ASM327600v1_genomic.fna.maker.output/GCA_003276005.1_ASM327600v1_genomic.fna_master_datastore_index.log
cp GCA_003276005.1_ASM327600v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276005.1_ASM327600v1_genomic.fa
rm GCA_003276005.1_ASM327600v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276015.1_ASM327601v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276015.1_ASM327601v1_genomic.fna.maker.output/GCA_003276015.1_ASM327601v1_genomic.fna_master_datastore_index.log
cp GCA_003276015.1_ASM327601v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276015.1_ASM327601v1_genomic.fa
rm GCA_003276015.1_ASM327601v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276075.1_ASM327607v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276075.1_ASM327607v1_genomic.fna.maker.output/GCA_003276075.1_ASM327607v1_genomic.fna_master_datastore_index.log
cp GCA_003276075.1_ASM327607v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276075.1_ASM327607v1_genomic.fa
rm GCA_003276075.1_ASM327607v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276085.1_ASM327608v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276085.1_ASM327608v1_genomic.fna.maker.output/GCA_003276085.1_ASM327608v1_genomic.fna_master_datastore_index.log
cp GCA_003276085.1_ASM327608v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276085.1_ASM327608v1_genomic.fa
rm GCA_003276085.1_ASM327608v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276135.1_ASM327613v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276135.1_ASM327613v1_genomic.fna.maker.output/GCA_003276135.1_ASM327613v1_genomic.fna_master_datastore_index.log
cp GCA_003276135.1_ASM327613v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276135.1_ASM327613v1_genomic.fa
rm GCA_003276135.1_ASM327613v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276145.1_ASM327614v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276145.1_ASM327614v1_genomic.fna.maker.output/GCA_003276145.1_ASM327614v1_genomic.fna_master_datastore_index.log
cp GCA_003276145.1_ASM327614v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276145.1_ASM327614v1_genomic.fa
rm GCA_003276145.1_ASM327614v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276185.1_ASM327618v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276185.1_ASM327618v1_genomic.fna.maker.output/GCA_003276185.1_ASM327618v1_genomic.fna_master_datastore_index.log
cp GCA_003276185.1_ASM327618v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276185.1_ASM327618v1_genomic.fa
rm GCA_003276185.1_ASM327618v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276205.1_ASM327620v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276205.1_ASM327620v1_genomic.fna.maker.output/GCA_003276205.1_ASM327620v1_genomic.fna_master_datastore_index.log
cp GCA_003276205.1_ASM327620v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276205.1_ASM327620v1_genomic.fa
rm GCA_003276205.1_ASM327620v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276235.1_ASM327623v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276235.1_ASM327623v1_genomic.fna.maker.output/GCA_003276235.1_ASM327623v1_genomic.fna_master_datastore_index.log
cp GCA_003276235.1_ASM327623v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276235.1_ASM327623v1_genomic.fa
rm GCA_003276235.1_ASM327623v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276245.1_ASM327624v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276245.1_ASM327624v1_genomic.fna.maker.output/GCA_003276245.1_ASM327624v1_genomic.fna_master_datastore_index.log
cp GCA_003276245.1_ASM327624v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276245.1_ASM327624v1_genomic.fa
rm GCA_003276245.1_ASM327624v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276295.1_ASM327629v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276295.1_ASM327629v1_genomic.fna.maker.output/GCA_003276295.1_ASM327629v1_genomic.fna_master_datastore_index.log
cp GCA_003276295.1_ASM327629v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276295.1_ASM327629v1_genomic.fa
rm GCA_003276295.1_ASM327629v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276325.1_ASM327632v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276325.1_ASM327632v1_genomic.fna.maker.output/GCA_003276325.1_ASM327632v1_genomic.fna_master_datastore_index.log
cp GCA_003276325.1_ASM327632v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276325.1_ASM327632v1_genomic.fa
rm GCA_003276325.1_ASM327632v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276365.1_ASM327636v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276365.1_ASM327636v1_genomic.fna.maker.output/GCA_003276365.1_ASM327636v1_genomic.fna_master_datastore_index.log
cp GCA_003276365.1_ASM327636v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276365.1_ASM327636v1_genomic.fa
rm GCA_003276365.1_ASM327636v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276375.1_ASM327637v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276375.1_ASM327637v1_genomic.fna.maker.output/GCA_003276375.1_ASM327637v1_genomic.fna_master_datastore_index.log
cp GCA_003276375.1_ASM327637v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276375.1_ASM327637v1_genomic.fa
rm GCA_003276375.1_ASM327637v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276395.1_ASM327639v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276395.1_ASM327639v1_genomic.fna.maker.output/GCA_003276395.1_ASM327639v1_genomic.fna_master_datastore_index.log
cp GCA_003276395.1_ASM327639v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276395.1_ASM327639v1_genomic.fa
rm GCA_003276395.1_ASM327639v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276405.1_ASM327640v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276405.1_ASM327640v1_genomic.fna.maker.output/GCA_003276405.1_ASM327640v1_genomic.fna_master_datastore_index.log
cp GCA_003276405.1_ASM327640v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276405.1_ASM327640v1_genomic.fa
rm GCA_003276405.1_ASM327640v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276485.1_ASM327648v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276485.1_ASM327648v1_genomic.fna.maker.output/GCA_003276485.1_ASM327648v1_genomic.fna_master_datastore_index.log
cp GCA_003276485.1_ASM327648v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276485.1_ASM327648v1_genomic.fa
rm GCA_003276485.1_ASM327648v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276495.1_ASM327649v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276495.1_ASM327649v1_genomic.fna.maker.output/GCA_003276495.1_ASM327649v1_genomic.fna_master_datastore_index.log
cp GCA_003276495.1_ASM327649v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276495.1_ASM327649v1_genomic.fa
rm GCA_003276495.1_ASM327649v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276505.1_ASM327650v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276505.1_ASM327650v1_genomic.fna.maker.output/GCA_003276505.1_ASM327650v1_genomic.fna_master_datastore_index.log
cp GCA_003276505.1_ASM327650v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276505.1_ASM327650v1_genomic.fa
rm GCA_003276505.1_ASM327650v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276515.1_ASM327651v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276515.1_ASM327651v1_genomic.fna.maker.output/GCA_003276515.1_ASM327651v1_genomic.fna_master_datastore_index.log
cp GCA_003276515.1_ASM327651v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276515.1_ASM327651v1_genomic.fa
rm GCA_003276515.1_ASM327651v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276605.1_ASM327660v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276605.1_ASM327660v1_genomic.fna.maker.output/GCA_003276605.1_ASM327660v1_genomic.fna_master_datastore_index.log
cp GCA_003276605.1_ASM327660v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276605.1_ASM327660v1_genomic.fa
rm GCA_003276605.1_ASM327660v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276615.1_ASM327661v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276615.1_ASM327661v1_genomic.fna.maker.output/GCA_003276615.1_ASM327661v1_genomic.fna_master_datastore_index.log
cp GCA_003276615.1_ASM327661v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276615.1_ASM327661v1_genomic.fa
rm GCA_003276615.1_ASM327661v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276625.1_ASM327662v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276625.1_ASM327662v1_genomic.fna.maker.output/GCA_003276625.1_ASM327662v1_genomic.fna_master_datastore_index.log
cp GCA_003276625.1_ASM327662v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276625.1_ASM327662v1_genomic.fa
rm GCA_003276625.1_ASM327662v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276635.1_ASM327663v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276635.1_ASM327663v1_genomic.fna.maker.output/GCA_003276635.1_ASM327663v1_genomic.fna_master_datastore_index.log
cp GCA_003276635.1_ASM327663v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276635.1_ASM327663v1_genomic.fa
rm GCA_003276635.1_ASM327663v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276725.1_ASM327672v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276725.1_ASM327672v1_genomic.fna.maker.output/GCA_003276725.1_ASM327672v1_genomic.fna_master_datastore_index.log
cp GCA_003276725.1_ASM327672v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276725.1_ASM327672v1_genomic.fa
rm GCA_003276725.1_ASM327672v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276735.1_ASM327673v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276735.1_ASM327673v1_genomic.fna.maker.output/GCA_003276735.1_ASM327673v1_genomic.fna_master_datastore_index.log
cp GCA_003276735.1_ASM327673v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276735.1_ASM327673v1_genomic.fa
rm GCA_003276735.1_ASM327673v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276745.1_ASM327674v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276745.1_ASM327674v1_genomic.fna.maker.output/GCA_003276745.1_ASM327674v1_genomic.fna_master_datastore_index.log
cp GCA_003276745.1_ASM327674v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276745.1_ASM327674v1_genomic.fa
rm GCA_003276745.1_ASM327674v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276755.1_ASM327675v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276755.1_ASM327675v1_genomic.fna.maker.output/GCA_003276755.1_ASM327675v1_genomic.fna_master_datastore_index.log
cp GCA_003276755.1_ASM327675v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276755.1_ASM327675v1_genomic.fa
rm GCA_003276755.1_ASM327675v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276845.1_ASM327684v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276845.1_ASM327684v1_genomic.fna.maker.output/GCA_003276845.1_ASM327684v1_genomic.fna_master_datastore_index.log
cp GCA_003276845.1_ASM327684v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276845.1_ASM327684v1_genomic.fa
rm GCA_003276845.1_ASM327684v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276855.1_ASM327685v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276855.1_ASM327685v1_genomic.fna.maker.output/GCA_003276855.1_ASM327685v1_genomic.fna_master_datastore_index.log
cp GCA_003276855.1_ASM327685v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276855.1_ASM327685v1_genomic.fa
rm GCA_003276855.1_ASM327685v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276865.1_ASM327686v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276865.1_ASM327686v1_genomic.fna.maker.output/GCA_003276865.1_ASM327686v1_genomic.fna_master_datastore_index.log
cp GCA_003276865.1_ASM327686v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276865.1_ASM327686v1_genomic.fa
rm GCA_003276865.1_ASM327686v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276875.1_ASM327687v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276875.1_ASM327687v1_genomic.fna.maker.output/GCA_003276875.1_ASM327687v1_genomic.fna_master_datastore_index.log
cp GCA_003276875.1_ASM327687v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276875.1_ASM327687v1_genomic.fa
rm GCA_003276875.1_ASM327687v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276965.1_ASM327696v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276965.1_ASM327696v1_genomic.fna.maker.output/GCA_003276965.1_ASM327696v1_genomic.fna_master_datastore_index.log
cp GCA_003276965.1_ASM327696v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276965.1_ASM327696v1_genomic.fa
rm GCA_003276965.1_ASM327696v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276975.1_ASM327697v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276975.1_ASM327697v1_genomic.fna.maker.output/GCA_003276975.1_ASM327697v1_genomic.fna_master_datastore_index.log
cp GCA_003276975.1_ASM327697v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276975.1_ASM327697v1_genomic.fa
rm GCA_003276975.1_ASM327697v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276985.1_ASM327698v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276985.1_ASM327698v1_genomic.fna.maker.output/GCA_003276985.1_ASM327698v1_genomic.fna_master_datastore_index.log
cp GCA_003276985.1_ASM327698v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276985.1_ASM327698v1_genomic.fa
rm GCA_003276985.1_ASM327698v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003276995.1_ASM327699v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003276995.1_ASM327699v1_genomic.fna.maker.output/GCA_003276995.1_ASM327699v1_genomic.fna_master_datastore_index.log
cp GCA_003276995.1_ASM327699v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003276995.1_ASM327699v1_genomic.fa
rm GCA_003276995.1_ASM327699v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277085.1_ASM327708v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277085.1_ASM327708v1_genomic.fna.maker.output/GCA_003277085.1_ASM327708v1_genomic.fna_master_datastore_index.log
cp GCA_003277085.1_ASM327708v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277085.1_ASM327708v1_genomic.fa
rm GCA_003277085.1_ASM327708v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277095.1_ASM327709v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277095.1_ASM327709v1_genomic.fna.maker.output/GCA_003277095.1_ASM327709v1_genomic.fna_master_datastore_index.log
cp GCA_003277095.1_ASM327709v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277095.1_ASM327709v1_genomic.fa
rm GCA_003277095.1_ASM327709v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277105.1_ASM327710v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277105.1_ASM327710v1_genomic.fna.maker.output/GCA_003277105.1_ASM327710v1_genomic.fna_master_datastore_index.log
cp GCA_003277105.1_ASM327710v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277105.1_ASM327710v1_genomic.fa
rm GCA_003277105.1_ASM327710v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277115.1_ASM327711v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277115.1_ASM327711v1_genomic.fna.maker.output/GCA_003277115.1_ASM327711v1_genomic.fna_master_datastore_index.log
cp GCA_003277115.1_ASM327711v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277115.1_ASM327711v1_genomic.fa
rm GCA_003277115.1_ASM327711v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277205.1_ASM327720v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277205.1_ASM327720v1_genomic.fna.maker.output/GCA_003277205.1_ASM327720v1_genomic.fna_master_datastore_index.log
cp GCA_003277205.1_ASM327720v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277205.1_ASM327720v1_genomic.fa
rm GCA_003277205.1_ASM327720v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277235.1_ASM327723v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277235.1_ASM327723v1_genomic.fna.maker.output/GCA_003277235.1_ASM327723v1_genomic.fna_master_datastore_index.log
cp GCA_003277235.1_ASM327723v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277235.1_ASM327723v1_genomic.fa
rm GCA_003277235.1_ASM327723v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277245.1_ASM327724v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277245.1_ASM327724v1_genomic.fna.maker.output/GCA_003277245.1_ASM327724v1_genomic.fna_master_datastore_index.log
cp GCA_003277245.1_ASM327724v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277245.1_ASM327724v1_genomic.fa
rm GCA_003277245.1_ASM327724v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277255.1_ASM327725v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277255.1_ASM327725v1_genomic.fna.maker.output/GCA_003277255.1_ASM327725v1_genomic.fna_master_datastore_index.log
cp GCA_003277255.1_ASM327725v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277255.1_ASM327725v1_genomic.fa
rm GCA_003277255.1_ASM327725v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277295.1_ASM327729v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277295.1_ASM327729v1_genomic.fna.maker.output/GCA_003277295.1_ASM327729v1_genomic.fna_master_datastore_index.log
cp GCA_003277295.1_ASM327729v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277295.1_ASM327729v1_genomic.fa
rm GCA_003277295.1_ASM327729v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277355.1_ASM327735v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277355.1_ASM327735v1_genomic.fna.maker.output/GCA_003277355.1_ASM327735v1_genomic.fna_master_datastore_index.log
cp GCA_003277355.1_ASM327735v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277355.1_ASM327735v1_genomic.fa
rm GCA_003277355.1_ASM327735v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277365.1_ASM327736v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277365.1_ASM327736v1_genomic.fna.maker.output/GCA_003277365.1_ASM327736v1_genomic.fna_master_datastore_index.log
cp GCA_003277365.1_ASM327736v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277365.1_ASM327736v1_genomic.fa
rm GCA_003277365.1_ASM327736v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277375.1_ASM327737v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277375.1_ASM327737v1_genomic.fna.maker.output/GCA_003277375.1_ASM327737v1_genomic.fna_master_datastore_index.log
cp GCA_003277375.1_ASM327737v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277375.1_ASM327737v1_genomic.fa
rm GCA_003277375.1_ASM327737v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277385.1_ASM327738v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277385.1_ASM327738v1_genomic.fna.maker.output/GCA_003277385.1_ASM327738v1_genomic.fna_master_datastore_index.log
cp GCA_003277385.1_ASM327738v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277385.1_ASM327738v1_genomic.fa
rm GCA_003277385.1_ASM327738v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277475.1_ASM327747v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277475.1_ASM327747v1_genomic.fna.maker.output/GCA_003277475.1_ASM327747v1_genomic.fna_master_datastore_index.log
cp GCA_003277475.1_ASM327747v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277475.1_ASM327747v1_genomic.fa
rm GCA_003277475.1_ASM327747v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277485.1_ASM327748v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277485.1_ASM327748v1_genomic.fna.maker.output/GCA_003277485.1_ASM327748v1_genomic.fna_master_datastore_index.log
cp GCA_003277485.1_ASM327748v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277485.1_ASM327748v1_genomic.fa
rm GCA_003277485.1_ASM327748v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277495.1_ASM327749v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277495.1_ASM327749v1_genomic.fna.maker.output/GCA_003277495.1_ASM327749v1_genomic.fna_master_datastore_index.log
cp GCA_003277495.1_ASM327749v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277495.1_ASM327749v1_genomic.fa
rm GCA_003277495.1_ASM327749v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277505.1_ASM327750v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277505.1_ASM327750v1_genomic.fna.maker.output/GCA_003277505.1_ASM327750v1_genomic.fna_master_datastore_index.log
cp GCA_003277505.1_ASM327750v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277505.1_ASM327750v1_genomic.fa
rm GCA_003277505.1_ASM327750v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277595.1_ASM327759v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277595.1_ASM327759v1_genomic.fna.maker.output/GCA_003277595.1_ASM327759v1_genomic.fna_master_datastore_index.log
cp GCA_003277595.1_ASM327759v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277595.1_ASM327759v1_genomic.fa
rm GCA_003277595.1_ASM327759v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277615.1_ASM327761v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277615.1_ASM327761v1_genomic.fna.maker.output/GCA_003277615.1_ASM327761v1_genomic.fna_master_datastore_index.log
cp GCA_003277615.1_ASM327761v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277615.1_ASM327761v1_genomic.fa
rm GCA_003277615.1_ASM327761v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277635.1_ASM327763v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277635.1_ASM327763v1_genomic.fna.maker.output/GCA_003277635.1_ASM327763v1_genomic.fna_master_datastore_index.log
cp GCA_003277635.1_ASM327763v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277635.1_ASM327763v1_genomic.fa
rm GCA_003277635.1_ASM327763v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277645.1_ASM327764v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277645.1_ASM327764v1_genomic.fna.maker.output/GCA_003277645.1_ASM327764v1_genomic.fna_master_datastore_index.log
cp GCA_003277645.1_ASM327764v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277645.1_ASM327764v1_genomic.fa
rm GCA_003277645.1_ASM327764v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277715.1_ASM327771v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277715.1_ASM327771v1_genomic.fna.maker.output/GCA_003277715.1_ASM327771v1_genomic.fna_master_datastore_index.log
cp GCA_003277715.1_ASM327771v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277715.1_ASM327771v1_genomic.fa
rm GCA_003277715.1_ASM327771v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277725.1_ASM327772v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277725.1_ASM327772v1_genomic.fna.maker.output/GCA_003277725.1_ASM327772v1_genomic.fna_master_datastore_index.log
cp GCA_003277725.1_ASM327772v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277725.1_ASM327772v1_genomic.fa
rm GCA_003277725.1_ASM327772v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277735.1_ASM327773v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277735.1_ASM327773v1_genomic.fna.maker.output/GCA_003277735.1_ASM327773v1_genomic.fna_master_datastore_index.log
cp GCA_003277735.1_ASM327773v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277735.1_ASM327773v1_genomic.fa
rm GCA_003277735.1_ASM327773v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277755.1_ASM327775v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277755.1_ASM327775v1_genomic.fna.maker.output/GCA_003277755.1_ASM327775v1_genomic.fna_master_datastore_index.log
cp GCA_003277755.1_ASM327775v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277755.1_ASM327775v1_genomic.fa
rm GCA_003277755.1_ASM327775v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277815.1_ASM327781v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277815.1_ASM327781v1_genomic.fna.maker.output/GCA_003277815.1_ASM327781v1_genomic.fna_master_datastore_index.log
cp GCA_003277815.1_ASM327781v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277815.1_ASM327781v1_genomic.fa
rm GCA_003277815.1_ASM327781v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277865.1_ASM327786v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277865.1_ASM327786v1_genomic.fna.maker.output/GCA_003277865.1_ASM327786v1_genomic.fna_master_datastore_index.log
cp GCA_003277865.1_ASM327786v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277865.1_ASM327786v1_genomic.fa
rm GCA_003277865.1_ASM327786v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003277875.1_ASM327787v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003277875.1_ASM327787v1_genomic.fna.maker.output/GCA_003277875.1_ASM327787v1_genomic.fna_master_datastore_index.log
cp GCA_003277875.1_ASM327787v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003277875.1_ASM327787v1_genomic.fa
rm GCA_003277875.1_ASM327787v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003326515.1_ASM332651v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003326515.1_ASM332651v1_genomic.fna.maker.output/GCA_003326515.1_ASM332651v1_genomic.fna_master_datastore_index.log
cp GCA_003326515.1_ASM332651v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003326515.1_ASM332651v1_genomic.fa
rm GCA_003326515.1_ASM332651v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003326525.1_ASM332652v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003326525.1_ASM332652v1_genomic.fna.maker.output/GCA_003326525.1_ASM332652v1_genomic.fna_master_datastore_index.log
cp GCA_003326525.1_ASM332652v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003326525.1_ASM332652v1_genomic.fa
rm GCA_003326525.1_ASM332652v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003326555.1_ASM332655v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003326555.1_ASM332655v1_genomic.fna.maker.output/GCA_003326555.1_ASM332655v1_genomic.fna_master_datastore_index.log
cp GCA_003326555.1_ASM332655v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003326555.1_ASM332655v1_genomic.fa
rm GCA_003326555.1_ASM332655v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_003326565.1_ASM332656v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_003326565.1_ASM332656v1_genomic.fna.maker.output/GCA_003326565.1_ASM332656v1_genomic.fna_master_datastore_index.log
cp GCA_003326565.1_ASM332656v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_003326565.1_ASM332656v1_genomic.fa
rm GCA_003326565.1_ASM332656v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864775.1_ASM1786477v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864775.1_ASM1786477v1_genomic.fna.maker.output/GCA_017864775.1_ASM1786477v1_genomic.fna_master_datastore_index.log
cp GCA_017864775.1_ASM1786477v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864775.1_ASM1786477v1_genomic.fa
rm GCA_017864775.1_ASM1786477v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864785.1_ASM1786478v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864785.1_ASM1786478v1_genomic.fna.maker.output/GCA_017864785.1_ASM1786478v1_genomic.fna_master_datastore_index.log
cp GCA_017864785.1_ASM1786478v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864785.1_ASM1786478v1_genomic.fa
rm GCA_017864785.1_ASM1786478v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864795.1_ASM1786479v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864795.1_ASM1786479v1_genomic.fna.maker.output/GCA_017864795.1_ASM1786479v1_genomic.fna_master_datastore_index.log
cp GCA_017864795.1_ASM1786479v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864795.1_ASM1786479v1_genomic.fa
rm GCA_017864795.1_ASM1786479v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864805.1_ASM1786480v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864805.1_ASM1786480v1_genomic.fna.maker.output/GCA_017864805.1_ASM1786480v1_genomic.fna_master_datastore_index.log
cp GCA_017864805.1_ASM1786480v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864805.1_ASM1786480v1_genomic.fa
rm GCA_017864805.1_ASM1786480v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864825.1_ASM1786482v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864825.1_ASM1786482v1_genomic.fna.maker.output/GCA_017864825.1_ASM1786482v1_genomic.fna_master_datastore_index.log
cp GCA_017864825.1_ASM1786482v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864825.1_ASM1786482v1_genomic.fa
rm GCA_017864825.1_ASM1786482v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864915.1_ASM1786491v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864915.1_ASM1786491v1_genomic.fna.maker.output/GCA_017864915.1_ASM1786491v1_genomic.fna_master_datastore_index.log
cp GCA_017864915.1_ASM1786491v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864915.1_ASM1786491v1_genomic.fa
rm GCA_017864915.1_ASM1786491v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864925.1_ASM1786492v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864925.1_ASM1786492v1_genomic.fna.maker.output/GCA_017864925.1_ASM1786492v1_genomic.fna_master_datastore_index.log
cp GCA_017864925.1_ASM1786492v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864925.1_ASM1786492v1_genomic.fa
rm GCA_017864925.1_ASM1786492v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864945.1_ASM1786494v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864945.1_ASM1786494v1_genomic.fna.maker.output/GCA_017864945.1_ASM1786494v1_genomic.fna_master_datastore_index.log
cp GCA_017864945.1_ASM1786494v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864945.1_ASM1786494v1_genomic.fa
rm GCA_017864945.1_ASM1786494v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017864975.1_ASM1786497v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017864975.1_ASM1786497v1_genomic.fna.maker.output/GCA_017864975.1_ASM1786497v1_genomic.fna_master_datastore_index.log
cp GCA_017864975.1_ASM1786497v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017864975.1_ASM1786497v1_genomic.fa
rm GCA_017864975.1_ASM1786497v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865125.1_ASM1786512v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865125.1_ASM1786512v1_genomic.fna.maker.output/GCA_017865125.1_ASM1786512v1_genomic.fna_master_datastore_index.log
cp GCA_017865125.1_ASM1786512v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865125.1_ASM1786512v1_genomic.fa
rm GCA_017865125.1_ASM1786512v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865225.1_ASM1786522v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865225.1_ASM1786522v1_genomic.fna.maker.output/GCA_017865225.1_ASM1786522v1_genomic.fna_master_datastore_index.log
cp GCA_017865225.1_ASM1786522v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865225.1_ASM1786522v1_genomic.fa
rm GCA_017865225.1_ASM1786522v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865245.1_ASM1786524v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865245.1_ASM1786524v1_genomic.fna.maker.output/GCA_017865245.1_ASM1786524v1_genomic.fna_master_datastore_index.log
cp GCA_017865245.1_ASM1786524v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865245.1_ASM1786524v1_genomic.fa
rm GCA_017865245.1_ASM1786524v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865275.1_ASM1786527v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865275.1_ASM1786527v1_genomic.fna.maker.output/GCA_017865275.1_ASM1786527v1_genomic.fna_master_datastore_index.log
cp GCA_017865275.1_ASM1786527v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865275.1_ASM1786527v1_genomic.fa
rm GCA_017865275.1_ASM1786527v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865285.1_ASM1786528v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865285.1_ASM1786528v1_genomic.fna.maker.output/GCA_017865285.1_ASM1786528v1_genomic.fna_master_datastore_index.log
cp GCA_017865285.1_ASM1786528v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865285.1_ASM1786528v1_genomic.fa
rm GCA_017865285.1_ASM1786528v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865375.1_ASM1786537v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865375.1_ASM1786537v1_genomic.fna.maker.output/GCA_017865375.1_ASM1786537v1_genomic.fna_master_datastore_index.log
cp GCA_017865375.1_ASM1786537v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865375.1_ASM1786537v1_genomic.fa
rm GCA_017865375.1_ASM1786537v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865415.1_ASM1786541v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865415.1_ASM1786541v1_genomic.fna.maker.output/GCA_017865415.1_ASM1786541v1_genomic.fna_master_datastore_index.log
cp GCA_017865415.1_ASM1786541v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865415.1_ASM1786541v1_genomic.fa
rm GCA_017865415.1_ASM1786541v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865475.1_ASM1786547v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865475.1_ASM1786547v1_genomic.fna.maker.output/GCA_017865475.1_ASM1786547v1_genomic.fna_master_datastore_index.log
cp GCA_017865475.1_ASM1786547v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865475.1_ASM1786547v1_genomic.fa
rm GCA_017865475.1_ASM1786547v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865645.1_ASM1786564v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865645.1_ASM1786564v1_genomic.fna.maker.output/GCA_017865645.1_ASM1786564v1_genomic.fna_master_datastore_index.log
cp GCA_017865645.1_ASM1786564v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865645.1_ASM1786564v1_genomic.fa
rm GCA_017865645.1_ASM1786564v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865655.1_ASM1786565v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865655.1_ASM1786565v1_genomic.fna.maker.output/GCA_017865655.1_ASM1786565v1_genomic.fna_master_datastore_index.log
cp GCA_017865655.1_ASM1786565v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865655.1_ASM1786565v1_genomic.fa
rm GCA_017865655.1_ASM1786565v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865695.1_ASM1786569v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865695.1_ASM1786569v1_genomic.fna.maker.output/GCA_017865695.1_ASM1786569v1_genomic.fna_master_datastore_index.log
cp GCA_017865695.1_ASM1786569v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865695.1_ASM1786569v1_genomic.fa
rm GCA_017865695.1_ASM1786569v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865715.1_ASM1786571v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865715.1_ASM1786571v1_genomic.fna.maker.output/GCA_017865715.1_ASM1786571v1_genomic.fna_master_datastore_index.log
cp GCA_017865715.1_ASM1786571v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865715.1_ASM1786571v1_genomic.fa
rm GCA_017865715.1_ASM1786571v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865745.1_ASM1786574v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865745.1_ASM1786574v1_genomic.fna.maker.output/GCA_017865745.1_ASM1786574v1_genomic.fna_master_datastore_index.log
cp GCA_017865745.1_ASM1786574v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865745.1_ASM1786574v1_genomic.fa
rm GCA_017865745.1_ASM1786574v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865875.1_ASM1786587v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865875.1_ASM1786587v1_genomic.fna.maker.output/GCA_017865875.1_ASM1786587v1_genomic.fna_master_datastore_index.log
cp GCA_017865875.1_ASM1786587v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865875.1_ASM1786587v1_genomic.fa
rm GCA_017865875.1_ASM1786587v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017865995.1_ASM1786599v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017865995.1_ASM1786599v1_genomic.fna.maker.output/GCA_017865995.1_ASM1786599v1_genomic.fna_master_datastore_index.log
cp GCA_017865995.1_ASM1786599v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017865995.1_ASM1786599v1_genomic.fa
rm GCA_017865995.1_ASM1786599v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866045.1_ASM1786604v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866045.1_ASM1786604v1_genomic.fna.maker.output/GCA_017866045.1_ASM1786604v1_genomic.fna_master_datastore_index.log
cp GCA_017866045.1_ASM1786604v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866045.1_ASM1786604v1_genomic.fa
rm GCA_017866045.1_ASM1786604v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866055.1_ASM1786605v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866055.1_ASM1786605v1_genomic.fna.maker.output/GCA_017866055.1_ASM1786605v1_genomic.fna_master_datastore_index.log
cp GCA_017866055.1_ASM1786605v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866055.1_ASM1786605v1_genomic.fa
rm GCA_017866055.1_ASM1786605v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866065.1_ASM1786606v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866065.1_ASM1786606v1_genomic.fna.maker.output/GCA_017866065.1_ASM1786606v1_genomic.fna_master_datastore_index.log
cp GCA_017866065.1_ASM1786606v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866065.1_ASM1786606v1_genomic.fa
rm GCA_017866065.1_ASM1786606v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866075.1_ASM1786607v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866075.1_ASM1786607v1_genomic.fna.maker.output/GCA_017866075.1_ASM1786607v1_genomic.fna_master_datastore_index.log
cp GCA_017866075.1_ASM1786607v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866075.1_ASM1786607v1_genomic.fa
rm GCA_017866075.1_ASM1786607v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866085.1_ASM1786608v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866085.1_ASM1786608v1_genomic.fna.maker.output/GCA_017866085.1_ASM1786608v1_genomic.fna_master_datastore_index.log
cp GCA_017866085.1_ASM1786608v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866085.1_ASM1786608v1_genomic.fa
rm GCA_017866085.1_ASM1786608v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866135.1_ASM1786613v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866135.1_ASM1786613v1_genomic.fna.maker.output/GCA_017866135.1_ASM1786613v1_genomic.fna_master_datastore_index.log
cp GCA_017866135.1_ASM1786613v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866135.1_ASM1786613v1_genomic.fa
rm GCA_017866135.1_ASM1786613v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866155.1_ASM1786615v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866155.1_ASM1786615v1_genomic.fna.maker.output/GCA_017866155.1_ASM1786615v1_genomic.fna_master_datastore_index.log
cp GCA_017866155.1_ASM1786615v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866155.1_ASM1786615v1_genomic.fa
rm GCA_017866155.1_ASM1786615v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866175.1_ASM1786617v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866175.1_ASM1786617v1_genomic.fna.maker.output/GCA_017866175.1_ASM1786617v1_genomic.fna_master_datastore_index.log
cp GCA_017866175.1_ASM1786617v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866175.1_ASM1786617v1_genomic.fa
rm GCA_017866175.1_ASM1786617v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866185.1_ASM1786618v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866185.1_ASM1786618v1_genomic.fna.maker.output/GCA_017866185.1_ASM1786618v1_genomic.fna_master_datastore_index.log
cp GCA_017866185.1_ASM1786618v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866185.1_ASM1786618v1_genomic.fa
rm GCA_017866185.1_ASM1786618v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866255.1_ASM1786625v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866255.1_ASM1786625v1_genomic.fna.maker.output/GCA_017866255.1_ASM1786625v1_genomic.fna_master_datastore_index.log
cp GCA_017866255.1_ASM1786625v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866255.1_ASM1786625v1_genomic.fa
rm GCA_017866255.1_ASM1786625v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866275.1_ASM1786627v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866275.1_ASM1786627v1_genomic.fna.maker.output/GCA_017866275.1_ASM1786627v1_genomic.fna_master_datastore_index.log
cp GCA_017866275.1_ASM1786627v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866275.1_ASM1786627v1_genomic.fa
rm GCA_017866275.1_ASM1786627v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866305.1_ASM1786630v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866305.1_ASM1786630v1_genomic.fna.maker.output/GCA_017866305.1_ASM1786630v1_genomic.fna_master_datastore_index.log
cp GCA_017866305.1_ASM1786630v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866305.1_ASM1786630v1_genomic.fa
rm GCA_017866305.1_ASM1786630v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866315.1_ASM1786631v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866315.1_ASM1786631v1_genomic.fna.maker.output/GCA_017866315.1_ASM1786631v1_genomic.fna_master_datastore_index.log
cp GCA_017866315.1_ASM1786631v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866315.1_ASM1786631v1_genomic.fa
rm GCA_017866315.1_ASM1786631v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866385.1_ASM1786638v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866385.1_ASM1786638v1_genomic.fna.maker.output/GCA_017866385.1_ASM1786638v1_genomic.fna_master_datastore_index.log
cp GCA_017866385.1_ASM1786638v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866385.1_ASM1786638v1_genomic.fa
rm GCA_017866385.1_ASM1786638v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866405.1_ASM1786640v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866405.1_ASM1786640v1_genomic.fna.maker.output/GCA_017866405.1_ASM1786640v1_genomic.fna_master_datastore_index.log
cp GCA_017866405.1_ASM1786640v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866405.1_ASM1786640v1_genomic.fa
rm GCA_017866405.1_ASM1786640v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866455.1_ASM1786645v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866455.1_ASM1786645v1_genomic.fna.maker.output/GCA_017866455.1_ASM1786645v1_genomic.fna_master_datastore_index.log
cp GCA_017866455.1_ASM1786645v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866455.1_ASM1786645v1_genomic.fa
rm GCA_017866455.1_ASM1786645v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866465.1_ASM1786646v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866465.1_ASM1786646v1_genomic.fna.maker.output/GCA_017866465.1_ASM1786646v1_genomic.fna_master_datastore_index.log
cp GCA_017866465.1_ASM1786646v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866465.1_ASM1786646v1_genomic.fa
rm GCA_017866465.1_ASM1786646v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866475.1_ASM1786647v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866475.1_ASM1786647v1_genomic.fna.maker.output/GCA_017866475.1_ASM1786647v1_genomic.fna_master_datastore_index.log
cp GCA_017866475.1_ASM1786647v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866475.1_ASM1786647v1_genomic.fa
rm GCA_017866475.1_ASM1786647v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866485.1_ASM1786648v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866485.1_ASM1786648v1_genomic.fna.maker.output/GCA_017866485.1_ASM1786648v1_genomic.fna_master_datastore_index.log
cp GCA_017866485.1_ASM1786648v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866485.1_ASM1786648v1_genomic.fa
rm GCA_017866485.1_ASM1786648v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866495.1_ASM1786649v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866495.1_ASM1786649v1_genomic.fna.maker.output/GCA_017866495.1_ASM1786649v1_genomic.fna_master_datastore_index.log
cp GCA_017866495.1_ASM1786649v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866495.1_ASM1786649v1_genomic.fa
rm GCA_017866495.1_ASM1786649v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866565.1_ASM1786656v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866565.1_ASM1786656v1_genomic.fna.maker.output/GCA_017866565.1_ASM1786656v1_genomic.fna_master_datastore_index.log
cp GCA_017866565.1_ASM1786656v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866565.1_ASM1786656v1_genomic.fa
rm GCA_017866565.1_ASM1786656v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866575.1_ASM1786657v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866575.1_ASM1786657v1_genomic.fna.maker.output/GCA_017866575.1_ASM1786657v1_genomic.fna_master_datastore_index.log
cp GCA_017866575.1_ASM1786657v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866575.1_ASM1786657v1_genomic.fa
rm GCA_017866575.1_ASM1786657v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866745.1_ASM1786674v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866745.1_ASM1786674v1_genomic.fna.maker.output/GCA_017866745.1_ASM1786674v1_genomic.fna_master_datastore_index.log
cp GCA_017866745.1_ASM1786674v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866745.1_ASM1786674v1_genomic.fa
rm GCA_017866745.1_ASM1786674v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866755.1_ASM1786675v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866755.1_ASM1786675v1_genomic.fna.maker.output/GCA_017866755.1_ASM1786675v1_genomic.fna_master_datastore_index.log
cp GCA_017866755.1_ASM1786675v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866755.1_ASM1786675v1_genomic.fa
rm GCA_017866755.1_ASM1786675v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866845.1_ASM1786684v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866845.1_ASM1786684v1_genomic.fna.maker.output/GCA_017866845.1_ASM1786684v1_genomic.fna_master_datastore_index.log
cp GCA_017866845.1_ASM1786684v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866845.1_ASM1786684v1_genomic.fa
rm GCA_017866845.1_ASM1786684v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866945.1_ASM1786694v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866945.1_ASM1786694v1_genomic.fna.maker.output/GCA_017866945.1_ASM1786694v1_genomic.fna_master_datastore_index.log
cp GCA_017866945.1_ASM1786694v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866945.1_ASM1786694v1_genomic.fa
rm GCA_017866945.1_ASM1786694v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866955.1_ASM1786695v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866955.1_ASM1786695v1_genomic.fna.maker.output/GCA_017866955.1_ASM1786695v1_genomic.fna_master_datastore_index.log
cp GCA_017866955.1_ASM1786695v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866955.1_ASM1786695v1_genomic.fa
rm GCA_017866955.1_ASM1786695v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017866995.1_ASM1786699v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017866995.1_ASM1786699v1_genomic.fna.maker.output/GCA_017866995.1_ASM1786699v1_genomic.fna_master_datastore_index.log
cp GCA_017866995.1_ASM1786699v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017866995.1_ASM1786699v1_genomic.fa
rm GCA_017866995.1_ASM1786699v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867005.1_ASM1786700v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867005.1_ASM1786700v1_genomic.fna.maker.output/GCA_017867005.1_ASM1786700v1_genomic.fna_master_datastore_index.log
cp GCA_017867005.1_ASM1786700v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867005.1_ASM1786700v1_genomic.fa
rm GCA_017867005.1_ASM1786700v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867015.1_ASM1786701v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867015.1_ASM1786701v1_genomic.fna.maker.output/GCA_017867015.1_ASM1786701v1_genomic.fna_master_datastore_index.log
cp GCA_017867015.1_ASM1786701v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867015.1_ASM1786701v1_genomic.fa
rm GCA_017867015.1_ASM1786701v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867205.1_ASM1786720v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867205.1_ASM1786720v1_genomic.fna.maker.output/GCA_017867205.1_ASM1786720v1_genomic.fna_master_datastore_index.log
cp GCA_017867205.1_ASM1786720v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867205.1_ASM1786720v1_genomic.fa
rm GCA_017867205.1_ASM1786720v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867255.1_ASM1786725v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867255.1_ASM1786725v1_genomic.fna.maker.output/GCA_017867255.1_ASM1786725v1_genomic.fna_master_datastore_index.log
cp GCA_017867255.1_ASM1786725v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867255.1_ASM1786725v1_genomic.fa
rm GCA_017867255.1_ASM1786725v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867285.1_ASM1786728v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867285.1_ASM1786728v1_genomic.fna.maker.output/GCA_017867285.1_ASM1786728v1_genomic.fna_master_datastore_index.log
cp GCA_017867285.1_ASM1786728v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867285.1_ASM1786728v1_genomic.fa
rm GCA_017867285.1_ASM1786728v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867295.1_ASM1786729v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867295.1_ASM1786729v1_genomic.fna.maker.output/GCA_017867295.1_ASM1786729v1_genomic.fna_master_datastore_index.log
cp GCA_017867295.1_ASM1786729v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867295.1_ASM1786729v1_genomic.fa
rm GCA_017867295.1_ASM1786729v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867305.1_ASM1786730v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867305.1_ASM1786730v1_genomic.fna.maker.output/GCA_017867305.1_ASM1786730v1_genomic.fna_master_datastore_index.log
cp GCA_017867305.1_ASM1786730v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867305.1_ASM1786730v1_genomic.fa
rm GCA_017867305.1_ASM1786730v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867515.1_ASM1786751v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867515.1_ASM1786751v1_genomic.fna.maker.output/GCA_017867515.1_ASM1786751v1_genomic.fna_master_datastore_index.log
cp GCA_017867515.1_ASM1786751v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867515.1_ASM1786751v1_genomic.fa
rm GCA_017867515.1_ASM1786751v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867525.1_ASM1786752v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867525.1_ASM1786752v1_genomic.fna.maker.output/GCA_017867525.1_ASM1786752v1_genomic.fna_master_datastore_index.log
cp GCA_017867525.1_ASM1786752v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867525.1_ASM1786752v1_genomic.fa
rm GCA_017867525.1_ASM1786752v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867555.1_ASM1786755v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867555.1_ASM1786755v1_genomic.fna.maker.output/GCA_017867555.1_ASM1786755v1_genomic.fna_master_datastore_index.log
cp GCA_017867555.1_ASM1786755v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867555.1_ASM1786755v1_genomic.fa
rm GCA_017867555.1_ASM1786755v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867585.1_ASM1786758v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867585.1_ASM1786758v1_genomic.fna.maker.output/GCA_017867585.1_ASM1786758v1_genomic.fna_master_datastore_index.log
cp GCA_017867585.1_ASM1786758v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867585.1_ASM1786758v1_genomic.fa
rm GCA_017867585.1_ASM1786758v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867595.1_ASM1786759v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867595.1_ASM1786759v1_genomic.fna.maker.output/GCA_017867595.1_ASM1786759v1_genomic.fna_master_datastore_index.log
cp GCA_017867595.1_ASM1786759v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867595.1_ASM1786759v1_genomic.fa
rm GCA_017867595.1_ASM1786759v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867735.1_ASM1786773v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867735.1_ASM1786773v1_genomic.fna.maker.output/GCA_017867735.1_ASM1786773v1_genomic.fna_master_datastore_index.log
cp GCA_017867735.1_ASM1786773v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867735.1_ASM1786773v1_genomic.fa
rm GCA_017867735.1_ASM1786773v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867755.1_ASM1786775v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867755.1_ASM1786775v1_genomic.fna.maker.output/GCA_017867755.1_ASM1786775v1_genomic.fna_master_datastore_index.log
cp GCA_017867755.1_ASM1786775v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867755.1_ASM1786775v1_genomic.fa
rm GCA_017867755.1_ASM1786775v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867835.1_ASM1786783v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867835.1_ASM1786783v1_genomic.fna.maker.output/GCA_017867835.1_ASM1786783v1_genomic.fna_master_datastore_index.log
cp GCA_017867835.1_ASM1786783v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867835.1_ASM1786783v1_genomic.fa
rm GCA_017867835.1_ASM1786783v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867885.1_ASM1786788v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867885.1_ASM1786788v1_genomic.fna.maker.output/GCA_017867885.1_ASM1786788v1_genomic.fna_master_datastore_index.log
cp GCA_017867885.1_ASM1786788v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867885.1_ASM1786788v1_genomic.fa
rm GCA_017867885.1_ASM1786788v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867915.1_ASM1786791v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867915.1_ASM1786791v1_genomic.fna.maker.output/GCA_017867915.1_ASM1786791v1_genomic.fna_master_datastore_index.log
cp GCA_017867915.1_ASM1786791v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867915.1_ASM1786791v1_genomic.fa
rm GCA_017867915.1_ASM1786791v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867955.1_ASM1786795v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867955.1_ASM1786795v1_genomic.fna.maker.output/GCA_017867955.1_ASM1786795v1_genomic.fna_master_datastore_index.log
cp GCA_017867955.1_ASM1786795v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867955.1_ASM1786795v1_genomic.fa
rm GCA_017867955.1_ASM1786795v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017867995.1_ASM1786799v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017867995.1_ASM1786799v1_genomic.fna.maker.output/GCA_017867995.1_ASM1786799v1_genomic.fna_master_datastore_index.log
cp GCA_017867995.1_ASM1786799v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017867995.1_ASM1786799v1_genomic.fa
rm GCA_017867995.1_ASM1786799v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868025.1_ASM1786802v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868025.1_ASM1786802v1_genomic.fna.maker.output/GCA_017868025.1_ASM1786802v1_genomic.fna_master_datastore_index.log
cp GCA_017868025.1_ASM1786802v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868025.1_ASM1786802v1_genomic.fa
rm GCA_017868025.1_ASM1786802v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868065.1_ASM1786806v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868065.1_ASM1786806v1_genomic.fna.maker.output/GCA_017868065.1_ASM1786806v1_genomic.fna_master_datastore_index.log
cp GCA_017868065.1_ASM1786806v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868065.1_ASM1786806v1_genomic.fa
rm GCA_017868065.1_ASM1786806v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868075.1_ASM1786807v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868075.1_ASM1786807v1_genomic.fna.maker.output/GCA_017868075.1_ASM1786807v1_genomic.fna_master_datastore_index.log
cp GCA_017868075.1_ASM1786807v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868075.1_ASM1786807v1_genomic.fa
rm GCA_017868075.1_ASM1786807v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868145.1_ASM1786814v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868145.1_ASM1786814v1_genomic.fna.maker.output/GCA_017868145.1_ASM1786814v1_genomic.fna_master_datastore_index.log
cp GCA_017868145.1_ASM1786814v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868145.1_ASM1786814v1_genomic.fa
rm GCA_017868145.1_ASM1786814v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868155.1_ASM1786815v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868155.1_ASM1786815v1_genomic.fna.maker.output/GCA_017868155.1_ASM1786815v1_genomic.fna_master_datastore_index.log
cp GCA_017868155.1_ASM1786815v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868155.1_ASM1786815v1_genomic.fa
rm GCA_017868155.1_ASM1786815v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868165.1_ASM1786816v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868165.1_ASM1786816v1_genomic.fna.maker.output/GCA_017868165.1_ASM1786816v1_genomic.fna_master_datastore_index.log
cp GCA_017868165.1_ASM1786816v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868165.1_ASM1786816v1_genomic.fa
rm GCA_017868165.1_ASM1786816v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868175.1_ASM1786817v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868175.1_ASM1786817v1_genomic.fna.maker.output/GCA_017868175.1_ASM1786817v1_genomic.fna_master_datastore_index.log
cp GCA_017868175.1_ASM1786817v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868175.1_ASM1786817v1_genomic.fa
rm GCA_017868175.1_ASM1786817v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868185.1_ASM1786818v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868185.1_ASM1786818v1_genomic.fna.maker.output/GCA_017868185.1_ASM1786818v1_genomic.fna_master_datastore_index.log
cp GCA_017868185.1_ASM1786818v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868185.1_ASM1786818v1_genomic.fa
rm GCA_017868185.1_ASM1786818v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868345.1_ASM1786834v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868345.1_ASM1786834v1_genomic.fna.maker.output/GCA_017868345.1_ASM1786834v1_genomic.fna_master_datastore_index.log
cp GCA_017868345.1_ASM1786834v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868345.1_ASM1786834v1_genomic.fa
rm GCA_017868345.1_ASM1786834v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868355.1_ASM1786835v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868355.1_ASM1786835v1_genomic.fna.maker.output/GCA_017868355.1_ASM1786835v1_genomic.fna_master_datastore_index.log
cp GCA_017868355.1_ASM1786835v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868355.1_ASM1786835v1_genomic.fa
rm GCA_017868355.1_ASM1786835v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868485.1_ASM1786848v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868485.1_ASM1786848v1_genomic.fna.maker.output/GCA_017868485.1_ASM1786848v1_genomic.fna_master_datastore_index.log
cp GCA_017868485.1_ASM1786848v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868485.1_ASM1786848v1_genomic.fa
rm GCA_017868485.1_ASM1786848v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868655.1_ASM1786865v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868655.1_ASM1786865v1_genomic.fna.maker.output/GCA_017868655.1_ASM1786865v1_genomic.fna_master_datastore_index.log
cp GCA_017868655.1_ASM1786865v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868655.1_ASM1786865v1_genomic.fa
rm GCA_017868655.1_ASM1786865v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868835.1_ASM1786883v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868835.1_ASM1786883v1_genomic.fna.maker.output/GCA_017868835.1_ASM1786883v1_genomic.fna_master_datastore_index.log
cp GCA_017868835.1_ASM1786883v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868835.1_ASM1786883v1_genomic.fa
rm GCA_017868835.1_ASM1786883v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868915.1_ASM1786891v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868915.1_ASM1786891v1_genomic.fna.maker.output/GCA_017868915.1_ASM1786891v1_genomic.fna_master_datastore_index.log
cp GCA_017868915.1_ASM1786891v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868915.1_ASM1786891v1_genomic.fa
rm GCA_017868915.1_ASM1786891v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868925.1_ASM1786892v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868925.1_ASM1786892v1_genomic.fna.maker.output/GCA_017868925.1_ASM1786892v1_genomic.fna_master_datastore_index.log
cp GCA_017868925.1_ASM1786892v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868925.1_ASM1786892v1_genomic.fa
rm GCA_017868925.1_ASM1786892v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017868975.1_ASM1786897v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017868975.1_ASM1786897v1_genomic.fna.maker.output/GCA_017868975.1_ASM1786897v1_genomic.fna_master_datastore_index.log
cp GCA_017868975.1_ASM1786897v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017868975.1_ASM1786897v1_genomic.fa
rm GCA_017868975.1_ASM1786897v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869095.1_ASM1786909v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869095.1_ASM1786909v1_genomic.fna.maker.output/GCA_017869095.1_ASM1786909v1_genomic.fna_master_datastore_index.log
cp GCA_017869095.1_ASM1786909v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869095.1_ASM1786909v1_genomic.fa
rm GCA_017869095.1_ASM1786909v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869105.1_ASM1786910v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869105.1_ASM1786910v1_genomic.fna.maker.output/GCA_017869105.1_ASM1786910v1_genomic.fna_master_datastore_index.log
cp GCA_017869105.1_ASM1786910v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869105.1_ASM1786910v1_genomic.fa
rm GCA_017869105.1_ASM1786910v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869165.1_ASM1786916v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869165.1_ASM1786916v1_genomic.fna.maker.output/GCA_017869165.1_ASM1786916v1_genomic.fna_master_datastore_index.log
cp GCA_017869165.1_ASM1786916v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869165.1_ASM1786916v1_genomic.fa
rm GCA_017869165.1_ASM1786916v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869305.1_ASM1786930v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869305.1_ASM1786930v1_genomic.fna.maker.output/GCA_017869305.1_ASM1786930v1_genomic.fna_master_datastore_index.log
cp GCA_017869305.1_ASM1786930v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869305.1_ASM1786930v1_genomic.fa
rm GCA_017869305.1_ASM1786930v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869315.1_ASM1786931v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869315.1_ASM1786931v1_genomic.fna.maker.output/GCA_017869315.1_ASM1786931v1_genomic.fna_master_datastore_index.log
cp GCA_017869315.1_ASM1786931v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869315.1_ASM1786931v1_genomic.fa
rm GCA_017869315.1_ASM1786931v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869325.1_ASM1786932v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869325.1_ASM1786932v1_genomic.fna.maker.output/GCA_017869325.1_ASM1786932v1_genomic.fna_master_datastore_index.log
cp GCA_017869325.1_ASM1786932v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869325.1_ASM1786932v1_genomic.fa
rm GCA_017869325.1_ASM1786932v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869335.1_ASM1786933v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869335.1_ASM1786933v1_genomic.fna.maker.output/GCA_017869335.1_ASM1786933v1_genomic.fna_master_datastore_index.log
cp GCA_017869335.1_ASM1786933v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869335.1_ASM1786933v1_genomic.fa
rm GCA_017869335.1_ASM1786933v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869405.1_ASM1786940v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869405.1_ASM1786940v1_genomic.fna.maker.output/GCA_017869405.1_ASM1786940v1_genomic.fna_master_datastore_index.log
cp GCA_017869405.1_ASM1786940v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869405.1_ASM1786940v1_genomic.fa
rm GCA_017869405.1_ASM1786940v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869445.1_ASM1786944v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869445.1_ASM1786944v1_genomic.fna.maker.output/GCA_017869445.1_ASM1786944v1_genomic.fna_master_datastore_index.log
cp GCA_017869445.1_ASM1786944v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869445.1_ASM1786944v1_genomic.fa
rm GCA_017869445.1_ASM1786944v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869495.1_ASM1786949v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869495.1_ASM1786949v1_genomic.fna.maker.output/GCA_017869495.1_ASM1786949v1_genomic.fna_master_datastore_index.log
cp GCA_017869495.1_ASM1786949v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869495.1_ASM1786949v1_genomic.fa
rm GCA_017869495.1_ASM1786949v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869545.1_ASM1786954v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869545.1_ASM1786954v1_genomic.fna.maker.output/GCA_017869545.1_ASM1786954v1_genomic.fna_master_datastore_index.log
cp GCA_017869545.1_ASM1786954v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869545.1_ASM1786954v1_genomic.fa
rm GCA_017869545.1_ASM1786954v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869555.1_ASM1786955v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869555.1_ASM1786955v1_genomic.fna.maker.output/GCA_017869555.1_ASM1786955v1_genomic.fna_master_datastore_index.log
cp GCA_017869555.1_ASM1786955v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869555.1_ASM1786955v1_genomic.fa
rm GCA_017869555.1_ASM1786955v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869615.1_ASM1786961v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869615.1_ASM1786961v1_genomic.fna.maker.output/GCA_017869615.1_ASM1786961v1_genomic.fna_master_datastore_index.log
cp GCA_017869615.1_ASM1786961v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869615.1_ASM1786961v1_genomic.fa
rm GCA_017869615.1_ASM1786961v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869625.1_ASM1786962v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869625.1_ASM1786962v1_genomic.fna.maker.output/GCA_017869625.1_ASM1786962v1_genomic.fna_master_datastore_index.log
cp GCA_017869625.1_ASM1786962v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869625.1_ASM1786962v1_genomic.fa
rm GCA_017869625.1_ASM1786962v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869635.1_ASM1786963v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869635.1_ASM1786963v1_genomic.fna.maker.output/GCA_017869635.1_ASM1786963v1_genomic.fna_master_datastore_index.log
cp GCA_017869635.1_ASM1786963v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869635.1_ASM1786963v1_genomic.fa
rm GCA_017869635.1_ASM1786963v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869715.1_ASM1786971v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869715.1_ASM1786971v1_genomic.fna.maker.output/GCA_017869715.1_ASM1786971v1_genomic.fna_master_datastore_index.log
cp GCA_017869715.1_ASM1786971v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869715.1_ASM1786971v1_genomic.fa
rm GCA_017869715.1_ASM1786971v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869735.1_ASM1786973v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869735.1_ASM1786973v1_genomic.fna.maker.output/GCA_017869735.1_ASM1786973v1_genomic.fna_master_datastore_index.log
cp GCA_017869735.1_ASM1786973v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869735.1_ASM1786973v1_genomic.fa
rm GCA_017869735.1_ASM1786973v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869745.1_ASM1786974v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869745.1_ASM1786974v1_genomic.fna.maker.output/GCA_017869745.1_ASM1786974v1_genomic.fna_master_datastore_index.log
cp GCA_017869745.1_ASM1786974v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869745.1_ASM1786974v1_genomic.fa
rm GCA_017869745.1_ASM1786974v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869765.1_ASM1786976v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869765.1_ASM1786976v1_genomic.fna.maker.output/GCA_017869765.1_ASM1786976v1_genomic.fna_master_datastore_index.log
cp GCA_017869765.1_ASM1786976v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869765.1_ASM1786976v1_genomic.fa
rm GCA_017869765.1_ASM1786976v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869795.1_ASM1786979v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869795.1_ASM1786979v1_genomic.fna.maker.output/GCA_017869795.1_ASM1786979v1_genomic.fna_master_datastore_index.log
cp GCA_017869795.1_ASM1786979v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869795.1_ASM1786979v1_genomic.fa
rm GCA_017869795.1_ASM1786979v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869965.1_ASM1786996v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869965.1_ASM1786996v1_genomic.fna.maker.output/GCA_017869965.1_ASM1786996v1_genomic.fna_master_datastore_index.log
cp GCA_017869965.1_ASM1786996v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869965.1_ASM1786996v1_genomic.fa
rm GCA_017869965.1_ASM1786996v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017869985.1_ASM1786998v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017869985.1_ASM1786998v1_genomic.fna.maker.output/GCA_017869985.1_ASM1786998v1_genomic.fna_master_datastore_index.log
cp GCA_017869985.1_ASM1786998v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017869985.1_ASM1786998v1_genomic.fa
rm GCA_017869985.1_ASM1786998v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870085.1_ASM1787008v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870085.1_ASM1787008v1_genomic.fna.maker.output/GCA_017870085.1_ASM1787008v1_genomic.fna_master_datastore_index.log
cp GCA_017870085.1_ASM1787008v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870085.1_ASM1787008v1_genomic.fa
rm GCA_017870085.1_ASM1787008v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870115.1_ASM1787011v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870115.1_ASM1787011v1_genomic.fna.maker.output/GCA_017870115.1_ASM1787011v1_genomic.fna_master_datastore_index.log
cp GCA_017870115.1_ASM1787011v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870115.1_ASM1787011v1_genomic.fa
rm GCA_017870115.1_ASM1787011v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870125.1_ASM1787012v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870125.1_ASM1787012v1_genomic.fna.maker.output/GCA_017870125.1_ASM1787012v1_genomic.fna_master_datastore_index.log
cp GCA_017870125.1_ASM1787012v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870125.1_ASM1787012v1_genomic.fa
rm GCA_017870125.1_ASM1787012v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870425.1_ASM1787042v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870425.1_ASM1787042v1_genomic.fna.maker.output/GCA_017870425.1_ASM1787042v1_genomic.fna_master_datastore_index.log
cp GCA_017870425.1_ASM1787042v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870425.1_ASM1787042v1_genomic.fa
rm GCA_017870425.1_ASM1787042v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870435.1_ASM1787043v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870435.1_ASM1787043v1_genomic.fna.maker.output/GCA_017870435.1_ASM1787043v1_genomic.fna_master_datastore_index.log
cp GCA_017870435.1_ASM1787043v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870435.1_ASM1787043v1_genomic.fa
rm GCA_017870435.1_ASM1787043v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870445.1_ASM1787044v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870445.1_ASM1787044v1_genomic.fna.maker.output/GCA_017870445.1_ASM1787044v1_genomic.fna_master_datastore_index.log
cp GCA_017870445.1_ASM1787044v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870445.1_ASM1787044v1_genomic.fa
rm GCA_017870445.1_ASM1787044v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870455.1_ASM1787045v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870455.1_ASM1787045v1_genomic.fna.maker.output/GCA_017870455.1_ASM1787045v1_genomic.fna_master_datastore_index.log
cp GCA_017870455.1_ASM1787045v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870455.1_ASM1787045v1_genomic.fa
rm GCA_017870455.1_ASM1787045v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870595.1_ASM1787059v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870595.1_ASM1787059v1_genomic.fna.maker.output/GCA_017870595.1_ASM1787059v1_genomic.fna_master_datastore_index.log
cp GCA_017870595.1_ASM1787059v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870595.1_ASM1787059v1_genomic.fa
rm GCA_017870595.1_ASM1787059v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870715.1_ASM1787071v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870715.1_ASM1787071v1_genomic.fna.maker.output/GCA_017870715.1_ASM1787071v1_genomic.fna_master_datastore_index.log
cp GCA_017870715.1_ASM1787071v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870715.1_ASM1787071v1_genomic.fa
rm GCA_017870715.1_ASM1787071v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870815.1_ASM1787081v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870815.1_ASM1787081v1_genomic.fna.maker.output/GCA_017870815.1_ASM1787081v1_genomic.fna_master_datastore_index.log
cp GCA_017870815.1_ASM1787081v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870815.1_ASM1787081v1_genomic.fa
rm GCA_017870815.1_ASM1787081v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870925.1_ASM1787092v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870925.1_ASM1787092v1_genomic.fna.maker.output/GCA_017870925.1_ASM1787092v1_genomic.fna_master_datastore_index.log
cp GCA_017870925.1_ASM1787092v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870925.1_ASM1787092v1_genomic.fa
rm GCA_017870925.1_ASM1787092v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870945.1_ASM1787094v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870945.1_ASM1787094v1_genomic.fna.maker.output/GCA_017870945.1_ASM1787094v1_genomic.fna_master_datastore_index.log
cp GCA_017870945.1_ASM1787094v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870945.1_ASM1787094v1_genomic.fa
rm GCA_017870945.1_ASM1787094v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017870955.1_ASM1787095v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017870955.1_ASM1787095v1_genomic.fna.maker.output/GCA_017870955.1_ASM1787095v1_genomic.fna_master_datastore_index.log
cp GCA_017870955.1_ASM1787095v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017870955.1_ASM1787095v1_genomic.fa
rm GCA_017870955.1_ASM1787095v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017871015.1_ASM1787101v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017871015.1_ASM1787101v1_genomic.fna.maker.output/GCA_017871015.1_ASM1787101v1_genomic.fna_master_datastore_index.log
cp GCA_017871015.1_ASM1787101v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017871015.1_ASM1787101v1_genomic.fa
rm GCA_017871015.1_ASM1787101v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017871025.1_ASM1787102v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017871025.1_ASM1787102v1_genomic.fna.maker.output/GCA_017871025.1_ASM1787102v1_genomic.fna_master_datastore_index.log
cp GCA_017871025.1_ASM1787102v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017871025.1_ASM1787102v1_genomic.fa
rm GCA_017871025.1_ASM1787102v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017871065.1_ASM1787106v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017871065.1_ASM1787106v1_genomic.fna.maker.output/GCA_017871065.1_ASM1787106v1_genomic.fna_master_datastore_index.log
cp GCA_017871065.1_ASM1787106v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017871065.1_ASM1787106v1_genomic.fa
rm GCA_017871065.1_ASM1787106v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_017871085.1_ASM1787108v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_017871085.1_ASM1787108v1_genomic.fna.maker.output/GCA_017871085.1_ASM1787108v1_genomic.fna_master_datastore_index.log
cp GCA_017871085.1_ASM1787108v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_017871085.1_ASM1787108v1_genomic.fa
rm GCA_017871085.1_ASM1787108v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393555.1_ASM1939355v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393555.1_ASM1939355v1_genomic.fna.maker.output/GCA_019393555.1_ASM1939355v1_genomic.fna_master_datastore_index.log
cp GCA_019393555.1_ASM1939355v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393555.1_ASM1939355v1_genomic.fa
rm GCA_019393555.1_ASM1939355v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393615.1_ASM1939361v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393615.1_ASM1939361v1_genomic.fna.maker.output/GCA_019393615.1_ASM1939361v1_genomic.fna_master_datastore_index.log
cp GCA_019393615.1_ASM1939361v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393615.1_ASM1939361v1_genomic.fa
rm GCA_019393615.1_ASM1939361v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393625.1_ASM1939362v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393625.1_ASM1939362v1_genomic.fna.maker.output/GCA_019393625.1_ASM1939362v1_genomic.fna_master_datastore_index.log
cp GCA_019393625.1_ASM1939362v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393625.1_ASM1939362v1_genomic.fa
rm GCA_019393625.1_ASM1939362v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393785.1_ASM1939378v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393785.1_ASM1939378v1_genomic.fna.maker.output/GCA_019393785.1_ASM1939378v1_genomic.fna_master_datastore_index.log
cp GCA_019393785.1_ASM1939378v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393785.1_ASM1939378v1_genomic.fa
rm GCA_019393785.1_ASM1939378v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393865.1_ASM1939386v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393865.1_ASM1939386v1_genomic.fna.maker.output/GCA_019393865.1_ASM1939386v1_genomic.fna_master_datastore_index.log
cp GCA_019393865.1_ASM1939386v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393865.1_ASM1939386v1_genomic.fa
rm GCA_019393865.1_ASM1939386v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393935.1_ASM1939393v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393935.1_ASM1939393v1_genomic.fna.maker.output/GCA_019393935.1_ASM1939393v1_genomic.fna_master_datastore_index.log
cp GCA_019393935.1_ASM1939393v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393935.1_ASM1939393v1_genomic.fa
rm GCA_019393935.1_ASM1939393v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393945.1_ASM1939394v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393945.1_ASM1939394v1_genomic.fna.maker.output/GCA_019393945.1_ASM1939394v1_genomic.fna_master_datastore_index.log
cp GCA_019393945.1_ASM1939394v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393945.1_ASM1939394v1_genomic.fa
rm GCA_019393945.1_ASM1939394v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393955.1_ASM1939395v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393955.1_ASM1939395v1_genomic.fna.maker.output/GCA_019393955.1_ASM1939395v1_genomic.fna_master_datastore_index.log
cp GCA_019393955.1_ASM1939395v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393955.1_ASM1939395v1_genomic.fa
rm GCA_019393955.1_ASM1939395v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019393975.1_ASM1939397v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019393975.1_ASM1939397v1_genomic.fna.maker.output/GCA_019393975.1_ASM1939397v1_genomic.fna_master_datastore_index.log
cp GCA_019393975.1_ASM1939397v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019393975.1_ASM1939397v1_genomic.fa
rm GCA_019393975.1_ASM1939397v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394065.1_ASM1939406v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394065.1_ASM1939406v1_genomic.fna.maker.output/GCA_019394065.1_ASM1939406v1_genomic.fna_master_datastore_index.log
cp GCA_019394065.1_ASM1939406v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394065.1_ASM1939406v1_genomic.fa
rm GCA_019394065.1_ASM1939406v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394075.1_ASM1939407v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394075.1_ASM1939407v1_genomic.fna.maker.output/GCA_019394075.1_ASM1939407v1_genomic.fna_master_datastore_index.log
cp GCA_019394075.1_ASM1939407v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394075.1_ASM1939407v1_genomic.fa
rm GCA_019394075.1_ASM1939407v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394085.1_ASM1939408v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394085.1_ASM1939408v1_genomic.fna.maker.output/GCA_019394085.1_ASM1939408v1_genomic.fna_master_datastore_index.log
cp GCA_019394085.1_ASM1939408v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394085.1_ASM1939408v1_genomic.fa
rm GCA_019394085.1_ASM1939408v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394095.1_ASM1939409v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394095.1_ASM1939409v1_genomic.fna.maker.output/GCA_019394095.1_ASM1939409v1_genomic.fna_master_datastore_index.log
cp GCA_019394095.1_ASM1939409v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394095.1_ASM1939409v1_genomic.fa
rm GCA_019394095.1_ASM1939409v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394105.1_ASM1939410v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394105.1_ASM1939410v1_genomic.fna.maker.output/GCA_019394105.1_ASM1939410v1_genomic.fna_master_datastore_index.log
cp GCA_019394105.1_ASM1939410v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394105.1_ASM1939410v1_genomic.fa
rm GCA_019394105.1_ASM1939410v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394265.1_ASM1939426v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394265.1_ASM1939426v1_genomic.fna.maker.output/GCA_019394265.1_ASM1939426v1_genomic.fna_master_datastore_index.log
cp GCA_019394265.1_ASM1939426v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394265.1_ASM1939426v1_genomic.fa
rm GCA_019394265.1_ASM1939426v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394275.1_ASM1939427v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394275.1_ASM1939427v1_genomic.fna.maker.output/GCA_019394275.1_ASM1939427v1_genomic.fna_master_datastore_index.log
cp GCA_019394275.1_ASM1939427v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394275.1_ASM1939427v1_genomic.fa
rm GCA_019394275.1_ASM1939427v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 5
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394305.1_ASM1939430v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394305.1_ASM1939430v1_genomic.fna.maker.output/GCA_019394305.1_ASM1939430v1_genomic.fna_master_datastore_index.log
cp GCA_019394305.1_ASM1939430v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394305.1_ASM1939430v1_genomic.fa
rm GCA_019394305.1_ASM1939430v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 6
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394315.1_ASM1939431v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394315.1_ASM1939431v1_genomic.fna.maker.output/GCA_019394315.1_ASM1939431v1_genomic.fna_master_datastore_index.log
cp GCA_019394315.1_ASM1939431v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394315.1_ASM1939431v1_genomic.fa
rm GCA_019394315.1_ASM1939431v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 7
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394325.1_ASM1939432v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394325.1_ASM1939432v1_genomic.fna.maker.output/GCA_019394325.1_ASM1939432v1_genomic.fna_master_datastore_index.log
cp GCA_019394325.1_ASM1939432v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394325.1_ASM1939432v1_genomic.fa
rm GCA_019394325.1_ASM1939432v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 8
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394485.1_ASM1939448v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394485.1_ASM1939448v1_genomic.fna.maker.output/GCA_019394485.1_ASM1939448v1_genomic.fna_master_datastore_index.log
cp GCA_019394485.1_ASM1939448v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394485.1_ASM1939448v1_genomic.fa
rm GCA_019394485.1_ASM1939448v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 9
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394495.1_ASM1939449v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394495.1_ASM1939449v1_genomic.fna.maker.output/GCA_019394495.1_ASM1939449v1_genomic.fna_master_datastore_index.log
cp GCA_019394495.1_ASM1939449v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394495.1_ASM1939449v1_genomic.fa
rm GCA_019394495.1_ASM1939449v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 10
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394525.1_ASM1939452v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394525.1_ASM1939452v1_genomic.fna.maker.output/GCA_019394525.1_ASM1939452v1_genomic.fna_master_datastore_index.log
cp GCA_019394525.1_ASM1939452v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394525.1_ASM1939452v1_genomic.fa
rm GCA_019394525.1_ASM1939452v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 11
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394535.1_ASM1939453v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394535.1_ASM1939453v1_genomic.fna.maker.output/GCA_019394535.1_ASM1939453v1_genomic.fna_master_datastore_index.log
cp GCA_019394535.1_ASM1939453v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394535.1_ASM1939453v1_genomic.fa
rm GCA_019394535.1_ASM1939453v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 12
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394585.1_ASM1939458v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394585.1_ASM1939458v1_genomic.fna.maker.output/GCA_019394585.1_ASM1939458v1_genomic.fna_master_datastore_index.log
cp GCA_019394585.1_ASM1939458v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394585.1_ASM1939458v1_genomic.fa
rm GCA_019394585.1_ASM1939458v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 13
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394645.1_ASM1939464v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394645.1_ASM1939464v1_genomic.fna.maker.output/GCA_019394645.1_ASM1939464v1_genomic.fna_master_datastore_index.log
cp GCA_019394645.1_ASM1939464v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394645.1_ASM1939464v1_genomic.fa
rm GCA_019394645.1_ASM1939464v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 14
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394655.1_ASM1939465v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394655.1_ASM1939465v1_genomic.fna.maker.output/GCA_019394655.1_ASM1939465v1_genomic.fna_master_datastore_index.log
cp GCA_019394655.1_ASM1939465v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394655.1_ASM1939465v1_genomic.fa
rm GCA_019394655.1_ASM1939465v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 15
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394665.1_ASM1939466v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394665.1_ASM1939466v1_genomic.fna.maker.output/GCA_019394665.1_ASM1939466v1_genomic.fna_master_datastore_index.log
cp GCA_019394665.1_ASM1939466v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394665.1_ASM1939466v1_genomic.fa
rm GCA_019394665.1_ASM1939466v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 16
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394675.1_ASM1939467v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394675.1_ASM1939467v1_genomic.fna.maker.output/GCA_019394675.1_ASM1939467v1_genomic.fna_master_datastore_index.log
cp GCA_019394675.1_ASM1939467v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394675.1_ASM1939467v1_genomic.fa
rm GCA_019394675.1_ASM1939467v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 17
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394725.1_ASM1939472v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394725.1_ASM1939472v1_genomic.fna.maker.output/GCA_019394725.1_ASM1939472v1_genomic.fna_master_datastore_index.log
cp GCA_019394725.1_ASM1939472v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394725.1_ASM1939472v1_genomic.fa
rm GCA_019394725.1_ASM1939472v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 18
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394775.1_ASM1939477v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394775.1_ASM1939477v1_genomic.fna.maker.output/GCA_019394775.1_ASM1939477v1_genomic.fna_master_datastore_index.log
cp GCA_019394775.1_ASM1939477v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394775.1_ASM1939477v1_genomic.fa
rm GCA_019394775.1_ASM1939477v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 19
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394785.1_ASM1939478v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394785.1_ASM1939478v1_genomic.fna.maker.output/GCA_019394785.1_ASM1939478v1_genomic.fna_master_datastore_index.log
cp GCA_019394785.1_ASM1939478v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394785.1_ASM1939478v1_genomic.fa
rm GCA_019394785.1_ASM1939478v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 20
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394805.1_ASM1939480v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394805.1_ASM1939480v1_genomic.fna.maker.output/GCA_019394805.1_ASM1939480v1_genomic.fna_master_datastore_index.log
cp GCA_019394805.1_ASM1939480v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394805.1_ASM1939480v1_genomic.fa
rm GCA_019394805.1_ASM1939480v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 1
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394815.1_ASM1939481v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394815.1_ASM1939481v1_genomic.fna.maker.output/GCA_019394815.1_ASM1939481v1_genomic.fna_master_datastore_index.log
cp GCA_019394815.1_ASM1939481v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394815.1_ASM1939481v1_genomic.fa
rm GCA_019394815.1_ASM1939481v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

wait
(cd 2
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394885.1_ASM1939488v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394885.1_ASM1939488v1_genomic.fna.maker.output/GCA_019394885.1_ASM1939488v1_genomic.fna_master_datastore_index.log
cp GCA_019394885.1_ASM1939488v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394885.1_ASM1939488v1_genomic.fa
rm GCA_019394885.1_ASM1939488v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 3
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394895.1_ASM1939489v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394895.1_ASM1939489v1_genomic.fna.maker.output/GCA_019394895.1_ASM1939489v1_genomic.fna_master_datastore_index.log
cp GCA_019394895.1_ASM1939489v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394895.1_ASM1939489v1_genomic.fa
rm GCA_019394895.1_ASM1939489v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

(cd 4
sed -i 's/masked_genome\/.*\.fna\.masked/masked_genome\/GCA_019394915.1_ASM1939491v1_genomic\.fna\.masked/g' maker_opts.ctl
maker -cpus 4 -fix_nucleotides 
fasta_merge -d GCA_019394915.1_ASM1939491v1_genomic.fna.maker.output/GCA_019394915.1_ASM1939491v1_genomic.fna_master_datastore_index.log
cp GCA_019394915.1_ASM1939491v1_genomic.fna.all.maker.proteins.fasta /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker/GCA_019394915.1_ASM1939491v1_genomic.fa
rm GCA_019394915.1_ASM1939491v1_genomic.fna.all*
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/maker ) &

