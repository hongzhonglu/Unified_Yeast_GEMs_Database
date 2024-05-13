#!/bin/bash

#SBATCH --job-name=mmseq_pan1900_vs_nr
#SBATCH --partition=64c512g
#SBATCH -N 2
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /dssg/home/acct-clslhz/clslhz/miniconda3/envs/mmseqs2

python blastp.py -i /dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/data/genome/pan1900_new_50_70_all9323_rename.fasta

