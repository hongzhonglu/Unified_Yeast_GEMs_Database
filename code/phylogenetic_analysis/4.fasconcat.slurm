#!/bin/bash

#SBATCH --job-name=fasconcat
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /dssg/home/acct-clslhz/clslhz/miniconda3/envs/hyena-d/iqtree

cd /dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/phylogenetic_analysis/fasconcat_input

perl FASconCAT-G_v1.06.1.pl -j -s