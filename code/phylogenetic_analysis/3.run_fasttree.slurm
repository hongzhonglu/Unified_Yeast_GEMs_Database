#!/bin/bash

#SBATCH --job-name=fast_tree
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load miniconda3
source activate /dssg/home/acct-clslhz/clslhz/miniconda3/envs/hyena-d/fasttree

cd /dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/phylogenetic_analysis

FastTree -lg  coreORFs_supermatrix_limit.fas > coreORFs_fasttree.tre