#!/bin/bash

#SBATCH --job-name=iqtree_sce
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

cd /dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database
# mkdir
if [ ! -d code/phylogenetic_analysis/output/coreORFs_tree ]; then
    mkdir code/phylogenetic_analysis/output/coreORFs_tree
fi

iqtree -s code/phylogenetic_analysis/iqtree_input --prefix code/phylogenetic_analysis/output/coreORFs_tree -B 1000 -T AUTO -m Q.yeast+G4
