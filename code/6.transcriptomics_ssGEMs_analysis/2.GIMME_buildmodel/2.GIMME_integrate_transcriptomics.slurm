#!/bin/bash
#SBATCH -J GIMME_run
#SBATCH -p 64c512g
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -n 10
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn

module load matlab/r2022a
cd /dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/7.transcriptomics_ssGEMs_analysis
matlab -r GIMME_run