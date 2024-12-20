# 1. use morphological mode to build tree by pan-genome existence/absence matrix
iqtree -s code/phylogenetic_analysis/output/pangenome_existence.fa -st MORPH -B 1000 -m MK+R9 --prefix code/phylogenetic_analysis/output/pan_existence_tree/pan_existence

# 2.Core genome concatenate-based tree
# set substitution model: Q.yeast+G4
iqtree -s code/phylogenetic_analysis/output/core_ORFs_mas_trimmed --prefix code/phylogenetic_analysis/output/coreORF_tree/coreORFs_conct_tree -B 1000 -T AUTO -m Q.yeast+G4

iqtree -s coreORFs_supermatrix_limit.fas --prefix code/phylogenetic_analysis/output/coreORF_tree/coreORFs_conct_tree -B 1000 -T AUTO -m Q.yeast+G4

# automatically find best model
iqtree -s code/phylogenetic_analysis/output/core_ORFs_mas_trimed --prefix code/phylogenetic_analysis/output/coreORF_tree/coreORFs_conct_tree -B 1000 -T AUTO


