# 1.pan-genome existence/absence matrix
mafft --thread 4 --auto --maxiterate 1000 pangenome_existence.fa > pangenome_existence_msa.fa
# trimming
trimal -in test -out aligned_trimed.fa -gappyout

# 2.Core geno