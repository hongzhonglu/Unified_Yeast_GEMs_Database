# blastp:cenpk1 vs pan1800_v2_50_70
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/cenpk1.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1800_50_70_v2.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/cenpk_vs_pan1800_50_70_v2.txt tmp --min-seq-id 0.6 -e 1e-5

# blastp:cenpk1 vs lgpangenome
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/cenpk1.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/lg_pangenome.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/cenpk_vs_lgpan_blastp.txt tmp --min-seq-id 0.6 -e 1e-5

# blastp:cenpk1 vs napangenome
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/cenpk1.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1011_v1.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/cenpk_vs_pan1011_blastp.txt tmp --min-seq-id 0.6 -e 1e-5

# blastp:s288c vs pan1800_v2_50_70
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/S288c_R64.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1800_50_70_v2.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/s288c_vs_pan1800_50_70_v2.txt tmp --min-seq-id 0.6 -e 1e-5


#blastp: pan1800_v2_50_70 vs napgenome
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1800_50_70_v2.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1011_v1.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/pan1800_v2_vs_pan1011_blastp.txt tmp --min-seq-id 0.6 -e 1e-5

#blastp: pan1800_v2_50_70 vs lgpangenome
mmseqs easy-search /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/pan1800_50_70_v2.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/data/genome/lg_pangenome.fasta /mnt/d/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database/code/3.pan-genome_construction/3.pan-genome_comparison/output/pan1800_v2_vs_lgpan_blastp.txt tmp --min-seq-id 0.6 -e 1e-5


