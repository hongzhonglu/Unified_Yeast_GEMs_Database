# mmseqs2 all vs all blastp
mmseqs easy-search data/genome/pan1900_new_50_70_all9323.fasta data/genome/pan1900_new_50_70_all9323.fasta code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_all9323_vs_all9323.txt tmp --min-seq-id 0.6 -e 1e-5 --format-output "query,target,pident,bits,qcov,tcov,qlen,tlen,evalue"

# pan1900_v4 all vs all blastp
mmseqs easy-search data/genome/pan1900_new_v4_50_70.fasta data/genome/pan1900_new_v4_50_70.fasta code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_v4_all_vs_all.txt tmp --min-seq-id 0.6 -e 1e-5 --format-output "query,target,pident,bits,qcov,tcov,qlen,tlen,evalue"

# pan1900_v4_50_70_v2 all vs all blastp
mmseqs easy-search code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2.fasta code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_new_v4_50_70_v2.fasta code/3.pan-genome_construction/1.new_pan-genome_pipeline/output/pan1900_v4_50_70_v2_all_vs_all.txt tmp --min-seq-id 0.6 -e 1e-5 --format-output "query,target,pident,bits,qcov,tcov,qlen,tlen,evalue"





