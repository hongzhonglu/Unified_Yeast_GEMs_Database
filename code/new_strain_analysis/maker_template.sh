strain="strain_name"
rm_output_dir="repeatmasker_output_dir"
genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome"
masked_genome_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/masked_genome"
allcds_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/allcds_maker"

RepeatMasker -species "saccharomyces cerevisiae" -xsmall -e ncbi -dir $rm_output_dir/$strain $genome_dir/$strain.fa
cp $rm_output_dir/$strain/$strain.fna.masked $masked_genome_dir/ 
sed -i 's/masked_genome\/.*\.fna\.masked/data\/strain\.fna\.masked/g' maker_opts.ctl
maker -fix_nucleotides
fasta_merge -d $strain.fna/$strain.fna_master_datastore_index.log
cp $strain.fna.all.maker.proteins.fasta $allcds_dir/$strain.fa