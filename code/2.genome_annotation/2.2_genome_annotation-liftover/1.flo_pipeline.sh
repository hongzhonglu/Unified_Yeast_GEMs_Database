# flo pipeline
# step1 : prepare reference transcripts gff file and assembly fasta file
less data/saccharomyces_cerevisiae_R64-2-1_20150113.gff|grep '^#' > data/R64.gff
vim R64.gff #to delete "###" "##Fasta"
less data/saccharomyces_cerevisiae_R64-2-1_20150113.gff|grep '^chr' >> data/R64.gff
# prepare the assembled genome fasta file
less data/saccharomyces_cerevisiae_R64-2-1_20150113.gff |grep -v '#'|grep -v '^chr'>data/R64.fa

#step2: clean the gff file——only leave the "mRNA" & "CDS" annotation
less data/R64.gff|grep '^#' > data/R64_cleaned_0.gff
less data/R64.gff|grep '^chr'|grep -E 'SGD\s(mRNA|CDS|gene)' >> data/R64_cleaned_0.gff
~/biosoft/flo/gff_remove_feats.rb gene data/R64_cleaned_0.gff > data/R64_cleaned.gff

#step3: rename the mRNA name
python3 code/genome_annotation-liftover/clean_parent.py

#step4:build flo scripts——3.build_flo_script.py

