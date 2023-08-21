'''combine gff files from MAKER and liftover methods'''

from Bio import SeqIO
import os
import numpy as np
import sys
sys.path.append("code/genome_annotaion-combining")
from 2.combine_gene_prediction import load_genes,genePreds

lift_gff_dir = 'data/genome/predicted_allcds/flo_gff_result/'
lift_prot_dir = 'data/genome/predicted_allcds/flo_proteomes/'
maker_gff_dir = 'data/genome/predicted_allcds/maker_gff_file/maker_GCA_019394805.gff/'
maker_prot_dir = 'data/genome/predicted_allcds/maker_proteomes/'
outdir = 'data/genome/predicted_allcds/combine_gff/'

strainlist=os.listdir(lift_gff_dir)
for strain in strainlist:
    strain_name=strain.replace(".gff","")
    lift_gff = os.path.join(lift_gff_dir,strain_name+".gff")
    lift_prot = os.path.join(lift_prot_dir,strain_name+".fasta")
    maker_gff = os.path.join(maker_gff_dir,strain_name+".gff")
    maker_prot = os.path.join(maker_prot_dir,strain_name+".fa")

    lift_genes = load_genes(lift_gff, lift_prot, 'lift')
    maker_genes = load_genes(maker_gff, maker_prot, 'maker')
    gp = genePreds(lift_genes,maker_genes)
    gp.combine()
    combined_genes=gp.combined_genes
    with open(outdir+"cb_"+strain_name+".gff","w") as f:
        for gene in combined_genes:
            gff_line=gene.gff
            gff_line=gff_line.replace("SGD","liftover")
            f.write(gff_line)






