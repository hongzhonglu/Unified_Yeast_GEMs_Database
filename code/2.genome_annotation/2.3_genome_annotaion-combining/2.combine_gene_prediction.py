# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：combine_gene_prediction.py
# 2022/6/29
# Combine the results from liftOver and MAKER
# (1) any sequences either from liftOver or MAKER with a length shorter than 30 would be discarded
# (2) sequences from liftOver with internal stop codon would be discarded
# (3) The priority was given to liftOver when the predictions from two methods overlapped.
# overlapped rules: coverage>0.6

from Bio import SeqIO
import os
import numpy as np

class Gene:
    def __init__(self,strain,gene_id,contig,start,end,seq,predictor,gff_record):
        self.strain = strain
        self.gene_id = gene_id
        self.contig = contig
        self.start = start
        self.end = end
        self.seq = seq
        self.predictor = predictor
        self.gff=gff_record

def load_genes(gff,fasta,typ):
    # load those genes without internal stop codon and longer than 30
    # gff: predicted gff file
    # fasta: protein sequences, including sequences shorter than 30 and with internal
    #        stop codon
    # typ: 'lift', 'aug'

    # if a transcript has two or more sequences
	#  1> discard all sequences shorter than 30
	#  2> discard all sequences with internal stop codon
	#  3> if there are still two or more sequences left, keep all of them
    strain = gff[:-4].split('/')[-1]
    print(strain)
    seqs = dict()
    for rec in SeqIO.parse(fasta,'fasta'):
        # discard sequences with internal stop condon or shorter than 30
        if len(rec.seq) <= 30: continue
        if '*' in rec.seq: continue
        seqs[rec.id] = seqs.get(rec.id, []) + [rec.seq]

    genes = list()


    # lift: find mRNA line and ID
    if typ=="lift":
        for line in open(gff):
            if line.startswith('#'): continue
            gff_record=line
            cont = line.strip().split('\t')


            if cont[2] != 'mRNA' and cont[2]!='transcript': continue
            if cont[2] == 'mRNA':
                predictor = 'liftOver'
                contig = cont[0]
                start, end = int(cont[3]), int(cont[4])

                for item in cont[8].split(';'):
                    if item.startswith('ID='): gene_id = item.replace('ID=','')
            try:
                seq = seqs[gene_id]
            except:
                print(gene_id, 'has no good sequence in fasta file')
                continue
            g = Gene(strain, gene_id, contig, start, end, seq, predictor,gff_record)
            genes.append(g)

    if typ=="maker":
        for line in open(gff):
            if line.startswith('#'): continue
            gff_record=line
            cont = line.strip().split('\t')
            try:
                cont[8]
            except:
                continue
            # if cont[1] not in ['maker','blastn','augustus','est']:continue
            if cont[2] == 'mRNA' or cont[2]=='match' or cont[2]=='expressed_sequence_match' or cont[2]=='protein_match':
                predictor = 'MAKER2'
                contig = cont[0]
                start, end = int(cont[3]), int(cont[4])
                gene_id = cont[8]
                gene_id=gene_id.split(";")[0]
                gene_id=gene_id.replace("ID=","")
                try:
                    seq = seqs[gene_id]
                except:
                    print(gene_id, 'has no good sequence in fasta file')
                    continue
                g = Gene(strain,gene_id,contig,start,end,seq,predictor,gff_record)
                genes.append(g)
    return genes


class genePreds:
    def __init__(self,lift_genes,maker_genes):
        self.lift_genes = lift_genes
        self.maker_genes = maker_genes
        self.strain = lift_genes[0].strain

    def combine(self):
        self.combined_genes = list(tuple(self.lift_genes))
        self.dual_genes = dict()
        for g_maker in self.maker_genes:
            for g_lift in self.lift_genes:
                if g_maker.contig != g_lift.contig: continue
                if np.min([g_maker.end, g_lift.end]) > np.max([g_maker.start,g_lift.start]):
                    g_maker_len=abs(g_maker.end-g_maker.start)
                    g_lift_len=abs(g_lift.end-g_lift.start)
                    cov=np.min([g_maker_len,g_lift_len]) / np.max([g_maker_len,g_lift_len])
                    if cov>0.6:
                        self.dual_genes[g_maker] = True
                        break

        for g_maker in self.maker_genes:
            if not self.dual_genes.get(g_maker,False):
                self.combined_genes.append(g_maker)

        self.report = '''>{}
        Total liftOver genes: {}
        Total MAKER2 genes: {}
        Overlapped genes: {}
        Combined genes: {}

        '''.format(self.strain, len(self.lift_genes),
                                len(self.maker_genes),
                                len(self.dual_genes),
                                len(self.combined_genes))
        print(self.report)
    def write_fasta(self,outname):
        # fasta format
        # >strain_name|gene_id
        fhand = open(outname, 'w')
        for g in self.combined_genes:
            if len(g.seq) < 2:
                fhand.write('>{}|{}\n{}\n'.format(g.strain,g.gene_id,g.seq[0]))
            else:
                print(g.strain, g.gene_id, 'has multi protein sequences')
                for i in range(len(g.gff)):
                    fhand.write('>{}|{}\n{}\n'.format(g.strain,g.gene_id + '_' + str(i+1),g.seq[i]))
        fhand.close()

    def write_gff(self,strain_name,outdir):
        with open(outdir + "cb_" + strain_name + "gff", "w") as f:
            for gene in self.combined_genes:
                    gff_line = gene.gff
                    gff_line = gff_line.replace("SGD", "liftover")
                    f.write(gff_line)


# combine maker & lift_over annotation result
os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database_from_13pro')
lift_gff_dir = r'E:\data\sce_ssGEMs\predicted_proteomes/flo_gff/'
lift_prot_dir = r'E:\data\sce_ssGEMs\predicted_proteomes/flo_proteomes/'
maker_gff_dir = r'E:\data\sce_ssGEMs\predicted_proteomes/maker_gff/'
maker_prot_dir = r'E:\data\sce_ssGEMs\predicted_proteomes\maker_proteomes/'

outdir = 'data/genome/predicted_allcds/combined_proteomes/'
# if not os.path.exists(outdir): os.mkdir(outdir)

report_file = 'data/genome/predicted_allcds/combine_proteomes_report.txt'

combined_list=os.listdir(outdir)
allstrain_list=os.listdir(maker_gff_dir)

with open(report_file, 'w') as f:
    for name in allstrain_list:
        # strain_id=name.replace("gff","fa")
        # if strain_id not in combined_list:
        #     i+=1
        #     print(i)
            lift_gff = os.path.join(lift_gff_dir,name[:-3]+'gff')
            lift_prot = os.path.join(lift_prot_dir,name[:-3]+'fasta')
            maker_gff = os.path.join(maker_gff_dir,name[:-3]+'gff')
            maker_prot = os.path.join(maker_prot_dir,name[:-3]+'fa')

            lift_genes = load_genes(lift_gff,lift_prot,'lift')
            maker_genes = load_genes(maker_gff,maker_prot,'maker')

            outfile = os.path.join(outdir,name[:-3]+'fa')
            gp = genePreds(lift_genes,maker_genes)
            gp.combine()
            f.write(gp.report)
            gp.write_fasta(outfile)
            # gp.write_gff(strain_name=name[:-3],outdir="data/genome/predicted_allcds/combine_gff/")   # write gff file


