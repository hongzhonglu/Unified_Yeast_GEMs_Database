# -*- coding: utf-8 -*-
# date : 2023/5/6 
# author : wangh
# file : 2.extract_nonref_genes.py
# project : Unified_Yeast_GEMs_Database
import os

import pandas as pd
import sys
from Bio import SeqIO
sys.path.append(r"D:\code\github\Unified_Yeast_GEMs_Database\code")
from mainFunction import get_gene_lens

def parse_blastp_result(blastp_file,blastp_dir,query,query_dir,subject,subject_dir):
    """
    parse the blastp result file and calculate COV for query and subject
    :param blastp_file: blastp result file
    :param query: query file name
    :param query_dir: query file directory
    :param subject: subject file name
    :param subject_dir: subject file directory
    :return: blastp_file with COV for query and subject
    """
    # load blastp result file
    df_blastp_file = pd.read_csv(blastp_dir+blastp_file, sep="\t", header=None, index_col=0)
    columns = ["subject", "identity", "alignment length", "mismatches", "gap opens", "q_start", "q_end", "s_start",
               "s_end", "evalue", "bit score"]
    df_blastp_file.columns = columns

    # get gene length for query and subject
    query_lens=get_gene_lens(query,in_folder=query_dir)
    subject_lens=get_gene_lens(subject,in_folder=subject_dir)
    query_lens.set_index("gene",inplace=True)
    subject_lens.set_index("gene",inplace=True)

    # map query lens to blastp_file and name the column as "query_lens"
    df_blastp_file=df_blastp_file.join(query_lens,how="left")
    df_blastp_file.rename(columns={"gene_length":"query_lens"},inplace=True)
    # map subject lens to blastp_file and name the column as "subject_lens"
    df_blastp_file=df_blastp_file.join(subject_lens,how="left",on="subject")
    df_blastp_file.rename(columns={"gene_length":"subject_lens"},inplace=True)

    # calculate COV for query and subject
    df_blastp_file["query_cov"]=(df_blastp_file["q_end"]-df_blastp_file["q_start"]+1)/df_blastp_file["query_lens"]
    df_blastp_file["subject_cov"]=(df_blastp_file["s_end"]-df_blastp_file["s_start"]+1)/df_blastp_file["subject_lens"]

    return df_blastp_file


def extract_nonref_genes(strainList,ref_genome,ref_dir,blastp_file_dir,proteomes_dir,cov_cutoff,pid_cutoff,output_dir,output_fasta,output_df):
    df_gene_count=pd.DataFrame(index=strainList,columns=["nonref_gene_numb","ref_gene_numb","s288c_gene_numb"])
    for strain in strainList:
        count=0
        blastp_file=parse_blastp_result(blastp_file=ref_genome.strip(".fasta")+"_vs_"+strain+"_blastp.txt",blastp_dir=blastp_file_dir,query=ref_genome,query_dir=ref_dir,subject=strain,subject_dir=proteomes_dir)
        # blastp_file=blastp_file[(blastp_file["query_cov"]>=cov_cutoff)&(blastp_file["subject_cov"]>=cov_cutoff)&(blastp_file["identity"]>=pid_cutoff)]
        blastp_file = blastp_file[(blastp_file["query_cov"] >= cov_cutoff) & (blastp_file["identity"] >= pid_cutoff)]

        hit_geneIDlist=list(set(blastp_file["subject"].tolist()))
        proteome=[gene for gene in SeqIO.parse(proteomes_dir+strain,"fasta")]
        with open(output_dir+output_fasta,'a+') as f:
            for gene in proteome:
                if gene.id not in hit_geneIDlist:
                    count+=1
                    f.write(">"+gene.id+"\n"+str(gene.seq)+"\n")
        df_gene_count.loc[strain,:]=[count,len(proteome)-count,len(set(blastp_file.index.tolist()))]
    df_gene_count.to_csv(output_dir+output_df)
    print("Done!")


if __name__ == '__main__':
    os.chdir(r"D:\code\github\Unified_Yeast_GEMs_Database")
    protomes_dir=r"data/genome/yeast_species/pep/"
    ref_genome="S288c_R64.fasta"
    ref_dir=r"data/genome/"
    blastp_file_dir=r"code/4.core_genome_analysis/compare_with_closed_yeasts/output/blastp_S288c_R64/"
    output_dir="code/4.core_genome_analysis/compare_with_closed_yeasts/output/"
    output_fasta_file="yeasts_s288c_nonref.fasta"
    output_df_file="yeasts_s288c_nonref_gene_count.csv"
    strainList=os.listdir(protomes_dir)
    extract_nonref_genes(strainList=strainList,
                         ref_genome=ref_genome,
                         ref_dir=ref_dir,
                         blastp_file_dir=blastp_file_dir,
                         proteomes_dir=protomes_dir,
                         cov_cutoff=0.5,
                         pid_cutoff=0.7,
                         output_dir=output_dir,
                         output_fasta=output_fasta_file,
                            output_df=output_df_file)
