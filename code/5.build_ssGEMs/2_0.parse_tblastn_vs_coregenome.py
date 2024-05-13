# -*- coding: utf-8 -*-

'''parse tblastn result for nature1011 pangenome core genes in panYeast,to confirm the existence of these core genes in each strain,
if exist, it will be add wo the ssGEM while the building ssGEM process.
'''
import pandas as pd
from cobra.io import read_sbml_model
import sys
import os
sys.path.append("code")
from mainFunction import get_gene_lens
from tqdm import tqdm

def parse_tblastn_result(tblastn_file,tblastn_dir,query,query_dir):
    """
    parse the tblastn result file and calculate COV for query and subject
    :param tblastn_file: tblastn result file
    :param query: query file name
    :param query_dir: query file directory
    :param subject: subject file name
    :param subject_dir: subject file directory
    :return: tblastn_file with COV for query and subject
    """
    # load tblastn result file
    tblastn_file = pd.read_csv(
        tblastn_dir+tblastn_file, sep="\t",
        index_col=0, header=None)
    columns = ["subject", "identity", "alignment length", "mismatches", "gap opens", "q_start", "q_end", "s_start",
               "s_end", "evalue", "bit score"]
    tblastn_file.columns = columns

    # get gene length for query and subject
    query_lens=get_gene_lens(query,in_folder=query_dir)
    query_lens.set_index("gene",inplace=True)

    # map query lens to tblastn_file and name the column as "query_lens"
    tblastn_file=tblastn_file.join(query_lens,how="left")
    tblastn_file.rename(columns={"gene_length":"query_lens"},inplace=True)

    # calculate COV for query and subject
    tblastn_file["query_cov"]=(tblastn_file["q_end"]-tblastn_file["q_start"]+1)/tblastn_file["query_lens"]
    return tblastn_file


def build_tblastn_geneMatrix(strainList,geneList,pangenome,tblastn_file_dir,cov_cutoff=0.6,pid_cutoff=0.8):
    """
    build gene matrix for each strain
    :param strainList: strain list
    :param geneList: gene list
    :param tblastn_file_dir: tblastn file directory
    :return: gene matrix for each strain
    """
    cnvMatrix=pd.DataFrame(index=geneList,columns=strainList)
    for strain in tqdm(strainList):
        tblastn_file=parse_tblastn_result(tblastn_file=pangenome.strip('.fasta')+"_vs_"+strain.strip('.fna')+"_tblastn.txt",
                                          tblastn_dir=tblastn_file_dir,
                                          query=pangenome,
                                          query_dir="data/genome/")
        tblastn_file=tblastn_file[tblastn_file.index.isin(geneList)]
        # keep more similiar genes with high coverage and high identity or incomplete genes with very high coverage
        tblastn_file=tblastn_file[((tblastn_file["query_cov"]>=cov_cutoff)&(tblastn_file["identity"]>=pid_cutoff))]
        # tblastn_file=tblastn_file[((tblastn_file["query_cov"]>=cov_cutoff)&(tblastn_file["identity"]>=pid_cutoff))|((tblastn_file["identity"]>=0.95)&(tblastn_file["query_cov"]>0.15))]
        tblastn_file['query']=tblastn_file.index
        tblastn_file.reset_index(inplace=True,drop=True)
        df_tblastn_group=tblastn_file.groupby("query")
        df_tblastn_gene_count=df_tblastn_group.count()
        cnvMatrix.loc[df_tblastn_gene_count.index,strain]=df_tblastn_gene_count['subject']
    cnvMatrix.fillna(0, inplace=True)
    cnvMatrix = cnvMatrix.astype(int)
    geneMatrix = cnvMatrix.copy()
    geneMatrix[geneMatrix > 0] = 1

    return geneMatrix,cnvMatrix


# load nacore100_panYeast_genes
df_nacore100=pd.read_excel("result/na1011_coregene.xlsx",index_col=0)
nacore100_geneList=df_nacore100[df_nacore100["core100"]==1].index.tolist()

panYeast_geneList=[g.id for g in read_sbml_model('model/panYeast_v4_5.xml').genes]
nacore100_panYeast_genes=[g for g in nacore100_geneList if g in panYeast_geneList]

etc_rxnList=['r_0226','r_0438','r_0439']
etc_geneList=[g.id for r in etc_rxnList for g in read_sbml_model('model/panYeast_v4_5.xml').reactions.get_by_id(r).genes]
etc_geneList=list(set(etc_geneList))

strainList=os.listdir(r"E:\data\sce_ssGEMs\genomes\assembled_sequence\1900_assembled_genome")
tblastn_dir=r"code/3.pan-genome_construction/2.build_geneMatrix/tblastn_geneMatrix/output/S288c_R64_tblastn_vs_1800/"
pangenome="S288c_R64.fasta"

all_strain_info=pd.read_excel("data/1897_strains_info.xlsx", index_col=0)
kept_strainList=all_strain_info[all_strain_info["remove"]==False]["genome_id"].tolist()
kept_strainList.append('s288c_R64')
kept_strainList=[s+".fna" for s in kept_strainList]
strainList=[s for s in strainList if s in kept_strainList]
# remove _Saccharomyces_cerevisiae_ in strain name to shorten the strain name
strainList=[s.replace("_Saccharomyces_cerevisiae_","_") for s in strainList]
# remove _genomic_
# strainList=[s.replace("_genomic","") for s in strainList]
tblastn_geneMatrix,tblastn_cnvMatrix=build_tblastn_geneMatrix(strainList=strainList,
                                                              geneList=nacore100_panYeast_genes,
                                                              pangenome=pangenome,
                                                                tblastn_file_dir=tblastn_dir,
                                                                cov_cutoff=0.5,
                                                                pid_cutoff=0.7)

etc_tblastn_geneMatrix,etc_tblastn_cnvMatrix=build_tblastn_geneMatrix(strainList=strainList,
                                                                geneList=etc_geneList,
                                                                pangenome=pangenome,
                                                                tblastn_file_dir=tblastn_dir,
                                                                cov_cutoff=0.5,
                                                                pid_cutoff=0.7)


# some genes in mitochondria are lost because of the worse quality for mitochondria genome sequencing.
# change all values to 1 if the index starts with 'Q'
etc_tblastn_geneMatrix.loc[etc_tblastn_geneMatrix.index.str.startswith('Q'),:]=1
etc_tblastn_cnvMatrix.loc[etc_tblastn_cnvMatrix.index.str.startswith('Q'),:]=1

tblastn_geneMatrix=pd.concat([tblastn_geneMatrix,etc_tblastn_geneMatrix],axis=0)
tblastn_cnvMatrix=pd.concat([tblastn_cnvMatrix,etc_tblastn_cnvMatrix],axis=0)

tblastn_cnvMatrix.columns=strainList
tblastn_geneMatrix.columns=strainList

# check does the index have duplicated genes
len(tblastn_geneMatrix.index.unique())
# remove duplicated genes
tblastn_geneMatrix=tblastn_geneMatrix[~tblastn_geneMatrix.index.duplicated(keep='first')]
tblastn_cnvMatrix=tblastn_cnvMatrix[~tblastn_cnvMatrix.index.duplicated(keep='first')]

tblastn_cnvMatrix.to_csv("code/5.build_ssGEMs/output/coregene_tblastn4_cnvMatrix.csv")
tblastn_geneMatrix.to_csv("code/5.build_ssGEMs/output/coregene_tblastn4_geneMatrix.csv")

check_geneList=['YFL058W','YGL012W','YHR020W','YNL247W','YPL160W','YDL141W','YML126C','YAR015W','YDR354W','YFR025C','YIL020C','YHR072W','YLL018C']
gene_ratio=tblastn_geneMatrix.sum(axis=1)/len(strainList)
df_check_gene_ratio=gene_ratio[check_geneList]
len(gene_ratio[gene_ratio>0.9])


tblastn_geneMatrix=pd.read_csv("code/5.build_ssGEMs/output/nacore100_tblastn2_geneMatrix.csv",index_col=0)
test_strains=["GCA_019394525.1_ASM1939452v1_genomic.fna","GCA_019394085.1_ASM1939408v1_genomic.fna","GCA_019394815.1_ASM1939481v1_genomic.fna"]
df_check=tblastn_geneMatrix.loc[:,test_strains]

# check the strain that lost 'YFL058W'
df_lost=tblastn_geneMatrix.loc['YFL058W',tblastn_geneMatrix.loc['YFL058W',:]==0]



