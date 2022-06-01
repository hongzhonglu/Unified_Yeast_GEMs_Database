# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：build_newstrain_geneMatrix.py
# 2022/5/27
'''build new ~600 strains gene presence/absence Matrix by lg_annotated genomes data on SJTU-HPC '''
import pandas as pd
import os
from Bio import SeqIO
from multiprocessing.pool import ThreadPool

def get_gene_lens(query, in_folder):
    from Bio import SeqIO
    file = '%s/%s.fa' % (in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []

    for record in records:
        out.append({'gene': record.name, 'gene_length': len(record.seq)})

    out = pd.DataFrame(out)
    return out


def get_bbh(query, subject, bbh_result_folder,seq_folder,subject_folder,cov):

    query_lengths = get_gene_lens(query, in_folder=seq_folder)
    subject_lengths = get_gene_lens(subject, in_folder=subject_folder)

    # Define the output file of this BLAST
    out_file="%s_vs_%s_parsed.csv"%(query,subject)
    already_bbh_results=os.listdir(bbh_result_folder)
    if out_file in already_bbh_results:
        print("%s has already got bbh parsed result"%query)
    else:
        # Combine the results of the protein BLAST into a dataframe
        print('parsing BBHs for', query, subject)
        cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd',
                    'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
        bbh = pd.read_csv('%s/%s_vs_%s.txt' % (bbh_result_folder, query, subject), sep='\t', names=cols)
        bbh = pd.merge(bbh, query_lengths)
        bbh['COV'] = bbh['alnLength'] / bbh['gene_length']

        bbh2 = pd.read_csv('%s/%s_vs_%s.txt' % (bbh_result_folder, subject, query), sep='\t', names=cols)
        bbh2 = pd.merge(bbh2, subject_lengths)
        bbh2['COV'] = bbh2['alnLength'] / bbh2['gene_length']
        out = pd.DataFrame()

        # Filter the genes based on coverage
        bbh = bbh[bbh.COV >= cov]
        bbh2 = bbh2[bbh2.COV >= cov]

        # Delineate the best hits from the BLAST
        for g in bbh.gene.unique():
            res = bbh[bbh.gene == g]
            if len(res) == 0:
                continue
            best_hit = res.loc[res.PID.idxmax()]
            best_gene = best_hit.subject
            res2 = bbh2[bbh2.gene == best_gene]
            if len(res2) == 0:
                continue
            best_hit2 = res2.loc[res2.PID.idxmax()]
            best_gene2 = best_hit2.subject
            if g == best_gene2:
                best_hit['BBH'] ="<=>"
            else:
                best_hit['BBH'] ="->"
            out = pd.concat([out, pd.DataFrame(best_hit).transpose()])

        # Save the final file to a designated CSV file
        out.to_csv(bbh_result_folder+"/"+out_file)

# combine bi-directional blast result
bbh_result_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/new600_bbh_resullt"
proteomes_dir="/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/new600_allproteomes_maker"
ref_id="pan1011_cds"
ref_dir='/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs'
strainlist=os.listdir("/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bbh/new600_blastdb")
strainlist.remove("pan1011_cds.dmnd")
# for strain in strainlist:
#     strain_name=strain.replace(".fa.dmnd","")
#     try:
#         get_bbh(query=strain_name,subject="pan1011_cds",bbh_result_folder=bbh_result_dir,seq_folder=proteomes_dir,subject_folder=ref_dir,cov=0.4)
#     except:
#         print("% can't get bbh result"%strain)
# multiprocessing to build lg_geneMatrix
# strainlist_v1=[i.replace(".fna","") for i in strainlist]
# pool_size = 4
# pool = ThreadPool(pool_size)  # 创建一个线程池
# pool.map(get_bbh(subject="pan1011_cds",bbh_result_folder=bbh_result_dir,seq_folder=proteomes_dir,subject_folder=ref_dir,cov=0.4), strainlist_v1)  # 往线程池中填线程
# pool.close()
# pool.join()

# Parse the BLAST Results into one Homology Matrix of the Reconstruction Gene
listGeneIDs=[]
records=SeqIO.parse("/lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/bio_libs/pan1011_cds.fa","fasta")
for gene in records:
    id=gene.id
    listGeneIDs.append(id)

ortho_matrix=pd.DataFrame(index=listGeneIDs,columns=strainlist)
geneIDs_matrix=pd.DataFrame(index=listGeneIDs,columns=strainlist)
error_strains=[]
for strain in strainlist:
    strain_name=strain.replace(".fa.dmnd","")
    try:
        bbh=pd.read_csv(bbh_result_dir+"/"+"%s_vs_%s_parsed.csv"%(strain_name,ref_id))
    except:
        print("cam't find %s file"%strain)
        error_strains.append(strain)
    else:
        bbh=bbh[bbh["BBH"]=="<=>"]
        listIDs = []
        listPID = []
        for r in ortho_matrix.iterrows():
            try:
                currentOrtholog = bbh[bbh['subject'] == r[0]].reset_index()
                listIDs.append(currentOrtholog.iloc[0]['subject'])
                listPID.append(currentOrtholog.iloc[0]['PID'])
                print(r[0])
            except:
                listIDs.append('None')
                listPID.append(0)
                # print(r[0])
        for col in ortho_matrix.columns:
            if col in strain:
                ortho_matrix[col] = listPID
                geneIDs_matrix[col] = listIDs

for column in ortho_matrix:
    ortho_matrix.loc[ortho_matrix[column]<=80.0,column]=0
    ortho_matrix.loc[ortho_matrix[column]>80.0,column]=1

erro_strainlist=pd.DataFrame()
erro_strainlist["error_strain"]=error_strains

ortho_matrix.to_csv('/lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/new600_geneMatrix.csv')
erro_strainlist.to_csv("/lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/new600_error_strains_for_geneGeneMatrix.csv")
# geneIDs_matrix.to_csv('test/lg1392_geneMatrix.csv')