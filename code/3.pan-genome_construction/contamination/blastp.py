# -*- coding: utf-8 -*-
# date : 2023/5/17 
# author : wangh
# file : blastp.py
# project : Unified_Yeast_GEMs_Database
import os
import sys
import time
from Bio import SeqIO
import multiprocess


def run_mmseqs_blastp(gene,blastp_output_dir):

    geneID=gene.id.replace("|","_")
    seq = str(gene.seq)
    if geneID+'.txt' in os.listdir(blastp_output_dir):
        print('This gene has been blasted:', geneID)
        return
    print('This is gene:', geneID)
    print('This is seq:', seq)
    with open('tmp/%s.fasta'%geneID, 'w+') as f:
        f.write('>' + geneID + '\n'+seq+'\n')
    if not os.path.exists(blastp_output_dir):
        os.system('mkdir %s' % blastp_output_dir)
    # Execute shell command line in Python using the os.system function
    cmd1 = 'mmseqs createdb tmp/%s.fasta tmp/%sDB' % (geneID, geneID)
    # use ncbi nr database to blastp
    cmd2='mmseqs search tmp/%sDB /dssg/home/acct-clslhz/clslhz/why_ssGEM/data/mmseqs2/database/ncbi_nr tmp/%s_resultDB tmp/%s_tmp --max-seqs 250 --min-seq-id 0.8 -s 1.0 --threads 120 --remove-tmp-files 1' % (geneID, geneID,geneID)
    cmd3 = 'mmseqs convertalis tmp/%sDB /dssg/home/acct-clslhz/clslhz/why_ssGEM/data/mmseqs2/database/ncbi_nr tmp/%s_resultDB %s/%s.txt --format-output \"query,target,pident,bits,qcov,tcov,qlen,tlen,evalue,taxid,taxname,taxlineage\"' % (geneID, geneID,blastp_output_dir, geneID)
    # use swissprot database to blastp
    # cmd2 = 'mmseqs search %sDB /dssg/home/acct-clslhz/clslhz/why_ssGEM/data/mmseqs2/database/swissprot %s_resultDB tmp --max-seqs 250 --min-seq-id 0.8 -s 4.0 --threads 120' % (geneID, geneID)
    # cmd3 = 'mmseqs convertalis %sDB /dssg/home/acct-clslhz/clslhz/why_ssGEM/data/mmseqs2/database/swissprot %s_resultDB %s/%s.txt --format-output \"query,target,pident,bits,qcov,tcov,qlen,tlen,evalue,taxid,taxname,taxlineage\"' % (geneID, geneID,blastp_output_dir, geneID)
    rmcmd1='rm tmp/%s.fasta' % geneID
    rmcmd2='rm tmp/%sDB*' % geneID
    rmcmd3='rm tmp/%s_resultDB*' % geneID
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(rmcmd1)
    os.system(rmcmd2)
    os.system(rmcmd3)

def main() :
    start_time = time.time()
    blastp_output_dir="./output/nr_blastp_files"
    fasta_file_index=sys.argv.index('-i')+1
    fasta_file = sys.argv[fasta_file_index]
    records=[record for record in SeqIO.parse(fasta_file,'fasta')]
    for gene in records:
        run_mmseqs_blastp(gene,blastp_output_dir)
        # create the thread pool include 10 threads
    # arg_list=[(gene,blastp_output_dir) for gene in records]
    # with multiprocess.Pool(processes=5) as pool:
    #     # 并行运算
    #     results = pool.starmap(run_mmseqs_blastp,arg_list)
    elapsed_time = time.time() - start_time

    print("Finished------------------")
    print("Running the Blast process: %.4fs" %(elapsed_time))

if __name__ == "__main__" :
    main()

