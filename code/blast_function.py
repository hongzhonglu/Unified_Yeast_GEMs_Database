# -*- coding: utf-8 -*-
# date : 2022/10/14 
# author : wangh
# file : blast_function.py
# project : Unified_Yeast_GEMs_Database_from_13pro
'''some function used to do blast processes:tblastn,blastp'''

def tblastn(query,subject,query_path,subject_path,work_dir,e_value=0.00001):
    ''' makeblastdb + tBlastn ;
    query:protein sequence ;
    subject:DNA sequence ;
    example: tblastn(query=AAD.fa,subject=AAR.fna,query_path=)
    tips:query file命名必须为<file>.fasta, subject file命名必须为 <file>.fna
    '''
    import os

    tblastn_output_dir = work_dir + "/tblastn_result"
    if not os.path.exists(tblastn_output_dir):
        os.mkdir(tblastn_output_dir)
    query_name = query.replace(".fasta", "")
    subject_name = subject.replace(".fna", "")
    all_blast_result=os.listdir(tblastn_output_dir)
    blast_out="%s_vs_%s.blast"%(query_name,subject_name)
    if blast_out in all_blast_result:
        print("%s already have the tblast result"%query_name)
    else:
        #make blastdb for subject sequence
        mkblastdb_output_dir=work_dir+"/blastdb"
        # check blastdb directory exists or not
        if not os.path.exists(mkblastdb_output_dir):
            os.mkdir(mkblastdb_output_dir)
        mkblastdb_cmd="makeblastdb -in %s/%s -dbtype nucl -out %s/%s"%(subject_path,subject,mkblastdb_output_dir,subject_name)
        # mkblastdb_cmd="D:/biosofts/diamond-windows/diamond.exe makedb --in %s/%s -d %s/%s"%(subject_path,subject,mkblastdb_output_dir,subject_name)
        print("making blastdb for %s"%subject)
        os.system(mkblastdb_cmd)

        # run tblastn
        # check output directory exist or not,if not, create it.
        tblastn_cmd="tblastn -query %s/%s -out %s/%s_vs_%s.blast -db %s/%s -outfmt 6 -evalue %s"%(query_path,query,tblastn_output_dir,query_name,subject_name,mkblastdb_output_dir,subject_name,e_value)
        # tblastn_cmd="D:/biosofts/diamond-windows/diamond.exe tblasn"
        print("running tblastn process for %s" % subject_name)
        os.system(tblastn_cmd)
        print("%s tblastn result has been saved to %s/%s.tblastn"%(query_name,work_dir,query_name))

# test_tblastn: √
# query_path="data/genome"
# subject_path="data/genome/1900_assembled_genome"
# tblastn(query="pan1011_v2.fasta",subject="AAA_6.re.fna",query_path=query_path,subject_path="data/genome/1900_assembled_genome",work_dir="result/blast/tblastn")


# blastp process
def blastp(query,subject,query_path,subject_path,work_dir,e_value=0.00001):
    ''' makeblastdb + Blastp ;
    query:protein sequence ;
    subject:DNA sequence ;
    example: tblastn(query=AAD.fa,subject=AAR.fasta,query_path=)
    tips:query file命名必须为<file>.fasta, subject file命名必须为 <file>.fasta
    '''
    import os
    #make blastdb for subject sequence
    if subject.endswith(".fasta"):
        subject_name=subject.replace(".fasta","")
    if subject.endswith(".fa"):
        subject_name=subject.replace(".fa","")
    mkblastdb_output_dir=work_dir+"/blastdb"
    # check blastdb directory exists or not,if not,create one
    if not os.path.exists(mkblastdb_output_dir):
        os.mkdir(mkblastdb_output_dir)
    mkblastdb_cmd="makeblastdb -in %s/%s -dbtype prot -out %s/%s"%(subject_path,subject,mkblastdb_output_dir,subject_name)
    print("making blastdb for %s"%subject)
    print(mkblastdb_cmd)
    os.system(mkblastdb_cmd)

    # run blastp
    query_name = query.replace(".fasta", "")
    blastp_output_dir=work_dir+"/blastp_result"
    # check output directory exist or not,if not, create it.
    if not os.path.exists(blastp_output_dir):
        os.mkdir(blastp_output_dir)
    blastp_cmd="blastp -query %s/%s -out %s/%s_vs_%s.blastp -db %s/%s -outfmt 6 -evalue %s"%(query_path,query,blastp_output_dir,query_name,subject_name,mkblastdb_output_dir,subject_name,e_value)
    print("running blastp process for %s" % query_name)
    print(blastp_cmd)
    os.system(blastp_cmd)
    print("%s blastp result has been saved to %s/%s_vs_%s.blastp"%(query_name,work_dir,query_name,subject_name))

# test blastp：×
# query_path="data/genome"
# subject_path="data/genome/predicted_allcds/combined_proteomes"
# blastp(query="pan1011_v2.fasta",subject="AAA_6.re.fasta",query_path=query_path,subject_path=subject_path,work_dir="result/blast/blastp")
# # error: BLAST Database error: No alias or index file found for protein database
# # [result\blast\blastp\blastdb\AAA_6.re] in search path [D:\code\github\Unified_Yeast_GEMs_Database_from_13pro;;]


def parese_blast_result(result_name,result_dir,PID=80,COV=0.6,best_hit=False):
    '''parse the -outfmt 6 format blast result to a Dataframe
    parameters:
        result_name:xxx.blast
        result_dir:directory svaing the blast output result
        PID: percent of identity
        COV: ['alnLength']/['queryEnd']
        '''
    import pandas as pd
    from mainFunction import get_gene_lens

    cols = [ 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd',
            'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    print("%s file is being parsed"%result_name)
    tblastn_result= pd.read_csv("%s/%s"%(result_dir,result_name),sep="\t",names=cols,index_col=0)

    # 整合基因长度信息
    query_length=get_gene_lens("pan1011_v3","data/genome/")
    query_length.set_index("gene",inplace=True)
    tblastn_result=tblastn_result.join(query_length)
    tblastn_result["COV"]=tblastn_result['alnLength']/tblastn_result['gene_length']
    # 根据PID,COV进行初步筛选
    tblastn_result = tblastn_result[(tblastn_result['PID'] > PID) & (tblastn_result['COV'] > COV)]
    tblastn_result["gene"]=tblastn_result.index
    if best_hit==False:
        return tblastn_result
    #取所有query的唯一比对结果
    elif best_hit:
        print("getting the best hit result of %s"%result_name)
        best_hit_tblastn_result = pd.DataFrame(columns=tblastn_result.columns)
        for gene in tblastn_result.gene.unique():
            res = tblastn_result[tblastn_result.gene == gene]
            res.sort_values("PID",ascending=False,inplace=True)
            best_hit = res.iloc[0,:]
            best_hit_tblastn_result = pd.concat([best_hit_tblastn_result, pd.DataFrame(best_hit).transpose()])
        return best_hit_tblastn_result


# test parse blast result：√
# blast_result=parese_blast_result(result_name="pan1011_v2.blast",result_dir="result/blast/tblastn/tblastn_result")
# len(blast_result.gene.unique())
# best_hit_result=parese_blast_result(result_name="pan1011_v2.blast",result_dir="result/blast/tblastn/tblastn_result",best_hit=True)
#
#
