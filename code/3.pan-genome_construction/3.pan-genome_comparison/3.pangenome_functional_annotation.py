# -*- coding: utf-8 -*-
# date : 2024/10/14 
# author : wangh
# file : 3.pangenome_functional_annotation.py
# project : Unified_Yeast_GEMs_Database
'''Do functional annotation for pan1807 by Uniprot and PLMsearch'''
import pandas as pd
import requests
from io import StringIO
from Bio import SeqIO

def fetch_uniprot_annotations(uniprot_ids):
    url = "https://rest.uniprot.org/uniprotkb/search"
    # Join multiple UniProt IDs with spaces or "+" for URL encoding
    query = " OR ".join(uniprot_ids)

    params = {
        'query': query,
        'fields': 'accession,id,gene_oln,protein_name,cc_catalytic_activity,ec,cc_function,cc_pathway,rhea,go_p,go_c,go,go_f,cc_subcellular_location,organism_name,go_id,cc_function,lineage',  # Fields you want to retrieve. https://www.uniprot.org/help/return_fields
        'format': 'tsv',  # You can choose 'tsv', 'xml', or 'json'
        'size': len(uniprot_ids)  # Limit the result size to the number of IDs you are querying
    }

    response = requests.get(url, params=params)

    if params['format'] == 'tsv':
        response = requests.get(url, params=params)

        if response.status_code == 200:
            # Convert the response content to a Pandas DataFrame
            tsv_data = StringIO(response.text)
            df = pd.read_csv(tsv_data, sep='\t')
            return df
        else:
            print(f"Error: {response.status_code}")
            return None

    elif params['format'] == 'json':
        if response.status_code == 200:
            data = response.json()['results']
            return data
        else:
            print(f"Error: {response.status_code}")
            return None


panID_list=[g.id for g in SeqIO.parse(r'data/genome/pan1800.fasta','fasta')]


# load s288c Uniprot annotation info
df_s288c=pd.read_excel(r'data/uniprotkb_s288c_annotation.xlsx')
df_s288c.index=df_s288c['Gene Names (ordered locus)']
s288c_unknown_geneIDlist=[i for i in panID_list if i not in df_s288c.index and 'sce' not in i]


# load plmsearch result for new ORFs
# thrshold=0.5
df_predict=pd.read_csv(r'code/3.pan-genome_construction/3.pan-genome_comparison/output/plmsearch_output.csv',index_col=0,sep='\t',header=None)
df_predict.columns=['uniprot_id','similarity']
# df_predict=df_predict[df_predict['similarity']>thrshold]
pre_geneIDlist=list(df_predict['uniprot_id'])

fetch_geneIDlist=s288c_unknown_geneIDlist+pre_geneIDlist
# split to 500 item each time:
batch_size=200
df_pre_annotation=pd.DataFrame()
for i in range(0,len(fetch_geneIDlist)//batch_size+1):
    if i < len(fetch_geneIDlist)//batch_size:
        df_pre_annotation_sub=fetch_uniprot_annotations(fetch_geneIDlist[i*batch_size:(i+1)*batch_size])
    else:
        df_pre_annotation_sub=fetch_uniprot_annotations(fetch_geneIDlist[i*batch_size:])
    if df_pre_annotation_sub is not None:
        df_pre_annotation=pd.concat([df_pre_annotation,df_pre_annotation_sub])
    else:
        print(f'{i*batch_size}:{(i+1)*batch_size} error')

df_pre_annotation.index=df_pre_annotation['Entry']
# remove duliacted index
df_pre_annotation=df_pre_annotation[~df_pre_annotation.index.duplicated(keep='first')]

# map to df_predict, according to uniprot_id
df_predict_annotation=pd.merge(df_predict,df_pre_annotation,left_on='uniprot_id',right_index=True,how='left')

# concatanate with s288c annotation
df_combined=pd.concat([df_s288c,df_predict_annotation],sort=False)
# add rows: 'YIL080W', 'YPR080W', 'YAR029W', 'YOR133W', 'YKR059W', 'YNL030W', 'YNL031C'
df_add=pd.DataFrame(index=['YIL080W', 'YPR080W', 'YAR029W', 'YOR133W', 'YKR059W', 'YNL030W', 'YNL031C'])
df_combined=pd.concat([df_combined,df_add],sort=False)

df_pangenome_annotation=df_combined.loc[panID_list]
# save result
df_pangenome_annotation.to_excel(r'data/genome/scepan1807_functional_annotation_uniprot.xlsx')