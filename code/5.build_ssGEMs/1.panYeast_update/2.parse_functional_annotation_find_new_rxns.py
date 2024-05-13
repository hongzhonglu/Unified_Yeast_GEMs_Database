# -*- coding: utf-8 -*-
# date : 2023/5/29 
# author : wangh
# file : 2.parse_functional_annotation_find_new_rxns.py
# project : Unified_Yeast_GEMs_Database
'''combine eggnog and kegg functional annotation result to find potential new reactions to update panYeast
'''
import pandas as pd
from Bio import SeqIO

#load s288c geneID list
s288c_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/S288c_R64.fasta","fasta")]
# load rename gene list
df_rename=pd.read_excel(r"code/5.build_ssGEMs/output/panYeast_rename_panID.xlsx",index_col=0)
already_geneList=s288c_geneIDlist+list(df_rename['panID'])

# check eggnog & KEGG annotation result
df_eggnog=pd.read_excel(r"data/genome/pan1900_new_v4_50_70_functional_annotations.xlsx",index_col=0)
df_eggnog_add=df_eggnog[~df_eggnog.index.isin(already_geneList)]

# load kegg annotation.and set koID as column
df_kegg=pd.read_excel(r"data/genome/pan1900_new_v4_50_70_KEGG_functional_annotations.xlsx",index_col=0)
df_kegg.dropna(inplace=True)
df_kegg_add=df_kegg[~df_kegg.index.isin(already_geneList)]

# combine eggnog and kegg annotation result to add geneList:
df_eggnog_add_ko=df_eggnog_add['KEGG_ko'].to_frame()
df_eggnog_add_ko.rename(columns={"KEGG_ko":"ko"},inplace=True)
# remove all 'ko:' in ko column
df_eggnog_add_ko['ko']=df_eggnog_add_ko['ko'].str.replace("ko:","")
df_eggnog_add_ko=df_eggnog_add_ko[~(df_eggnog_add_ko['ko']=="-")]
df_eggnog_add_ko.dropna(inplace=True)
#combine df_eggnog_add_ko and df_kegg_add
df_eggnog_kegg_add=df_eggnog_add_ko.copy()
for id in df_kegg_add.index:
    if id not in df_eggnog_kegg_add.index:
        df_eggnog_kegg_add.loc[id]=df_kegg_add.loc[id]
    else:
        kegg_result=df_kegg_add.loc[id].values[0]
        eggnog_result=df_eggnog_add_ko.loc[id].values[0].replace("ko:","")
        if kegg_result==eggnog_result:
            eggnog_result=df_eggnog_add_ko.loc[id].values[0].replace("ko:","")
            df_eggnog_kegg_add.loc[id]=eggnog_result
        else:
            df_eggnog_kegg_add.loc[id]=kegg_result


# use biopython to get reaction acording to koID
from Bio.KEGG import REST


link=REST.kegg_link("ko","reaction").read()
koList=[]
rnList=[]
for item in link.split("\n"):
    if item.startswith("rn"):
        ko=item.split("\t")[1].strip("ko:")
        rn=item.split("\t")[0].strip("rn:")
        koList.append(ko)
        rnList.append(rn)
df_reaction=pd.DataFrame({"ko":koList,"rn":rnList})


add_rxn_list=[]
geneIDlist=[]
for id in df_eggnog_kegg_add.index:
    koIDlist=df_eggnog_kegg_add.loc[id].values[0].split(",")
    for koID in koIDlist:
        if koID in df_reaction['ko'].tolist():
            add_rxn_list.append(df_reaction[df_reaction['ko']==koID]['rn'].values[0])
            geneIDlist.append(id)
        else:
            pass
df_add_rxn=pd.DataFrame({"geneID":geneIDlist,"rxn":add_rxn_list})

# check reaction information in KEGG according to reaction ID
from Bio.KEGG import REST
expressionList=[]
rxnNameList=[]
for id in df_add_rxn['rxn'].tolist():
    #sleep 1 second per request
    # time.sleep(2)
    try:
        record=REST.kegg_get(id).read()
    except:
        print(id+' not found')
        expressionList.append("NaN")
        rxnNameList.append("NaN")
        # break this loop and continue next loop
        continue
    rxnName = 'NaN'
    expression = 'NaN'
    for item in record.split("\n"):
        if item.startswith("NAME"):
            rxnName=item.split("NAME")[1].strip()
        elif item.startswith("DEFINITION"):
            expression=item.split("DEFINITION")[1].strip()
        else:
            continue
    rxnNameList.append(rxnName)
    expressionList.append(expression)

df_add_rxn['rxnName']=rxnNameList
df_add_rxn['equation']=expressionList

# remove weather incomplete or collapsed genes exist in df_add_rxn
df_add_rxn=pd.read_excel("code/5.build_ssGEMs/output/panYeast_add_gene&reaction.xlsx")
pangenomeList=[g.id for g in SeqIO.parse("data/genome/pan1900_new_v4_50_70_filter.fasta","fasta")]
df_add_rxn=df_add_rxn[df_add_rxn['geneID'].isin(pangenomeList)]

# save result: potential genes and reactions to add
df_add_rxn.to_excel("code/5.build_ssGEMs/output/panYeast_add_gene&reaction.xlsx",index=False)

