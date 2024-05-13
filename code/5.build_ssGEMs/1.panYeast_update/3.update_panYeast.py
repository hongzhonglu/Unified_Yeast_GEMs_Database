# -*- coding: utf-8 -*-
# date : 2023/5/30 
# author : wangh
# file : 3.update_panYeast.py
# project : Unified_Yeast_GEMs_Database
'''update the panYeast according pan1900_v4 sequence and pan-genome functional annotation
1. update GPR reactions
2. remove other pan1011 genes
3. add new reactions
'''
from cobra.io import read_sbml_model,write_sbml_model
import pandas as pd
from Bio import SeqIO
from cobra.manipulation import remove_genes
from cobra.core.gene import GPR

panYeast=read_sbml_model("model/panYeast_v3.xml")
tmp_model=panYeast.copy()
tmp_model.id="panYeast_v4.5"

# 1. update GPR reactions: update 98 rxns' GPR
df_update_gpr=pd.read_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_rxns.xlsx",index_col=0)
# reset index
df_update_gpr=df_update_gpr.set_index('rxnID')

# 1.1 update pan1011 gene related GPR in model
for rxnID in df_update_gpr.index.tolist():
    rxn=tmp_model.reactions.get_by_id(rxnID)
    new_gpr_string=df_update_gpr.loc[rxnID,'new_gpr']
    gpr=GPR.from_string(new_gpr_string)
    rxn.gpr=gpr
    rxn.gpr.update_genes()
    print(rxnID+" : "+rxn.gene_reaction_rule)

tmp_model.slim_optimize()  # 1454 genes, 4042 rxns(add 82 new genes)

# 1.2 update s288c genes related GPR in model: change 4 rxns' GPR
import re
df_update_gpr_s288c=pd.read_excel("code/5.build_ssGEMs/output/panYeast_update_gpr_s288c.xlsx")
c_existpanID=0
c_newpanID=0
for i in df_update_gpr_s288c.index:
    rxnID=df_update_gpr_s288c.loc[i,'rxnID']
    old_gpr_string=tmp_model.reactions.get_by_id(rxnID).gene_reaction_rule
    panID=df_update_gpr_s288c.loc[i,'panID']
    s288cID=df_update_gpr_s288c.loc[i,'s288cID']
    if panID in old_gpr_string:
        c_existpanID+=1
        continue
    else:
        c_newpanID+=1
        if 'and' in old_gpr_string:
            # set complex pattern to find all elments like:"(A and B)" or "(C and D)"
            pattern=r'\(.*\)'  # 以(开头，以)结尾，中间任意字符，非贪婪模式
            complexes=re.findall(pattern,old_gpr_string)
            complexes_new=complexes.copy()
            for complex in complexes:
                if s288cID in complex:
                    new_complex=complex.replace(s288cID,panID)
                    complexes_new.append(new_complex)
            new_gpr_string=" or ".join(complexes_new)
        else:
            new_gpr_string=old_gpr_string+" or "+panID
    new_gpr=GPR.from_string(new_gpr_string)
    tmp_model.reactions.get_by_id(rxnID).gpr=new_gpr
    tmp_model.reactions.get_by_id(rxnID).gpr.update_genes()
    print(old_gpr_string+"***************>"+tmp_model.reactions.get_by_id(rxnID).gene_reaction_rule)
tmp_model.slim_optimize()  # 1454 genes, 4042 rxns(update 7 rxns GPR)


tmp_model2=tmp_model.copy()
# 2.remove other pan1011 non-s288c genes: remove 221 pan1011 non-s288c genes from panYeast
s288c_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/S288c_R64.fasta","fasta")]
pan1011_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/pan1011_v1.fasta","fasta")]
pan1900_geneIDlist=[i.id for i in SeqIO.parse(r"data/genome/pan1800_50_70_v2.fasta","fasta")]
pan1011_nons288c=list(set(pan1011_geneIDlist)-set(s288c_geneIDlist))
pan1011_non1900=list(set(pan1011_geneIDlist)-set(pan1900_geneIDlist))

remove_geneList=[gene.id for gene in tmp_model.genes if gene.id in pan1011_nons288c]
remove_genes(tmp_model2,remove_geneList,remove_reactions=True)
tmp_model2.slim_optimize()  # 1233 genes, 4023 rxns(221 pan1011 genes removed)

# save GPR change result: total 639 rxns' GPR have been changed
rxnIDlist=list()
rxnNamelist=list()
old_gprlist=list()
new_gprlist=list()
for rxn in tmp_model2.reactions:
    rxnID=rxn.id
    gpr=rxn.gene_reaction_rule
    if panYeast.reactions.get_by_id(rxnID).gene_reaction_rule!=gpr:
        old_gpr=panYeast.reactions.get_by_id(rxnID).gene_reaction_rule
        print(rxnID+" : "+old_gpr+"--->"+gpr)
        rxnIDlist.append(rxnID)
        rxnNamelist.append(rxn.name)
        old_gprlist.append(old_gpr)
        new_gprlist.append(gpr)

df_update_gpr_result=pd.DataFrame({'rxnID':rxnIDlist,'rxnName':rxnNamelist,'old_gpr':old_gprlist,'new_gpr':new_gprlist})
df_update_gpr_result.to_excel("code/5.build_ssGEMs/output/panYeast_v4.5_update_gpr_result.xlsx")
# write_sbml_model(tmp_model2,"model/panYeast_v4.xml")

write_sbml_model(tmp_model2,"model/panYeast_v4_5.xml")



