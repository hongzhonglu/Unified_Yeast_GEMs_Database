# -*- coding: utf-8 -*-
# date : 2024/2/26 
# author : wangh

import pandas as pd
import cobra

# load data
df_human=pd.read_excel('code/human_diary_analysis/output/human.xlsx',index_col=1,sheet_name='Sheet1')
df_human=df_human[df_human['Category']=='KEGG_PATHWAY']

df_check=df_human.iloc[:13,:]

model=cobra.io.read_sbml_model(r'model/yeast-GEM.xml')
model_geneList=[gene.id for gene in model.genes]

df_check_model=pd.DataFrame(index=df_check.index,columns=['Genes','Model_geneNumb','Lost_geneNumb'])
for pathway in df_check.index:
    geneList=df_check.loc[pathway,'Genes'].split(',')
    geneList=[gene.strip() for gene in geneList]
    model_geneNumb=len(set(geneList).intersection(model_geneList))
    lost_geneNumb=len(set(geneList).difference(model_geneList))
    df_check_model.loc[pathway,:]=[geneList,model_geneNumb,lost_geneNumb]

