# -*- coding: utf-8 -*-
# date : 2023/6/15 
# author : wangh
# file : 4.build_annotation_for_itol.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd


# load aniMatrx to get all strainIDlist
aniMatrix=pd.read_csv('code/1.genome_collection/outputs/sce1900_aniMatrix.csv',index_col=0)
strainIDlist=aniMatrix.index.tolist()

all_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)
# use genome_id column as index and choose type column as column in df_straintype
df_straintype=all_strain_info[['genome_id','sub_type']]
df_straintype.set_index('genome_id',inplace=True)

# change all 0 value as unknown
df_straintype.replace(0,'unknown',inplace=True)
df_straintype=df_straintype[df_straintype.index.isin(strainIDlist)]

# set different color for different strain type
# total strain type: Wine Wild Fermentation Human Beer unknown Sake Bakery Bioethanol Soil Industrial Distillery Dairy Insect Lab Cider
# color list: #ff0000 #ffa500 #ffff00 #008000 #00ffff #0000ff #ee82ee #800000 #808000 #008080 #000080 #4b0082 #ff69b4 #00ff7f #ffa07a #a9a9a9
df_straintype['color']=''
df_straintype.loc[df_straintype['sub_type']=='Wine','color']='#ff0000'
df_straintype.loc[df_straintype['sub_type']=='Wild','color']='#ffa500'
df_straintype.loc[df_straintype['sub_type']=='Fermentation','color']='#ffff00'
df_straintype.loc[df_straintype['sub_type']=='Human','color']='#008000'
df_straintype.loc[df_straintype['sub_type']=='Beer','color']='#00ffff'
df_straintype.loc[df_straintype['sub_type']=='unknown','color']='#0000ff'
df_straintype.loc[df_straintype['sub_type']=='Sake','color']='#ee82ee'
df_straintype.loc[df_straintype['sub_type']=='Bakery','color']='#800000'
df_straintype.loc[df_straintype['sub_type']=='Bioethanol','color']='#808000'
df_straintype.loc[df_straintype['sub_type']=='Soil','color']='#008080'
df_straintype.loc[df_straintype['sub_type']=='Industrial','color']='#000080'
df_straintype.loc[df_straintype['sub_type']=='Distillery','color']='#4b0082'
df_straintype.loc[df_straintype['sub_type']=='Dairy','color']='#ff69b4'
df_straintype.loc[df_straintype['sub_type']=='Insect','color']='#00ff7f'
df_straintype.loc[df_straintype['sub_type']=='Lab','color']='#ffa07a'
df_straintype.loc[df_straintype['sub_type']=='Cider','color']='#a9a9a9'

# change the order of columns
df_straintype=df_straintype[['color','sub_type']]
# write df_straintype behind the code/1.genome_collection/outputs/sce1900_strain_type_color_strip.txt and use space as separator
df_straintype.to_csv('code/1.genome_collection/outputs/sce1900_strain_type_color_strip.txt',sep=' ')

# add content in sce1900_strain_type_color_strip0.txt to the head of sce1900_strain_type_color_strip.txt
with open('code/1.genome_collection/outputs/sce1900_strain_type_color_strip0.txt','r') as f:
    content=f.read()
    with open('code/1.genome_collection/outputs/sce1900_strain_type_color_strip.txt','r+') as f1:
        original=f1.read()
        f1.seek(0)   # move the cursor to the head of file
        f1.write(content+'\n'+original)



