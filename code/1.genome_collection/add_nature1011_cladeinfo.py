# -*- coding: utf-8 -*-
# date : 2023/6/27 
# author : wangh
# file : add_nature1011_cladeinfo.py
# project : Unified_Yeast_GEMs_Database
'''add nature1011 calde information into all strain information'''

import pandas as pd
import os

# load data
all_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)

# load nature1011 clade information
nature1011_strain_info=pd.read_excel('data/strain_information/nature1011_strain_info.xls',sheet_name='Table S1',skiprows=3)

# all strain info set ssGEM as index
all_strain_info=all_strain_info.set_index('ssGEM')

# nature1011 strain info set ssGEM as index
nature1011_strain_info=nature1011_strain_info.set_index('Standardized name')
# remove SACE_ prefix
nature1011_strain_info.index=nature1011_strain_info.index.str.replace('SACE_','')

# check does all nature1011 strains in all_strain_info
len(set(nature1011_strain_info.index)-set(all_strain_info.index))  # 0

# add nature1011 clade information into all strain information
all_strain_info['nature_clade']=all_strain_info.index.map(lambda x:nature1011_strain_info.loc[x,'Clades'] if x in nature1011_strain_info.index else 'nan')

# save data
all_strain_info.to_excel('data/1897_strains_info.xlsx')