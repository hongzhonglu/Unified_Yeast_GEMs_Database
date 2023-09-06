# -*- coding: utf-8 -*-
# date : 2023/7/17 
# author : wangh
# file : check_growth_data.py
# project : Unified_Yeast_GEMs_Database
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# check liquid culture growth data
df_growth1=pd.read_csv('data/transcriptomics/969_liquid_culture.csv',index_col=0)


# check the t_mid_norm distribution
fig,axes=plt.subplots(1,2,figsize=(10,5))
axis_fontdict={'family':'Times New Roman','size':12,'weight':'bold'}
df_growth1['t_mid_norm'].hist(bins=100,ax=axes[0])
axes[0].set_xlabel('t_mid_norm',fontdict=axis_fontdict)
axes[0].set_ylabel('strain count',fontdict=axis_fontdict)
axes[0].set_title('969 strains 96-well liquid cultivation',fontsize=15,weight='bold')


# colony size evaluation
df_growth_0=pd.read_csv('data/transcriptomics/raw_969_colony_size.txt',sep='\t')

# only keep Stardardized_name as index and YPD_40h as column
df_growth2=df_growth_0[['Standardized_name','YPD_40h']]

# remove

# check the distribution of YPD_40h
df_growth2['YPD_40h'].hist(bins=100,ax=axes[1])
axes[1].set_xlabel('Colony size in YPD at 40h',fontdict=axis_fontdict)
axes[1].set_ylabel('strain count',fontdict=axis_fontdict)
axes[1].set_title('969 strains colony size in YPD at 40h',fontsize=15,weight='bold')
plt.tight_layout()
plt.show()

# normalized the df_growth1 and df_growth2
df_growth1_norm=df_growth1.copy()
df_growth2_norm=df_growth2.copy()

# normalized the df_growth1 by min-max normalization
df_growth1_norm['t_mid_norm']=(df_growth1_norm['t_mid_norm']-df_growth1_norm['t_mid_norm'].min())/(df_growth1_norm['t_mid_norm'].max()-df_growth1_norm['t_mid_norm'].min())

# normalized the df_growth2 by min-max normalization
df_growth2_norm['YPD_40h']=(df_growth2_norm['YPD_40h']-df_growth2_norm['YPD_40h'].min())/(df_growth2_norm['YPD_40h'].max()-df_growth2_norm['YPD_40h'].min())


# plot the kde of normalized data of df_growth1 and df_growth2
fig,ax=plt.subplots(figsize=(10,5))
sns.set()
axis_fontdict={'family':'Times New Roman','size':12,'weight':'bold'}
sns.kdeplot(df_growth1_norm['t_mid_norm'],ax=ax,label='96-well SC liquid culture',shade=True)
sns.kdeplot(df_growth2_norm['YPD_40h'],ax=ax,label='YPD_40h colony size',shade=True)
ax.set_xlabel('Normalized value',fontdict=axis_fontdict)
ax.set_ylabel('Density',fontdict=axis_fontdict)
ax.set_title('Normalized distribution of 969 strains growth data',fontsize=15,weight='bold')
plt.tight_layout()
plt.legend()
plt.show()


# map strain ID in grwoth data to genome ID
# load strain info
df_all_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)

# map strain ID in grwoth data to genome ID
#only keep strain that index include in df_all_strain_info
# remove SACE_ in index
df_growth1.index=df_growth1.index.str.replace('SACE_','')
df_growth1=df_growth1.loc[df_growth1.index.isin(df_all_strain_info.index.tolist()),:]
df_growth1['genome_id']=df_growth1.index.map(lambda x:df_all_strain_info.loc[x,'genome_id'])

# set Standardized_name as index for df_growth2
df_growth2.set_index('Standardized_name',inplace=True)
#remove SACE_ in index
df_growth2.index=df_growth2.index.str.replace('SACE_','')
df_growth2=df_growth2.loc[df_growth2.index.isin(df_all_strain_info.index.tolist()),:]
df_growth2['genome_id']=df_growth2.index.map(lambda x:df_all_strain_info.loc[x,'genome_id'])

# combine df_growth1 and df_growth2
df_growth=pd.concat([df_growth1,df_growth2],axis=1)
# merge duplicated genome_id columns
df_growth=df_growth.loc[:,~df_growth.columns.duplicated()]

# rename YPD_40 to YPD_colony_size
df_growth.rename(columns={'YPD_40h':'YPD_colony_size'},inplace=True)
# save df_growth
df_growth.to_csv('data/transcriptomics/969growth_data_t_mid&colonysize.csv')
