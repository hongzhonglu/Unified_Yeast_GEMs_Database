# -*- coding: utf-8 -*-
# date : 2024/3/5 
# author : wangh
import omicverse as ov
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import json

def de_analysis(data,group1,group2,fc_threshold=2,pval_threshold=0.05,logp_max=10):
    '''
    This function is used to do differential analysis for fluxes
    :param df_expression: pd.DataFrame, expression data
    :param group1: list, group1 sample list
    :param group2: list, group2 sample list
    :param method: str, method of differential expression analysis
    :return: pd.DataFrame, result of differential expression analysis
    '''
    # calculate the difference of fluxes: bioethanol vs wt
    dds=ov.bulk.pyDEG(data)
    # calculate the different flux reaction
    result=dds.deg_analysis(group1,group2,method='DEseq2')
    # set the threshold
    dds.foldchange_set(fc_threshold=fc_threshold,
                       pval_threshold=pval_threshold,
                       logp_max=logp_max)

    return result


def run_pathway_enrichment_analysis(data,group1,group2,background,pathway_dict,fc_threshold=2,pval_threshold=0.05,logp_max=10):
    '''
    This function is used to run pathway enrichment analysis
    :param sig_genes: list, significant reaction list
    :param pathway_dict: dict, pathway dict
    :param background: list, background gene list
    :return: pd.DataFrame, result of pathway enrichment analysis
    '''
    # check if the group1 and group2 are in the data
    group1=[i for i in group1 if i in data.columns]
    group2=[i for i in group2 if i in data.columns]
    # run differential expression analysis
    # result,dds=run_de_analysis(data,group1,group2,method='ttest')
    result=de_analysis(data,group1,group2,
                       fc_threshold=fc_threshold,
                       pval_threshold=pval_threshold,
                       logp_max=logp_max)
    # pathway enrichment analysis
    sig_up_genes=result[result['sig']=='up'].index.tolist()
    sig_dn_genes=result[result['sig']=='down'].index.tolist()
    print('sig up:',len(sig_up_genes))
    print('sig down:',len(sig_dn_genes))

    if len(sig_up_genes)==0:
        sig_up_enr=None
    else:
        sig_up_enr=ov.bulk.geneset_enrichment(gene_list=sig_up_genes,
                                        pathways_dict=pathway_dict,
                                        pvalue_type='auto',
                                        background=background,
                                       description='None')
    if len(sig_dn_genes)==0:
        sig_dn_enr=None
    else:
        sig_dn_enr=ov.bulk.geneset_enrichment(gene_list=sig_dn_genes,
                                        pathways_dict=pathway_dict,
                                        pvalue_type='auto',
                                        background=background,
                                        description='None')

    return sig_up_enr,sig_dn_enr,result



def enr_barplot(enr,title,colunm='Adjusted P-value'):
    # ignore the upper or lower case, check if up exists in title
    if 'up' in title.lower():
        color='red'
    elif 'down' in title.lower():
        color='skyblue'
    # only keep the significant pathways(adj_pvalue<0.05)
    enr=enr[enr[colunm]<0.05]
    # sort the pathways by pvalue
    enr=enr.sort_values(by=colunm,ascending=False)
    # plot the barh plot,set Term as y axis and -log10(pvalue) as x axis
    if len(enr)==0:
        print('No significant pathways')
        return
    n=len(enr)
    fig,ax=plt.subplots(1,1,figsize=(6,5))
    # set style as white
    sns.set_style('white')
    y=enr['Term']
    x=-np.log10(enr[colunm])
    # set the bar width as 0.4
    ax.barh(y,x,color=color,height=0.4)
    ax.set_xlabel(f'-log10({colunm})')
    ax.set_title(title,weight='bold')
    # rotate the y axis label
    plt.yticks(rotation=30,fontsize=10)
    plt.tight_layout()
    plt.show()


# load expression data
df_expression=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_countMatrix.csv',index_col=0)

# load KEGG pathway dict
kegg_pathway_dict=json.load(open('data/sce_kegg_geneset.json','r'))
all_kegg_genes=[gene for genes in kegg_pathway_dict.values() for gene in genes]
backgroud=df_expression.index.intersection(all_kegg_genes).tolist()

#load strainlist
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()
human_strainList=df_strain_info[(df_strain_info['nature_clade']=='10. French Guiana human ')&(df_strain_info['type']=='Human')].index.tolist()
dairy_strainList=df_strain_info[(df_strain_info['nature_clade']=='5. French dairy ')&(df_strain_info['type']=='Fermentation')].index.tolist()

wt_strainList=[strain for strain in wt_strainList if strain in df_expression.columns]
bioethanol_strainList=[strain for strain in bioethanol_strainList if strain in df_expression.columns]
human_strainList=[strain for strain in human_strainList if strain in df_expression.columns]
dairy_strainList=[strain for strain in dairy_strainList if strain in df_expression.columns]

# set all Nan values to 0
df_expression=df_expression.fillna(0)
# set as integer
df_expression=df_expression.astype(int)

# bioethanol vs wt
bioethanol_sig_up_enr,bioethanol_sig_dn_enr,bioethanol_result=run_pathway_enrichment_analysis(data=df_expression,
                                                                                                group1=bioethanol_strainList,
                                                                                                group2=wt_strainList,
                                                                                                background=backgroud,
                                                                                                pathway_dict=kegg_pathway_dict,
                                                                                                fc_threshold=2)

enr_barplot(bioethanol_sig_up_enr,'bioethanol significant up')
enr_barplot(bioethanol_sig_dn_enr,'bioethanol significant down')

# human vs wt
human_sig_up_enr,human_sig_dn_enr,human_result=run_pathway_enrichment_analysis(data=df_expression,
                                                                                group1=human_strainList,
                                                                                group2=wt_strainList,
                                                                                background=backgroud,
                                                                                pathway_dict=kegg_pathway_dict,
                                                                                fc_threshold=2)
enr_barplot(human_sig_up_enr,'human significant up')
enr_barplot(human_sig_dn_enr,'human significant down')


# dairy vs wt
dairy_sig_up_enr,dairy_sig_dn_enr,dairy_result=run_pathway_enrichment_analysis(data=df_expression,
                                                                                group1=dairy_strainList,
                                                                                group2=wt_strainList,
                                                                                background=backgroud,
                                                                                pathway_dict=kegg_pathway_dict,
                                                                                fc_threshold=2)


enr_barplot(dairy_sig_up_enr,'dairy significant up')
enr_barplot(dairy_sig_dn_enr,'dairy significant down')

# save result
bioethanol_sig_dn_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_bioethanol_sig_dn_enr.csv')
dairy_sig_dn_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_dairy_sig_dn_enr.csv')
human_sig_dn_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_human_sig_dn_enr.csv')
bioethanol_sig_up_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_bioethanol_sig_up_enr.csv')
dairy_sig_up_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_dairy_sig_up_enr.csv')
human_sig_up_enr.to_csv(r'code/7.human_diary_ethanol_analysis/output/expression_human_sig_up_enr.csv')

