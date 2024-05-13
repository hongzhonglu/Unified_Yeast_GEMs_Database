'''Different flux analysis of the model'''
import omicverse as ov
import pandas as pd
from scipy import stats
import numpy as np
# import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
os.getcwd()



def de_analysis(data,group1,group2,pval_threshold=0.05,logp_max=10):
    '''
    This function is used to do differential analysis for fluxes
    :param df_expression: pd.DataFrame, expression data
    :param group1: list, group1 sample list
    :param group2: list, group2 sample list
    :param method: str, method of differential expression analysis
    :return: pd.DataFrame, result of differential expression analysis
    '''
    # calculate the difference of fluxes: bioethanol vs wt
    # dds=ov.bulk.pyDEG(data)
    # # normalize the fluxes
    # dds.normalize()
    # # calculate the different flux reaction
    # result=dds.deg_analysis(group1,group2,method='ttest')
    # # set the threshold
    # dds.foldchange_set(fc_threshold=fc_threshold,
    #                    pval_threshold=pval_threshold,
    #                    logp_max=logp_max)
    #calculate the mean flux and p-value
    result=pd.DataFrame(index=data.index,columns=['mean1','mean2','pvalue','sig'])
    for rxn in data.index:
        group1_mean=data.loc[rxn,group1].mean()
        group2_mean=data.loc[rxn,group2].mean()
        group1_rxn_fluxes=data.loc[rxn,group1]
        group2_rxn_fluxes=data.loc[rxn,group2]
        ttest,pval=stats.ttest_ind(group1_rxn_fluxes,group2_rxn_fluxes)
        if pval>pval_threshold:
            sig='Normal'
        else:
            if group1_mean>group2_mean:
                sig='up'
            else:
                sig='down'
        result.loc[rxn,:]=[group1_mean,group2_mean,pval,sig]


    print('sig up:',result[result['sig']=='up'].shape[0])
    print('sig down:',result[result['sig']=='down'].shape[0])
    return result


def run_pathway_enrichment_analysis(data,group1,group2,background,pathway_dict):
    '''
    This function is used to run pathway enrichment analysis
    :param sig_rxns: list, significant reaction list
    :param pathway_dict: dict, pathway dict
    :param background: list, background gene list
    :return: pd.DataFrame, result of pathway enrichment analysis
    '''
    # check if the group1 and group2 are in the data
    group1=[i for i in group1 if i in data.columns]
    group2=[i for i in group2 if i in data.columns]
    # run differential expression analysis
    # result,dds=run_de_analysis(data,group1,group2,method='ttest')
    result=de_analysis(data,group1,group2)
    # pathway enrichment analysis
    sig_up_rxns=result[result['sig']=='up'].index.tolist()
    sig_dn_rxns=result[result['sig']=='down'].index.tolist()

    if len(sig_up_rxns)==0:
        sig_up_enr=None
    else:
        sig_up_enr=ov.bulk.geneset_enrichment(gene_list=sig_up_rxns,
                                        pathways_dict=pathway_dict,
                                        pvalue_type='auto',
                                        background=background,
                                       description='None')
        sig_up_enr=sig_up_enr[~sig_up_enr['Term'].str.contains('Transport|Exchange')]

    if len(sig_dn_rxns)==0:
        sig_dn_enr=None
    else:
        sig_dn_enr=ov.bulk.geneset_enrichment(gene_list=sig_dn_rxns,
                                        pathways_dict=pathway_dict,
                                        pvalue_type='auto',
                                        background=background,
                                        description='None')

        # remove pathway include Transcport or Exchange
        sig_dn_enr=sig_dn_enr[~sig_dn_enr['Term'].str.contains('Transport|Exchange')]

    return sig_up_enr,sig_dn_enr,result


def pre_process_fluxes(data):
    # remove columns with 0 value in the row named 'growth'
    data=data.loc[:,(data.loc['r_2111',:]!=0)]
    #each column divide the row named r_1714 and multiply 100
    data = data.div(data.loc['r_1714', :], axis=1) * -100
    # remove columns with all NaN values
    data=data.dropna(axis=0,how='all')
    # set all Nan values to 0
    data=data.fillna(0)
    return data


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


#load strainlist
df_strain_info=pd.read_excel('data/1897_strains_info.xlsx',index_col=0)
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()
human_strainList=df_strain_info[(df_strain_info['nature_clade']=='10. French Guiana human ')&(df_strain_info['type']=='Human')].index.tolist()
dairy_strainList=df_strain_info[(df_strain_info['nature_clade']=='5. French dairy ')&(df_strain_info['type']=='Fermentation')].index.tolist()

# # only keep the strains in available_strains
# wt_strainList=[i for i in wt_strainList if i in available_strains]
# bioethanol_strainList=[i for i in bioethanol_strainList if i in available_strains]
# human_strainList=[i for i in human_strainList if i in available_strains]
# dairy_strainList=[i for i in dairy_strainList if i in available_strains]

# load pathway dictionary
import cobra
model=cobra.io.read_sbml_model('model/yeast-GEM.xml')
pathway_dict=dict()
for group in model.groups:
    if len(group.members)>3:
        pathway_dict[group.name]=[i.id for i in group.members]

all_rxnList=[i.id for i in model.reactions]

data_dir=r'code/6.transcriptomics_ssGEMs_analysis/RIPTiDe_integrate_transcriptome/output'



# file='humandairy_riptide_growth0.3_ethanolyield0.15_mean_flux.csv'
file='humandairy_riptide_growth0.7_ethanolyield0.18_mean_flux.csv'

df_fluxes = pd.read_csv(os.path.join(data_dir, file), index_col=0)
# pre-process the fluxes
df_fluxes = pre_process_fluxes(df_fluxes)
# pathway enrichment analysis
# bioethanol vs wt
bioethanol_sig_up_enr, bioethanol_sig_dn_enr, bioethanol_flux_result = run_pathway_enrichment_analysis(data=df_fluxes,
                                                                                        group1=bioethanol_strainList,
                                                                                        group2=wt_strainList,
                                                                                        background=all_rxnList,
                                                                                        pathway_dict=pathway_dict)


bioethanol_flux_result.loc['r_1761']

enr_barplot(bioethanol_sig_up_enr,'bioethanol significant up')
enr_barplot(bioethanol_sig_dn_enr,'bioethanol significant down')


# human vs wt
human_sig_up_enr, human_sig_dn_enr,human_flux_result = run_pathway_enrichment_analysis(data=df_fluxes,
                                                                                group1=human_strainList,
                                                                                group2=wt_strainList,
                                                                                background=all_rxnList,
                                                                                pathway_dict=pathway_dict)
enr_barplot(human_sig_up_enr,'human significant up')
enr_barplot(human_sig_dn_enr,'human significant down')

human_flux_result.loc['r_1761']


# dairy vs wt
dairy_sig_up_enr, dairy_sig_dn_enr,dairy_flux_result = run_pathway_enrichment_analysis(data=df_fluxes,
                                                                                group1=dairy_strainList,
                                                                                group2=wt_strainList,
                                                                                background=all_rxnList,
                                                                                pathway_dict=pathway_dict)

enr_barplot(dairy_sig_up_enr,'dairy significant up')
enr_barplot(dairy_sig_dn_enr,'dairy significant down')

dairy_flux_result.loc['r_1761']


# save enrichment result
bioethanol_sig_up_enr.to_csv('code/7.human_diary_ethanol_analysis/output/bioethanol_flux_sig_up_enr.csv')
bioethanol_sig_dn_enr.to_csv('code/7.human_diary_ethanol_analysis/output/bioethanol_flux_sig_dn_enr.csv')
dairy_sig_dn_enr.to_csv('code/7.human_diary_ethanol_analysis/output/dairy_flux_sig_dn_enr.csv')
dairy_sig_up_enr.to_csv('code/7.human_diary_ethanol_analysis/output/dairy_flux_sig_up_enr.csv')
human_sig_up_enr.to_csv('code/7.human_diary_ethanol_analysis/output/human_flux_sig_up_enr.csv')
human_sig_dn_enr.to_csv('code/7.human_diary_ethanol_analysis/output/human_flux_sig_dn_enr.csv')


# combine flux results
flux_result=pd.DataFrame(index=bioethanol_flux_result.index,columns=['bioethanol','human','dairy','wt'])
flux_pvalue=pd.DataFrame(index=bioethanol_flux_result.index,columns=['bioethanol','human','dairy','wt'])
flux_result['bioethanol']=bioethanol_flux_result['mean1']
flux_result['human']=human_flux_result['mean1']
flux_result['dairy']=dairy_flux_result['mean1']
flux_result['wt']=bioethanol_flux_result['mean2']
flux_pvalue['bioethanol']=bioethanol_flux_result['pvalue']
flux_pvalue['human']=human_flux_result['pvalue']
flux_pvalue['dairy']=dairy_flux_result['pvalue']
flux_pvalue['wt']=1

op_rxnlist=['r_0773','r_0439','r_0438','r_0226']
flux_result.loc[op_rxnlist]
flux_pvalue.loc[op_rxnlist]

# save the flux result
flux_result.to_csv('code/7.human_diary_ethanol_analysis/output/riptide_mean_flux_result.csv')
flux_pvalue.to_csv('code/7.human_diary_ethanol_analysis/output/riptide_mean_flux_pvalue.csv')

