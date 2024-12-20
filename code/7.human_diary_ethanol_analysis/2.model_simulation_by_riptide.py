# -*- coding: utf-8 -*-
# date : 2024/2/29 
# author : wangh
'''
This script is used to model simulation by RIPTiDe for human diary ethanol wildtype strains.
1. Integrate growth data by set relative growth as the fraction
2. Integrate RNAseq data
'''
import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database/code')
from cobra.io import read_sbml_model,write_sbml_model
from model_modifications import set_SCmedium
from cobra.flux_analysis import pfba,gapfill
import math
import riptide


def add_ethanol_soft_constraint(model,ethanol_yield_lb,ethanol='r_1761',gluc_uptake='r_1714'):
    ethanol_mw=46    # g/mol
    glucose_mw=180  # g/mol

    ethanol_constraint = model.problem.Constraint(model.reactions.get_by_id(ethanol).flux_expression * ethanol_mw +
                                                  model.reactions.get_by_id(gluc_uptake).flux_expression * glucose_mw * ethanol_yield_lb,
                                                  lb=0)
    model.add_cons_vars(ethanol_constraint)
    return model


# set working directory
os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database')

# set output directory
output_dir=r'code/7.human_diary_ethanol_analysis/output/'
# check if the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ethanol_lb=5
ethanol_yield_lb=0.18 # g/g glucose
growth_scale=0.7

def get_relative_growth(pow=None):
    # load growth data
    df_exp_data = pd.read_csv('data/transcriptomics/combined_969growth_data.csv', index_col=0)
    # remove rows with genome_id is NaN
    df_exp_data = df_exp_data[df_exp_data['genome_id'].notnull()]
    # set genome_id as index
    df_exp_data = df_exp_data.set_index('genome_id')
    exp_growth = df_exp_data['t_mid_norm']
    exp_growth=1/exp_growth

    # find max growth rate index
    # option 1: select the max growth rate as reference
    # max_growth_index = exp_growth.idxmax()
    # print('The fastest growing strain is %s' % max_growth_index)
    # # transform the growth rate to relative growth rate
    # relative_growth = exp_growth / exp_growth[max_growth_index]

    # option 2: select the mean of the top 5% growth rate as reference
    max_growth_index=exp_growth.nlargest(int(len(exp_growth)*0.05)).index.tolist()
    max_growth=exp_growth.loc[max_growth_index].mean()
    relative_growth=exp_growth/max_growth
    # if relative_growth>1, set it as 1
    relative_growth[relative_growth>1]=1

    if pow:
        # normalize the relative growth rate more close to 1 by power transformation
        relative_growth=relative_growth.apply(lambda x: math.pow(x,pow))

    return relative_growth,max_growth_index


def prepare_expression_for_riptide(df_expressions,strainName):
    # save strain_exp to a tsv file
    strain_expression = df_expressions[strainName]
    strain_expression.to_csv(f'growth{growth_scale}_ethanolyield{ethanol_yield_lb}_strain_expression.tsv', sep='\t', header=False)
    transcript_abundances = riptide.read_transcription_file(f'growth{growth_scale}_ethanolyield{ethanol_yield_lb}_strain_expression.tsv', header=False)

    return transcript_abundances


# load growth data
relative_growth,max_growth_index=get_relative_growth(pow=growth_scale)

# load RNA-seq data
df_expressions=pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.csv',index_col=0)

# load strain list
ssGEM_dir = 'model/ssGEMs'
df_strain_info = pd.read_excel('data/1897_strains_info.xlsx', index_col=0)
wildtypelist = ['14. CHNIII ', '20. CHN V ', '15. CHNII ', '17. Taiwanese ', '24. Asian islands ',
                '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList = df_strain_info[
    (df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type'] == 'Wild')].index.tolist()
bioethanol_strainList = df_strain_info[(df_strain_info['nature_clade'] == '3. Brazilian bioethanol ') & (
        df_strain_info['type'] == 'Industry')].index.tolist()
human_strainList = df_strain_info[(df_strain_info['nature_clade'] == '10. French Guiana human ') & (
        df_strain_info['type'] == 'Human')].index.tolist()
diary_strainList = df_strain_info[(df_strain_info['nature_clade'] == '5. French dairy ') & (
        df_strain_info['type'] == 'Fermentation')].index.tolist()

strainList = wt_strainList + bioethanol_strainList + human_strainList + diary_strainList
# only keep the strains with available ssGEM
ssGEM_List = [f.replace('.xml', '') for f in os.listdir(ssGEM_dir) if f.endswith('.xml')]
strainList = [s for s in strainList if s in ssGEM_List]
strainList = list(set(strainList).intersection(df_expressions.columns))

wt_mean_growth = relative_growth[relative_growth.index.isin(wt_strainList)].mean()
bioethanol_mean_growth = relative_growth[relative_growth.index.isin(bioethanol_strainList)].mean()
human_mean_growth = relative_growth[relative_growth.index.isin(human_strainList)].mean()
diary_mean_growth = relative_growth[relative_growth.index.isin(diary_strainList)].mean()

# load reference model
ref_model = read_sbml_model('model/panYeast.xml')
ref_model = set_SCmedium(ref_model)
max_growth = ref_model.slim_optimize()
rxnList = [rxn.id for rxn in ref_model.reactions]

# simulate each strain
df_mean_fluxes=pd.DataFrame(index=rxnList,columns=strainList)
df_median_fluxes=pd.DataFrame(index=rxnList,columns=strainList)
df_std_fluxes=pd.DataFrame(index=rxnList,columns=strainList)
concordances=dict()
for strainName in strainList:
    print('simulate strain:',strainName)
    # load model
    model = read_sbml_model(os.path.join(ssGEM_dir, strainName + '.xml'))

    # set SC medium
    model=set_SCmedium(model)

    # check if the model need to be gapfilled
    pre_growth = model.slim_optimize()

    if pre_growth < max_growth * 0.95:
        print('The model %s need to be gapfilled from %s to %s!' % (strainName, pre_growth, max_growth*0.95))
        # gapfilling
        try:
            gapfill_solution = gapfill(model, universal=ref_model, lower_bound=max_growth*0.95,
                                       demand_reactions=False)
            # add gapfilling reactions to the model
            for rxn in gapfill_solution[0]:
                model.add_reactions([rxn])
                # print('Add reaction %s to the model!' % rxn.id)
        except:
            print('Gapfilling failed!')
            continue

    # prepare transcriptome data
    transcript_abundances=prepare_expression_for_riptide(df_expressions,strainName)

    # get relative growth rate
    try:
        relative_growth_rate = relative_growth[strainName]
    except:
        if strainName in wt_strainList:
            relative_growth_rate = wt_mean_growth
        elif strainName in bioethanol_strainList:
            relative_growth_rate = bioethanol_mean_growth
        elif strainName in human_strainList:
            relative_growth_rate = human_mean_growth
        elif strainName in diary_strainList:
            relative_growth_rate = diary_mean_growth

    # set objective function
    model.objective = 'r_2111'

    # set ethanol product not 0: In batch culture, we hypothesis that at least 10 mmol/gDW/h ethanol produced when glucose uptake is 20
    # ethanol = 'r_1761'
    # model.reactions.get_by_id(ethanol).lower_bound = ethanol_lb
    # Set ethanol soft constraint make sure the ethanol product is activated. We hypothesis that all strains' ethanol yield are more than 0.25 g/g glucose
    model= add_ethanol_soft_constraint(model,ethanol_yield_lb=ethanol_yield_lb)

    # run RIPTiDe
    riptide_object = riptide.contextualize(model=model, transcriptome=transcript_abundances, fraction=relative_growth_rate,
                                           silent=True)
    # do not consider the growth data
    # riptide_object = riptide.contextualize(model=model, transcriptome=transcript_abundances, fraction=0.7,
    #                                          silent=True)

    print('concordance:',riptide_object.concordance)
    concordances[strainName]=riptide_object.concordance
    flux_samples=riptide_object.flux_samples
    mean_fluxes=flux_samples.mean()
    median_fluxes=flux_samples.median()
    std_fluxes=flux_samples.std()
    df_mean_fluxes[strainName]=mean_fluxes
    df_median_fluxes[strainName]=median_fluxes
    df_std_fluxes[strainName]=std_fluxes
    print('Ethanol production:',mean_fluxes['r_1761'],median_fluxes['r_1761'])

df_concordance=pd.DataFrame.from_dict(concordances,orient='index',columns=['p','r'])

# check the correlation between predicted growth and relative growth
mean_growth_pred=df_mean_fluxes.loc['r_2111',:].copy()
mediam_growth_pred=df_median_fluxes.loc['r_2111',:].copy()

df_growth=pd.DataFrame({'mean_growth_pred':mean_growth_pred,'mediam_growth_pred':mediam_growth_pred,'relative_growth':relative_growth})
# remove the strains with NaN value
df_growth=df_growth.dropna()
# set as float
df_growth=df_growth.astype(float)


# calculate the correlation
correlation_mean=df_growth['mean_growth_pred'].corr(df_growth['relative_growth'])
correlation_median=df_growth['mediam_growth_pred'].corr(df_growth['relative_growth'])
print('correlation between predicted growth and relative growth (mean):',correlation_mean)
print('correlation between predicted growth and relative growth (median):',correlation_median)


# # # save result
df_median_fluxes.to_csv(os.path.join(output_dir,f'humandairy_riptide_growth{growth_scale}_ethanolyield{ethanol_yield_lb}_median_flux.csv'))
df_mean_fluxes.to_csv(os.path.join(output_dir,f'humandairy_riptide_growth{growth_scale}_ethanolyield{ethanol_yield_lb}_mean_flux.csv'))
df_concordance.to_csv(os.path.join(output_dir,f'humandairy_riptide_growth{growth_scale}_ethanolyield{ethanol_yield_lb}_concordance.csv'))
#
# ethanol='r_1761'
# df_mean_fluxes.loc[ethanol,:].describe()
# df_median_fluxes.loc[ethanol,:].describe()

