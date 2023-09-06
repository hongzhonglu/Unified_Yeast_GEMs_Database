# -*- coding: utf-8 -*-
# date : 2023/7/6 
# author : wangh
# file : 2.new_method_integrate_transcriptome.py
# project : Unified_Yeast_GEMs_Database
'''integrate transcriptome data into ssGEM by constrain the bounds of reactions according to the expression level of genes.
* rules:
    1. if the gene is not expressed, set the bounds of the corresponding reaction as 0
'''
import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database\code')
from cobra.io import read_sbml_model,write_sbml_model
from model_modifications import set_SCmedium
from cobra.flux_analysis import pfba,gapfill,geometric_fba
from cobra.sampling import sample
import numpy as np
import math
import tqdm


def tpm_to_foldchange(tpmMatrix):
    df_expression_log2=tpmMatrix.applymap(lambda x:math.log2(x+1))
    # reference: mean value for rows
    df_expression_log2_mean=df_expression_log2.apply(lambda x: x[x!=0].mean(),axis=1)
    df_expression_foldchange=df_expression_log2.apply(lambda x: x/df_expression_log2_mean,axis=0)
    # fill inf as 10, and nan as 0
    df_expression_foldchange.replace([np.inf, -np.inf], 10,inplace=True)
    df_expression_foldchange.fillna(0,inplace=True)
    return df_expression_foldchange


def calculate_bounds(fva_ub,fva_lb,ref_flux,foldchange):
    '''rules:
        1. if the gene is not expressed, set the bounds of the corresponding reaction as 0
        2. if the gene is expressed, set the bounds of the corresponding reaction as the foldchange of the gene
    '''
    flux = ref_flux * foldchange
    if ref_flux>0:
        if flux>fva_lb and flux<fva_ub:
            if foldchange>1:
                new_lb=flux
                new_ub=fva_ub
            elif foldchange==1:
                new_lb=fva_lb
                new_ub=fva_ub
            elif foldchange<1:
                new_lb=fva_lb
                new_ub=flux
        else:
            new_lb=fva_lb
            new_ub=fva_ub
    elif ref_flux==0:
        new_lb=fva_lb
        new_ub=fva_ub
    elif ref_flux<0:
        if flux>fva_lb and flux<fva_ub:
            if foldchange>1:
                new_lb=fva_lb
                new_ub=flux
            elif foldchange==1:
                new_lb=fva_lb
                new_ub=fva_ub
            elif foldchange<1:
                new_lb=flux
                new_ub=fva_ub
        else:
            new_lb=fva_lb
            new_ub=fva_ub
    return new_lb,new_ub


def scaling_foldchange(rxns_fcMatrix,scaling_factor):
    '''scaling the foldchange of reactions according to the scaling factor'''
    mean_fc=rxns_fcMatrix.mean(axis=1)
    std_fc=rxns_fcMatrix.std(axis=1)

    # standardization: (x-mean)/std
    standardized_fcMatrix=rxns_fcMatrix.apply(lambda x: (x-mean_fc)/std_fc,axis=0)
    scaled_data=standardized_fcMatrix*scaling_factor

    # rescaling: x*std+mean
    rxns_fcMatrix_scaled=scaled_data.apply(lambda x: x*std_fc+mean_fc,axis=0)

    # fill nan as 1
    rxns_fcMatrix_scaled.fillna(1,inplace=True)
    return rxns_fcMatrix_scaled


def build_tissGEM(model, df_fva_bounds, ref_fluxes, rxns_foldchange, ssGEM_dir):
    '''rules:
        1. if the gene is not expressed, set the bounds of the corresponding reaction as 0
        2. if the gene is expressed, set the bounds of the corresponding reaction as the foldchange of the gene
        '''
    # only bound the reactions with more than 0.001 for absolute flux
    rxnIDlist = [rxn.id for rxn in model.reactions if abs(ref_fluxes[rxn.id]) > 0.001]
    print(len(rxnIDlist),'reactions need to be constrained')
    for rxnID in rxnIDlist:
        fva_lb= float(df_fva_bounds.loc[rxnID, 'minimum'])
        fva_ub= float(df_fva_bounds.loc[rxnID, 'maximum'])
        ref_flux = ref_fluxes[rxnID]
        foldchange = rxns_foldchange[rxnID]
        # ignore the reactions with foldchange==0
        if foldchange==0:
            model.reactions.get_by_id(rxnID).bounds = (fva_lb, fva_ub)
            continue

        # add new bounds
        new_lb, new_ub = calculate_bounds(fva_ub=fva_ub,
                                          fva_lb=fva_lb,
                                          ref_flux=ref_flux,
                                          foldchange=foldchange)
        # only keep 4 digits after the decimal point
        new_lb = round(new_lb, 4)
        new_ub = round(new_ub, 4)
        try:
            if new_lb>=0:
                model.reactions.get_by_id(rxnID).bounds = (new_lb*0.9, new_ub*1.1)
            elif new_ub<=0:
                model.reactions.get_by_id(rxnID).bounds = (new_lb*1.1, new_ub*0.9)
            else:
                model.reactions.get_by_id(rxnID).bounds = (new_lb, new_ub)
        except:
            model.reactions.get_by_id(rxnID).bounds = (fva_lb, fva_ub)
            # print()
        # if new bounds make the model infeasible, ignore the reaction new bounds
        # try:
        #     gr = model.slim_optimize()
        #     if gr < 0.04:
        #         model.reactions.get_by_id(rxnID).bounds = (fva_lb, fva_ub)
        #     if not gr:
        #         model.reactions.get_by_id(rxnID).bounds = (fva_lb, fva_ub)
        # except:
        #     model.reactions.get_by_id(rxnID).bounds = (fva_lb, fva_ub)

    return model


def ssGEM_sampling(model,min_growth,process,n):
    # set the objective function as biomass
    model.objective='r_2111'  # set the objective function as biomass
    # set the lower bound of biomass as the minimum growth rate
    model.reactions.get_by_id('r_2111').lower_bound=min_growth
    sample_sol=sample(model=model,
                      n=n,
                      method='optgp',
                      thinning=10,
                      processes=process)
    mean_fluxes=sample_sol.mean(axis=0)
    # calculate the deviation of the fluxes
    std_fluxes=sample_sol.std(axis=0)
    # calculate the coefficient of variation
    cv_fluxes=std_fluxes/mean_fluxes

    return mean_fluxes,std_fluxes,cv_fluxes




if __name__ == '__main__':

    # set work dir
    os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database_from_13pro\Unified_Yeast_GEMs_Database')

    reference=read_sbml_model('model/panYeast_v4_5.xml')
    reference=set_SCmedium(reference)
    ref_growth=reference.slim_optimize()
    ref_solution=pfba(reference,fraction_of_optimum=0.90)
    ref_fluxes=ref_solution.fluxes

    # load the transcriptome data
    rxn_tpmMatrix=pd.read_csv('code/7.transcriptomics_ssGEMs_analysis/output/sce969_rxn_tpmMatrix.csv',index_col=0)

    # convert tpm to foldchange
    rxn_foldchange=tpm_to_foldchange(rxn_tpmMatrix)

    # scaling the foldchange data to make it more closed to 1 , and keep the originial distribution
    scaling_factor=0.2
    scaled_rxn_foldchange=scaling_foldchange(rxns_fcMatrix=rxn_foldchange,
                                             scaling_factor=scaling_factor)

    # load fva bounds
    df_fva_bounds=pd.read_excel('code/7.transcriptomics_ssGEMs_analysis/output/panYeast_fva_result.xlsx',index_col=0)

    ssGEM_dir='model/pan1800_ssGEMs'

    r_0711=reference.reactions.get_by_id('r_0711')
    # r_0486=reference.reactions.get_by_id('r_0486')
    rxnList=list(df_fva_bounds.index)
    strainList=list(rxn_foldchange.columns)

    df_flux=pd.DataFrame(index=['growth']+rxnList)
    df_std_flux=pd.DataFrame(index=rxnList)
    df_cv_flux=pd.DataFrame(index=rxnList)

    for strain in tqdm.tqdm(strainList):

        # load model
        ssGEM_name = strain + '.xml'
        if not os.path.exists(os.path.join(ssGEM_dir, ssGEM_name)):
            print('The ssGEM of {} is not exist!'.format(strain))
            continue
        model = read_sbml_model(os.path.join(ssGEM_dir, ssGEM_name))
        model = set_SCmedium(model)
        gr=model.slim_optimize()
        print('The growth rate of %s is %f' %(strain,gr))

        # check predicted growth rate, and do gapfilling if the model grow unnormally
        if gr<0.95*ref_growth:
            print('The growth rate of %s is %f, and need to be gapfilled!' % (strain,gr))
            # gapfilling
            try:
                gapfill_solution=gapfill(model,universal=reference,lower_bound=0.95*ref_growth,demand_reactions=False)
                # add gapfilling reactions to the model
                for rxn in gapfill_solution[0]:
                    model.add_reactions([rxn])
                    print('Add reaction %s to the model!' % rxn.id)
            except:
                print('Gapfilling failed!')

        # extract the transcriptome data of the strain
        rxns_foldchange=scaled_rxn_foldchange[strain]
        model=build_tissGEM(model=model,
                            df_fva_bounds=df_fva_bounds,
                            ref_fluxes=ref_fluxes,
                            rxns_foldchange=rxns_foldchange,
                            ssGEM_dir=ssGEM_dir)

        if model is not None:
            model=set_SCmedium(model)

            # supply essential reaction: r_0711
            try:
                model.reactions.get_by_id('r_0711').bounds = (0, 1000)
            except:
                print('%s doesn\'t have r_0711' % strain)
                model.add_reactions([r_0711])
                model.reactions.get_by_id('r_0711').bounds = (0, 1000)

            model.objective = 'r_2111'  # set the objective function as biomass

            try:
                # run FBA
                # solution = model.optimize()
    #             # pfba
    #             # solution = pfba(model, fraction_of_optimum=0.90)
    #             fluxes = solution.fluxes
    #             # use sample to get the fluxes
                fluxes,std_fluxes,cv_fluxes = ssGEM_sampling(model=model,
                                           n=100,
                                           process=10,
                                           min_growth=0.07)
                gr = fluxes['r_2111']
            except:
                print('%s doesn\'t have solution' % strain)
                fluxes = pd.Series(0, index=rxnList)
                std_fluxes = pd.Series(0, index=rxnList)
                cv_fluxes = pd.Series(0, index=rxnList)
                gr = 0
    #
            fluxes['growth'] = gr
            # set 'growth as the
            df_flux.loc[:, strain] = fluxes
            df_std_flux.loc[:, strain] = std_fluxes
            df_cv_flux.loc[:, strain] = cv_fluxes
        else:
            print('%s doesn\'t have ssGEM' % strain)
    #
    #
    # # check how many 0 in the growth
    df_flux.loc['growth'].value_counts()


    # # save result
    df_flux.to_csv('code/7.transcriptomics_ssGEMs_analysis/output/PNAS_v3_tissGEMs_sample_flux.csv')
    df_cv_flux.to_csv('code/7.transcriptomics_ssGEMs_analysis/output/PNAS_v3_tissGEMs_sample_cv_flux.csv')
    df_std_flux.to_csv('code/7.transcriptomics_ssGEMs_analysis/output/PNAS_v3_tissGEMs_sample_std_flux.csv')

