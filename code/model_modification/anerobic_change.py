# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：anerobic_change.py
# 2022/4/26

def anaerobicsimulation(model):
    '''Note that the reaction coefficients can't be directly changed to 0 by present code,so the function need
    some parameters tuning for specific yeast model
    reference:https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/anaerobicModel.m'''
    # 切换anaerobic growth medium
    from cobra.io import read_sbml_model
    import os
    import pandas as pd

    anerobic_model = model.copy()

    # 1. change: Refit GAM and NGAM to exp. data, change biomass composition
    GAM=30.49
    P=0.461     #change protein fraction——unfinished
    NGAM=0
    NGAM_id = 'r_4046'
    GAM_id = 'r_4041'
    anerobic_model.reactions.get_by_id(NGAM_id).bounds = NGAM,NGAM
    biomass_rxn=anerobic_model.reactions.get_by_id('r_4041')
    biomass_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_0434[c]'): 55.4 - 30.49,
                             anerobic_model.metabolites.get_by_id('s_0803[c]'): 55.4 - 30.49,
                             anerobic_model.metabolites.get_by_id('s_0394[c]'): -55.4 + 30.49,
                             anerobic_model.metabolites.get_by_id('s_0794[c]'): -55.4 + 30.49,
                             anerobic_model.metabolites.get_by_id('s_1322[c]'): -55.4 + 30.49})

    # 2. Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
    cofactor_rxn = anerobic_model.reactions.get_by_id('r_4598')
    # yeastGEM8.5 co-factors reaction coeffiecients
    # cofactor_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_3714[c]'):9.99999997475243e-07,
    #                                 anerobic_model.metabolites.get_by_id('s_1198[c]'):0.00264999992214143,
    #                                 anerobic_model.metabolites.get_by_id('s_1203[c]'):0.000150000007124618 ,
    #                                 anerobic_model.metabolites.get_by_id('s_1207[c]'):0.000569999974686652,
    #                                 anerobic_model.metabolites.get_by_id('s_1212[c]'):0.00270000007003546,
    #                                 anerobic_model.metabolites.get_by_id('s_0529[c]'):0.000190000006114133 })

    # panYeast co-factors reaction coefficients
    cofactor_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_3714[c]'): 9.99999997475243e-07,
                                  anerobic_model.metabolites.get_by_id('s_1198[c]'): 0.00100000004749745,
                                  anerobic_model.metabolites.get_by_id('s_1203[c]'): 0.000788461999036372,
                                  anerobic_model.metabolites.get_by_id('s_1207[c]'): 6.54000032227486e-05,
                                  anerobic_model.metabolites.get_by_id('s_1212[c]'): 7.6900003477931e-05,
                                  anerobic_model.metabolites.get_by_id('s_0529[c]'): 0.000190000006114133})

    # 3. change: Changes media to anaerobic————no O2；add sterol and fatty acid exchanges
    change_modium = anerobic_model.medium
    change_modium['r_1992'] = 0
    # 改变培养基成分，允许甾醇和脂肪酸的交换反应
    change_modium_lists = ['r_1757', 'r_1915', 'r_1994', 'r_2106', 'r_2134', 'r_2137', 'r_2189']
    for met in change_modium_lists:
        change_modium[met] = 1000.0
    anerobic_model.medium = change_modium
        # print(test_model.reactions.get_by_id(met).bounds)

    # 4. Blocked pathways for proper glycerol production.
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    anerobic_model.reactions.get_by_id('r_0713').lower_bound = 0    #Mithocondria
    anerobic_model.reactions.get_by_id('r_0714').lower_bound = 0    #Cytoplasm
    # Block glycerol dehydroginase (only acts in microaerobic conditions)
    anerobic_model.reactions.get_by_id('r_0487').upper_bound=0
    # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    anerobic_model.reactions.get_by_id('r_0472').upper_bound = 0

    return anerobic_model