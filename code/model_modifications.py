# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：model_modifications.py
# 2022/5/2
'''functions for modifying the model:content——
1. anaerobic_simulation
2. scale_biomass
3. rescale_pseudorxn
'''


def anaerobic_simulation(model):
    '''Note that the reaction coefficients can't be directly changed to 0 by present code,so the function need
    some parameters tuning for specific yeast model
    reference:https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/anaerobicModel.m'''
    # swich to anaerobic growth state

    anaerobic_model = model.copy()

    # 1. change: Refit GAM and NGAM to exp. data, change biomass composition
    GAM = 30.49
    P = 0.461  # change protein fraction——unfinished
    NGAM = 0
    NGAM_id = 'r_4046'
    GAM_id = 'r_4041'    # biomass pseudoreaction
    anaerobic_model=scale_biomass(anaerobic_model,'protein',P,'carbohydrate')   #change the protein fraction in biomass
    anaerobic_model.reactions.get_by_id(NGAM_id).bounds = NGAM, NGAM            #change non-growth associated maintainance(NGAM)
    biomass_rxn = anaerobic_model.reactions.get_by_id('r_4041')                 #change GAM
    biomass_rxn.add_metabolites({'s_0434[c]': -GAM,
                                 's_0803[c]': -GAM,
                                 's_0394[c]': GAM,
                                 's_0794[c]': GAM,
                                 's_1322[c]': GAM,
                                 }, combine=False)

    # 2. Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
    cofactor_rxn = anaerobic_model.reactions.get_by_id('r_4598')
    cofactor_rxn.add_metabolites({'s_3714[c]': 0,   # heme a
                                  's_1198[c]': 0,   # NAD
                                  's_1203[c]': 0,   # NADH
                                  's_1207[c]': 0,   # NADP+
                                  's_1212[c]': 0,   # NADPH
                                  's_0529[c]': 0},  # coenzyme A
                                 combine=False)

    # 3. change: Changes media to anaerobic————no O2；add sterol and fatty acid exchanges
    change_modium = anaerobic_model.medium
    change_modium['r_1992'] = 0
    # change medium:
    change_modium_lists = ['r_1757', 'r_1915', 'r_1994', 'r_2106', 'r_2134', 'r_2137', 'r_2189']
    for met in change_modium_lists:
        change_modium[met] = 1000.0
    anaerobic_model.medium = change_modium

    # 4. Blocked pathways for proper glycerol production.
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    anaerobic_model.reactions.get_by_id('r_0713').lower_bound = 0  # Mithocondria
    anaerobic_model.reactions.get_by_id('r_0714').lower_bound = 0  # Cytoplasm
    # Block glycerol dehydroginase (only acts in microaerobic conditions)
    anaerobic_model.reactions.get_by_id('r_0487').upper_bound = 0
    # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    anaerobic_model.reactions.get_by_id('r_0472').upper_bound = 0

    return anaerobic_model


def scale_biomass(model,component,new_value,balance_out):
    '''Scales the biomass composition
    model: metabolic model in COBRA format
    component: name of the component to rescale (e.g. "protein")——
        'carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor'
    new_value: new total fraction for target component
    balance_out: if chosen, the name of another component with which the model will be balanced out so that
                    the total mass remains = 1 g/gDW——
        'carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor'
    disp_output: if output from sumBioMass should be displayed (default = true)

    Usage: model = scale_biomass(model,component,new_value,balance_out,dispOutput)
    original_values get from SysBioChalmers/yeast-GEM repo:
    P -> 0.46 g/gDW
    C -> 0.38067 g/gDW
    R -> 0.061 g/gDW
    D -> 0.0037021 g/gDW
    L -> 0.087299 g/gDW
    I -> 0.0024815 g/gDW
    F -> 0.0048478 g/gDW
    X -> 1 gDW/gDW
    Growth = 0.083748 1/h
    reference:https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/scaleBioMass.m
    '''

    P=0.46
    C=0.38
    R=0.061
    D=0.0037
    L=0.0873
    I=0.00248
    F=0.00485
    content_Cap = [C, P, L, R, D, I, F]
    content_all = ['carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor']
    for i in range(len(content_all)):
        if content_all[i]==component:
            pos=i

    if pos:
        old_value=content_Cap[pos]
        f = new_value / old_value
        model=rescale_pseudorxn(model,component,f)
        if balance_out in content_all:
            for i in range(len(content_all)):
                if content_all[i] == component:
                    pos = i
            balance_value=content_Cap[pos]
            f=(balance_value - (new_value - old_value)) / balance_value
            model = rescale_pseudorxn(model, balance_out, f)
        else:
            print('input balance_out argument is wrong!╯︿╰')
    else:
        print('input component argument is wrong!╯︿╰')

    return model


def rescale_pseudorxn(model,met_name,f):
    '''Rescales a specific pseudoreaction by a given value
    model: the YeastGEM
    met_name: name of the component to rescale(eg. 'protein')
    f: fraction to use for rescaling

    usage: model=scaleBioMass(model,met_name,f)'''

    if met_name=='lipid':
        model = rescale_pseudorxn(model, 'lipid backbone', f)
        model = rescale_pseudorxn(model, 'lipid chain', f)
    else:
        rxn_name=met_name+' pseudoreaction'
        for reaction in model.reactions:
            if reaction.name==rxn_name:
                target_pseudorxn=reaction
                break

        for met,old_value in target_pseudorxn.metabolites.items():
            if old_value!=1:
                change=f-1
                target_pseudorxn.add_metabolites({met: old_value * change})

        return model


def set_SCmedium(model):
    # 'r_1663'; ... % bicarbonate exchange
    #  'r_4062'; ... % lipid backbone exchange
    #  'r_4064'};    % lipid chain exchange

    # block some exchange reactions
    blockedExchanges = ['r_1663', 'r_4062', 'r_4064']
    for rxn in blockedExchanges:
        try:
            model.reactions.get_by_id(rxn).bounds = (0, 0)
        except:
            pass
    minimal_medium={'r_1654': 1000, # ammonium exchange
                'r_1992': 1000, # oxygen exchange
                'r_2005': 1000, # phosphate exchange
                'r_2060': 1000, # sulphate exchange
                'r_1861': 1000, # iron exchange, for test of expanded biomass def
                'r_1832': 1000, # hydrogen exchange
                'r_2100': 1000, # water exchange
                'r_4593': 1000, # chloride exchange
                'r_4595': 1000, # Mn(2+) exchange
                'r_4596': 1000, # Zn(2+) exchange
                'r_4597': 1000, # Mg(2+) exchange
                'r_2049': 1000, # sodium exchange
                'r_4594': 1000, # Cu(2+) exchange
                'r_4600': 1000, # Ca(2+) exchange
                'r_2020': 1000,# potassium exchange
                'r_1714': 20    # glucose exchange
                    }
    # add 20 ammino acids
    ammino_acids={'r_1810':0.01,
                  'r_1873':0.01,
                  'r_1879':0.01,
                  'r_1880':0.01,
                  'r_1881':0.01,
                  'r_1883':0.01,
                  'r_1889':0.01,
                  'r_1891':0.01,
                  'r_1893':0.01,
                  'r_1897':0.01,
                  'r_1899':0.01,
                  'r_1900':0.01,
                  'r_1902':0.01,
                  'r_1903':0.01,
                  'r_1904':0.01,
                  'r_1906':0.01,
                  'r_1911':0.01,
                  'r_1912':0.01,
                  'r_1913':0.01,
                  'r_1914':0.01
                  }

    # merge dictionaries
    sc_medium = minimal_medium.copy()
    sc_medium.update(ammino_acids)
    model.medium = sc_medium

    return model


def set_minimal_medium(model):
    # 'r_1663'; ... % bicarbonate exchange
    #  'r_4062'; ... % lipid backbone exchange
    #  'r_4064'};    % lipid chain exchange

    # block some exchange reactions
    blockedExchanges = ['r_1663', 'r_4062', 'r_4064']
    for rxn in blockedExchanges:
        try:
            model.reactions.get_by_id(rxn).bounds = (0, 0)
        except:
            pass
    minimal_medium={'r_1654': 1000, # ammonium exchange
                'r_1992': 1000, # oxygen exchange
                'r_2005': 1000, # phosphate exchange
                'r_2060': 1000, # sulphate exchange
                'r_1861': 1000, # iron exchange, for test of expanded biomass def
                'r_1832': 1000, # hydrogen exchange
                'r_2100': 1000, # water exchange
                'r_4593': 1000, # chloride exchange
                'r_4595': 1000, # Mn(2+) exchange
                'r_4596': 1000, # Zn(2+) exchange
                'r_4597': 1000, # Mg(2+) exchange
                'r_2049': 1000, # sodium exchange
                'r_4594': 1000, # Cu(2+) exchange
                'r_4600': 1000, # Ca(2+) exchange
                'r_2020': 1000,# potassium exchange
                'r_1714': 20    # glucose exchange
                    }
    model.medium = minimal_medium
    return model
