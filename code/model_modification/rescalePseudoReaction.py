# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：rescalePseudoReaction.py
# 2022/4/26

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


