# -*- coding: utf-8 -*-
# date : 2024/11/20 
'''Find all essential reaction for growth in synthetic complete medium.'''
import sys
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database/code')
import cobra
from cobra.flux_analysis import single_reaction_deletion
from model_modifications import set_SCmedium
import tqdm

model=cobra.io.read_sbml_model(r'model/panYeast.xml')
model=set_SCmedium(model)

growth=model.slim_optimize()

essential_rxnlist=list()
for rxn in tqdm.tqdm(model.reactions):
    id=rxn.id
    with model:
        model.reactions.get_by_id(id).knock_out()
        try:
            mutant_gr=model.slim_optimize()
        except:
            mutant_gr=0
        if mutant_gr < growth*0.5:
            essential_rxnlist.append(id)

# save result as json
import json
with open(r'code/6.transcriptomics_ssGEMs_analysis/2.tinit_buildmodel/output/essential_rxnlist.json','w') as f:
    json.dump(essential_rxnlist,f)
