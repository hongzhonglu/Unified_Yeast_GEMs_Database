# -*- coding: utf-8 -*-
# date : 2024/11/20
import os
import tqdm
import pandas as pd
import cobra
import re
from cobra.flux_analysis import gapfill
import sys

from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper,ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
from troppo.tasks.task_io import ExcelTaskIO
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database/code')
from model_modifications import set_SCmedium

def run_init(model,expression_data,essential_rxnlist):

    model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999)
    data_map = expression_data.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                      and_func=min, or_func=sum)
    def score_apply(data_map,protected_core):
        '''Ignore reaction with no gene realtionship'''
        # return {k:0  if v is None else v for k, v in reaction_map_scores.items()}
        maxv = max([k for k in data_map.get_scores().values() if k is not None])
        scores = {k: (k/maxv if v < 0 else v) if v is not None else 0 for k, v in data_map.get_scores().items()}
        scores.update({x: 100*max(scores.values()) for x in protected_core})
        return scores

    # set protected rxnlist: essential rxns and exchange rxn
    exchange_rxnlist=[rxn.id for rxn in model.exchanges]
    protected_core=essential_rxnlist+exchange_rxnlist
    # protected_core=essential_rxnlist
    scores=score_apply(data_map,protected_core=protected_core)
    with model:
        gr=model.slim_optimize()
        gr=round(gr,2)
        # release rxn bounds
        for rxn in model.reactions:
            model.reactions.get_by_id(rxn.id).bounds = (-1000, 1000)

        # set growth constraint
        model.reactions.get_by_id('r_2111').bounds = 0.5*gr, 1000

        model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999)

    tinit_result = model_wrapper.run_from_omics(
        omics_data=scores, algorithm='tinit',integration_strategy='continuous')

    toremove = [id for id, value in tinit_result.items() if value == False]
    print(f'{len(toremove)} reactions are removed!')

    mutant=model.copy()
    mutant.remove_reactions(toremove, remove_orphans=True)
    return mutant

input_dir='model/ssGEMs/'
output_dir='code/6.transcriptomics_ssGEMs_analysis/output/tinit_ssGEMs/'
panModel=cobra.io.read_sbml_model('model/panYeast.xml')
panModel=set_SCmedium(panModel)
ref_gr=panModel.slim_optimize()

# load essential reactions
import json
with open(r'code/6.transcriptomics_ssGEMs_analysis/2.tinit_buildmodel/output/essential_rxnlist.json','r') as f:
    essential_rxnlist=json.load(f)

# 2. load expression data
# expression_data = pd.read_csv('code/6.transcriptomics_ssGEMs_analysis/2.tinit_buildmodel/example/CCLE_breast_cancer_preprocessed.csv', index_col=0)
expression_data = pd.read_csv(r'code/6.transcriptomics_ssGEMs_analysis/output/sce969_transcriptome_tpmMatrix.csv',index_col=0).T
omics_container = TabularReader(path_or_df=expression_data, nomenclature='sce_ssGEMs',
                                omics_type='transcriptomics').to_containers()


ssGEM_list=os.listdir(input_dir)
df_ssGEMs_size=pd.DataFrame(columns=['gene_number','reaction_number','metabolite_number'])
for sample_expression in tqdm.tqdm(omics_container):
    id=sample_expression.condition+'.xml'
    # Check if the model has been processed
    if os.path.exists(os.path.join(output_dir,id)):
        continue
    if id in ssGEM_list:
        model = cobra.io.read_sbml_model(os.path.join(input_dir,id))
        # add essential rxn
        model_rxnlist=[x.id for x in model.reactions]
        toadd=[panModel.reactions.get_by_id(id) for id in essential_rxnlist if id not in model_rxnlist]
        if len(toadd)>0:
            model.add_reactions(toadd)
        model=set_SCmedium(model)
        gr=model.slim_optimize()
        if gr < 0.5*ref_gr:
            try:
                gapfill_solution = gapfill(model, universal=panModel, lower_bound=ref_gr * 0.5,
                                           demand_reactions=False)
                # add gapfilling reactions to the model
                for rxn in gapfill_solution[0]:
                    model.add_reactions([rxn])
                    # print('Add reaction %s to the model!' % rxn.id)
            except:
                print('Gapfilling failed!')
                continue
        mutant=run_init(model=model,expression_data=sample_expression,essential_rxnlist=essential_rxnlist)
        print(f'{id} have {len(mutant.reactions)} reactions, {len(mutant.genes)} genes')
        df_ssGEMs_size.loc[sample_expression.condition]=[len(mutant.genes),len(mutant.reactions),len(mutant.metabolites)]
        cobra.io.write_sbml_model(mutant,os.path.join(output_dir,id))

# save result
df_ssGEMs_size.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/tinit_ssGEMs_size.csv')