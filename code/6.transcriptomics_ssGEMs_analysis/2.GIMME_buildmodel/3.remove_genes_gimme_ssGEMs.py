'''Shrink the ssGEMs from GIMME method by removing unsed genes and metabolites.'''

from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import prune_unused_metabolites, remove_genes
import os
import pandas as pd
import tqdm

# set work directory
os.chdir(r'D:\github\Unified_Yeast_GEMs_Database')

input_model_dir = 'code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs'
output_model_dir = 'code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_shrinked'
strainList = os.listdir(input_model_dir)

df_gimme_ssGEMs_size=pd.DataFrame()
gene_numbs=[]
rxn_numbs=[]
met_numbs=[]
strain_names=[]

for strain in tqdm.tqdm(strainList):
    # check if the model have existed in the output directory
    if os.path.exists(os.path.join(output_model_dir, strain)):
        print('The ssGEM of {} have existed in the output directory!'.format(strain))
        continue
    # load model
    model = read_sbml_model(os.path.join(input_model_dir, strain))

    ori_gene_num = len(model.genes)
    ori_rxn_num = len(model.reactions)
    ori_met_num = len(model.metabolites)

    # remove unused genes
    to_remove = []
    for g in model.genes:
        if len(g.reactions) == 0:
            to_remove.append(g.id)
    remove_genes(model, to_remove)

    # remove unused metabolites
    prune_unused_metabolites(model)

    new_gene_num = len(model.genes)
    new_rxn_num = len(model.reactions)
    new_met_num = len(model.metabolites)

    # print('For strain %s', strain)
    # print('rxn number: %d -> %d' % (ori_rxn_num, new_rxn_num))
    # print('gene number: %d -> %d' % (ori_gene_num, new_gene_num))
    # print('metabolite number: %d -> %d' % (ori_met_num, new_met_num))

    strainname = strain.replace('.re.xml','')

    strain_names.append(strainname)
    gene_numbs.append(new_gene_num)
    rxn_numbs.append(new_rxn_num)
    met_numbs.append(new_met_num)

    model.id= strainname
    # save the model
    write_sbml_model(model, os.path.join(output_model_dir, strain))

df_gimme_ssGEMs_size['strain'] = strain_names
df_gimme_ssGEMs_size['gene_number'] = gene_numbs
df_gimme_ssGEMs_size['reaction_number'] = rxn_numbs
df_gimme_ssGEMs_size['metabolite_number'] = met_numbs
df_gimme_ssGEMs_size.to_csv('code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_size.csv', index=False)


