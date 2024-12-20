import pandas as pd
import os
import sys
sys.path.append(r'/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code')
from cobra.io import read_sbml_model,write_sbml_model
from cobra.flux_analysis import gapfill
from model_modifications import set_SCmedium
import shutil


# set work directory
os.chdir('/dssg/home/acct-clslhz/clslhz/why_ssGEM/Unified_Yeast_GEMs_Database/code/7.transcriptomics_ssGEMs_analysis')

# load the transcriptome data
rxn_expMatrix = pd.read_csv('output/sce969_rxn_expressionMatrix_normalized.csv',
                            index_col=0)

ssGEM_dir = '../../model/ssGEMs'
output_dir = 'output/gapfilled_ssGEMs'

ref_model= read_sbml_model('../../model/panYeast.xml')
ref_model= set_SCmedium(ref_model)
ref_growth = ref_model.slim_optimize()

strainList = list(rxn_expMatrix.columns)

# check if output directory exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for strain in strainList:

    # check if the model have existed in the output directory
    if os.path.exists(os.path.join(output_dir, strain + '.xml')):
        print('The ssGEM of {} have existed in the output directory!'.format(strain))
        continue
    # load model
    ssGEM_name = strain + '.xml'
    if not os.path.exists(os.path.join(ssGEM_dir, ssGEM_name)):
        print('The ssGEM of {} is not exist!'.format(strain))
        continue
    model = read_sbml_model(os.path.join(ssGEM_dir, ssGEM_name))
    model.solver.configuration.timeout = 100
    # set medium
    model = set_SCmedium(model)
    gr = model.slim_optimize()
    # print('The growth rate of %s is %f' %(strain,gr))

    # check predicted growth rate, and do gapfilling if the model grow unnormally
    if gr < 0.95 * ref_growth:
        print('The growth rate of %s is %f, and need to be gapfilled!' % (strain, gr))
        # gapfilling
        try:
            gapfill_solution = gapfill(model, universal=ref_model, lower_bound=0.95 * ref_growth,
                                       demand_reactions=False)
            # add gapfilling reactions to the model
            for rxn in gapfill_solution[0]:
                model.add_reactions([rxn])
                print('Add reaction %s to the model!' % rxn.id)

            # save the gapfilled model
            write_sbml_model(model, os.path.join(output_dir, ssGEM_name))

        except:
            print('Gapfilling failed!')
            continue

    else:
        # copy the model to the output directory
        shutil.copy(os.path.join(ssGEM_dir, ssGEM_name), os.path.join(output_dir, ssGEM_name))

