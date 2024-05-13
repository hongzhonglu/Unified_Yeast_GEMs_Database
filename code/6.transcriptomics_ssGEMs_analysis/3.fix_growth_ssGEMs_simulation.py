import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\Unified_Yeast_GEMs_Database/code')
from cobra.io import read_sbml_model,write_sbml_model
from model_modifications import set_SCmedium
from cobra.flux_analysis import pfba,gapfill
from cobra.sampling import sample
import tqdm
import cobra
import os
import math

# set work dir
os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database')

# set configuration of cobra
configuration= cobra.Configuration()
configuration.solver = 'gurobi'

growth_scale=0.5

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
def ssGEM_sampling(model,process,n):
    sample_sol=sample(model=model,
                      n=n,
                      method='achr',
                      thinning=10,
                      processes=process)
    mean_fluxes=sample_sol.mean(axis=0)

    return mean_fluxes


if __name__ == '__main__':

    # set work dir
    os.chdir(r'D:\code\github\Unified_Yeast_GEMs_Database')

    ref_model=read_sbml_model('model/panYeast.xml')
    model=set_SCmedium(ref_model)
    # ref_growth=ref_fluxes['r_2111']
    ref_growth=ref_model.slim_optimize()

    rxnList=[r.id for r in ref_model.reactions]

    # load the relative growth rate data
    relative_growth,max_growth_index=get_relative_growth(pow=growth_scale)

    ssGEM_dir='code/6.transcriptomics_ssGEMs_analysis/2.GIMME_buildmodel/output/gimme_ssGEMs_shrinked'

    strainList=os.listdir(ssGEM_dir)
    strainList=[strain for strain in strainList if strain.replace('.xml','') in relative_growth.index]

    df_flux=pd.DataFrame(index=rxnList)

    for strain in tqdm.tqdm(strainList):
        # load model
        strainName=strain.replace('.xml','')
        if not os.path.exists(os.path.join(ssGEM_dir, strain)):
            print('The ssGEM of {} is not exist!'.format(strain))
            continue
        model = read_sbml_model(os.path.join(ssGEM_dir, strain))

        # 1. set medium
        model = set_SCmedium(model)

        # 2. fix growth rate according to the experimental data
        gr_max=0.448
        gr_strain=relative_growth[strainName]*gr_max

        # check if the strain need to be gapfilled
        pre_gr=model.slim_optimize()
        if pre_gr<gr_strain:
            print('The growth rate of %s is %f, and need to be gapfilled!' % (strain,pre_gr))
            # gapfilling
            try:
                gapfill_solution=gapfill(model,universal=ref_model,lower_bound=gr_strain,demand_reactions=False)
                # add gapfilling reactions to the model
                for rxn in gapfill_solution[0]:
                    model.add_reactions([rxn])
                    # print('Add reaction %s to the model!' % rxn.id)
            except:
                print('Gapfilling failed!')

        tol_ratio=0.01
        model.reactions.get_by_id('r_2111').bounds=gr_strain*(1-tol_ratio),gr_strain


        # 3. pFBA simulation: minimize glucose uptake rate
        # objective: minimize glucose uptake rate
        model.objective = 'r_1714'  # glucose exchange
        model.objective_direction = 'max'
        # open the glucose exchange reaction
        model.reactions.get_by_id('r_1714').bounds = (-1000, 0)
        try:
            # run FBA
            solution = model.optimize()
            # pfba
            # solution = pfba(model, fraction_of_optimum=1.05)
            fluxes = solution.fluxes
            # sample
            # fluxes=ssGEM_sampling(model,process=1,n=20)

        except:
            print('%s doesn\'t have solution' % strain)
            fluxes = pd.Series(0, index=rxnList)

        df_flux[strainName]=fluxes

    # fill nan with 0
    df_flux.fillna(0,inplace=True)
    # save fluxes
    # df_flux.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/fix_growth_gimme_pfba_fluxes_shrinked.csv')
    df_flux.to_csv(f'code/6.transcriptomics_ssGEMs_analysis/output/fix_growth{growth_scale}_gimme_fba_fluxes_shrinked.csv')
    # df_flux.to_csv('code/6.transcriptomics_ssGEMs_analysis/output/fix_growth_gimme_sample_fluxes.csv')


