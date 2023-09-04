# -*- coding: utf-8 -*-
# date : 2023/6/27 
# author : wangh
# file : add_fva_bound.py
# project : Unified_Yeast_GEMs_Databas
import sys
sys.path.append('code')
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis,pfba
import pandas as pd
from model_modifications import set_SCmedium

def set_max_sum_of_fluxes(model,max_sum_of_fluxes):
    '''set the maximum bound of sum of absolute fluxes'''
    coefficients = dict()
    for rxn in model.reactions:
        coefficients[rxn.forward_variable] = 1.
        coefficients[rxn.reverse_variable] = 1.
    constraint = model.problem.Constraint(0, lb=0, ub=max_sum_of_fluxes, name='sum_of_fluxes_constraint')
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

panYeast=read_sbml_model('model/panYeast_v4_5.xml')
panYeast=set_SCmedium(panYeast)

# # calculate the sum of absolute fluxes
# sum_of_fluxes = pfba(panYeast).objective_value
# # max_sum_of_fluxes = sum_of_fluxes * 3
#
# # constrain the sum of absolute fluxes
# panYeast=set_max_sum_of_fluxes(panYeast,max_sum_of_fluxes)

rxnList=[rxn for rxn in panYeast.reactions]
fva_sol=flux_variability_analysis(panYeast,rxnList,fraction_of_optimum=0.2)
fva_sol.to_excel("code/7.transcriptomics_ssGEMs_analysis/output/panYeast_fva_result.xlsx")