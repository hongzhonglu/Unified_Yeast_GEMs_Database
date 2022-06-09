# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：read_model.py
# 2022/4/29

from cobra.io import read_sbml_model
from cobra import *

model=read_sbml_model('auxotroph_predict/c_saudii_YCFA.xml')
model.objective.expression
obj=model.reactions.get_by_id('EX_cpd11416_c0')
obj_met=model.metabolites.get_by_id('cpd11416_c0')

model.slim_optimize()
# test_auxofind
import sys
import os

sys.path.append('auxotroph_predict/AuxoFind')
from AuxoFind import *

dir=os.getcwd()

CEGs=get_CEGs(model, growth_medium = [])

cobra_model=model

