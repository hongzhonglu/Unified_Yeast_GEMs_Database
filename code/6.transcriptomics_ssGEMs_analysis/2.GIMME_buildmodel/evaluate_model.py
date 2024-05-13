import cobra

model = cobra.io.read_sbml_model('code/6.transcriptomics_ssGEMs_analysis/GIMME_buildmodel/output/gapfilled_ssGEMs/AGK_1.re.xml')

model_gimme= cobra.io.read_sbml_model('code/6.transcriptomics_ssGEMs_analysis/GIMME_buildmodel/output/gimme_ssGEMs/YBK.re.xml')
model_gimme.slim_optimize()