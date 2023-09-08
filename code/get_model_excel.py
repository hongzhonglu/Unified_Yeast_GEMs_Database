import cobra
import pandas as pd
import os
import re

model = cobra.io.read_sbml_model(r'model/panYeast_v3.xml')

# write reaction infor to excel
col = ['ID', 'Name', 'Equation','GPR']
sheet_rxn = pd.DataFrame(index=list(range(len(model.reactions))), columns=col)

rxnid = [r.id for r in model.reactions]
rxnName = [n.name for n in model.reactions]
tmprxnEquation = [e.reaction for e in model.reactions]
rxnGene = [g.gene_reaction_rule for g in model.reactions]


# equation处理
def get_eq(model, equation):
    ss = re.compile(r's_\d{4}\[[a-z]\]')
    sss = ss.findall(equation)
    for s in sss:
        equation = equation.replace(s, model.metabolites.get_by_id(s).name)
    return equation

rxnEquation = [get_eq(model, eq) for eq in tmprxnEquation]

sheet_rxn.loc[:, 'ID'] = rxnid
sheet_rxn.loc[:, 'Name'] = rxnName
sheet_rxn.loc[:, 'Equation'] = rxnEquation
sheet_rxn.loc[:, 'GPR'] = rxnGene

# write mets infor to excel
colmet = ['ID', 'Name', 'Formula', 'compartment','kegg', 'chebi']
sheet_met = pd.DataFrame(index=list(range(len(model.metabolites))), columns=colmet)

metID = [m.id for m in model.metabolites]
metName = [mn.name for mn in model.metabolites]
metFormula = [mf.formula for mf in model.metabolites]
metcom = [mcom.compartment for mcom in model.metabolites]
metkg = []
for mk in model.metabolites:
    if 'kegg.compound' in mk.annotation:
        metkg.append(mk.annotation['kegg.compound'])
    else:
        metkg.append('')

metchebi = []
for mk in model.metabolites:
    if 'chebi' in mk.annotation:
        metchebi.append(mk.annotation['chebi'])
    else:
        metchebi.append('')

sheet_met.loc[:, 'ID'] = metID
sheet_met.loc[:, 'Name'] = metName
sheet_met.loc[:, 'Formula'] = metFormula
sheet_met.loc[:, 'kegg'] = metkg
sheet_met.loc[:, 'chebi'] = metchebi
sheet_met.loc[:, 'compartment'] = metcom

# write genes info to excel
colgene = ['ID', 'Name','rxnID']
sheet_gene = pd.DataFrame(index=list(range(len(model.genes))), columns=colgene)

geneIDs=[g.id for g in model.genes]
geneNames=[g.name for g in model.genes]
geneRxnID=[]
for g in model.genes:
    rxnList=g.reactions
    geneRxn=';'.join([r.id for r in rxnList])
    geneRxnID.append(geneRxn)


sheet_gene["ID"]=geneIDs
sheet_gene["Name"]=geneNames
sheet_gene["rxnID"]=geneRxnID


# 直接用to_excel每次都会新建一个文件，要先创建一个ExcelWriter才能在一个excel里面写多个sheet
writer = pd.ExcelWriter(r'model/panYeast_v3.xlsx')

sheet_rxn.to_excel(writer,
                   sheet_name='reactions',
                   index=False)
sheet_met.to_excel(writer,
                   sheet_name='metabolites',
                   index=False)
sheet_gene.to_excel(writer,sheet_name="genes",index=False)
writer.save()

print('finish')
