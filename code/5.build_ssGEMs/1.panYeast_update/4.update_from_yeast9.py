# -*- coding: utf-8 -*-
# date : 2024/10/28 
'''Update panYeast according to yeast9.
yeast9 vs yeast8: 29 new genes, 202 new reactions, and 139 new metabolites.'''
import pandas as pd
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite, Gene
from cobra.core.gene import GPR

def gpr_to_list(gpr_string):
    '''Convert gene_reaction_rule string to a list'''
    gpr_list=gpr_string.split(' or ')
    return gpr_list

yeast9=read_sbml_model(r'model/yeast-GEM.xml')
panmodel=read_sbml_model(r'model/panYeast.xml')

panmodel_new=panmodel.copy()

# 1.Find new mets need to update
yeast9_to_metsID_mapdict=dict()
new_metsList=[]
panmodel_metsIDlist=[m.id for m in panmodel_new.metabolites]
panmodel_metsID_re_list=[i.split('[')[0] for i in panmodel_metsIDlist]
error_list=[]
for met in yeast9.metabolites:
    id=met.id
    name=met.name
    formula=met.formula
    if id not in panmodel_metsID_re_list:
        new_metsList.append(id)
        # print(id,name)
    else:
        index=panmodel_metsID_re_list.index(id)
        panmodel_name=panmodel_new.metabolites.get_by_id(panmodel_metsIDlist[index]).name
        pan_formula=panmodel.metabolites.get_by_id(panmodel_metsIDlist[index]).formula
        if (name != panmodel_name.split(' [')[0]) & (formula != pan_formula):
            error_list.append(id)
            # print('panmodel name is :',panmodel_name)
            # print('yeast9 name is :',name)
        # if formula != pan_formula:
        #     error_list.append(id)
        #     print('panmodel formula is :',pan_formula)
        #     print('yeast9 formula is :',formula)
        else:
            yeast9_to_metsID_mapdict[id]=panmodel_metsIDlist[index]
toadd_metsIDList=new_metsList+error_list
i=0
toadd_mets_info=dict()
for id in toadd_metsIDList:
    yeast9_id=id
    met=yeast9.metabolites.get_by_id(yeast9_id)
    compartment=met.compartment
    name=met.name
    formula=met.formula
    if id in error_list:
        pan_id='r_'+str(4332+i)
        i+=1
    else:
        pan_id=id.split('[')[0]
    pan_id=pan_id+'['+compartment+']'
    toadd_mets_info[id]=[id,pan_id,name,formula,compartment]

df_toadd_mets_info=pd.DataFrame.from_dict(toadd_mets_info,orient='index',columns=['yeast9_id','pan_id','name','formula','compartment'])
print(f'add {len(df_toadd_mets_info)} new mets')
for id in df_toadd_mets_info.index:
    panID=df_toadd_mets_info.loc[id,'pan_id']
    yeast9_to_metsID_mapdict[id]=panID

# add 123 new mets to panmodel by add a pseudo reaction(it will be removed later by panmodel_new.remove_reactions(reactions=['add_new_mets'],remove_orphans=True))
metsList=[]
add_mets_Reaction=Reaction(id='add_new_mets',name='add new mets to panmodel')
for id in df_toadd_mets_info.index:
    annotation=yeast9.metabolites.get_by_id(id).annotation
    met=Metabolite(id=df_toadd_mets_info.loc[id,'pan_id'],
                   name=df_toadd_mets_info.loc[id,'name'],
                   formula=df_toadd_mets_info.loc[id,'formula'],
                   compartment=df_toadd_mets_info.loc[id,'compartment'])
    met.annotation=annotation
    met.notes={'yeast9_id':df_toadd_mets_info.loc[id,'yeast9_id']}
    metsList.append(met)
    add_mets_Reaction.add_metabolites({met:1})
panmodel_new.add_reactions([add_mets_Reaction])



# 2. Find new genes to add
yeast9_genesIDlist=[g.id for g in yeast9.genes]
panmodel_genesIDlist=[g.id for g in panmodel_new.genes]
new_geneIDList=[i for i in yeast9_genesIDlist if i not in panmodel_genesIDlist]


# 3. Find new reactions to add
new_rxnsList=[]
gpr_update_rxnsList=[]
conflict_rxnsList=[]
panmodel_rxnsIDlist=[r.id for r in panmodel_new.reactions]
for rxn in yeast9.reactions:
    id=rxn.id
    name=rxn.name
    gpr_string=rxn.gene_reaction_rule
    gpr_list=gpr_to_list(gpr_string)
    if id not in panmodel_rxnsIDlist:
        new_rxnsList.append(id)
    else:
        # check if different reaction with the same id
        pan_name=panmodel_new.reactions.get_by_id(id).name
        if name != pan_name:
            conflict_rxnsList.append(id)
            print(id,name)
            print(pan_name)
        else:
            # check if the same reaction need to update gpr
            if len(rxn.genes)==0:
                continue
            pan_gpr_string=panmodel_new.reactions.get_by_id(id).gene_reaction_rule
            pan_gpr_list=gpr_to_list(pan_gpr_string)

            if len(set(gpr_list).difference(set(pan_gpr_list)))>0:
                gpr_update_rxnsList.append(id)
                # print(id,gpr_string)
                # print(pan_gpr_string)

# remove first 14 conflict rxns
conflict_rxnsList=conflict_rxnsList[14:]
toadd_rxnsIDList=new_rxnsList+conflict_rxnsList

print(f'{len(toadd_rxnsIDList)} new rxns need to add, {len(gpr_update_rxnsList)} exist rxn gpr need update.')
# 149 reactions need to be add, and 66 existed reactions's GPR need to be updated

# add 149 new rxns
to_add_rxns_info=dict()
i=0
for id in toadd_rxnsIDList:
    if id in conflict_rxnsList:
        pan_id='r_'+str(4780+i)
        i+=1
    else:
        pan_id=id
    to_add_rxns_info[pan_id]=dict()
    rxn=yeast9.reactions.get_by_id(id)
    name=rxn.name
    gpr_string=rxn.gene_reaction_rule
    annotation=rxn.annotation
    bounds=rxn.bounds
    subsystem=rxn.subsystem
    mets=dict()
    for key,value in rxn.metabolites.items():
        yeast9_metid=key.id
        pan_metid=yeast9_to_metsID_mapdict[yeast9_metid]
        mets[pan_metid]=value

    to_add_rxns_info[pan_id]['yeast9_id']=id
    to_add_rxns_info[pan_id]['name']=name
    to_add_rxns_info[pan_id]['gpr']=gpr_string
    to_add_rxns_info[pan_id]['annotation']=annotation
    to_add_rxns_info[pan_id]['bounds']=bounds
    to_add_rxns_info[pan_id]['subsystem']=subsystem
    to_add_rxns_info[pan_id]['metabolites']=mets


old_rxnnum=len(panmodel_new.reactions)
for id in to_add_rxns_info.keys():
    rxn=Reaction(id=id,name=to_add_rxns_info[id]['name'])
    rxn.gene_reaction_rule=to_add_rxns_info[id]['gpr']
    rxn.annotation=to_add_rxns_info[id]['annotation']
    rxn.bounds=to_add_rxns_info[id]['bounds']
    rxn.subsystem=to_add_rxns_info[id]['subsystem']
    for metid in to_add_rxns_info[id]['metabolites'].keys():
        met=panmodel_new.metabolites.get_by_id(metid)
        rxn.add_metabolites({met:to_add_rxns_info[id]['metabolites'][metid]})
    rxn.notes={'yeast9_id':to_add_rxns_info[id]['yeast9_id']}
    panmodel_new.add_reactions([rxn])
add_rxnnum=len(panmodel_new.reactions)
print(f'add {add_rxnnum-old_rxnnum} reactions')

# panmodel_new.slim_optimize()

# update 66 GPRs
old_genenum=len(panmodel_new.genes)
i=0
for id in gpr_update_rxnsList:
    rxn=panmodel_new.reactions.get_by_id(id)
    old_gpr_string=rxn.gene_reaction_rule
    new_gpr_string=yeast9.reactions.get_by_id(id).gene_reaction_rule
    # if no gene in old gpr, replace with new gpr
    if len(rxn.genes)==0:
        # rxn.gene_reaction_rule=new_gpr_string
        gpr=GPR.from_string(new_gpr_string)
        rxn.gpr=gpr
        rxn.gpr.update_genes()
        if rxn.gene_reaction_rule!=old_gpr_string:
            i+=1
            print(f'Update {id} GPR: from {old_gpr_string} to {rxn.gene_reaction_rule}')

    # if there are genes in old gpr, add new genes to gpr
    else:
        old_gpr_list=gpr_to_list(old_gpr_string)
        new_gpr_list=gpr_to_list(new_gpr_string)
        to_add_gpr_list=list(set(new_gpr_list).difference(set(old_gpr_list)))
        new_gpr_list=old_gpr_list+to_add_gpr_list
        if len(new_gpr_list)>1:
            # add ( & ) to each element, if it doesn't have and 'and' exists in the element
            new_gpr_list_modified=[]
            for element in new_gpr_list:
                if 'and' in element and '(' not in element:
                    element='('+element+')'
                new_gpr_list_modified.append(element)
            new_gpr_list=new_gpr_list_modified
        new_gpr_string=' or '.join(new_gpr_list)
        gpr=GPR.from_string(new_gpr_string)
        rxn.gpr=gpr
        rxn.gpr.update_genes()
        if rxn.gene_reaction_rule!=old_gpr_string:
            i+=1
            print(f'Update {id} GPR: from {old_gpr_string} to {rxn.gene_reaction_rule}')
new_genenum=len(panmodel_new.genes)
print(f'{i} rxn gpr were updated, and {new_genenum-old_genenum} new genes were added')

# 4. Summary: add 121 new mets, 23 new genes, 149 new reaction and update 66 reactions' GPR
# remove add_new_mets pseudo reaction
panmodel_new.remove_reactions(reactions=['add_new_mets'],remove_orphans=True)
# update new genes information
for geneID in new_geneIDList:
    try:
        gene=panmodel_new.genes.get_by_id(geneID)
        gene.annotation=yeast9.genes.get_by_id(geneID).annotation
    except:
        print(f'Error: {geneID} gene not found')

print(f'old panmodel: {len(panmodel.metabolites)} metabolites, {len(panmodel.reactions)} reactions, {len(panmodel.genes)} genes')
print(f'updated panmodel: {len(panmodel_new.metabolites)} metabolites, {len(panmodel_new.reactions)} reactions, {len(panmodel_new.genes)} genes')

# save updated panmodel
from cobra.io import write_sbml_model
write_sbml_model(panmodel_new,'model/panYeast.xml')
