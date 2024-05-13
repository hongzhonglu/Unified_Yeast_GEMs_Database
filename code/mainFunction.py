# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-


# Import packages
import re
import numpy as np
import pandas as pd
import os    ##for directory
import sys
import pprint



'''general function for easy use of python'''
def splitAndCombine(gene, rxn, sep0, moveDuplicate=False):
    ## one rxn has several genes, this function was used to splite the genes
    ## used for the dataframe data

    gene = gene.fillna('NA')  # fill the NaN with 'NA'
    gene0 = gene.tolist()
    rxn0 = rxn.tolist()
    s1 = list()
    s2 = list()
    for i in range(len(gene0)):
        s1 = s1 + [rxn0[i]] * len(gene0[i].split(sep0))
        s2 = s2 + gene0[i].split(sep0)
    df0 = pd.DataFrame({'V1': s1,
                        'V2': s2}
                       )
    if moveDuplicate == True:
        df00 = df0.drop_duplicates()
    else:
        df00 = df0
    return df00


def getSimilarTarget(rxn_yeast0,rxn_newGPR0,ss):
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    rxn_yeast1 = np.array(rxn_yeast0)  # np.ndarray()
    rxn_yeast2 = rxn_yeast1.tolist()
    rxn_yeast3 = pd.Series((v[0] for v in rxn_yeast2))
    rxn_newGPR1 = np.array(rxn_newGPR0)  # np.ndarray()
    rxn_newGPR2 = rxn_newGPR1.tolist()
    rxn_newGPR3 = pd.Series((v[0] for v in rxn_newGPR2))
    similarTarget = [None] * ss
    for i in range(ss):
        similarTarget[i] = process.extract(rxn_newGPR3[i], rxn_yeast3, limit=2)

    return similarTarget
'''
#example
newMet = pd.read_excel('new metabolite for check.xlsx')
newMet0 = newMet[['name_unify']]
gemMet = pd.read_excel('unique metabolite in yeastGEM.xlsx')
gemMet0 = gemMet[['Description_simple']]
ss0 = len(newMet0)
similarTarget0 = getSimilarTarget(gemMet0,newMet0,ss=ss0)
'''


def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result
'''
w=['a','b','c']
v=[1,2,3]
s=[3,1,2,4]
singleMapping(w,v,s,dataframe=False)
'''

def multiMapping (description, item1, item2, dataframe=True, sep=";", removeDuplicates=True):
    """get multiple description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    #used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    result = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index0 = [description[index] for index in range(len(item1)) if item1[index] == item2[i]]
            if removeDuplicates:
                index1 = pd.unique(index0).tolist()
            else:
                index1 = index0
            result[i] = sep.join(str(e) for e in index1) #string cat
        else:
            result[i] = None
    return result

'''
# example data to test all the above function
df1 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2}
                   )

df2 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'D' : np.random.randn(12)})


df2['C'] = singleMapping(df1['C'], df1['A'], df2['A'])
df2['C'] = multiMapping(df1['C'], df1['A'], df2['A'])
'''


def updateOneColumn(df1, df2, key0, value0):
    """
    using dataframe df2 to update the df1

    :param df1:
    :param df2:
    :param key0: the common column name, a string, used for the mapping
    :param value0: the column in df2 used to update the df1
    :return:
    example
    df10 = pd.DataFrame({'A': ['a', 'b', 'c'],
                 'B': ['x', 'y', 'z']})

    df20 = pd.DataFrame({'A':['c','b'],
                       'B': ['e', 'd']})
    updateOneColumn(df10,df20,key0='A',value0='B')
    """
    df10 = df1.copy()
    df11 = df1.copy()
    df10[value0] = multiMapping(df2[value0], df2[key0], df10[key0])
    for i, x in df10.iterrows():
        print(x[value0])
        if x[value0] is None:
            df11[value0][i] = df11[value0][i]
        else:
            df11[value0][i] = df10[value0][i]
    return df11[value0]






def RemoveDuplicated(s1):
    """
    example:
    s1=['a // a', 'b // a', None, 'non']

    """
    s2=list()
    for x in s1:
        print(x)
        if x =='non':
            s2.append('')
        elif x is None:
            s2.append('')
        else:
            if "//" in x:
                s0= x.split(' // ')
                s0 = [x.strip() for x in s0]
                s01= list(set(s0))
                if len(s01)==1:
                    s2.append(s01[0])
                else:
                    s2.append(' // '.join(s01))
            else:
                s2.append(x)
    return s2



def nz(value):

    '''
    Convert None to string else return value.
    '''

    if value == None:
        return 'none'
    return value



def AutoUpdate(description1, para1, description2, para2):
    # using the description1 in para1 to update the description2 in para2
    description1 = description1.tolist()
    para1 = para1.tolist()
    description2 = description2.tolist()
    para2 = para2.tolist()
    ss = [None]*len(para2)
    for i in range(len(para2)):
       if para2[i] in para1:
          ss[i] = para1.index(para2[i])
       else:
          ss[i] = None

    for i in range(len(para2)):
        if ss[i] != None:
            description2[i] = description1[ss[i]]
        else:
            description2[i] = description2[i]

    return description2



'''
# example data to test the followed function
df1 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2}
                   )

df2 = df1.iloc[[1,2]]
df2['C'] = ['good','good']
df1['C'] = AutoUpdate(df2['C'],df2['A'],df1['C'],df1['A'])
'''


def calculateFrequency(list0, item0):
    '''
    This function is used to calculate the frequency occured in a list and turn the frequency list into a dataframe
    :param list0:  ['a','b','a']
    :param item0:
    :return: a dataframe with two columns
    '''
    summary = pd.Series(list0).value_counts()
    summary = summary.to_frame(name='number')
    summary.index.name = item0
    summary.reset_index(inplace=True)
    return summary






"""function for model part"""
from cobra.manipulation import remove_genes
def getStrainGEMrxn(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the rxn list for each new reaction3

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    rxn = []
    for x in newModel.reactions:
        rxn.append(x.id)
    return rxn

def getStrainGEM(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the rxn list for each new reaction3

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    return newModel

# n=0
# for i in gene_exist:
#     if i==1.0:n+=1
# print(n)





def getRemoveGeneList(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the gene list removed from each strain specific model

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    gene = []
    for x in newModel.genes:
        gene.append(x.id)

    gene_remove_from_model = list(set(templateGene)-set(gene))
    return  gene_remove_from_model




def updateGPR(gpr0, nameMapping):
    '''
    This function is used to update the gpr reaction only with 'or' relation. It is used to replace the old gene name using
    the new gene name. Also it did not remove the duplicated value.
    :param: gpr0
    :nameMapping: a dataframe contains the mapping relation between the old and new gene name, has two columns-'geneID', 'panID'
    :return: gpr with the replaced new gene name
    '''
    #this function is mainly used to update the gene relation with 'or'
    s1 = gpr0
    s2 = s1.split(' ')
    s3 = singleMapping(nameMapping['panID'].tolist(),nameMapping['geneID'].tolist(),s2, dataframe=False)
    for i, x in enumerate(s3):
        if x is None:
            s3[i]=s2[i]
        else:
            s3[i] = s3[i]
    s4 = ' '.join(s3)
    return s4



def getCompartment(rxn):
    """
    This function is used to obtain the compartment information from reaction of yeastGEM
    :param rxn:  example acetyl-CoA[m] + L-glutamate[m]  -> coenzyme A[m] + H+[m] + N-acetyl-L-glutamate[m]'
    :return:
    """
    cp1 = ['[c]','[ce]','[e]','[er]','[erm]','[g]','[gm]','[lp]','[m]','[mm]','[n]','[p]','[v]','[vm]']
    cp2 = ['cytoplasm','cell envelope','extracellular','endoplasmic reticulum','endoplasmic reticulum membrane','Golgi','Golgi membrane','lipid particle',
             'mitochondrion','mitochondrial membrane','nucleus','peroxisome','vacuole','vacuolar membrane']


    cp = [None]*len(cp1)
    for i in range(len(cp1)):
       if cp1[i] in rxn:
         cp[i] = cp2[i]
       else:
          cp[i] = None
    cp1 = [x for i,x in enumerate(cp) if x is not None]
    cp0 = ';'.join(str(e) for e in cp1)
    return cp0



def getCommonCompartment(c1,c2, sep0=";"):
    '''this function could get the common part between string c1 and c2
    for example, c1="a;b", c2="a;c" '''
    if c1 is None:
        c10 = 'NONE'
    else:
        c10 = c1.split(sep0)
        c10 = [x.strip() for x in c10]
    if c2 is None:
        c20 = 'NONE'
    else:
        c20 = c2.split(sep0)
        c20 = [x.strip() for x in c20]
    c3 = list(set(c10).intersection(c20))
    c4 = sep0.join(str(e) for e in c3)
    return c4


def getRXNgeneMapping(rxn0, gpr0):
    '''this function is used to split the GPR;
    input, for example rxn0=['r1','g2']
    gpr0=['a or c','a and b']
    output, each rxn related with each gene'''
    s1 = rxn0
    s2 = gpr0
    s2 = s2.str.replace('and','@')
    s2 = s2.str.replace('or','@')
    s2 = s2.str.replace('\\( ','')
    s2 = s2.str.replace('\\(\\( ','')
    s2 = s2.str.replace('\\(', '')
    s2 = s2.str.replace('\\(\\(', '')
    s2 = s2.str.replace(' \\)','')
    s2 = s2.str.replace(' \\)\\) ','')
    s2 = s2.str.replace('\\)', '')
    s2 = s2.str.replace('\\)\\) ', '')
    s3 = splitAndCombine(s2,s1,sep0="@")
    s3['V2'] = s3['V2'].str.strip()
    s3.columns = ['rxnID', 'gene']
    return s3

def getRXNmetaboliteMapping(rxn0, met0):
    '''this function is used to split the equation of metabolites; used to produce the dataframe format of GEM using
    cobrapy
    input, for example rxn0=['r1','g2']
    gpr0=['a => c','a => b']
    output, each rxn related with each gene'''
    met_annotation = pd.read_excel('../data/met_panYeast_v3.xlsx')
    s1 = rxn0
    s2 = met0
    s3 = splitAndCombine(s2,s1,sep0=" ")
    s3['V2'] = s3['V2'].str.strip()
    s3.columns = ['rxnID', 'met']
    s3['met_name'] = singleMapping(met_annotation['description'],met_annotation['m_name'],s3['met'])
    for i, x in s3.iterrows():
        if s3['met_name'][i] is None:
            s3['met_name'][i] = s3['met'][i]
        else:
            s3['met_name'][i] = s3['met_name'][i]
    return s3


def correctSomeWrongFormat(model0):
  """
  This function is used to correct some wrong format when read yeastGEM model from cobratoolbox
  """
  # Correct metabolite ids:
  for met in model0.metabolites:
    met.id = met.id.replace('__91__', '_')
    met.id = met.id.replace('__93__', '')
  #for reaction in model0.reactions:
  #    reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace('__45__','-')
  for gene in model0.genes:
      gene.id = gene.id.replace('__45__', '-')

  return model0


def produceMetaboliteList(model0):
  #produce the dataframe for the metabolites from yeastGEM
  met_list =[None]*len(model0.metabolites)
  met_dataframe = pd.DataFrame({'m_name':met_list,
                              'description':met_list,
                              'formula':met_list})

  for i, met in enumerate(model0.metabolites):
      print(i)
      met_dataframe['m_name'][i] = met.id
      met_dataframe['description'][i] = met.name
      met_dataframe['formula'][i] = met.formula
  #s2 = met_dataframe['m_name'].str.split('_', expand=True)
  #met_dataframe['description'] = met_dataframe['description'].str.replace('\s\[', '@')
  #s3 = met_dataframe['description'].str.split('@', expand=True)
  #met_dataframe['description'] = s3.iloc[:, 0] + '[' + s2.iloc[:, 2] + ']'
  return met_dataframe


def produceGeneList(model0):
  #produce the gene list from GEM
  genelist = []
  for i in model0.genes:
      print(i)
      genelist.append(i.id)
  return genelist


def produceRxnList(model0):
  #produce the dataframe for the rxn from yeastGEM
  reaction_list =[None]*len(model0.reactions)
  gem_dataframe = pd.DataFrame({'name':reaction_list,
                              'equation':reaction_list,
                              'GPR':reaction_list,
                              'rxnID':reaction_list,
                              'formula':reaction_list
                              })

  for i, reaction in enumerate(model0.reactions):
      print(i)
      gem_dataframe['name'][i] = reaction.name
      gem_dataframe['equation'][i] = reaction.reaction
      gem_dataframe['GPR'][i] = reaction.gene_reaction_rule
      gem_dataframe['rxnID'][i] = reaction.id
  gem_dataframe['ID'] = ['R'+ str(i) for i in range(0, len(model0.reactions))]
  gem_dataframe['GPR'] = gem_dataframe['GPR'].str.replace('__45__', '-')
  #replace the metabolite name in gem_dataframe
  s0 = getRXNmetaboliteMapping(rxn0=gem_dataframe['rxnID'], met0=gem_dataframe['equation'])
  # map
  gem_dataframe['formula'] = multiMapping(description=s0['met_name'],item1=s0['rxnID'],item2=gem_dataframe['rxnID'],removeDuplicates=False)
  gem_dataframe['formula'] = gem_dataframe['formula'].str.replace(";", " ")
  return gem_dataframe



def exchange(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem


def exchange_ecYeast(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem




















#SLIME rxn
def SLIME(rxnName, subsystem):
    """
    if the rxnName contains the SLIME, classify the reaction into SLIME reaction
    """
    for i,x in enumerate(rxnName):
        if 'SLIME' in x:
            subsystem[i] = 'SLIME reaction'
            print(subsystem[i])
        else:
            subsystem[i] = subsystem[i]
    return subsystem



def transport(s1, subsysem):
    """
    this function is used to define the transport reaction
    #example
     s1 =['2-methylbutyl acetate [cytoplasm] --> 2-methylbutyl acetate [extracellular]', 'H+ [extracellular] + phosphoenolpyruvate [extracellular] <=> H+ [cytoplasm] + phosphoenolpyruvate [cytoplasm]']
     subsysem = ['a','b']

    :param s1:
    :param subsysem:
    :return:
    """
    for i, x0 in enumerate(s1):
        x1 = re.findall(r"\[([A-Za-z0-9_\s]+)\]", x0)
        x0 = x0.replace('(','[')
        x0 = x0.replace(')',']')
        x2 = re.sub(r"\[([A-Za-z0-9_\s+]+)\]", '', x0)

        if "<=>" in x2:
            x3 = x2.split("<=>")
        elif "<->" in x2: #bigg database format
            x3 = x2.split("<->")
        else:
            x3 = x2.split("-->")
        x3 = [x.strip() for x in x3]
        x1=pd.unique(x1).tolist() #remove the duplicated
        if '+' in x3[0]:
            x30=x3[0].split('+')
        else:
            x30=x3[0]
        x30=[x.strip() for x in x30]
        x30 = [x for x in x30 if x != '']
        if '+' in x3[1]:
            x31 = x3[1].split('+')
        else:
            x31=x3[1]
        x31 = [x.strip() for x in x31]
        x31 = [x for x in x31 if x != '']

        if set(x30) == set(x31):
            subsysem[i] ='Transport' + '['+', '.join(x1)+']'
            print(subsysem[i])
        elif set(x30)-set(['ATP','H2O']) == set(x31) - set(['ADP','phosphate','H']):
            subsysem[i] = 'Transport' + '[' + ', '.join(x1) + ']'
            print(subsysem[i])
        else:
            subsysem[i] = subsysem[i]
    return subsysem



def findRemoveRxnBasedOnGene(rxnRemovedGene, rxnAllGene):
    '''this function is used to remove rxn based on the removed gene list
    if the all genes in a reaction were in the removed gene list, then this reaction was removed'''
    #x0 = gem_dataframe['removed_gene'].tolist()
    #y0 = gem_dataframe['all_gene'].tolist()
    x0=rxnRemovedGene.tolist()
    y0=rxnAllGene.tolist()
    removed_rxn = list()
    for x,y in zip(x0,y0):
        if x is None:
            removed_rxn.append('NO')
        else:
            if len(x) ==len(y):
                removed_rxn.append('YES')
            else:
                removed_rxn.append('NO')
    return removed_rxn



def saveExcel(infile, outfile):
    '''
    function to save the dataframe into xlsx format
    :param infile:
    :param outfile:
    :return:
    '''
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer,'Sheet1')
    writer.save()


def find(rxn_name_list, specific_string, equal=False):
    '''
    function to find the index of element in a list contains a specific string, refer the function
    from matlab
    :param rxn_name_list:
    :param specific_string:
    :param equal: if true will return the indexes in rxn_name_list which is equal to specific_string
    :return: index of the element from the list contains the specific string
    example:
    rxn_name_list = ['a_b','e_b','c','d']
    specific_string = '_b'
    '''
    s0 = rxn_name_list
    if equal ==False:
        index =[i for i,x in enumerate(s0) if specific_string in x]
    else:
        index = [i for i, x in enumerate(s0) if specific_string == x]
    return index


def sampling_to_test(fileDir, tarDir, sample_ratio):
    '''sample a part number of the ssGEMs for building analysing pipeline
    fileDir: path to the files storage dict
    tarDir: path to the samples storage dict
    sample_ratio: haw many samples used to test
    usage: sampling_to_test(fileDir='result/ssGEMs/',tarDir='result/test/ssGEMs/',sample_ratio=0.1)
    '''
    import os, random, shutil
    pathDir = os.listdir(fileDir)
    filenumber = len(pathDir)
    picknumber = int(filenumber * sample_ratio)
    sample = random.sample(pathDir, picknumber)
    print(sample)
    for name in sample:
        shutil.copy(fileDir + name, tarDir + name)
    return


def anaerobicsimulation(model):
    '''Note that the reaction coefficients can't be directly changed to 0 by present code,so the function need
    some parameters tuning for specific yeast model
    reference:https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/anaerobicModel.m'''
    # 切换anaerobic growth medium
    from cobra.io import read_sbml_model
    import os
    import pandas as pd

    anerobic_model = model.copy()

    # 1. change: Refit GAM and NGAM to exp. data, change biomass composition
    GAM=30.49
    P=0.461     #change protein fraction——unfinished
    NGAM=0
    NGAM_id = 'r_4046'
    GAM_id = 'r_4041'
    anerobic_model.reactions.get_by_id(NGAM_id).bounds = NGAM,NGAM
    biomass_rxn=anerobic_model.reactions.get_by_id('r_4041')
    biomass_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_0434[c]'): 55.4 - 30.49,
                             anerobic_model.metabolites.get_by_id('s_0803[c]'): 55.4 - 30.49,
                             anerobic_model.metabolites.get_by_id('s_0394[c]'): -55.4 + 30.49,
                             anerobic_model.metabolites.get_by_id('s_0794[c]'): -55.4 + 30.49,
                             anerobic_model.metabolites.get_by_id('s_1322[c]'): -55.4 + 30.49})

    # 2. Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
    cofactor_rxn = anerobic_model.reactions.get_by_id('r_4598')
    # yeastGEM8.5 co-factors reaction coeffiecients
    # cofactor_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_3714[c]'):9.99999997475243e-07,
    #                                 anerobic_model.metabolites.get_by_id('s_1198[c]'):0.00264999992214143,
    #                                 anerobic_model.metabolites.get_by_id('s_1203[c]'):0.000150000007124618 ,
    #                                 anerobic_model.metabolites.get_by_id('s_1207[c]'):0.000569999974686652,
    #                                 anerobic_model.metabolites.get_by_id('s_1212[c]'):0.00270000007003546,
    #                                 anerobic_model.metabolites.get_by_id('s_0529[c]'):0.000190000006114133 })

    # panYeast co-factors reaction coefficients
    cofactor_rxn.add_metabolites({anerobic_model.metabolites.get_by_id('s_3714[c]'): 9.99999997475243e-07,
                                  anerobic_model.metabolites.get_by_id('s_1198[c]'): 0.00100000004749745,
                                  anerobic_model.metabolites.get_by_id('s_1203[c]'): 0.000788461999036372,
                                  anerobic_model.metabolites.get_by_id('s_1207[c]'): 6.54000032227486e-05,
                                  anerobic_model.metabolites.get_by_id('s_1212[c]'): 7.6900003477931e-05,
                                  anerobic_model.metabolites.get_by_id('s_0529[c]'): 0.000190000006114133})

    # 3. change: Changes media to anaerobic————no O2；add sterol and fatty acid exchanges
    change_modium = anerobic_model.medium
    change_modium['r_1992'] = 0
    # 改变培养基成分，允许甾醇和脂肪酸的交换反应
    change_modium_lists = ['r_1757', 'r_1915', 'r_1994', 'r_2106', 'r_2134', 'r_2137', 'r_2189']
    for met in change_modium_lists:
        change_modium[met] = 1000.0
    anerobic_model.medium = change_modium
        # print(test_model.reactions.get_by_id(met).bounds)

    # 4. Blocked pathways for proper glycerol production.
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    anerobic_model.reactions.get_by_id('r_0713').lower_bound = 0    #Mithocondria
    anerobic_model.reactions.get_by_id('r_0714').lower_bound = 0    #Cytoplasm
    # Block glycerol dehydroginase (only acts in microaerobic conditions)
    anerobic_model.reactions.get_by_id('r_0487').upper_bound=0
    # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    anerobic_model.reactions.get_by_id('r_0472').upper_bound = 0

    return anerobic_model


def scale_biomass(model,component,new_value,balance_out):
    '''Scales the biomass composition
    model: metabolic model in COBRA format
    component: name of the component to rescale (e.g. "protein")——
        'carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor'
    new_value: new total fraction for target component
    balance_out: if chosen, the name of another component with which the model will be balanced out so that
                    the total mass remains = 1 g/gDW——
        'carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor'
    disp_output: if output from sumBioMass should be displayed (default = true)

    Usage: model = scale_biomass(model,component,new_value,balance_out,dispOutput)
    original_values get from SysBioChalmers/yeast-GEM repo:
    P -> 0.46 g/gDW
    C -> 0.38067 g/gDW
    R -> 0.061 g/gDW
    D -> 0.0037021 g/gDW
    L -> 0.087299 g/gDW
    I -> 0.0024815 g/gDW
    F -> 0.0048478 g/gDW
    X -> 1 gDW/gDW
    Growth = 0.083748 1/h
    reference:https://github.com/SysBioChalmers/yeast-GEM/blob/main/code/otherChanges/scaleBioMass.m
    '''

    P=0.46
    C=0.38
    R=0.061
    D=0.0037
    L=0.0873
    I=0.00248
    F=0.00485
    content_Cap = [C, P, L, R, D, I, F]
    content_all = ['carbohydrate', 'protein', 'lipid', 'RNA', 'DNA', 'ion', 'cofactor']
    for i in range(len(content_all)):
        if content_all[i]==component:
            pos=i

    if pos:
        old_value=content_Cap[pos]
        f = new_value / old_value
        model=rescale_pseudorxn(model,component,f)
        if balance_out in content_all:
            for i in range(len(content_all)):
                if content_all[i] == component:
                    pos = i
            balance_value=content_Cap[pos]
            f=(balance_value - (new_value - old_value)) / balance_value
            model = rescale_pseudorxn(model, balance_out, f)
        else:
            print('input balance_out argument is wrong!╯︿╰')
    else:
        print('input component argument is wrong!╯︿╰')

    return model


def rescale_pseudorxn(model,met_name,f):
    '''Rescales a specific pseudoreaction by a given value
    model: the YeastGEM
    met_name: name of the component to rescale(eg. 'protein')
    f: fraction to use for rescaling

    usage: model=scaleBioMass(model,met_name,f)'''

    if met_name=='lipid':
        model = rescale_pseudorxn(model, 'lipid backbone', f)
        model = rescale_pseudorxn(model, 'lipid chain', f)
    else:
        rxn_name=met_name+' pseudoreaction'
        for reaction in model.reactions:
            if reaction.name==rxn_name:
                target_pseudorxn=reaction
                break

        for met,old_value in target_pseudorxn.metabolites.items():
            if old_value!=1:
                change=f-1
                target_pseudorxn.add_metabolites({met: old_value * change})

        return model


def rename_ssGEM_id(id=False,ids_list=False):
    '''modify ssGEM id, make it correspond to the all_strain_infomation'''
    if id:
        print("modify one ssGEM id,and reture a id str")
        new_id=id.strip(".xml")
        if new_id.startswith("GCA"):
            new_id="_".join(new_id.split("_")[:2])
        else:
            new_id.replace(".re","")
            new_id.split("_")[0]
        return  new_id
    if ids_list:
        print("modeify a list of ssGEM's id,and reture a list of standardize id")
        new_ids_list=[]
        for id in ids_list:
            new_id = id.strip(".xml")
            if new_id.startswith("GCA"):
                new_id = "_".join(new_id.split("_")[:2])
            else:
                new_id=new_id.replace(".re", "")
                new_id=new_id.split("_")[0]
            new_ids_list.append(new_id)
        return new_ids_list


def check_strains_info(ssGEM_ids_list):
    '''strain_id:the strain you are interested in.
    usage: info=check_info('AAA')/check_info('AAA.xml')
    '''
    import pandas as pd

    ssGEM_ids_list_v2 = rename_ssGEM_id(id=False, ids_list=ssGEM_ids_list)
    all_strains_info=pd.read_excel('data/strain_information/1897_strains_info.xlsx',index_col=0)

    strains_info=all_strains_info.loc[all_strains_info["ssGEM"].isin(ssGEM_ids_list_v2)]
    return strains_info


def make_blast_db(seq, folder, db_type='prot'):
    '''construct blast database.
    *parameters:
    seq: query sequence used to build blast_db
    in_folder:path to query sequence
    out_folder: the path used to store db file
    db_type: 'prot'——protein sequence； or 'nucl' ——nucleotide sequence
    diamond_path:path to run diamond
    *example：
    make_blast_db(ref_id, folder='~/why/S.cerevisiae_ssGEMs_auto-construction/data/strain_genome_diamond_db', db_type='prot',diamond_path='~/why')
    '''
    import os
    from glob import glob

    # out_file = '%s/%s.dmnd' % (out_folder, seq)
    # files = glob('%s/*.dmnd' % out_folder)

    out_file = '%s.fa.pin'%seq
    files = os.listdir(folder)

    if out_file in files:
        print(seq, 'already has a blast db')
        return
    if db_type == 'nucl':
        ext = 'fna'
    else:
        ext = 'fa'

    # cmd_line = '%s/./diamond makedb --in %s/%s -d %s/%s.dmnd'%(diamond_path,in_folder, seq,out_folder,seq)   #diamond命令
    cmd_line = 'makeblastdb -in %s/%s.%s -dbtype %s ' % (folder,seq,ext, db_type)                  #ncbi blast命令
    print('making blast db with following command line...')
    print(cmd_line)
    os.system(cmd_line)


def trans_genome(strain,out_folder,in_folder):
    '''translate genome cds sequence into protein sequence'''
    import os
    from Bio import SeqIO
    in_file='%s/%s.fa'%(in_folder,strain)
    out_file='%s/%s.fa'%(out_folder,strain)
    all_strain=os.listdir(in_folder)
    print(all_strain)
    strain_name=strain+'.fa'
    if strain_name in all_strain:
        print('translating %s'%strain)
    else:
        print('can not find %s'%strain)
        print(in_file)
        return
    records=SeqIO.parse(in_file,'fasta')

    with open(out_file,'w') as output:
        for record in records:
            seq=record.seq.translate()
            seqid=record.id
            output.write('>%s\n%s\n'%(seqid,seq))
    return


# define a function to run BLASTp
def run_blastp(seq, db, in_folder='prots', out_folder='bbh', out=None, outfmt=6, evalue=0.001, threads=1):
    import os
    from glob import glob
    if out == None:
        out = '%s/%s_vs_%s.txt' % (out_folder, seq, db)
        print(out)

    files = glob('%s/*.txt' % out_folder)
    if out in files:
        print(seq, 'already blasted')
        return

    print('blasting %s vs %s' % (seq, db))

    db = '%s/%s.fa' % (in_folder, db)
    seq = '%s/%s.fa' % (in_folder, seq)
    cmd_line = 'blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
               % (db, seq, out, evalue, outfmt, threads)

    print('running blastp with following command line...')
    print(cmd_line)
    os.system(cmd_line)
    return out


# define a function to get sequence length
def get_gene_lens(query, in_folder='test/build_geneMatrix/all_cds'):
    from Bio import SeqIO
    file = '%s/%s' % (in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []

    for record in records:
        out.append({'gene': record.name, 'gene_length': len(record.seq)})

    out = pd.DataFrame(out)
    return out


# define a function to get Bi-Directional BLASTp Best Hits
def get_bbh(query, subject, in_folder='bbh'):
    from glob import glob

    out_file = '%s\\%s_vs_%s_parsed.csv' % (in_folder, query, subject)
    files = glob('%s\\*_parsed.csv' % in_folder)
    if out_file in files:
        print('%s has already done the bbh'%query)
        return
    else:
        # Utilize the defined protein BLAST function
        run_blastp(query, subject,in_folder='test/build_geneMatrix/blast_db',out_folder=in_folder)
        run_blastp(subject, query,in_folder='test/build_geneMatrix/blast_db',out_folder=in_folder)

        query_lengths = get_gene_lens(query, in_folder='test/build_geneMatrix/blast_db')
        subject_lengths = get_gene_lens(subject, in_folder='test/build_geneMatrix/blast_db')

        # Define the output file of this BLAST

        # Combine the results of the protein BLAST into a dataframe
        print('parsing BBHs for', query, subject)
        cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd',
                'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
        bbh = pd.read_csv('%s/%s_vs_%s.txt' % (in_folder, query, subject), sep='\t', names=cols)
        bbh = pd.merge(bbh, query_lengths)
        bbh['COV'] = bbh['alnLength'] / bbh['gene_length']

        bbh2 = pd.read_csv('%s/%s_vs_%s.txt' % (in_folder, subject, query), sep='\t', names=cols)
        bbh2 = pd.merge(bbh2, subject_lengths)
        bbh2['COV'] = bbh2['alnLength'] / bbh2['gene_length']
        out = pd.DataFrame()

        # Filter the genes based on coverage
        bbh = bbh[bbh.COV >= 0.4]
        bbh2 = bbh2[bbh2.COV >= 0.4]

        # Delineate the best hits from the BLAST
        for g in bbh.gene.unique():
            res = bbh[bbh.gene == g]
            if len(res) == 0:
                continue
            best_hit = res.loc[res.PID.idxmax()]
            best_gene = best_hit.subject
            res2 = bbh2[bbh2.gene == best_gene]
            if len(res2) == 0:
                continue
            best_hit2 = res2.loc[res2.PID.idxmax()]
            best_gene2 = best_hit2.subject
            if g == best_gene2:
                best_hit['BBH'] = '<=>'
            else:
                best_hit['BBH'] = '->'
            out = pd.concat([out, pd.DataFrame(best_hit).transpose()])

        # Save the final file to a designated CSV file
        out.to_csv(out_file)