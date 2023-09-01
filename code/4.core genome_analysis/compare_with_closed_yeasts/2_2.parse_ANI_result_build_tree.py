# -*- coding: utf-8 -*-
# date : 2023/6/12 
# author : wangh
# file : 2_2.parse_ANI_result_build_tree.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import os
import tqdm

os.chdir(r"D:/code/github/Unified_Yeast_GEMs_Database_from_13pro/Unified_Yeast_GEMs_Database")

def parse_fastani_result(file_name,file_dir):
    '''parse the ANI triangle matrix result from fastANI, and return a full matrix in dataframe
    '''
    triangle_matrix = []
    with open(file_dir+file_name, 'r') as f:
        for line in f.readlines():
            if "/mnt/e/data/343_yeast/0_332yeast_genomes/332_genome_assemblies/332_genome_assemblies/" not in line:
                continue
            line=line.strip('\n')
            line=line.replace("/mnt/e/data/343_yeast/0_332yeast_genomes/332_genome_assemblies/332_genome_assemblies/","")
            ani_list=line.split('\t')+[100]
            triangle_matrix.append(ani_list)

    # 在triangle_matrix最前端加上一行
    strainList=[i[0] for i in triangle_matrix]

    # 去除所有列表第一个位置元素
    triangle_matrix=[i[1:] for i in triangle_matrix]

    n=len(triangle_matrix)
    full_matrix = [[0] * n for i in range(n)]
    for i in range(n):
        for j in range(i+1):
            full_matrix[i][j] = triangle_matrix[i][j]
            full_matrix[j][i] = triangle_matrix[i][j]

    # build dataframe
    df = pd.DataFrame(full_matrix,columns=strainList,index=strainList)
    # replace all "NA" to 0 in df
    df=df.replace("NA",0)
    # set all values into float
    df=df.astype(float)
    return df

df_7yeasts_aniMatrix=parse_fastani_result(r"7yeasts_ANI.matrix",r"code/4.pan-genome_analysis/compare_with_closed_yeasts/output/")

# convert aniMatrix to distanceMatrix
df_7yeasts_distanceMatrix=100-df_7yeasts_aniMatrix

# use distanceMatrix to build tree
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
import matplotlib.pyplot as plt
import seaborn as sns

# calculate the linkage matrix
linkage_matrix = linkage(df_7yeasts_distanceMatrix, 'ward')
tree=to_tree(linkage_matrix)

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

newick=get_newick(tree, tree.dist, df_7yeasts_distanceMatrix.index.tolist())

# write tree to newick file
with open("code/4.pan-genome_analysis/compare_with_closed_yeasts/output/7yeasts_tree.newick","w") as f:
    f.write(newick)


