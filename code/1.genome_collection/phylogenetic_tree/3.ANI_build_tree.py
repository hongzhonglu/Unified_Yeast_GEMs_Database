# -*- coding: utf-8 -*-
# date : 2023/6/14 
# author : wangh
# file : 4_4.ANI_build_tree.py
# project : Unified_Yeast_GEMs_Database
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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


# load aniMatrix
aniMatrix=pd.read_csv("code/1.genome_collection/outputs/sce1900_aniMatrix.csv",index_col=0)
distanceMatrix=100-aniMatrix
# calculate the linkage matrix
linkage_matrix = linkage(distanceMatrix, 'ward')
tree=to_tree(linkage_matrix)

newick=get_newick(tree, tree.dist,distanceMatrix.index.tolist())

# write tree to newick file
with open("code/1.genome_collection/outputs/sce1900_tree.newick","w") as f:
    f.write(newick)