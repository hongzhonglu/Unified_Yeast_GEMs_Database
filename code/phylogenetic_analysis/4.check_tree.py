# -*- coding: utf-8 -*-
# date : 2024/10/16 
# author : wangh
# file : 4.check_tree.py
# project : Unified_Yeast_GEMs_Database
from ete3 import Tree
import pandas as pd
import numpy as np

species_tree=Tree('code/phylogenetic_analysis/output/coreORFs_fasttree.tre'
                  ,quoted_node_names=False, format=2)    # set format: 0:flexible with support values; 1:flexible with internal node names


# count how many leaf nodes
len(species_tree)

# extract the support values for all internal nodes
internal_nodes_support={}
leaf_nodes_support={}
i=0
for node in species_tree.traverse():
    if not node.is_leaf():
        internal_nodes_support[i]=node.support
        node.name=i
        i+=1
    else:
        leaf_nodes_support[node.name]=node.support

df_internal_nodes_support=pd.DataFrame(internal_nodes_support.items(), columns=['node','support'])
df_leaf_nodes_support=pd.DataFrame(leaf_nodes_support.items(), columns=['node','support'])

df_internal_nodes_support['support'].describe()
# plot the distribution
from matplotlib import pyplot as plt
df_internal_nodes_support['support'].plot.hist(bins=100)
plt.show()

# threshold
for threshold in np.linspace(0.9,0.99,10):
    print(df_internal_nodes_support[df_internal_nodes_support['support']>threshold].shape[0])

