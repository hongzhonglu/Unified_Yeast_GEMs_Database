# -*- coding: utf-8 -*-
# date : 2024/10/23 
# author : wangh
# file : plot_network.py
# project : Unified_Yeast_GEMs_Database
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


G = nx.Graph()


G = nx.complete_graph(10)
pos = nx.shell_layout(G)
nx.draw(G, pos=pos)  # Draw the original graph
# Draw a subgraph, reusing the same node positions
nx.draw(G.subgraph([0, 1, 2]), pos=pos, node_color="red")


G = nx.karate_club_graph()
nx.draw_circular(G, with_labels=True)
plt.show()

G = nx.path_graph(10)
shells = [[0,7,8,9], [1, 2, 3],[4,5,6]]
nx.draw_shell(G, nlist=shells)
plt.show()