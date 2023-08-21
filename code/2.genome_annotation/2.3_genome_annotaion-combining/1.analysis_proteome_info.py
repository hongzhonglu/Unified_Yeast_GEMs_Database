# -*- coding: utf-8 -*-
# date : 2022/9/23 
# author : wangh
# file : 1.analysis_proteome_info.py
# project : Unified_Yeast_GEMs_Database_from_13pro


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

infile = r'data\genome\predicted_allcds\combine_proteomes_report.txt'

genomes = list()
data = dict()
for line in open(infile):
    line = line.strip()
    if len(line) < 1: continue
    if line.startswith('>'):
        genome = line.strip()[1:]
        genomes.append(genome)
        continue
    num = line.split(':')[-1].strip()
    num = int(num)
    if 'liftOver' in line:
        data['liftOver'] = data.get('liftOver', []) + [num]
        continue

    if 'MAKER2' in line:
        data['MAKER'] = data.get('MAKER', []) + [num]
        continue

    if 'Overlapped' in line:
        data['Overlapped'] = data.get('Overlapped', []) + [num]
        continue

    if 'Combined' in line:
        data['Combined'] = data.get('Combined', []) + [num]
        continue
df = pd.DataFrame(data=data, index=genomes)
df.to_csv('data/genome/predicted_allcds/combine_reports.csv')

