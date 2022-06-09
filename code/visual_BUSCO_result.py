# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：visual_BUSCO_result.py
# 2022/6/8
'''visualize the BUSCO results to evaluate the genome annotation process'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lg1392_busco_result=pd.read_csv("result/lg1392_busco_scores.txt")
new600_busco_resullt=pd.read_csv("result/new600_busco_scores.txt")


# lg1392_genomes annotated by AUGUSTUS+liftover
plt.figure(figsize=(13,2.5))
plt.suptitle("lg1392_genome_completeness_check",fontsize=16)
plt.subplot(1,4,1)
plt.ylim((0,100))
sns.histplot(lg1392_busco_result['C(%)'],color="darkorange",stat="percent",binrange=(80,100),binwidth=0.25)
plt.xlabel('Complete BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')


plt.subplot(1,4,2)
plt.xlim((0,10))
plt.ylim((0,100))
sns.histplot(lg1392_busco_result['F(%)'],color="darkorange",stat="percent",binrange=(0,20),binwidth=0.2)
plt.xlabel('Fragmented BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')

plt.subplot(1,4,3)
plt.xlim((0,20))
plt.ylim((0,100))
sns.histplot(lg1392_busco_result['M(%)'],color="darkorange",stat="percent",binrange=(0,40),binwidth=0.5)
plt.xlabel('Missing BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')

plt.subplot(1,4,4)
plt.xlim((0,30))
plt.ylim((0,100))
sns.histplot(lg1392_busco_result['D(%)'],color="darkorange",stat="percent",binrange=(0,30),binwidth=0.5)
plt.xlabel('Complete and duplicated BUSCOs(%)')
plt.ylabel('Percent of strains(%)')

plt.tight_layout()
# plt.savefig('result/figures/lg_busco_scores.png')
plt.show()


# new600_genomes annotated by MAKER2
plt.figure(figsize=(13,2.5))
plt.suptitle("new600_genome_completeness_check",fontsize=16)
plt.subplot(1,4,1)
plt.ylim((0,100))
sns.histplot(new600_busco_resullt['C(%)'],color='deepskyblue',stat="percent",binrange=(80,100),binwidth=0.25)
plt.xlabel('Complete BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')

plt.subplot(1,4,2)
plt.xlim((0,10))
plt.ylim((0,100))
sns.histplot(new600_busco_resullt['F(%)'],color='deepskyblue',stat="percent",binrange=(0,20),binwidth=0.2)
plt.xlabel('Fragmented BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')

plt.subplot(1,4,3)
plt.xlim((0,20))
plt.ylim((0,100))
sns.histplot(new600_busco_resullt['M(%)'],color='deepskyblue',stat="percent",binrange=(0,40),binwidth=0.5)
plt.xlabel('Missing BUSCOs (%)')
plt.ylabel('Percent of strains(%)')
# plt.yscale('log')

plt.subplot(1,4,4)
plt.xlim((0,30))
plt.ylim((0,100))
sns.histplot(new600_busco_resullt['D(%)'],color='deepskyblue',stat="percent",binrange=(0,30),binwidth=0.5)
plt.xlabel('Complete and duplicated BUSCOs(%)')
plt.ylabel('Percent of strains(%)')

plt.tight_layout()
# plt.savefig('result/figures/lg_busco_scores.png')
plt.show()