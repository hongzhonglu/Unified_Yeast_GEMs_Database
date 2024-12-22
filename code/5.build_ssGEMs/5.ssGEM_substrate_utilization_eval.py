# -*- coding: utf-8 -*-
# date : 2023/5/31 
# author : wangh
# file : 5.1_ssGEM_substrate_utilization_eval.py
# project : Unified_Yeast_GEMs_Database
'''simulate 3 ssGEMs different substrate utilization with abundant experiment data to prove the quality of ssGEMs
GCA_019394525.1,GCA_019394085.1,GCA_019394815.1
 '''
import pandas as pd
from cobra.io import read_sbml_model,write_sbml_model
import os
import random
import cobra
from cobra.flux_analysis import gapfill


panmodel=read_sbml_model("model/panYeast.xml")
# test_strains
ssGEM_dir="model/ssGEMs/"
test_strains=["GCA_019394525.1_ASM1939452v1_genomic.xml","GCA_019394085.1_ASM1939408v1_genomic.xml","GCA_019394815.1_ASM1939481v1_genomic.xml"]
for strain in test_strains:
    model=read_sbml_model(ssGEM_dir+strain)
    try:
        sol=gapfill(model,panmodel,lower_bound=0.04)
        print(sol)
        # model.add_reactions(sol[0])
    except:
        print('%s has trouble in gapfill' %strain)
    print(model.slim_optimize())

# according to gapfilling, r_0061,r_1838 need to be added to the 3 ssGEMs
add_rxnidList=['r_0061','r_1838']
add_rxnList=[panmodel.reactions.get_by_id(id) for id in add_rxnidList]

# 61种不同碳源
carbon_sources_0=pd.read_excel('data/panYeast_exchangeRXNs.xlsx',sheet_name='carbon_sources')
carbon_sources=carbon_sources_0['carbon sources'].values.tolist()
carbob_IDlist=carbon_sources_0['panID'].values.tolist()


# predict different carbon source utilization capacity
growth_simulation=pd.DataFrame(index=carbon_sources)
model=read_sbml_model("model/yeast-GEM.xml")
growth_values=[]
print(model.slim_optimize())
for carbon in carbob_IDlist:
        with model:
            medium = model.medium
            medium['r_1714'] = 0.0  # 移除原培养基中的glucose
            medium[carbon]=1.0
            model.medium=medium
            growth_value=model.slim_optimize()
            growth_values.append(growth_value)
growth_simulation["yeastGEM"]=growth_values

# predict 3 ssGEMs
for strain in test_strains:
    model=read_sbml_model(ssGEM_dir+strain)
    model.add_reactions(add_rxnList)
    growth_values=[]
    # model = anerobicsimulation(model)
    print(model.slim_optimize())
    for carbon in carbob_IDlist:
        with model:
            medium = model.medium
            medium['r_1714'] = 0.0  # 移除原培养基中的glucose
            medium[carbon]=1.0
            model.medium=medium
            growth_value=model.slim_optimize()
            growth_values.append(growth_value)
    growth_simulation[strain]=growth_values
strain_id_list=[i.strip(".xml") for i in test_strains]
strain_id_list.insert(0,"yeastGEM")
growth_simulation.columns=strain_id_list

# check test ssGEM size
# for strain in test_strains:
#     print(strain)
#     model=read_sbml_model(ssGEM_dir+strain)
#     print("rxn_numb:"+str(len(model.reactions)))
#     print("gene_numb:"+str(len(model.genes)))



# load giga36_data
giga36_pheno_data=pd.read_excel("data/strain_information/Giga_36_supplemental_files/Dataset_S5.xlsx",sheet_name="Carbon sources")
giga36_pheno_data.drop(labels=0,inplace=True)

# build experiment growthMatrix
columns=["CEN.PK","s288c","ethanol_red"]
for column in columns:
    giga36_pheno_data.loc[giga36_pheno_data[column]<=1,column]=0
    giga36_pheno_data.loc[giga36_pheno_data[column]>1,column]=1

growthMatrix=giga36_pheno_data[["Carbon source","CEN.PK","s288c","ethanol_red"]]
growthMatrix=growthMatrix[growthMatrix["Carbon source"].isin(carbon_sources)]
growthMatrix.set_index(keys="Carbon source",inplace=True)
# build simulation growthMatrix
simul_growthMatrix=growth_simulation
simu_columns=["yeastGEM","ethanol_red","s288c","CEN.PK"]
simul_growthMatrix.columns=simu_columns
for column in simu_columns:
    simul_growthMatrix.loc[simul_growthMatrix[column]<=0.01,column]=0
    simul_growthMatrix.loc[simul_growthMatrix[column]>0.01,column]=1
# simul_growthMatrix.fillna(0,inplace=True)
# simul_growthMatrix.loc[simul_growthMatrix['yeastGEM']>0.01,'yeastGEM']=1
# simul_growthMatrix.loc[simul_growthMatrix['yeastGEM']<=0.01,'yeastGEM']=0

# save simulation growthMatrix and experiment growthMatrix into one excel
writer=pd.ExcelWriter("result/model_simulation/61substrate_prediction.xlsx")
simul_growthMatrix.to_excel(writer,sheet_name="simulation")
growthMatrix.to_excel(writer,sheet_name="experiment")
writer.save()

# caculate the prediction accuracy
def cal_accuracy(exp_data,sim_data,exp_column,sim_column):
    exp_gro_list=exp_data[exp_data[exp_column]==1].index.tolist()
    exp_ungro_list=exp_data[exp_data[exp_column]==0].index.tolist()
    sim_gro_list=sim_data[sim_data[sim_column]==1].index.tolist()
    sim_ungro_list=sim_data[sim_data[sim_column]==0].index.tolist()
    TF=len(set(exp_gro_list).intersection(set(sim_gro_list)))
    FP=len(set(exp_ungro_list).intersection(set(sim_gro_list)))
    FN=len(set(exp_gro_list).intersection(set(sim_ungro_list)))
    TN=len(set(exp_ungro_list).intersection(set(sim_ungro_list)))
    accuracy_result=[TF,FP,FN,TN]
    return accuracy_result

cols=["TF","FP","FN","TN"]
accuracy_result={}
for strain in columns:
    ssGEM_accuracy=cal_accuracy(exp_data=growthMatrix,sim_data=simul_growthMatrix,exp_column=strain,sim_column=strain)
    accuracy_result[strain]=ssGEM_accuracy

yeast9_accuracy=cal_accuracy(exp_data=growthMatrix,sim_data=simul_growthMatrix,exp_column="s288c",sim_column="yeastGEM")
accuracy_result["yeastGEM"]=yeast9_accuracy

df_accuracy_result=pd.DataFrame.from_dict(accuracy_result,orient="index",columns=cols)
# remove s288c row
df_accuracy_result.drop(labels="s288c",axis=0,inplace=True)

# ssGEM predit results: combine 0 and 1 rows
df_accuracy_result.loc["ssGEMs(CEN.PK&ethanol red)",:]=df_accuracy_result.loc["CEN.PK",:]+df_accuracy_result.loc["ethanol_red",:]
df_accuracy_result.drop(labels=["CEN.PK","ethanol_red"],axis=0,inplace=True)
# yeastGEM predict results
df_accuracy_result.iloc[2,]

df_accuracy_result.astype(int).to_excel("result/ssGEM_simulation/ssGEM_substrates_prediction_result.xlsx")
df_accuracy_result=pd.read_excel("result/ssGEM_simulation/ssGEM_substrates_prediction_result.xlsx",index_col=0)
# calculate F1 score,accuracy,precision,sensitivity,specificity
def cal_F1score(TP,FP,FN,TN):
    precision=TP/(TP+FP)
    sensitivity=TP/(TP+FN)
    specificity=TN/(TN+FP)
    accuracy=(TP+TN)/(TP+TN+FP+FN)
    F1_score=2*precision*sensitivity/(precision+sensitivity)
    return F1_score,accuracy,precision,sensitivity,specificity

ssGEM_score=cal_F1score(TP=df_accuracy_result.loc["ssGEMs(CEN.PK&ethanol red)","TF"],
                        FP=df_accuracy_result.loc["ssGEMs(CEN.PK&ethanol red)","FP"],
                        FN=df_accuracy_result.loc["ssGEMs(CEN.PK&ethanol red)","FN"],
                        TN=df_accuracy_result.loc["ssGEMs(CEN.PK&ethanol red)","TN"])
cenpk_score=cal_F1score(TP=df_accuracy_result.loc["CEN.PK","TF"],
                        FP=df_accuracy_result.loc["CEN.PK","FP"],
                        FN=df_accuracy_result.loc["CEN.PK","FN"],
                        TN=df_accuracy_result.loc["CEN.PK","TN"])
ethanolred_score=cal_F1score(TP=df_accuracy_result.loc["ethanol_red","TF"],
                        FP=df_accuracy_result.loc["ethanol_red","FP"],
                        FN=df_accuracy_result.loc["ethanol_red","FN"],
                        TN=df_accuracy_result.loc["ethanol_red","TN"])
yeastGEM_score=cal_F1score(TP=df_accuracy_result.loc["yeastGEM","TF"],
                        FP=df_accuracy_result.loc["yeastGEM","FP"],
                        FN=df_accuracy_result.loc["yeastGEM","FN"],
                        TN=df_accuracy_result.loc["yeastGEM","TN"])
yeastGEM_score=cal_F1score(TP=yeast9_accuracy[0],FP=yeast9_accuracy[1],FN=yeast9_accuracy[2],TN=yeast9_accuracy[3])
df_score=pd.DataFrame([ssGEM_score,cenpk_score,ethanolred_score,yeastGEM_score],index=["ssGEMs(CEN.PK&ethanol red)","CEN_PK","Ethanol red","yeastGEM"],columns=["F1 score","Accuracy","Precision","Sensitivity","Specificity"])

# save df_score in another sheet of ssGEM_substrates_prediction_result.xlsx
with pd.ExcelWriter("result/ssGEM_simulation/ssGEM_substrates_prediction_result.xlsx",mode="a") as writer:
    df_score.to_excel(writer,sheet_name="F1_score")


# plot radar chart
import matplotlib.pyplot as plt
from math import pi
import seaborn as sns

sns.set_style("whitegrid")
# number of variable
categories=df_score.columns.tolist()
N = len(categories)

# What will be the angle of each axis in the plot? (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]
# Initialise the spider plot
fig,ax=plt.subplots(figsize=(8,6),subplot_kw=dict(polar=True))

ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

# Draw one axe per variable + add labels, set the label position
plt.xticks(angles[:-1], categories, size=12,fontweight="bold")
# Draw ylabels
ax.set_rlabel_position(0)
plt.yticks([0.25,0.50,0.75,1.00], size=10)
plt.ylim(0, 1)

# 在不同角度分别使用text添加[0.25,0.50,0.75,1.00]标签
rotate_dict={angles[1]:-72,angles[2]:35,angles[3]:-35,angles[4]:-288}
for i in [0.25,0.50,0.75]:
    for angle in angles[1:5]:
        ax.text(angle, i, str(i), ha="center",fontsize=10,rotation=rotate_dict[angle])

# Ind1
values = df_score.loc["yeastGEM"].values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=2, linestyle='solid', label="YeastGEM")
ax.fill(angles, values, 'red', alpha=0.1)

# Ind2
values = df_score.loc["ssGEMs(CEN.PK&ethanol red)"].values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=2, linestyle='solid', label="ssGEMs this study")
ax.fill(angles, values, 'blue', alpha=0.1)

# set legend for ax
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
# set title for ax
plt.title("ssGEMs substrates utilization prediction result",fontdict=dict(fontsize=15,fontweight="bold"))
plt.tight_layout()
# Show the graph
plt.show()

