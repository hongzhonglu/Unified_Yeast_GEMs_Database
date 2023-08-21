# -*- coding: utf-8 -*-
# date : 2023/6/12 
# author : wangh
# file : 4_3.merge_split_aniMatrix.py
# project : Unified_Yeast_GEMs_Database
import pandas as pd
import os
import tqdm

def parse_fastani_result(file_name,file_dir):
    '''parse the ANI triangle matrix result from fastANI, and return a full matrix in dataframe
    '''
    triangle_matrix = []
    with open(file_dir+file_name, 'r') as f:
        for line in f.readlines():
            if "/dssg/home/acct-clslhz/clslhz/why_ssGEM/data/assembled_genome/" not in line:
                continue
            line=line.strip('\n')
            line=line.replace("/dssg/home/acct-clslhz/clslhz/why_ssGEM/data/assembled_genome/","")
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
    # choose those columns which have 0 values, and change all values in these columns to 0
    df.loc[:,(df==0).any(axis=0)]=0
    # set all values into float
    df=df.astype(float)
    return df


def merge_split_aniMatrix(fastani_result_dir,fastani_resultList):
    '''merge all split fastANI result into one matrix
    '''
    df=parse_fastani_result(fastani_resultList[0],fastani_result_dir)
    aniMatrix=df.copy()
    for file in tqdm.tqdm(fastani_resultList[1:]):
        df2=parse_fastani_result(file,fastani_result_dir)
        aniMatrix=aniMatrix+df2
    return aniMatrix


fastani_result_dir=r"code/1.genome_collection/outputs/fastANI_split_result/"
matrix_fileList=[i for i in os.listdir(fastani_result_dir) if i.endswith(".matrix")]
aniMatrix=merge_split_aniMatrix(fastani_result_dir=fastani_result_dir,
                                fastani_resultList=matrix_fileList)


# check how many values > 100
print(aniMatrix[aniMatrix<95].count().sum())

# change index and columns name by remove all .fna
aniMatrix.index=[i.replace(".fna","") for i in aniMatrix.index]
aniMatrix.columns=[i.replace(".fna","") for i in aniMatrix.columns]

# save the result
aniMatrix.to_csv("code/1.genome_collection/outputs/sce1900_aniMatrix.csv")


