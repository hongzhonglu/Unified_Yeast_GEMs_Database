# -*- coding: utf-8 -*-
# date : 2024/11/16 
# author : wangh
# file : chrmosome_ID_convert.py
# project : Unified_Yeast_GEMs_Database
import re

# 定义一个映射字典，将原有的染色体 ID 替换为新的格式
chromosome_mapping = {
    'chrI': 'chromosome1',
    'chrII': 'chromosome2',
    'chrIII': 'chromosome3',
    'chrIV': 'chromosome4',
    'chrV': 'chromosome5',
    'chrVI': 'chromosome6',
    'chrVII': 'chromosome7',
    'chrVIII': 'chromosome8',
    'chrIX': 'chromosome9',
    'chrX': 'chromosome10',
    'chrXI': 'chromosome11',
    'chrXII': 'chromosome12',
    'chrXIII': 'chromosome13',
    'chrXIV': 'chromosome14',
    'chrXV': 'chromosome15',
    'chrXVI': 'chromosome16',
}

# 读取原始 GFF 文件并进行替换
def replace_chromosome_ids(input_gff, output_gff):
    with open(input_gff, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    with open(output_gff, 'w', encoding='utf-8') as outfile:
        for line in lines:
            # 通过正则表达式查找并替换染色体 ID
            for old_id, new_id in chromosome_mapping.items():
                line = re.sub(r'\b' + re.escape(old_id) + r'\b', new_id, line)
            outfile.write(line)

# 使用示例
input_gff_file = 'data/genome/s288c_R64-5-1.gff'
output_gff_file = 'data/genome/s288c_R64-5-1.re.gff'
replace_chromosome_ids(input_gff_file, output_gff_file)