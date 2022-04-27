# -*- coding: utf-8 -*-
# Unified_Yeast_GEMs_Database
# why
# file：sample_ssGEMs.py
# 2022/4/27

'''randomly sample 10% ssGEMs to build the analysis pipline '''

import sys
sys.path.append('code')
from mainFunction import *

fileDir='result/1011_ssGEMs/'
tarDir='test/ssGEMs/'
sample_ratia=0.1

sampling_to_test(fileDir=fileDir, tarDir=tarDir, sample_ratio=sample_ratia)