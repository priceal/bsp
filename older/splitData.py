#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""



import pandas as pd
from sklearn.model_selection import train_test_split
import os

'''
###############################################################################
###############################################################################
###############################################################################
'''

# input data
inputFile = 'F5-8-17_gold.csv'
inputDir = 'data'

# filter bining site lengths
siteFilter = ( 1, 20 )      # format: (min,max)

# filter sequence lengths
seqFilter = ( 100, 2000 )      # format: (min,max)

# define splits - leave 3rd value 0 for train/test only
split = ( 0.8, 0.2, 0 )   # format: (train,validate,test)

# output file names:
outputFiles = ( 'train.csv', 'test.csv' , None )


####################################################################
dataDf = pd.read_csv(os.path.join(inputDir,inputFile))

# filter sequence and site lengths
dataDf = dataDf[ dataDf['sequence'].str.len() >= seqFilter[0] ]
dataDf = dataDf[ dataDf['sequence'].str.len() <= seqFilter[1] ]
dataDf = dataDf[ dataDf['site'].str.len() >= siteFilter[0]]
dataDf = dataDf[ dataDf['site'].str.len() <= siteFilter[1] ]

if split[2] == 0:
    a, b = train_test_split(dataDf, test_size=split[1])
    a.to_csv(outputFiles[0], index=False)
    b.to_csv(outputFiles[1], index=False)

else:
    a, temp = train_test_split(dataDf, test_size=split[1]+split[2])
    b, c = train_test_split( temp, test_size=(split[2]/(split[1]+split[2])))


        
        