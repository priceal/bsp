#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""
import pandas as pd
import os
'''
###############################################################################
###############################################################################
###############################################################################
'''

dataFile = 'goldAndProto.csv'
dataDir = '/home/allen/projects/DATA/bsp'

###############################################################################
data=pd.read_csv(os.path.join(dataDir,dataFile),usecols=(1,2,3))

seqLengths =[ len(s) for s in data['sequence'] ]
siteLengths =[ len(s) for s in data['site'] ]

data['seqlen']=pd.Series(seqLengths)
data['sitelen']=pd.Series(siteLengths)

data[['sitelen','seqlen']].hist(bins=20)
print(data.describe())


