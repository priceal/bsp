#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""
import pandas as pd
import numpy as np
import torch
import bsp_utils as bsp
import os
'''
###############################################################################
###############################################################################
###############################################################################
'''

dataFile = 'goldAndProto.csv'
dataDir = '/home/allen/projects/DATA/bsp'

siteVocab = ' ^ACGTUWSMKRYBDHVN'
aaVocab = " ARNDCEQGHILKMFPSTWYV"

cropLength = 500
###############################################
data=pd.read_csv(os.path.join(dataDir,dataFile),index_col=0)
maxSiteLength =  max( [ len(s) for s in data['site'] ] )
maxSeqLength = max( [ len(s) for s in data['sequence'] ] )

siteList=[]
sequenceList=[]
for i in data.index:
    sequenceList.append(
        bsp.encode(
            data.at[i,'sequence'],
            vocab=aaVocab,
            length=cropLength)
        )
    siteList.append(
        bsp.encode(
            data.at[i,'site'],
            vocab=siteVocab,
            length=cropLength)
        )
sequenceTensor=torch.tensor(np.array(sequenceList))
siteTensor=torch.tensor(np.array(siteList))    

data['seqlen']=pd.Series([ len(s) for s in data['sequence'] ])
