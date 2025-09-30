#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""
import pandas as pd
import numpy as np
import torch
import bsp_utils as bsp

dataFile = 'data/goldAndProto.csv'

siteVocab = ' ^ACGTUWSMKRYBDHVN'
aaVocab = " ARNDCEQGHILKMFPSTWYV"

###############################################
data=pd.read_csv(dataFile)
maxSiteLength =  max( [ len(s) for s in data['site'] ] )
maxSeqLength = max( [ len(s) for s in data['sequence'] ] )

siteList=[]
sequenceList=[]
for i in data.index:
    sequenceList.append(
        bsp.encode(
            data.at[i,'sequence'],
            vocab=aaVocab,
            length=maxSeqLength)
        )
    siteList.append(
        bsp.encode(
            data.at[i,'site'],
            vocab=siteVocab,
            length=maxSiteLength)
        )
sequenceTensor=torch.tensor(np.array(sequenceList))
siteTensor=torch.tensor(np.array(siteList))    
