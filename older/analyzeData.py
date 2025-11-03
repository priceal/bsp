#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
'''
###############################################################################
###############################################################################
###############################################################################
'''

dataFile = 'F5-8-17_gold.csv'
dataDir = '.'

siteVocab = ' ACGTNUWSMKRYBDHV'
aaVocab = " ARNDCEQGHILKMFPSTWYV"
  
sequenceCrop = 500
siteCrop = 15
###############################################
data=pd.read_csv(os.path.join(dataDir,dataFile))
siteLength =  [ len(s) for s in data['site'] ] 
seqLength = [ len(s) for s in data['sequence'] ] 

data['site length']=pd.Series(siteLength)
data['seq length']=pd.Series(seqLength)

print('\n',data.describe())

# create data frame with sequence lengths, print stats
seqCounts = [0]*len(aaVocab)
siteCounts = [0]*len(siteVocab)
for idx in data.index:
  
    # count and increment character use
    for i,c in enumerate(aaVocab):
        seqCounts[i] += data.loc[idx]['sequence'].count(c)
    for i,c in enumerate(siteVocab):
        siteCounts[i] += data.loc[idx]['site'].count(c)

data[['site length','seq length']].hist()

# now print characters and plot use
print( '\ncharacter set:', list(aaVocab) )
plt.figure(2)
plt.title('sequence character use')
plt.bar(range(len(seqCounts)),seqCounts,tick_label=list(aaVocab)) 

# now print characters and plot use
print( '\ncharacter set:', list(siteVocab) )
plt.figure(3)
plt.title('site character use')
plt.bar(range(len(siteCounts)),siteCounts,tick_label=list(siteVocab)) 

