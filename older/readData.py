#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""
import pandas as pd
import os
from Bio import SeqIO
'''
###############################################################################
###############################################################################
###############################################################################
'''
siteFile = 'proto_clean.csv'
sequenceFile = 'Type_II_restriction_enzymes_Gold_Standards_Protein.txt'
saveFile = 'test.csv'
dataDir = '/home/allen/projects/DATA/bsp'

delimitChar = '\t'

###############################################################################
dataDf=pd.read_csv(os.path.join(dataDir,siteFile),
                   delimiter=delimitChar,
                   names=['RE','site'])
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),'fasta')

siteFileSet = set(dataDf['RE'])
dataDf['sequence'] = ['']*len(dataDf)

for rec in record:
    print('\nlooking at', rec.name)
    if rec.name in siteFileSet:
        print(' ...there!')
        idx=dataDf.RE[dataDf.RE == rec.name].index.tolist()[0]
        dataDf['sequence'].at[idx] = rec.seq

dataComplete = dataDf[ dataDf['sequence'] != '' ]
print(len(dataComplete))

dataComplete.to_csv(os.path.join(dataDir,saveFile))

