#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clean site file
"""
import os
import pandas as pd
'''
###############################################################################
###############################################################################
###############################################################################
'''
# data file name and directory
dataDict =    { 'file': 'F29.txt',
                'cols': (0,1),
                'sep' : ',' 
               }
dataDir = '/home/allen/projects/DATA/bsp' 

# name for new file created
saveName = 'testfile.csv'    # use None if not saving

# options
lowerCase = True    # start with False
nameStripChars = '+'         # start with ''
siteStripChars = '|()-1234567890/'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read, add length column and print stats
print("now reading",dataDict['file'])
data=pd.read_csv(os.path.join(dataDir,dataDict['file']),
                 usecols=dataDict['cols'],
                 sep=dataDict['sep'],
                 names=['RE','site'])

if lowerCase:
    data[ 'RE' ] = data[ 'RE' ].apply(str.lower)
    data[ 'site' ] = data[ 'site' ].apply(str.lower)
   
for c in nameStripChars:
    data[ 'RE' ] = data[ 'RE' ].apply( str.replace, args=(c,'') )

for c in siteStripChars:
    data[ 'site' ] = data[ 'site' ].apply( str.replace, args=(c,'') )

print(data.describe())
if saveName:
    data.to_csv(saveName, index=False)
    
    
    
    
    