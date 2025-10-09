#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cleans site data file and saves it in two column CSV format.
specify which columns to input and separator (delimiter) using
key words 'cols' and 'sep'

options include 
 (1) making RE names or site sequences upper or lower case.
 (2) removing or stripping characters from RE names or site sequences
 (3) the name for the saved file
 
for initial inspection of data file, set all options to None

"""
import os
import pandas as pd
'''
###############################################################################
###############################################################################
###############################################################################
'''
# particulars of input data (required)
dataDict =    { 'file': 'itype2.txt',
                'cols': (0,2),
                'sep' : '\t' 
               }
dataDir = '/home/allen/projects/DATA/bsp' 

# option: chose str.upper or str.lower for case
# use None to not alter case
nameCase = str.lower
siteCase = str.upper

# option: chars to remove --- set to '' or None to not strip any
nameStripChars = ''
siteStripChars = "^/1234567890()-,"
'''
nameStripChars = '+'
siteStripChars = '|()-1234567890/'
'''

# option: name for new file created --- use None if not saving
saveName = 'data/itype2.csv'    
#saveName = None

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

print("now reading",dataDict['file'])
data=pd.read_csv(os.path.join(dataDir,dataDict['file']),
                 usecols=dataDict['cols'],
                 sep=dataDict['sep'],
                 names=['RE','site'],
                 engine='python')

if nameCase:
    data[ 'RE' ] = data[ 'RE' ].apply(nameCase)

if siteCase:
    data[ 'site' ] = data[ 'site' ].apply(siteCase)

if nameStripChars:
    for c in nameStripChars:
        data[ 'RE' ] = data[ 'RE' ].apply( str.replace, args=(c,'') )

if siteStripChars:
    for c in siteStripChars:
        data[ 'site' ] = data[ 'site' ].apply( str.replace, args=(c,'') )

if saveName:
    data.to_csv(saveName, index=False, header=False)
    
print(data.describe())    
    
    
    