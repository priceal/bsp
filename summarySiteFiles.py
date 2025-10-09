#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read in a list of site data files and create a list of all unique RE names
that occur at least once in all the files. Data files should all be in same
directory. site data file format should be two columns separated by a 
delimiter:

    RE Name (1)     site sequence (2)

will create a data frame and save if requested that contains :

1) the unique RE name 
2) the site sequences 
3) list of files contained in
 
for now, will remove cut site symbol from strings, so character list is 
assumed to be 

    'ACGTNBDHKMRSUVWY'

Will need to remove non-recognized characters and convert to upper case

note: should check for consistency. If not, what? perhaps save in a separate 
dateframe that records all sequences and file names for inspection

finally, will print out:
    number of total entries
    number of unique entries
    number of unique entries with inconsistent site sequences
    site length stats
    histogram of site lengths
    histogram of character use in site strings
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
'''
###############################################################################
###############################################################################
###############################################################################
'''

# data file name and directory
fileList = [ 'gcg.csv', 'itype2.csv', 'strider.csv'
                ]

dataDir = '/home/allen/projects/bsp/data' # '.'


###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read, add length column and print stats
dataList = []
for i,file in enumerate(fileList):
    print("now reading",file)
    dataList.append(pd.read_csv(os.path.join(dataDir,file),
                     usecols=(0,1),
                     sep=',',
                     names=['RE','site'])
                    )
    siteLengths =[ len(s) for s in dataList[i]['site'] ]
    dataList[i]['length']=pd.Series(siteLengths)
    
for fn,frame in zip(fileList,dataList):
    print('\n'+fn)
    print(frame.describe())
    reSet =  set(frame['RE'])
    print( 'unique REs:', )
    
    uniqNames = uniqNames.union(  )

# create set of unique names


