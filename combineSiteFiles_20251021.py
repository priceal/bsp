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
    
    
notes: 
    
    1) format 15 has multiple entries for promiscuous REs, eschewing use of
    ambiguity code

    2) format 2 has flanking Ns and lists reverse complement for some
    
    3) format 13 & 15 has reverse complements and leading/trailing Ns
    
    4) format 33 has leading spaces

    5) formats 17 & 23 appear to have some entries reverse complement of 
    each other
    
    6) true inconsistencies seem limited to lpnpi,hpyum037x
    
    7) formats 17 and 23 appear to have reverse complements of each
    other for some entries
"""
import os
import pandas as pd
from Bio.Seq import reverse_complement
'''
###############################################################################
###############################################################################
###############################################################################
'''

# data file names and directory

fileList = [   'Format2_C.csv',
             'Format5_C.csv',
             'Format8_C.csv',
            'Format13_C.csv',
            'Format17_C.csv',
            'Format23_C.csv',
            'Format29_C.csv',
             'Format33_C.csv',
             'Format37_C.csv',
]

excludeList = ['lpnpi','hpyum037x']  # actual inconsistencies
stripChars = ' N'
dataDir = 'kylie'

saveName = 'sites_combined_20251021.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read, add length column and print stats
dataList=[]
for i,file in enumerate(fileList):
    print("now reading",file)
    dataList.append(
                     pd.read_csv(os.path.join(dataDir,file),
                     usecols=(0,1),
                     sep=',',
                     names=['RE','site'])
                    )
    for c in stripChars:    # strip prelim/terminal Ns
        dataList[i]['site'] = dataList[i]['site'].map( lambda x: x.strip(c) )
    
    # add columns for site length and source file
    dataList[i]['length']=pd.Series(
                                    [ len(s) for s in dataList[i]['site'] ]
                                    )
    dataList[i]['source']=pd.Series( [file]*len(dataList[i]) )

# print stats on all data sets    
for fn,frame in zip(fileList,dataList):
    print('\n'+fn)
    print(frame.describe())

# now concatenate
dataAll = pd.concat( dataList, ignore_index=True )

# remove duplicates
# first pass use subset= ['RE', 'site'] to capture multiple site entries
# for same RE. After eliminating all inconsistent sites, you can drop
# duplicates for subset = ['RE'], since duplicates are either same or
# reverse compliment
#dataAll.drop_duplicates(subset=['RE','site'],keep='first',inplace=True)
dataAll.drop_duplicates(subset=['RE'],keep='first',inplace=True)

# remove exclude list
if excludeList:
    for name in excludeList:
        dataAll = dataAll[ dataAll['RE'] != name ]

# group by RE name and count how many have multi entries
dataGroup = dataAll.groupby(by='RE')
groupCount = dataGroup.apply(len)
groupsMulti = groupCount[ groupCount >1 ]

multiGroups = dataGroup

# print out REs with inconsistent sites
rcList = []
otherList = []
for k in groupsMulti.index:
    names = list( dataGroup.get_group(k)['site'] )
    if len(names) == 1:
        continue
    elif len(names) == 2:
        print(k,'2 entries --- ', end=' ')
        if names[0] == reverse_complement(names[1]):
            print( 'reverse complement')
            rcList.append(k)
        else:
            print( 'other')
            otherList.append(k)
    elif len(names) > 2:
        print(k, len(names),'entries ??')
        otherList.append(k)
        


'''
for k in groupsMulti.index:
    names = set( dataGroup.get_group(k)['site'] )
    if len(names)>1:
        print(dataGroup.get_group(k))
'''
if saveName:
    dataAll.to_csv(saveName,header=True,index=False)

