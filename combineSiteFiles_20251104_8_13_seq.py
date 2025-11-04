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

    2) format 2 has flanking Ns and lists reverse complement for some
    
    3) format 13 & 15 has reverse complements and leading/trailing Ns
    
    6) true inconsistencies seem limited to lpnpi,hpyum037x
 
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

fileList = [  'sites_combined_20251103_8_13.csv',
            'All_Type_II_restriction_enzyme_genes_Protein_sites.csv'
            ]


excludeList = ['lpnpi','hpyum037x']  # actual inconsistencies
#excludeList = []
stripChars = ' N'
dataDir = 'data'

saveName = None # 'sites_combined_20251103_8_13.csv' 

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
    # strip leading/trailing whitespaces from RE names!!!
    dataList[i]['RE'] = dataList[i]['RE'].map( lambda x: x.strip() )

    # strip leading/trailing whitespaces from sites
    dataList[i]['site'] = dataList[i]['site'].map( lambda x: x.strip() )

    # now, if there are special chars, strip them from site
    for c in stripChars:    
        dataList[i]['site'] = dataList[i]['site'].map( lambda x: x.strip(c) )
        
    # add column for source    
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
#dataAll.drop_duplicates(subset=['RE'],keep='first',inplace=True)

# remove exclude list
if excludeList:
    for name in excludeList:
        dataAll = dataAll[ dataAll['RE'] != name ]

# group by RE name and count how many have multi entries
dataGroup = dataAll.groupby(by='RE')
groupCount = dataGroup.apply(len)
groupsMulti = groupCount[ groupCount >1 ]
groupCount.hist(bins=20)

# print out REs with inconsistent sites
rcList = []
otherList = []
for k in groupsMulti.index:
    names = list( dataGroup.get_group(k)['site'] )
    if len(names) == 1:
        print('error: not a multiple!') 
        continue
    elif len(names) == 2:
        if names[0] == reverse_complement(names[1]):
            #print( k,'2 entries --- reverse complement')
            rcList.append(k)
        else:
            print( k,'2 entries --- other')
            otherList.append(k)
    elif len(names) > 2:
        print(k, len(names), 'entries ??')
        otherList.append(k)
        
print('total entries: ', len(dataAll) )
print( 'RE names with multiple entries:', len(groupsMulti))

print('sources:')
for s in set( dataAll.source ):
    print( s, ':', len( dataAll[ dataAll.source==s ] ) )
        


'''
for k in groupsMulti.index:
    names = set( dataGroup.get_group(k)['site'] )
    if len(names)>1:
        print(dataGroup.get_group(k))
'''
if saveName:
    dataAll.to_csv(saveName,header=True,index=False)

