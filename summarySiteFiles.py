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
2) the site sequence (if all in agreement)
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
siteFileList = ['proto.csv',
                ]
dataDir = '/home/allen/projects/DATA/bsp' # '.'

# character separating columns---what if not same in all files?
delimitChar = ' '   #   ','  or '\t'. what if

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read, add length column and print stats
for fileName in siteFileList:
    print("now reading",fileName)
    data=pd.read_csv(os.path.join(dataDir,fileName),
                     delimiter=delimitChar,
                     names=['RE','site'])
    siteLengths =[ len(s) for s in data['site'] ]
    data['length']=pd.Series(siteLengths)
    print(data.describe())

# create site length histogram
data[['length']].hist(bins=20)

# now create character use list, print
sitesConcat = ''.join(data.site)
charSet = set( sitesConcat ).difference({'A','C','G','T','N'})
charList = list(charSet)
charList.sort()
charList = ['A','C','G','T','N'] + charList
print( '\ncharacter set:', charList )

# count character use and plot
charCounts = []
for c in charList:
    charCounts.append( sitesConcat.count(c) )
plt.figure(2)
plt.title('character use')
plt.bar(range(len(charCounts)),charCounts,tick_label=charList) 
