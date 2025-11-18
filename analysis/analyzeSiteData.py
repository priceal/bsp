#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze a data file containing cognate site data.
format should be two columns separated by a delimiter
    RE Name (1)     site sequence (2)
    
output:
    number of entries
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
#siteFile = 'All_Type_II_restriction_enzyme_genes_Protein_sites.csv'

# 20251021 file is created from formats 2 (GCG), 8 (itype2) and 13 (Bionet)
# and contains 7913 entries: 980 (GCG), 3253 (itype2) and 3680 (Bionet)  
siteFile = 'train.csv' # formats

dataDir = '.' 

# character separating columns
delimitChar = ','   #   ','  or '\t'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read, add length column and print stats
data=pd.read_csv(os.path.join(dataDir,siteFile),
                   delimiter=delimitChar)
#                   names=['RE','site'])
siteLengths =[ len(s) for s in data['site'] ]
data['length']=pd.Series(siteLengths)
print('\n',data.describe())

# create site length histogram
data[['length']].hist(bins=20)

# now create character use list, print
sitesConcat = ''.join(data.site)
charList = list( set( sitesConcat ).difference({'A','C','G','T','N'}) )
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
