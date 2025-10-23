#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze a data file containing protein sequences
format should be FASTA
    
output:
    number of entries
    site length stats
    histogram of site lengths
    histogram of character use in site strings
"""
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
'''
###############################################################################
###############################################################################
###############################################################################
'''
# data file name and directory
sequenceFile = 'All_Type_II_restriction_enzyme_genes_Protein.txt'
dataDir = '/home/allen/projects/DATA/bsp'

# option to print sequence names as processed
printNames = True  
reportCycle = 1000   # only print every reportCycle entries

# file format
fileFormat = 'fasta-pearson'  # if error message, try 'fasta-pearson'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),fileFormat)

# create data frame with sequence lengths, print stats
reNames = []
seqLengths = []
sites=[]
siteLengths =[]
charCounts = [0]*20
step = 0
chars = "ARNDCEQGHILKMFPSTWYV"
for rec in record:
    reNames.append(rec.name.lower())
    seqLengths.append( len(rec.seq) )
    split = rec.description.split()
    if len(split)==5:
        sites.append(split[1])
        siteLengths.append( len(split[1]) )
    elif len(split)==4 and split[-1] == 'aa':
        sites.append(split[1])
        siteLengths.append( len(split[1]) )
    else:
        sites.append('x')
        siteLengths.append( 'x' )
        print('      ',rec.description)
    # count and increment character use
    for i,c in enumerate(chars):
        charCounts[i] += rec.seq.count(c)
        
    # print out RE name if needed
    if step % reportCycle == 0:
        if printNames:
            print(step,rec.name)
            
    step += 1

# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': reNames, 'seqLength': seqLengths, 
                        'site':sites, 'siteLength': siteLengths  } )
print('\nall sequences\n',dataDf.describe())

# create siteDf and 
siteDf = dataDf[ dataDf['site']!='x']
print('\nall sequences with sites\n',siteDf.describe())








