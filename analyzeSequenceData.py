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
sequenceFile = 'Type_II_restriction_enzymes_Gold_Standards_Protein.txt'
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
charCounts = [0]*20
step = 0
chars = "ARNDCEQGHILKMFPSTWYV"
for rec in record:
    reNames.append(rec.name.lower())
    seqLengths.append( len(rec.seq) )
  
    # count and increment character use
    for i,c in enumerate(chars):
        charCounts[i] += rec.seq.count(c)
        
    # print out RE name if needed
    if step % reportCycle == 0:
        if printNames:
            print(step,rec.name)
            
    step += 1

# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': reNames, 'length': seqLengths } )
print('\n\n',dataDf.describe())
dataDf[['length']].hist(bins=20)

# now print characters and plot use
print( '\ncharacter set:', list(chars) )
plt.figure(2)
plt.title('character use')
plt.bar(range(len(charCounts)),charCounts,tick_label=list(chars)) 
