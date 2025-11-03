#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze a data file containing protein sequences and possible sites
format should be FASTA
    
output for (1) all sequences, and (2) all sequences with site data:
    statistics on lengths
    histogram of lengths
    
option for saving data in csv format as either 

    (1) site data file:   2 columns:  RE    SITE
    (2) combined data:    3 columns:  RE    SITE    SEQUENCE
    
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
# data file name and directory for input sequence data
sequenceFile='All_REBASE_Gold_Standards_Protein.txt'
#sequenceFile = 'All_Type_II_restriction_enzyme_genes_Protein.txt'
dataDir = '/home/allen/projects/DATA/bsp'

# input file format
fileFormat = 'fasta-pearson'  #if error message, try 'fasta-pearson' or 'fasta'

# option to print sequence names as processed
printNames = True  
reportCycle = 500   # only print every reportCycle entries

# output file names --- 'None' if not saving reformated data
siteFileOutput = None
combinedFileOutput = None

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),fileFormat)

# create data frame with sequence/site data, print stats
reNames = []
sites=[]
sequences = []
seqLengths = []
seqCharCounts = [0]*20
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
