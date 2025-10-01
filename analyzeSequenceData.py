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
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
'''
###############################################################################
###############################################################################
###############################################################################
'''

# data file name and directory
sequenceFile = 'All_Type_II_restriction_enzyme_genes_Protein_head100000.txt'
dataDir = '/home/allen/projects/DATA/bsp'

printNames = False

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),'fasta-pearson')

# create data frame with sequence lengths, print stats
reNames = []
#sequences = ''
seqLengths = []
charCounts = np.zeros(20)
vocab = 
for c in charList:
    charCounts.append( sequences.count(c) )
for rec in record:
    if printNames:
        print(rec.name,end=' ')
    reNames.append(rec.name)
#    sequences = sequences + rec.seq
    seqLengths.append( len(rec.seq) )
    for c in charList:
        charCounts[.append( sequences.count(c) )
dataDict = { 'RE': reNames, 'length': seqLengths }
dataDf = pd.DataFrame( dataDict )
print(dataDf.describe())

# create site length histogram
dataDf[['length']].hist(bins=20)

# now create character use list, print
charList = list("ARNDCEQGHILKMFPSTWYV")
print( 'character set:', charList )

# count character use and plot
charCounts = []
for c in charList:
    charCounts.append( sequences.count(c) )
plt.figure(2)
plt.title('character use')
plt.bar(range(len(charCounts)),charCounts,tick_label=charList) 
