#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse a dsequence ata file containing protein sequences and site data into
two data frames
input format should be FASTA
outputs are csv files, one a site file, the other a combined site/sequence
file
    
output:
    number of entries
    site length stats
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

# file format
fileFormat = 'fasta'  # if error message, try 'fasta-pearson'

# option to print sequence names as processed
printNames = True  
reportCycle = 500   # only print every reportCycle entries

# output file names
siteFileOutput = 'data/All_Type_II_restriction_enzyme_genes_Protein_sites.csv'
combinedFileOutput = 'All_Type_II_restriction_enzyme_genes_Protein_combined.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),fileFormat)

# create data frame with sequence lengths, print stats
reNames = []
sites=[]
sequences = []

chars = "ARNDCEQGHILKMFPSTWYV"
for i,rec in enumerate(record):
    reNames.append(rec.name)
    sequences.append( str(rec.seq) )
    split = rec.description.split()
    if len(split)==5:
        sites.append(split[1])
    elif len(split)==4 and split[-1] == 'aa':
        sites.append(split[1])
    else:
        sites.append('x')
        
    # print out RE name if needed
    if i % reportCycle == 0:
        if printNames:
            print(i,rec.name)

# create dataframe, print stats 
dataDf = pd.DataFrame( { 'RE': reNames, 
                        'site':sites,
                        'sequence': sequences} )
print('\nall sequences\n',dataDf.describe())

# create siteDf and 
siteDf = dataDf[ dataDf['site']!='x']
print('\nall sequences with sites\n',siteDf.describe())

if siteFileOutput:
    siteDf[ ['RE','site' ] ].to_csv(siteFileOutput,index=False)
if combinedFileOutput:
    siteDf.to_csv(combinedFileOutput,index=False)






