#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze a data file containing protein sequences and possible sites
format should be FASTA
    
output for (1) all sequences, and (2) all sequences with site data:
    statistics on lengths
    histogram of lengths (only for (2))
    
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
#sequenceFile='All_REBASE_Gold_Standards_Protein.txt'
sequenceFile = 'All_Type_II_restriction_enzyme_genes_Protein.txt'
dataDir = '/home/allen/projects/DATA/bsp'

# input file format
fileFormat = 'fasta'  #if error message, try 'fasta-pearson' or 'fasta'

# option to print sequence names as processed
printNames = True  
reportCycle = 1000   # only print every reportCycle entries

# output file names --- 'None' if not saving reformated data
siteFileOutput = 'data/All_Type_II_restriction_enzyme_genes_Protein_sites.csv'
combinedFileOutput = 'data/All_Type_II_restriction_enzyme_genes_Protein_combined.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),fileFormat)

# create data frame with sequence/site data, print stats
reNames = []
sites = []
sequences = []
seqLengths = []
for i,rec in enumerate(record):
    reNames.append(rec.name.lower())
    sequences.append( str(rec.seq) )
    seqLengths.append( len(rec.seq) )
    
    # parse out the site 
    split = rec.description.split()
    if len(split)==5:
        sites.append(split[1])
    elif len(split)==4 and split[-1] == 'aa':
        sites.append(split[1])
    else:
        sites.append('x') # missing site
            
    # print out RE name if needed
    if i % reportCycle == 0:
        if printNames:
            print(i,rec.name)
            
# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': reNames, 
                        'site': sites,
                        'sequence': sequences,
                        'length': seqLengths } )

print('\nall sequences\n',dataDf.describe())
print( 'unique sequences:', len(set(dataDf.sequence)) )
      
'''
plt.figure(1)
dataDf['length'].hist(bins=20)
'''
print('\nsequences with sites\n',dataDf[ dataDf.site!='x' ].describe())
dataDf[ dataDf.site!='x' ]['length'].hist(bins=20)
print( 'unique sequences:', len(set(dataDf[ dataDf.site!='x' ].sequence)) )
print( 'unique sites:', len(set(dataDf[ dataDf.site!='x' ].site)) )

if siteFileOutput:
    dataDf[ dataDf.site!='x' ][ ['RE','site' ] ].to_csv(siteFileOutput,index=False)
if combinedFileOutput:
    dataDf[ dataDf.site!='x' ][ ['RE','site','sequence' ] ].to_csv(combinedFileOutput,index=False)





