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
# source data file name and directory for sequence, format
sequenceFile = 'Type_II_restriction_enzymes_Gold_Standards_Protein.txt'
sequenceDir = '/home/allen/projects/DATA/bsp'
sequenceFormat = 'fasta'  # if error message, try 'fasta-pearson'

# source data for site data
siteFile = 'F5-8-17.csv'
siteDir = '.'

# name for saved data file
saveFileName = 'F5-8-17_gold.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in site file, create set of RE names
siteDf = pd.read_csv(os.path.join(siteDir,siteFile))
siteSet = set( siteDf['RE'] )

# read in sequence file
record = SeqIO.parse(os.path.join(sequenceDir,sequenceFile),sequenceFormat)

# create data frame with sequence lengths, print stats
sequences = []
names=[]
sites=[]
count=0
for rec in record:
    count+=1
    if rec.name.lower() in siteSet:
        sequences.append( str(rec.seq) )
        names.append( rec.name )
        sites.append( 
            siteDf[ siteDf['RE']==rec.name.lower() ]['site'].iloc[0] 
                    )
# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': names, 'site': sites, 'sequence': sequences} )
print(dataDf.describe() )
dataDf.to_csv(saveFileName,header=True, index=False)

