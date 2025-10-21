#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
combine a site data file and a sequence data file.
site data file must be in csv format with one column 'RE' and the
second one 'site'
sequence data must be in fasta format with name of enzyme leading in title 
line. for example:

>AatII   GACGTC  345 aa

Code will find intersection of RE names and create combined file with only
those

also prints out stats on resulting dataframe
  
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
sequenceFile = 'All_Type_II_restriction_enzyme_genes_Protein.txt'
sequenceDir = '/home/allen/projects/DATA/bsp'
sequenceFormat = 'fasta-pearson'  # if error message, try 'fasta-pearson'

# source data for site data
siteFile = 'sites_combined_20251021.csv'
siteDir = 'data'

# name for saved data file. 'None' to not save output
saveFileName = None # 'F5-8-17_gold.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in site file, create set of RE names
siteDf = pd.read_csv(os.path.join(siteDir,siteFile))
siteSet = set( siteDf['RE'] )

# read in sequence file
record = SeqIO.parse(os.path.join(sequenceDir,sequenceFile),sequenceFormat)

# create data frame with sequence lengths, print stats
seqNames=[]
sequences = []
names=[]
sites=[]
count = 0
for rec in record:
    count += 1
    seqNames.append( rec.name.lower() )
    if rec.name.lower() in siteSet:
        print( rec.name , end=' ')
        sequences.append( str(rec.seq) )
        names.append( rec.name )
        sites.append( 
            siteDf[ siteDf['RE']==rec.name.lower() ]['site'].iloc[0] 
                    )
# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': names, 'site': sites, 'sequence': sequences} )
print('\nprocessed',count,'sequences' )
print(dataDf.describe() )
if saveFileName:
    dataDf.to_csv(saveFileName,header=True, index=False)

seqSet=set(seqNames)
remain=siteSet.difference(seqSet)
for name in remain:
    for n in seqNames:
        if name in n:
            print(name,n)




