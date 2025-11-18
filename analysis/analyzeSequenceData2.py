#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze a data file containing protein sequences and possible sites
format should be FASTA, with headers to each sequence of following format:
    

>REBASE:M.Aac9709I	EnzType:Type II methyltransferase	RecSeq:GATC	...
GenBank:UFSG01000001	SeqLength:284	Locus:NCTC9709_01044	...
ProteinId:SSY84327.1

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
sequenceFile = 'protein_seqs.txt'
dataDir = '/home/allen/projects/DATA/bsp'

# input file format
fileFormat = 'fasta-pearson'  #if error message, try 'fasta-pearson' or 'fasta'

# option to print sequence names as processed
printNames = True  
reportCycle = 5000   # only print every reportCycle entries

# set to True to filter out putative sequences
filterPut = True

# output file names --- 'None' if not saving reformated data
siteFileOutput = None # 'data/All_Type_II_restriction_enzyme_genes_Protein_sites.csv'
combinedFileOutput =  'data/protein_seqs_cleaned_type.csv'

###############################################################################
################ DOT NOT CHANGE ANYTHING UNDER THIS SEPARATOR #################
###############################################################################

# read in sequence file
record = SeqIO.parse(os.path.join(dataDir,sequenceFile),fileFormat)

# create data frame with sequence/site data, print stats
names = []
sites = []
sequences = []
lengths = []
enzymeTypes = []
for i, rec in enumerate(record):
    headerDict = dict(
                    [s.split(':',maxsplit=1) for s in rec.description.split('\t')]
                    )
    sequences.append( str(rec.seq).strip(' <>') ) # strip to be safe
    lengths.append( len(str(rec.seq).strip(' <>') ) )
    
    if 'REBASE' in headerDict.keys():
        names.append( headerDict['REBASE'].strip() ) # strip to be safe
    else:
        names.append( '?' ) 
        
    if 'RecSeq' in headerDict.keys():
        sites.append( headerDict['RecSeq'].strip(' N') )
    else:
        sites.append('x')
        
    if 'EnzType' in headerDict.keys():
        enzymeTypes.append( headerDict['EnzType'].strip() )
    else:
        enzymeTypes.append('no type')
        
    # print out RE name if needed
    if i % reportCycle == 0:
        if printNames:
            print(i,rec.name)
    
# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': names, 
                        'type': enzymeTypes,
                        'site':sites,
                        'sequence': sequences,
                        'length': lengths } )

print('\nall sequences\n',dataDf.describe())
print( 'unique sequences:', len(set(dataDf.sequence)) )
      
# filter out putative sequences if requested
if filterPut:
    dataDf = dataDf[ dataDf['RE'].map(
        lambda x: not x.split()[0].endswith('P') 
        ) ] # some have more than one name 
    
    print('\nnon-putative sequences\n',dataDf.describe())
    print( 'unique non-putative sequences:', len(set(dataDf.sequence)) )
    
print('\nsequences with sites\n',dataDf[ dataDf.site!='x' ].describe())
dataDf[ dataDf.site!='x' ]['length'].hist(bins=20)
print( 'unique sequences:', len(set(dataDf[ dataDf.site!='x' ].sequence)) )
print( 'unique sites:', len(set(dataDf[ dataDf.site!='x' ].site)) )

if siteFileOutput:
    dataDf[ dataDf.site!='x' ][ ['RE','site' ] ].to_csv(siteFileOutput,index=False, header=False)
if combinedFileOutput:
    dataDf[ dataDf.site!='x' ][ ['RE','type','site','sequence' ] ].to_csv(combinedFileOutput,index=False)





