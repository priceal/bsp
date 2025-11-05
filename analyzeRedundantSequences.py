#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 12:31:21 2025

@author: allen
"""
import os
from Bio import SeqIO
import pandas as pd

'''
###############################################################################
###############################################################################
###############################################################################
'''
# source data 
sourceFile = 'All_Type_II_restriction_enzyme_genes_Protein_nonP_site.fasta'
sourceDir = '../DATA/bsp'

# name for saved data file. 'None' to not save output
saveFileName = None #'Type_II_methyltransferase_genes_Protein_nonP_site_unique.fasta'
saveDir = '../DATA/bsp'

################################################################################
        
# ---- Read input FASTA ----
records = SeqIO.parse( os.path.join( sourceDir, sourceFile) , "fasta-pearson")

# create data frame with RE name, site and sequence
sequences = []
names=[]
sites=[]
for rec in records:
    names.append( rec.name.lower().strip() )
    sequences.append( str(rec.seq).upper().strip() )
    
    # parse out the site 
    split = rec.description.split()
    if len(split)==5:
        sites.append(split[1].upper().strip())
    elif len(split)==4 and split[-1] == 'aa':
        sites.append(split[1].upper().strip())
    else:
        print('missing site', rec.description)
        sites.append('x') # missing site


# create dataframe, print stats and plot length histogram
dataDf = pd.DataFrame( { 'RE': names, 'site': sites, 'sequence': sequences} )


group = dataDf.groupby('sequence')
groupCount = group.apply(len)
groupsMulti = groupCount[ groupCount >1 ]
groupCount.hist(bins=20)
whoops=[]
for s in groupsMulti.index:
    temp = dataDf[ dataDf.sequence == s]
    testList = list(temp.site)
    if len(set(testList)) != 1:
        print('    whoops!')
        print( temp )
        whoops.append(temp)
        
 



# ---- Write output FASTA ----
if saveFileName:
    SeqIO.write(unique_records, os.path.join(saveDir,saveFileName), "fasta")
