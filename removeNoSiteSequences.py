#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 12:31:21 2025

@author: allen
"""
import os
from Bio import SeqIO
'''
###############################################################################
###############################################################################
###############################################################################
'''
# source data 
sourceFile = 'All_Type_II_restriction_enzyme_genes_Protein_nonP.fasta'
sourceDir = '../DATA/bsp'

# name for saved data file. 'None' to not save output
saveFileName = 'All_Type_II_restriction_enzyme_genes_Protein_nonP_site.fasta'
saveDir = '../DATA/bsp'

################################################################################

# ---- Read input FASTA ----
records = SeqIO.parse( os.path.join( sourceDir, sourceFile) , "fasta-pearson")

# ---- Remove sequences with no site data----
siteRecords = []
for rec in records:
    
    # parse out the site 
    split = rec.description.split()
    if (len(split)==3) or (len(split)>5):
        continue
    elif (len(split)==4) and (split[-1].strip()=='fragment'):
        continue
    else:
        siteRecords.append(rec)

# ---- Write output FASTA ----
if saveFileName:
    SeqIO.write(siteRecords, os.path.join(saveDir,saveFileName), "fasta")
