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
sourceFile = 'Type_II_methyltransferase_genes_Protein.txt'
sourceDir = '../DATA/bsp'

# name for saved data file. 'None' to not save output
saveFileName = 'Type_II_methyltransferase_genes_Protein_nonP.fasta'
saveDir = '../DATA/bsp'

################################################################################

# ---- Read input FASTA ----
records = SeqIO.parse( os.path.join( sourceDir, sourceFile) , "fasta-pearson")

# ---- Remove putative sequences ----
nonPRecords = []
nonP = []
P = []

for rec in records:
    if not str(rec.name).endswith('P'):
        nonPRecords.append(rec)
        nonP.append(rec.name )
    else:
        P.append( rec.name)

# ---- Write output FASTA ----
SeqIO.write(nonPRecords, os.path.join(saveDir,saveFileName), "fasta")
