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
sourceFile = 'Type_II_methyltransferase_genes_Protein_nonP_site.fasta'
sourceDir = '../DATA/bsp'

# name for saved data file. 'None' to not save output
saveFileName = 'Type_II_methyltransferase_genes_Protein_nonP_site_unique.fasta'
saveDir = '../DATA/bsp'

################################################################################

# ---- Read input FASTA ----
records = SeqIO.parse( os.path.join( sourceDir, sourceFile) , "fasta-pearson")

# ---- Remove duplicate sequences ----
unique_records = []
seen_sequences = set()

for i,record in enumerate(records):
    seq_str = str(record.seq).upper().strip()  # normalize case and spaces
    if seq_str not in seen_sequences:
        seen_sequences.add(seq_str)
        unique_records.append(record)

# ---- Write output FASTA ----
SeqIO.write(unique_records, os.path.join(saveDir,saveFileName), "fasta")
