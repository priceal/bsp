#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
c
  
"""
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
'''
###############################################################################
###############################################################################
###############################################################################
'''
# source data 
sourceFile = 'All_Type_II_restriction_enzyme_genes_Protein_combined.csv'
sourceDir = 'data'

# name for saved data file. 'None' to not save output
saveFileName = 'data/All_Type_II_restriction_enzyme_genes_Protein_combined.fasta'

################################################################################

sourceDf = pd.read_csv( os.path.join( sourceDir, sourceFile) )

# ---- Convert to SeqRecord objects ----
records = []
for _, row in sourceDf.iterrows():
    seq = Seq(str(row["sequence"]).strip())
    records.append(
                SeqRecord(seq, id=str(row["RE"]), description="")
        
        )

# ---- Write to FASTA ----
SeqIO.write(records, saveFileName, "fasta")




