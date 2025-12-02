#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 13:59:15 2025

@author: allen
"""

import pandas as pd
from sklearn.model_selection import train_test_split

dataDf = pd.read_csv('data/protein_seqs_cleaned_6types_reps.csv')

trainDf, reserveDf  = train_test_split( dataDf , test_size=0.15 )


trainDf.to_csv( 'protein_seqs_cleaned_6types_reps_train.csv' , index=False )
reserveDf.to_csv( 'protein_seqs_cleaned_6types_reps_reserve.csv' , index=False )
