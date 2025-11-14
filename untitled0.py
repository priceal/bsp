#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 13:59:15 2025

@author: allen
"""

import pandas as pd
data0 = pd.read_csv('notebook/data0.csv')
data1 = pd.read_csv('notebook/data1.csv')
data2 = pd.read_csv('notebook/data2.csv')


res = pd.read_csv('data/sites_combined_20251104_8_13.csv',names=['RE','site','source'])


resneb = set(res.RE.map( lambda x: x.lower()))
res0 = set(data0.RE.map( lambda x: x.lower()))
res1 = set(data1.RE.map( lambda x: x.lower()))
res2 = set(data2.RE.map( lambda x: x.lower()))

temp=None
for i, row in res.iterrows():
    temp=data0[ data0.RE.map( lambda x: (x.lower().count(row.RE)>0) and (x.lower()!=row.RE) ) ]
    if len(temp)>0:
        print('\n'+row.RE+' in NEB site data')
        print( temp )