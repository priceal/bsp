#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to train the classifier model.  

Script will load a training set and perform a given number of 
epochs of stochastic gradient descent. During optimization it tracks the 
training loss. At the end of optimization it will plot 
these losses versus iteration of training.

Filtering of input sequences can be done with sequence 
limits parameters. Accepted data entries can also be 
cropped using the crop parameters. 
    
"""

# import libraries -----------------------------------------------------
import os
import numpy as np
import torch
from torch.utils.data import DataLoader
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import bsp_utils as bsp
from ClassifierModel_20251202 import classModel

# set learning parameters --------------------------------------------------
numBatches = 1 # if non-zero, ignore batchSize and set to N/numBatches
batchSize = 0 # only use if numBatches = 0
numberEpochs = 500
learningRate = 0.1
reportCycle = 10
refine = False     # creates new model if False
                
# files to load and optional file directory -----------------------------
# can leave undefined '' or '.'
inputTrain = 'data/train.csv'

# data parameters for screening sequence length and cropping-----------------
# format ( min, max, crop )
sequenceLimits = ( 1, 2000, 502 )  

###########################################################################
###########################################################################
###########################################################################

# symbol usage, uses space as padding -------------------------------------
aaVocab = " ARNDCEQGHILKMFPSTWYV"

# load data -------------------------------------------
xTrain, yTrain = bsp.classDataReader(inputTrain, 
                                seq = sequenceLimits
                              )
yTrain.swapaxes_(1, 2)   # put in correct order for CNN
dataTrain = bsp.seqDataset(xTrain, yTrain ) # needed for batches

# print data/batch stats ------------------------------------
print("DATA SET")
print("{:<20} {:<15} {:<15}".format('DATA', 'ENTRIES', 'LENGTH'))
rows = ['training data']
ds = [xTrain]
for r, d in zip(rows, ds):
    a, b = d.shape
    print(f"{r:<20} {a:<15} {b:<15}")

###########################################################################

# determine batch size and number of batches -----------------------    
if numBatches > 0:
    batchSize = int(len(xTrain)/numBatches)
else:
    numBatches = int(len(xTrain)/batchSize)

print('number of batches:', numBatches)
print('size of batches:', batchSize)
dataloader = DataLoader(dataTrain, batch_size=batchSize, shuffle=True)

###########################################################################
# create model ----------------------------------------------------
if not refine:     # if refining pre-existing, don't create new model
    model = classModel()
print('\nMODEL ')
print("{0:20} {1:20}".format("MODULES", "PARAMETERS"))
total_params = 0
for name, parameter in model.named_parameters():
    if not parameter.requires_grad:
        continue
    params = parameter.numel()
    print("{0:20} {1:<20}".format(name, params))
    total_params += params
print("{0:20} {1:<20}".format("TOTAL", total_params))

###########################################################################
# run cycles of optimization ----------------------------------------

optimizer = torch.optim.SGD(model.parameters(), lr=learningRate)
'''<
define an approprirate loss function here for a classification model with 
6 classes
>'''
print('\nOPTIMIZATION')
print('{:10} {:10} {:10} {:10}'.format('epoch','batch','loss-train') )
stepCount = 0
trainLosses = []
for i in range(numberEpochs):
    for j, batch in enumerate(dataloader):
        
        xx, yy = batch[0], batch[1]
        prediction = model(xx)
'''<
        loss = ...use loss function here to calculate correct loss 
>'''      
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # print out if report cycle done
        if stepCount % reportCycle == 0:
            
            print(f"{i:<10} {j:<10} {loss.item():<10.5}")
            trainLosses.append(loss.item())        
        stepCount += 1
        
plt.figure(1)
plt.plot(trainLosses, '.k', label='Training Loss')
plt.legend()
plt.xlabel('training iteration')
plt.ylabel('loss')
plt.show()

###########################################################################
