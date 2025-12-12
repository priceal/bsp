#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to train the classifier model.  

Script will load a training set and perform a given number of 
epochs of stochastic gradient descent. During optimization it tracks the 
training loss. At the end of optimization it will plot 
these losses versus iteration of training.

Filtering of input sequences can be done with sequence limits parameters. 
Accepted data entries can also be cropped using the crop parameters. 

Created: December 2024
"""

# import libraries -----------------------------------------------------
'''  USE IF ERROR
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'  # this is dangerous, but works!
'''

import torch
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt
import bsp_utils as bsp
from ClassifierModel_20251202 import classModel

# set learning parameters --------------------------------------------------
batchSize = 64         
numberEpochs = 1
learningRate = 0.05
reportCycle = 1
refine = False  # creates new model if False
                
# files to load and optional file directory -----------------------------
inputTrain = 'protein_seqs_cleaned_6types_reps_train.csv'

# data parameters for screening sequence length and cropping-----------------
# format ( min, max, crop )
sequenceLimits = ( 1, 2000, 1708 )  # architecture 1 needs input length 1708

################################################################################
################################################################################
################################################################################

# symbol usage, uses space as padding -------------------------------------
aaVocab = " ARNDCEQGHILKMFPSTWYV"

# load data -------------------------------------------
# datareader should return encoded sequence and class index (no onhot vectors)
xTrain, yTrain = bsp.classDataReader(inputTrain, 
                                seq = sequenceLimits
                              )
dataTrain = bsp.seqDataset(xTrain, yTrain ) # needed for batches

# print data/batch stats ------------------------------------
print("DATA SET")
print("{:<20} {:<15} {:<15}".format('DATA', 'ENTRIES', 'LENGTH'))
rows = ['training data']
ds = [xTrain]
for r, d in zip(rows, ds):
    a, b = d.shape
    print(f"{r:<20} {a:<15} {b:<15}")

################################################################################
# determine batch size and number of batches -----------------------    
if len(xTrain)%batchSize == 0:
    print( f'{len(xTrain)/batchSize} batches of size {batchSize}' )
else:
    print(f'{int(len(xTrain)/batchSize)} batches of size {batchSize}' )
    print(f'1 batch of size {len(xTrain)%batchSize}' )
dataloader = DataLoader(dataTrain, batch_size=batchSize, shuffle=True)

################################################################################

# create model ----------------------------------------------------
if not refine:     # if refining pre-existing, don't create new model
    model = classModel()
    trainLosses = []        # clear records if not refining
    iterationCount = []
    stepCount = 0
print('\nMODEL ')
print("{0:20} {1:20}".format("MODULES", "PARAMETERS"))
total_params = 0
for name, parameter in model.named_parameters():
    if not parameter.requires_grad:
        continue
    params = parameter.numel()
    print("{0:20} {1:>20}".format(name, params))
    total_params += params
print("{0:20} {1:>20}".format("TOTAL", total_params))

################################################################################

# run cycles of optimization ----------------------------------------
optimizer = torch.optim.SGD(model.parameters(), lr=learningRate)

# CrossEntropyLoss for multi-class classification (6 classes)
# Note: CrossEntropyLoss expects raw logits, so softmax should be removed from
# model during training. Alternatively, use NLLLoss with log of softmax output.
# Here we use CrossEntropyLoss and assume class indices.
lossFunction = torch.nn.CrossEntropyLoss()

print('\nOPTIMIZATION')
print('{:10} {:10} {:10} {:10}'.format('iteration','epoch','batch','loss-train') )

for i in range(numberEpochs):
    for j, batch in enumerate(dataloader):
        
        xx, yy = batch[0], batch[1]
        prediction = model(xx)
        
        # assumiung target yy are class indices and prediction are raw logits
        # of dim 6
        loss = lossFunction(prediction, yy)
        
        # cycle optimizer
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # print out if report cycle done
        if stepCount % reportCycle == 0:
            
            print(f"{stepCount:<10} {i:<10} {j:<10} {loss.item():<10.5}")
            iterationCount.append(stepCount)
            trainLosses.append(loss.item())        
        stepCount += 1

################################################################################

# plot training loss curve
plt.figure(0)
plt.cla()
plt.plot(iterationCount,trainLosses, '.k', label='Training Loss')
plt.legend()
plt.xlabel('training iteration')
plt.ylabel('loss')
plt.show()

################################################################################
