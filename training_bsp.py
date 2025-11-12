#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to train the BSP model.  

Script will load a training and a test set and perform a given number of 
epochs of stochastic gradient descent. During optimization it tracks the 
training loss and the test loss. At the end of optimization it will plot 
these losses versus iteration of training and also calculates some 
evaluation metrics including the recall and precision per class and 
plots confusion matrices for both the training and test sets. 

Filtering of input sequences and sites can be done with sequence 
limits and site limits variables. Accepted data entries can also be 
cropped using the crop parameters. 

features included:
    1) torch optimization tools      
    2) scikit learn metrics             
    3) test/train split                 
    4) batch optimization            
    
"""

# import libraries -----------------------------------------------------
import os
import numpy as np
import torch
from torch.utils.data import DataLoader
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import bsp_utils as bsp
from model_20251014 import bspModel

# set learning parameters --------------------------------------------------
numBatches = 1 # if non-zero, ignore batchSize and set to N/numBatches
batchSize = 0 # only use if numBatches = 0
numberEpochs = 1000
learningRate = 0.1
reportCycle = 10
refine = False     # creates new model if False
                
# files to load and optional file directory -----------------------------
# can leave undefined '' or '.'
inputTrain = 'train.csv'
inputTest = 'test.csv'
fileDirectory = '.'

# data parameters for screening sequence length and cropping-----------------
# format ( min, max, crop )
sequenceLimits = ( 1, 2000, 502 )  # crop must be 502 for model_20251014
siteLimits = ( 1, 20, 6 )   # note: crop must be 6 for model_20251014 !!!!

###########################################################################
###########################################################################
###########################################################################

# symbol usage, uses space as padding -------------------------------------
siteVocab = ' ACGTNUWSMKRYBDHV'
aaVocab = " ARNDCEQGHILKMFPSTWYV"

# load data -------------------------------------------
xTrain, yTrain = bsp.dataReader(os.path.join(fileDirectory, inputTrain), 
                                site=siteLimits, 
                                seq = sequenceLimits
                              )
yTrain.swapaxes_(1, 2)   # put in correct order for CNN
xTest, yTest = bsp.dataReader(os.path.join(fileDirectory, inputTest), 
                                site=siteLimits, 
                                seq = sequenceLimits
                              )
yTest.swapaxes_(1, 2)   # put in correct order for CNN
dataTrain = bsp.seqDataset(xTrain, yTrain ) # needed for batches

# print data/batch stats ------------------------------------
print("DATA SET")
print("{:<20} {:<15} {:<15}".format('DATA', 'ENTRIES', 'LENGTH'))
rows = ['training data', 'test data']
ds = [xTrain, xTest]
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
    model = bspModel()
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
#lossFunc = torch.nn.CrossEntropyLoss( weight=weights,ignore_index=-1)
print('\nOPTIMIZATION')
print('{:10} {:10} {:10} {:10}'.format('epoch','batch','loss-train','loss-test') )
stepCount = 0
trainLosses = []
testLosses = []
for i in range(numberEpochs):
    for j, batch in enumerate(dataloader):
        
        xx, yy = batch[0], batch[1]
#        yyMask = yy.sum(axis=1).unsqueeze(1)
        prediction = model(xx)
        
        lossTerms = -yy*torch.log(prediction)
        loss = lossTerms.sum()/yy.shape.numel() # normalize by num of site bps
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # print out if report cycle done
        if stepCount % reportCycle == 0:
            
            # calc test loss
            testPrediction = model( xTest ) 
            testLossTerms = -yTest*torch.log(testPrediction)
            testLoss = testLossTerms.sum()/yTest.shape.numel() # normalize by num of AAs

            print(f"{i:<10} {j:<10} {loss.item():<10.5} {testLoss.item():<10.5}")
            trainLosses.append(loss.item())
            testLosses.append(testLoss.item())
        
        stepCount += 1
        
plt.figure(1)
plt.plot(trainLosses, '.k', label='Training Loss')
plt.plot(testLosses, '.r', label='Test Loss')
plt.legend()
plt.xlabel('training iteration')
plt.ylabel('loss')
plt.show()

###########################################################################
# metrics -------------------------------------------------------
# must convert probability-logits to one-hots ---
# convert max logit value to 1, others 0
print('\nFINAL METRICS')
titles = ['training', 'test']
xSets = [ xTrain, xTest ]
ySets = [ yTrain, yTest ]
for t,xs,ys in zip(titles,xSets,ySets):

    # problem: whenever zero-padding is encountered, argmax returns class 0!
    # that is checked against some 'random' prediction !!
    # makes it look worse than it is!
    print('\n'+t+' set performance')
    mask = (ys.sum(axis=1)).flatten()
    yCheck = np.argmax(ys.detach().numpy(), axis=1).flatten()
    prediction = model(xs)
    pCheck = np.argmax(prediction.detach().numpy(), axis=1).flatten()
    cm = confusion_matrix(yCheck, pCheck, sample_weight=mask,
                          labels=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=siteVocab)
    disp.plot()
    plt.title(t+' data')
    plt.show()
    recall = np.diagonal(cm)/cm.sum(axis=1)
    precision = np.diagonal(cm)/cm.sum(axis=0)
    print('{:10} {:10} {:10}'.format('class', 'recall', 'precision'))
    for n, r, p in zip(siteVocab, recall, precision):
        print(f'{n:<10} {r:<10.4} {p:<10.4}')
        

