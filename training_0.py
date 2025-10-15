#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 10:24:09 2025
@author: allen



features to include:
    1) employ torch optimization tools      
    2) use scikit learn metrics             
    3) split test/train                     
    4) employ batch optimization            
    5) drop out back propagation
    6) use weights to correct imbalance     
    7) add embedding                        
    8) padding to handle variable length
    9) add attention layer
    
"""

import os
import numpy as np
import torch
from torch.utils.data import DataLoader
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.utils.class_weight import compute_class_weight
import matplotlib.pyplot as plt

import bsp_utils as bsp
from model_20251014 import bspModel

'''
###############################################################################
############################# main ############################################
###############################################################################
'''



# learning parameters
numBatches = 1 # if non-zero, ignore batchSize and set to N/numBatches
batchSize = 0 # only use if numBatches = 0
numberEpochs = 200
learningRate = 0.1

reportCycle = 1

refine = True  # creates new model if False
weights = None    # None: unweighted. 
                    # (WH, WE, WC): use fixed weights
                    # 'calc' : calculated weights to use
                
# file to load and optional file directory---can leave undefined '' or '.'
inputTrain = 'F5-8-17_gold.csv'
fileDirectory = '.'

# data parameters
sequenceLimits = (100,600,502)  # screen data for seq lengths in this interval
siteLimits = (4,10,6)  # screen data for seq lengths in this interval

siteVocab = 'NACGTUWSMKRYBDHV'    # uses N as padding as well as internal N
aaVocab = " ARNDCEQGHILKMFPSTWYV" 
###########################################################################
###########################################################################
###########################################################################

# load data -------------------------------------------
xTrain, yTrain, df = bsp.dataReader(os.path.join(fileDirectory, inputTrain), 
                                site=siteLimits, 
                                seq = sequenceLimits
                              )
yTrain.swapaxes_(1, 2)   # put in correct order for CNN
dataTrain = bsp.seqDataset(xTrain, yTrain ) # needed for batches

###########################################################################
# print data/batch stats ------------------------------------
print("DATA SET")
print("{:<20} {:<15} {:<15}".format('DATA', 'ENTRIES', 'LENGTH'))
rows = ['training data']
ds = [xTrain]
for r, d in zip(rows, ds):
    a, b = d.shape
    print(f"{r:<20} {a:<15} {b:<15}")

###########################################################################
# create weights for classes--should broadcast correctly in loss calc
'''
yMask = yTrain.sum(axis=1).detach().numpy()
yClasses = np.argmax(yTrain.detach().numpy(), axis=1)
yAdjusted = yClasses + yMask
uniqueClasses, numClasses = np.unique(yAdjusted, return_counts=True)
uniqueClasses = np.delete( uniqueClasses, np.where(uniqueClasses==0))
numClasses = numClasses[1:]
if not weights:
    weights=[1.0,1.0,1.0]
elif weights=='calc':
    weights = numClasses.sum()/numClasses/3 # dims=(3)
# add dim in place to get dims = (3,1) for broadcasting
weights = torch.tensor(weights).unsqueeze_(1)
'''

# determine batch size and number of batches -----------------------    
if numBatches > 0:
    batchSize = int(len(xTrain)/numBatches)
else:
    numBatches = int(len(xTrain)/batchSize)

'''
print('{:<10} {:<10} {:<10} {:<10} {:<10}'.format('index','label','count','fraction','weight') )
for i,tl in enumerate(targetLabels):
    print('{:<10} {:<10} {:<10} {:<10.4} {:<10.4}'.format(
        i,tl,int(numClasses[i]),numClasses[i]/numClasses.sum(), float(weights[i]) 
        ) )
'''
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
plt.figure(1)
optimizer = torch.optim.SGD(model.parameters(), lr=learningRate)
#lossFunc = torch.nn.CrossEntropyLoss( weight=weights,ignore_index=-1)
print('\nOPTIMIZATION')
print('{:10} {:10} {:10} {:10}'.format('epoch','batch','loss-train','loss-test') )
stepCount = 0
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
            print(f"{i:<10} {j:<10} {loss.item():<10.5} ")
            plt.plot([stepCount], [loss.item()], '.k')
        
        stepCount += 1

plt.show()

###########################################################################
# metrics -------------------------------------------------------
# must convert probability-logits to one-hots ---
# convert max logit value to 1, others 0
print('\nFINAL METRICS')
titles = ['training']
xSets = [ xTrain ]
ySets = [ yTrain ]
for t,xs,ys in zip(titles,xSets,ySets):

    # problem: whenever zero-padding is encountered, argmax returns class 0!
    # that is checked against some 'random' prediction !!
    # makes it look worse than it is!
    print(t+' set performance')
    mask = (ys.sum(axis=1)).flatten()
    yCheck = np.argmax(ys.detach().numpy(), axis=1).flatten()
    prediction = model(xs)
    pCheck = np.argmax(prediction.detach().numpy(), axis=1).flatten()
    cm = confusion_matrix(yCheck, pCheck, sample_weight=mask)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=siteVocab)
    disp.plot()
    recall = np.diagonal(cm)/cm.sum(axis=1)
    precision = np.diagonal(cm)/cm.sum(axis=0)
    print('{:10} {:10} {:10}'.format('class', 'recall', 'precision'))
    for n, r, p in zip(siteVocab, recall, precision):
        print(f'{n:<10} {r:<10.4} {p:<10.4}')
        

