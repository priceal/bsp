#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to train the classifier model with k-fold cross-validation.

Script will load a training set and perform k-fold cross-validation
to evaluate model performance across multiple train/validation splits.
This enables hyperparameter tuning by comparing validation losses and
accuracies across different parameter settings.

Created: December 2024
"""

# import libraries -----------------------------------------------------
import numpy as np
import torch
from torch.utils.data import DataLoader, Subset
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import bsp_utils as bsp
from ClassifierModel_20251202 import classModel

# set learning parameters --------------------------------------------------
batchSize = 64        
numberEpochs = 1
learningRate = 0.1
reportCycle = 10

# cross-validation parameters ----------------------------------------------
numFolds = 5          # number of folds for cross-validation
randomSeed = 42       # for reproducibility
                
# files to load and optional file directory -----------------------------
inputTrain = 'protein_seqs_cleaned_6types_reps_train.csv'

# data parameters for screening sequence length and cropping-----------------
sequenceLimits = ( 1, 2000, 1708 )  

################################################################################
################################################################################
################################################################################

# load data -------------------------------------------
# datareader should return encoded sequence and class index (no onhot vectors)
xTrain, yTrain = bsp.classDataReader(inputTrain, seq=sequenceLimits )
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
    
# we set up dataLoader AFTER the train/val split
#dataloader = DataLoader(dataTrain, batch_size=batchSize, shuffle=True)
################################################################################

# k-fold cross-validation setup using sklearn KFold
kfold = KFold(n_splits=numFolds, shuffle=True, random_state=randomSeed)

# storage for cross-validation results
foldTrainLosses = []
foldValLosses = []
foldValAccuracies = []

print(f'\nCROSS-VALIDATION ({numFolds} folds)')
print('=' * 60)

# iterate over each fold. trainIdx and valIdx will be lists of indices, with
# (n-1)/n of indices in trainIdx and 1/n in valIdx, where n = number of folds
for fold, (trainIdx, valIdx) in enumerate(kfold.split(xTrain)):
    
    print(f'\nFOLD {fold + 1}/{numFolds}')
    print('-' * 40)
    
    # create subset datasets for this fold - these can be treated as training
    # and val sets for this fold calc
    trainSubset = Subset(dataTrain, trainIdx)
    valSubset = Subset(dataTrain, valIdx)
    
    # create dataloaders for train and validation sets
    trainLoader = DataLoader(trainSubset, batch_size=batchSize, shuffle=True)
    valLoader = DataLoader(valSubset, batch_size=batchSize, shuffle=False)
    
    # create fresh model for each fold
    model = classModel()
    optimizer = torch.optim.SGD(model.parameters(), lr=learningRate)
    lossFunction = torch.nn.CrossEntropyLoss()
    
    # storage for this fold's training losses
    trainLosses = []
    
    # training loop for this fold
    print('{:10} {:10} {:10}'.format('epoch', 'batch', 'loss-train'))
    stepCount = 0
    for i in range(numberEpochs):
        model.train()
        for j, batch in enumerate(trainLoader):
            xx, yy = batch[0], batch[1]
            prediction = model(xx)
            
        # assumiung target yy are class indices and prediction are raw logits
        # of dim 6
        loss = lossFunction(prediction, yy)
            
        # cycle optimizer
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # record and report training loss
        if stepCount % reportCycle == 0:
            print(f"{i:<10} {j:<10} {loss.item():<10.5}")
            trainLosses.append(loss.item())
        stepCount += 1
    
    # validation phase for this fold
    model.eval()
    with torch.no_grad():
        
        xx, yy = valSubset[:]   # returns 2-tuple of val x's and y's
        prediction = model(xx)
        valLoss = lossFunction(prediction, yy)
        
        # calculate accuracy
        predicted_classes = torch.argmax(prediction, dim=1)
        correct = (predicted_classes == yy).sum().item()
        total = yy.size(0)

    # compute average validation loss and accuracy for this fold
    valAccuracy = correct / total
    
    # store fold results
    foldTrainLosses.append(trainLosses)
    foldValLosses.append(valLoss)
    foldValAccuracies.append(valAccuracy)
    
    print(f'\nFold {fold + 1} Results:')
    print(f'  Validation Loss: {valLoss:.5f}')
    print(f'  Validation Accuracy: {valAccuracy:.4f} ({correct}/{total})')

################################################################################

# summarize cross-validation results
print('\n' + '=' * 60)
print('CROSS-VALIDATION SUMMARY')
print('=' * 60)
print(f'{"Fold":<10} {"Val Loss":<15} {"Val Accuracy":<15}')
for fold in range(numFolds):
    print(f'{fold + 1:<10} {foldValLosses[fold]:<15.5f} {foldValAccuracies[fold]:<15.4f}')

print('-' * 40)
print(f'{"Mean":<10} {np.mean(foldValLosses):<15.5f} {np.mean(foldValAccuracies):<15.4f}')
print(f'{"Std":<10} {np.std(foldValLosses):<15.5f} {np.std(foldValAccuracies):<15.4f}')

# plot training losses for all folds
plt.figure(1, figsize=(10, 6))
for fold in range(numFolds):
    plt.plot(foldTrainLosses[fold], '.', label=f'Fold {fold + 1}', alpha=0.7)
plt.legend()
plt.xlabel('Training Iteration')
plt.ylabel('Loss')
plt.title('Training Loss Across Folds')
plt.show()

# plot validation results summary
plt.figure(2, figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.bar(range(1, numFolds + 1), foldValLosses)
plt.axhline(y=np.mean(foldValLosses), color='r', linestyle='--', label='Mean')
plt.xlabel('Fold')
plt.ylabel('Validation Loss')
plt.title('Validation Loss by Fold')
plt.legend()

plt.subplot(1, 2, 2)
plt.bar(range(1, numFolds + 1), foldValAccuracies)
plt.axhline(y=np.mean(foldValAccuracies), color='r', linestyle='--', label='Mean')
plt.xlabel('Fold')
plt.ylabel('Validation Accuracy')
plt.title('Validation Accuracy by Fold')
plt.legend()
plt.tight_layout()
plt.show()

################################################################################