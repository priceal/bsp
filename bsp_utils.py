#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 10:24:09 2025
@author: allen

utilities for ssp
    
"""
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

'''
###############################################################################
######################### functions/classes ###################################
###############################################################################
'''
def oneHot(string, vocab="ARNDCEQGHILKMFPSTWYV"):
    '''
    create the one-hot array for the string. Any un-recognized character 
    (or space) will return as a zero vector.

    Args:
        string (TYPE): input sequence
        vocab (TYPE, optional): symbol list. order defines encoding. 
        Defaults to "ARNDCEQGHILKMFPSTWYV".

    Returns:
        array: NxM where N is length of input string and M is length of vocab 
        string

    '''
    result = []
    for c in string:
        code = np.zeros(len(vocab))
        if c not in vocab:
            result.append(code)
            continue
        code[ vocab.find(c) ] = 1
        result.append(code)
        
    return np.array(result)

#######################################################################
def encode(string, vocab=" ARNDCEQGHILKMFPSTWYV", length=100):
    '''
    create numeric encoding for the string. 

    Args:
        string (TYPE): input sequence
        vocab (TYPE, optional): symbol list. order defines encoding. 
        Defaults to " ARNDCEQGHILKMFPSTWYV", which assumes padding is space
        and assigns to index 0

    Returns:
        array: N where N is length of input string

    '''
    string = f'{string[:length]:<{length}}'   
    result = []
    for c in string:
        if c not in vocab:
            result.append(0)
            continue
        result.append( vocab.find(c) )
        
    return np.array(result)


#######################################################################
def dataReader(filePath, siteCrop=15, sequenceCrop=500, 
               siteVocab = ' ACGTNUWSMKRYBDHV', 
               aaVocab = " ARNDCEQGHILKMFPSTWYV" ):
    '''
    reads text file of protein site, sequences data.
    format of file should be RE, site, sequence, with first line the header
    
    Args:
        filePath (string): input file path.
        siteCrop (TYPE, optional): DESCRIPTION. Defaults to 15.
        sequenceCrop (TYPE, optional): DESCRIPTION. Defaults to 500.

    Returns:
        tensor, tensor: first is data in shape (N,sequenceCrop), the second
        is target in shape (N, siteCrop)
    '''
    dataDf = pd.read_csv( filePath )
    
    siteList=[]
    sequenceList=[]
    for i in dataDf.index:
        sequenceList.append(
                encode(
                dataDf.at[i,'sequence'],
                vocab=aaVocab,
                length=sequenceCrop)
            )
        siteList.append(
                encode(
                dataDf.at[i,'site'],
                vocab=siteVocab,
                length=siteCrop)
            )

    return torch.tensor(np.array(sequenceList), dtype=torch.int, requires_grad=False), \
        torch.tensor(np.array(siteList), dtype=torch.int, requires_grad=False)

###############################################################################
class seqDataset(Dataset):
    """  """

    def __init__(self, x, y):
        '''

        Parameters
        ----------
        
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.x = x
        self.y = y

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):

        return self.x[idx], self.y[idx]
