#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 10:24:09 2025
@author: allen

utilities for bsp
    
"""
import numpy as np
import pandas as pd

import torch
from torch.utils.data import Dataset



classDict = {'Type II methyltransferase': 0,
                'Type III methyltransferase': 1, 
                 'Type I methyltransferase': 2,
                 'Type I specificity subunit': 3, 
                'Type IIG restriction enzyme/methyltransferase': 4, 
                'Type II restriction enzyme': 5
                }
                


'''
###############################################################################
######################### functions/classes ###################################
###############################################################################
'''
def oneHot(string, vocab=" ARNDCEQGHILKMFPSTWYV", length=100):
    '''
    create the one-hot array for the string. Any un-recognized character 
    (or space) will return as a zero vector.

    Args:
        string (TYPE): input sequence, any spaces will be replaced with 
        padding char
        vocab (TYPE, optional): symbol list. order defines encoding. 
        Defaults to " ARNDCEQGHILKMFPSTWYV". expacts first character to be
        padding..

    Returns:
        array: NxM where N is length of input string and M is length of vocab 
        string

    '''
    # crop and left justify, extending to length if needed, replace all 
    # spaces with padding char
    string = f'{string[:length]:<{length}}'.replace(' ',vocab[0])
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
        string (str): input sequence, any spaces will be replaced with 
        padding char - assumed first char in vocab string
        vocab (TYPE, optional): symbol list. order defines encoding. 
        Defaults to " ARNDCEQGHILKMFPSTWYV". expacts first character to be
        padding..

    Returns:
        array: N where N is length of input string

    '''
    # crop and left justify, extending to length if needed, replace all 
    # spaces with padding char
    string = f'{string[:length]:<{length}}'.replace(' ',vocab[0])

    result = []
    for c in string:
        if c not in vocab:
            result.append(0)
            continue
        result.append( vocab.find(c) )
        
    return np.array(result)


#######################################################################
def oneHotClass(string, cdict=classDict ):
    '''
    create numeric encoding for the string. 

    Args:
        string (str): input sequence, any spaces will be replaced with 
        padding char - assumed first char in vocab string
        vocab (TYPE, optional): symbol list. order defines encoding. 
        Defaults to " ARNDCEQGHILKMFPSTWYV". expacts first character to be
        padding..

    Returns:
        array: N where N is length of input string

    '''
   
    oneHot = np.zeros(6,dtype=int)
    try:
        oneHot[ cdict[string] ] = 1
    except:
        print( string )
        print( cdict )
    return oneHot

#######################################################################
def dataReader(filePath, site=(0,15,15), seq=(0,500,500), 
               siteVocab = " NACGTUWSMKRYBDHV", 
               aaVocab = " ARNDCEQGHILKMFPSTWYV" ):
    '''
    reads text file of protein site, sequences data.
    format of file should be RE, site, sequence, with first line the header
    
    Args:
        filePath (string): input file path.
        site : (min,max,crop), optional
            filters inputs between min/max and crops. The default is (0,15,15).
        seq : (min,max,crop), optional
            filters inputs between min/max and crops. The default is (0,15,15).
        siteVocab : TYPE, optional
            for encoding. The default is " NACGTUWSMKRYBDHV".
        aaVocab : TYPE, optional
            for encoding. The default is " ARNDCEQGHILKMFPSTWYV".

    Returns:
        tensor, tensor: first is data in shape (N,seqCrop), the second
        is target in shape (N, siteCrop, len(siteVocab) )
    '''
    try:
        dataDf = pd.read_csv( filePath )
    except:
        print('error reading data frame:', filePath)
        
    siteMin, siteMax, siteCrop = site
    seqMin, seqMax, seqCrop = seq
    
    # filter on site/seq lengths
    dataDf = dataDf[ dataDf.site.str.len() >= siteMin ]
    dataDf = dataDf[ dataDf.site.str.len() <= siteMax ]
    dataDf = dataDf[ dataDf.sequence.str.len() >= seqMin ]
    dataDf = dataDf[ dataDf.sequence.str.len() <= seqMax ]

    siteList=[]
    sequenceList=[]
    for i in dataDf.index:
        # the aa sequence is encoding numerically to char position in aaVocab
        sequenceList.append(
                encode(
                dataDf.at[i,'sequence'],
                vocab=aaVocab,
                length=seqCrop)
            )
        # the site dna sequence is one-hot encoded using siteVocab
        siteList.append(
                oneHot(
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





#######################################################################
def classDataReader(filePath, seq=(0,500,500), 
                    aaVocab = " ARNDCEQGHILKMFPSTWYV" ):
    '''
    reads csv file of protein type/sequence data.
    format of file should be RE, site, sequence, with first line the header
    
    Args:
        filePath (string): input file path.
        seq : (min,max,crop), optional
            filters inputs between min/max and crops. The default is (0,500,500).
        aaVocab : TYPE, optional
            for encoding. The default is " ARNDCEQGHILKMFPSTWYV".

    Returns:
        tensor, tensor: first is data in shape (N,seqCrop), the second
        is target in shape (N, 6)
    '''
    try:
        dataDf = pd.read_csv( filePath )
    except:
        print('error reading data frame:', filePath)
        
    seqMin, seqMax, seqCrop = seq
    
    # filter on site/seq lengths
    dataDf = dataDf[ dataDf.sequence.str.len() >= seqMin ]
    dataDf = dataDf[ dataDf.sequence.str.len() <= seqMax ]

    classList=[]
    sequenceList=[]
    for i in dataDf.index:
        # the aa sequence is encoding numerically to char position in aaVocab
        sequenceList.append(
                encode(
                dataDf.at[i,'sequence'],
                vocab=aaVocab,
                length=seqCrop)
            )
        # the site dna sequence is one-hot encoded using siteVocab
        classList.append(
                oneHotClass( 
                dataDf.at[i,'type'],
                cdict = classDict )
            )

    return torch.tensor(np.array(sequenceList), dtype=torch.int, requires_grad=False), \
        torch.tensor(np.array(classList), dtype=torch.int, requires_grad=False)
