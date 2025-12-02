#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""

import torch

###############################################################################
class classModel(torch.nn.Module):
    '''
    A model that takes a sequence and classifies into one of 6 categories:
        
        ['Type II methyltransferase', 
        'Type III methyltransferase', 
         'Type I methyltransferase', 
         'Type I specificity subunit', 
        'Type IIG restriction enzyme/methyltransferase', 
        'Type II restriction enzyme'] 
        encoded 0 to 5 as onehot vectors
 
    the aa sequence is encoded 0 to 20, with 0 a space.

    aaVocab, 21 chars
    Out[4]: ' ARNDCEQGHILKMFPSTWYV'







    '''

    def __init__(self):
        super(classModel, self).__init__()

        # encoding and embedding parameters
        self.aaCodeSize = 21      # aa code size
        self.padidx = 0           # change if needed
        self.ed = 5               # embedding dimension for aa sequence

        # embedding transformation, input tensor must have dims (*), output
        # will have (*,self.ed). since CNN uses (*,channels,length), must swap last
        # two dims afterward. Has feature that can fix padding index so that embedding
        # coeffs for that are not updated---converts to 0.
        self.embedding = torch.nn.Embedding( 
                                                self.aaCodeSize, 
                                                self.ed, 
                                                padding_idx=self.padidx 
                                            )

        # convolutional layers, in format 
        #(in_channels, out_channels, kernel_size, stride, padding, activation)
        # each convolution w/o padding reduces size by kernel-1
        # so minimum size must be greater than SUM(kernel-1)
        paramsCL = (
                         ( self.ed, 4, 11 ),
                         ( 4, 4, 11 ),
                         ( 4, 4, 11 ),
                         ( 4, 4, 11 )
                    )
       
        # layer parameters = (input chans, output chans, kernel)
        convLayers= []; batchLayers = []; 
        for inch, outch, ksize in paramsCL:
            convLayers.append( torch.nn.Conv1d(in_channels=inch,
                                           out_channels=outch,
                                           kernel_size=ksize,
                                           stride=1,
                                           bias=False,
                                           padding='valid' )
                              )
            batchLayers.append( torch.nn.BatchNorm1d(num_features=outch) )

        self.conv = torch.nn.ModuleList( convLayers )
        self.batch = torch.nn.ModuleList( batchLayers )
        
        self.pooling = torch.nn.MaxPool1d(2, stride=2)
        self.relu = torch.nn.ReLU()
        
        self.outLayer = torch.nn.Conv1d(in_channels=4,
                                       out_channels=self.siteCodeSize,
                                       kernel_size=11,
                                       stride=1,
                                       bias=False,
                                       padding='valid' )
        self.outBatch = torch.nn.BatchNorm1d(num_features=self.siteCodeSize)
        self.outActiv = torch.nn.Softmax( dim = 1 )

    ###########################################################################
    def forward(self, x ):
        '''
        Args:
            x (TYPE): data batch, with shape (N,length)
            mask (TYPE): binary mask, with shape (N,1,length) to matmult with
            any x of shape (N,channels,length)

        Returns:
            x (TYPE): DESCRIPTION.
        '''

        # since CNN uses (*, channels, length), must swap last two dims after 
        # embedding operation which produces (*, length, embedding_dims)
        x = self.embedding(x)
        x = torch.transpose(x,-1,-2)
        
        # convolutional layers --- note, if manually padded on C-term, the
        # layers after first hidden will have c-term residues influenced by
        # a non-zero padding, which will propagate through layers!
        # might want to zero out the padding after each layer in future
        for cl, bl in zip( self.conv, self.batch ):
            x = cl(x)
            x = bl(x)
            x = self.pooling(x)
            x = self.relu(x)
            
        # final output layer
        x = self.outLayer(x)
        x = self.outBatch(x)
        x = self.pooling(x)
        x = self.outActiv(x)

        return x 
