#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    
"""

import torch

###############################################################################
class bspModel(torch.nn.Module):
    '''
    CNN that takes a numerical coding of sequence as input and outputs a one-hot
    rep of the secondary structure prediction. 
    input is expected to be (N,length), where each element is an index from 0
    to 20. 0 is assumed padding symbol.
    each layer is a 1d conv.
    followed by an activation.
    
    this is just a straw man model or null model which will be a CNN that 
    downsizes the input sequence of length 500 down to output site of 15
    
    architecture: a few layers of CNN/pooling, followed by a few fully 
    connected layers

    aaVocab, 21 chars
    Out[4]: ' ARNDCEQGHILKMFPSTWYV'

    siteVocab, 17 chars
    Out[5]: ' ACGTNUWSMKRYBDHV'
    
    
    kernel	11
	
    layer	
    input	502
    cnn	    492
    pool	246
    cnn	    236
    pool	118
    cnn	    108
    pool	54
    cnn	    44
    pool	22
    cnn	    12
    pool	6

    '''

    def __init__(self):
        super(bspModel, self).__init__()

        # encoding and embedding parameters
        self.siteCodeSize = 17    # binding site code
        self.aaCodeSize = 21      # aa code
        self.padidx = 0           # change this for another padding index
        self.ed = 5               # embedding dimension

        # convolutional layers, in format 
        #(in_channels, out_channels, kernel_size, stride, padding, activation)
        # each convolution w/o padding reduces size by kernel-1
        # so minimum size must be greater than SUM(kernel-1)
        paramsCL = (
                         ( self.ed, 4, 11, ),
                         ( 4, 4, 11 ),
                         ( 4, 4, 11 ),
                         ( 4, 4, 11 )
                    )
       
        # embedding transformation, input tensor must have dims (*), output
        # will have (*,self.ed). since CNN uses (*,channels,length), must swap last
        # two dims afterward. Has feature that can fix padding index so that embedding
        # coeffs for that are not updated---converts to 0.
        self.embedding = torch.nn.Embedding( 
                                                self.aaCodeSize, 
                                                self.ed, 
                                                padding_idx=self.padidx 
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
        
        self.outLayer = torch.nn.Conv1d(in_channels=4,
                                       out_channels=4,
                                       kernel_size=11,
                                       stride=1,
                                       bias=False,
                                       padding='valid' )
        self.outBatch = torch.nn.BatchNorm1d(num_features=4)
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
        # might need to zero out the padding after each layer!
        for cl, bl in zip( self.conv, self.batch ):
            x = cl(x)
            x = bl(x)
            x = torch.nn.MaxPool1d(2, stride=2)
            x = torch.nn.ReLU(x)
 
        x = self.outLayer(x)
        x = self.outBatch(x)
        x = torch.nn.MaxPool1d(2, stride=2)
        x = self.outActiv(x)


        return x 
