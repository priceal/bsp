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


    # CrossEntropyLoss for multi-class classification (6 classes)
    # Note: CrossEntropyLoss expects raw logits, so softmax should be removed from
    # model during training. Alternatively, use NLLLoss with log of softmax output.
    # Here we use CrossEntropyLoss and convert one-hot targets to class indices.




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
                         ( self.ed, 9, 427 ),
                         ( 9, 13, 160 ),
                         ( 13, 17, 60 )
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
        
        self.lastConv = torch.nn.Conv1d(in_channels=17,
                                       out_channels=21,
                                       kernel_size=91,
                                       stride=1,
                                       bias=False,
                                       padding='valid' )
        
        self.lastBatch = torch.nn.BatchNorm1d(num_features=21)

        
# INSERT DEFINITION OF FULLY CONNECTED LATER HERE

        self.outBatch = torch.nn.BatchNorm1d(num_features=self.siteCodeSize)
        self.outActiv = torch.nn.Softmax( dim = 1 )

    ###########################################################################
    def forward(self, x ):
      

        # since CNN uses (*, channels, length), must swap last two dims after 
        # embedding operation which produces (*, length, embedding_dims)
        x = self.embedding(x)
        x = torch.transpose(x,-1,-2)
        
        # convolutional layers --- 
        for cl, bl in zip( self.conv, self.batch ):
            x = cl(x)
            x = bl(x)
            x = self.pooling(x)
            x = self.relu(x)
            
        # final cnn layer
        x = self.lastConv(x)
        x = self.lastBatch(x)
        x = self.relu(x)
        
        # fc layer
        # INSERT CORRECT APPLICATON OF FULLY CONNECTED LAYER HERE

        return x 
