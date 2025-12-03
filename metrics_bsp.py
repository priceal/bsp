#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to evaluate the BSP model.  

      
    
"""

# import libraries -----------------------------------------------------
# should be run right after training_bsp.py and should not need to 


###########################################################################
###########################################################################
###########################################################################

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
    
    print("total accuracy:",np.diagonal(cm).sum()/cm.sum())
 #   print("total precision:",np.diagonal(cm).sum()/cm.sum())
    recall = np.diagonal(cm)/cm.sum(axis=1)
    precision = np.diagonal(cm)/cm.sum(axis=0)
    print('{:10} {:10} {:10}'.format('class', 'recall', 'precision'))
    for n, r, p in zip(siteVocab, recall, precision):
        print(f'{n:<10} {r:<10.4} {p:<10.4}')
        

