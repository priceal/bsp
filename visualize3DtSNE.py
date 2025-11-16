#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 13:59:15 2025

@author: allen
"""
import numpy as np
import matplotlib.pyplot as plt

tsne_results= np.load('notebook/tsne-2.npz')['arr_0']
#unique_sites = np.load('notebook/uniqueSites-0.npy',allow_pickle=True)['arr_0']
#lengths = np.array([len(seq) for seq in unique_sites])
lengths = None

# Create 3D plot
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

if lengths:
    scatter = ax.scatter(tsne_results[:, 0], tsne_results[:, 1], tsne_results[:, 2],
                     c=lengths, cmap='viridis', alpha=0.6, s=10,
                     edgecolors='w', linewidth=0.5)
else:
    scatter = ax.scatter(tsne_results[:, 0], tsne_results[:, 1], tsne_results[:, 2],
                     alpha=0.6, s=10,
                     edgecolors='w', linewidth=0.5)

    

# Add colorbar
plt.colorbar(scatter, ax=ax, pad=0.1)

# Set labels and title
ax.set_xlabel('t-SNE Component 1')
ax.set_ylabel('t-SNE Component 2')
ax.set_zlabel('t-SNE Component 3')
ax.set_title('3D t-SNE Visualization of DNA Binding Sites (colored by length)', fontsize=10)



