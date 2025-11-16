#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 13:59:15 2025

@author: allen

prompt:
'I'd like to explore another type of visualization of this data. 
concentrate on a case without classes for coloring. 
imagine an animation of the 3D t-SNE data, such that frame 1 
of the animation shows components 1 and 2 at a fixed value of 
component 3 (say the minimum value of component 3). Each 
successive frame of the animation shows successive values of 
component 3 until the last frame is reached, which shows the 
2D visualization of components 1 and 2 at the fixed value of 
component 3 equal to the maximum value of component 3.
'

Claude OPUS 4.1
'I'll create an animation that slices through the 3D t-SNE data 
along the third component, showing 2D cross-sections from minimum 
to maximum values of component 3. This creates a "traveling through 
depth" effect.'


Key Features:

Main Animation:

Left panel shows the 2D slice (components 1 & 2) at the current value of component 3
Right panel shows the full 3D scatter plot with a red plane indicating the current slice position
Points are colored and sized based on their distance from the slice plane (brighter/larger = closer to the exact plane)


Slice Mechanics:

Starts at the minimum value of component 3
Progressively moves through 50 frames to the maximum value
Uses a "slice thickness" to capture points near each plane (not just exact values)
Points fade in/out smoothly as the slice passes through them


Visual Enhancements:

Color gradient shows proximity to the slice plane (brighter = closer)
Point sizes vary with distance (larger = closer)
Text overlay shows current component 3 value and point count
3D reference rotates slightly during animation for better depth perception


Static Snapshots:

Additionally creates 6 static frames showing key positions through the data
Useful for publication or detailed analysis


Customization Options:

n_frames: Adjust for smoother/faster animation (currently 50)
slice_thickness: Control how many points appear in each frame (thicker = more points visible)
interval: Animation speed in milliseconds per frame (currently 100ms)
Uncomment the save section to export as a GIF file

This visualization is particularly useful for understanding:

How data clusters are distributed along the third dimension
Whether there are distinct layers or continuous transitions in your data
The density distribution through the depth of your 3D embedding

The animation effectively "CT scans" through your 3D t-SNE result, revealing internal structure that might not be visible in a static 3D plot.


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation, PillowWriter
from sklearn.manifold import TSNE
from scipy.interpolate import griddata


tsne_results= np.load('notebook/tsne.npz')['arr_0']
unique_sites = np.load('notebook/uniqueSites.npz',allow_pickle=True)['arr_0']
lengths = np.array([len(seq) for seq in unique_sites])


# Extract components
comp1 = tsne_results[:, 0]
comp2 = tsne_results[:, 1]
comp3 = tsne_results[:, 2]

# Get min and max values for component 3
comp3_min = comp3.min()
comp3_max = comp3.max()

# Define number of frames for animation
n_frames = 50
comp3_values = np.linspace(comp3_min, comp3_max, n_frames)

# Define slice thickness (how many points to include in each frame)
# Larger values = more points visible per frame
slice_thickness = (comp3_max - comp3_min) / 20

# Create figure and axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Initialize scatter plots
scatter1 = ax1.scatter([], [], alpha=0.6, s=50, edgecolors='w', linewidth=0.5)
scatter2 = ax2.scatter([], [], alpha=0.6, s=30)

# Set up the main plot (left)
ax1.set_xlim(comp1.min() - 1, comp1.max() + 1)
ax1.set_ylim(comp2.min() - 1, comp2.max() + 1)
ax1.set_xlabel('t-SNE Component 1', fontsize=12)
ax1.set_ylabel('t-SNE Component 2', fontsize=12)
ax1.set_title('2D Slice Through 3D t-SNE', fontsize=14)
ax1.grid(True, alpha=0.3)

# Set up the position indicator plot (right)
ax2_3d = fig.add_subplot(122, projection='3d')
ax2_3d.scatter(comp1, comp2, comp3, alpha=0.1, s=1, c='gray')
ax2_3d.set_xlabel('Comp 1', fontsize=10)
ax2_3d.set_ylabel('Comp 2', fontsize=10)
ax2_3d.set_zlabel('Comp 3', fontsize=10)
ax2_3d.set_title('Current Slice Position', fontsize=12)

# Initialize the slice plane in 3D plot
xx, yy = np.meshgrid(np.linspace(comp1.min(), comp1.max(), 10),
                      np.linspace(comp2.min(), comp2.max(), 10))
zz = np.ones_like(xx) * comp3_min
plane = ax2_3d.plot_surface(xx, yy, zz, alpha=0.3, color='red')

# Text for current component 3 value
text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes, 
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Color map for points based on their distance from the slice
cmap = plt.cm.viridis

def init():
    """Initialize animation"""
    scatter1.set_offsets(np.empty((0, 2)))
    return scatter1, text, plane

def animate(frame):
    """Animation function for each frame"""
    current_comp3 = comp3_values[frame]
    
    # Find points within the slice thickness
    mask = np.abs(comp3 - current_comp3) <= slice_thickness
    points_in_slice = np.sum(mask)
    
    if points_in_slice > 0:
        # Get points in current slice
        x_slice = comp1[mask]
        y_slice = comp2[mask]
        
        # Calculate distance from the exact slice plane for coloring
        distances = np.abs(comp3[mask] - current_comp3)
        colors = 1 - (distances / slice_thickness)  # Normalize to [0, 1]
        
        # Update main scatter plot
        scatter1.set_offsets(np.c_[x_slice, y_slice])
        scatter1.set_array(colors)
        scatter1.set_cmap(cmap)
        scatter1.set_clim(0, 1)
        
        # Calculate point sizes based on distance (closer = larger)
        sizes = 50 * (1 - distances / slice_thickness) + 20
        scatter1.set_sizes(sizes)
    else:
        scatter1.set_offsets(np.empty((0, 2)))
    
    # Update the slice plane position in 3D plot
    global plane
    plane.remove()
    zz_new = np.ones_like(xx) * current_comp3
    plane = ax2_3d.plot_surface(xx, yy, zz_new, alpha=0.3, color='red')
    
    # Update text
    text.set_text(f'Component 3: {current_comp3:.2f}\n'
                  f'Points in slice: {points_in_slice}\n'
                  f'Frame: {frame+1}/{n_frames}')
    
    # Rotate 3D view slightly for better visualization
    ax2_3d.view_init(elev=20, azim=45 + frame * 2)
    
    return scatter1, text, plane

# Create animation
anim = FuncAnimation(fig, animate, init_func=init, frames=n_frames, 
                     interval=100, blit=False, repeat=True)

# Display the animation
plt.tight_layout()
plt.show()

# Save animation as GIF (optional)
# Uncomment the following lines to save the animation
# writer = PillowWriter(fps=10)
# anim.save('tsne_slice_animation.gif', writer=writer)
# print("Animation saved as 'tsne_slice_animation.gif'")

# Alternative: Create static snapshots at key positions
fig2, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

# Select 6 key positions through the data
key_positions = np.linspace(comp3_min, comp3_max, 6)

for idx, (ax, position) in enumerate(zip(axes, key_positions)):
    # Find points near this position
    mask = np.abs(comp3 - position) <= slice_thickness
    
    if np.sum(mask) > 0:
        x_slice = comp1[mask]
        y_slice = comp2[mask]
        distances = np.abs(comp3[mask] - position)
        colors = 1 - (distances / slice_thickness)
        
        scatter = ax.scatter(x_slice, y_slice, c=colors, cmap='viridis',
                            s=50, alpha=0.6, edgecolors='w', linewidth=0.5,
                            vmin=0, vmax=1)
        
    ax.set_xlim(comp1.min() - 1, comp1.max() + 1)
    ax.set_ylim(comp2.min() - 1, comp2.max() + 1)
    ax.set_xlabel('Component 1', fontsize=9)
    ax.set_ylabel('Component 2', fontsize=9)
    ax.set_title(f'Slice at Component 3 = {position:.2f}', fontsize=10)
    ax.grid(True, alpha=0.3)

plt.suptitle('Static Snapshots: Slicing Through 3D t-SNE Data', fontsize=14)
plt.tight_layout()
plt.show()











