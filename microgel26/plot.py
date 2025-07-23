#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import Rbf

# Plot style definitions
plt.rc('axes', labelsize=12) 
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

# Create 3 subplots (stacked vertically)
fig, axes = plt.subplots(3, figsize=(10, 12), sharex=True)
plt.subplots_adjust(hspace=0.4)

#------------------------------------------
# Read data file
y1 = np.genfromtxt('datos.csv', delimiter=',', skip_header=1)

#------------------------------------------
# Extract data
x1 = y1[:, 0]

# Inertia moments (normalized)
traza = y1[:, 4] + y1[:, 5] + y1[:, 6]
Ix = y1[:, 4] / traza
Iy = y1[:, 5] / traza
Iz = y1[:, 6] / traza

# Medium bond length
lbond = y1[:, 7]

# Density
densidad = y1[:, 8]

#------------------------------------------
# Plot inertia moments
axes[0].plot(x1, Ix, 'g', label='Ix/tr(I)')
axes[0].plot(x1, Iy, 'b', label='Iy/tr(I)')
axes[0].plot(x1, Iz, 'k', label='Iz/tr(I)')
axes[0].set_ylim(0,1)
axes[0].set_ylabel('Normalized Inertia Moments')
axes[0].legend()
axes[0].grid(True, which='both', linewidth=1.5, linestyle='--')
axes[0].tick_params(width=2)
for spine in ['top', 'bottom', 'left', 'right']:
    axes[0].spines[spine].set_linewidth(2)

# Plot medium bond length
axes[1].plot(x1, lbond, 'r')
axes[1].set_ylabel('Mean Bond Length')
axes[1].grid(True, which='both', linewidth=1.5, linestyle='--')
axes[1].tick_params(width=2)
for spine in ['top', 'bottom', 'left', 'right']:
    axes[1].spines[spine].set_linewidth(2)

# Plot density
axes[2].plot(x1, densidad, 'c')
axes[2].set_ylabel('Density')
axes[2].set_xlabel('Cycle')
axes[2].grid(True, which='both', linewidth=1.5, linestyle='--')
axes[2].tick_params(width=2)
for spine in ['top', 'bottom', 'left', 'right']:
    axes[2].spines[spine].set_linewidth(2)

# Show plot
plt.show()
