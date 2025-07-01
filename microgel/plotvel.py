#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import Rbf

# Plot definitions
plt.rc('axes', labelsize=30) 
plt.rc('ytick', labelsize=24)
plt.rc('xtick', labelsize=24)

# subplots
figure, axis = plt.subplots(1) 

#------------------------------------------
# Read data file
y1 = np.genfromtxt('datos.csv',delimiter=',',skip_header=1)

#------------------------------------------
# subplots data

x1 = y1[:,0]
dx = np.diff(y1[:,0])
vel= y1[:,1]
# dvel = np.diff(y1[:,1]) / np.diff(y1[:,0])
# dvel2 = np.diff(dvel) / np.diff(y1[1:,0])

#------------------------------------------
# Plotting the Graph

# Vx
axis.plot(x1,vel,'g')
# axis.plot(x1[1:],dvel2,'r')
# axis.plot(x1[2:],dvel2,'r')

# Plot definitions
axis.set_ylabel('Vx')
axis.set_xlabel('Cycle')
# axis.set_ylim(-1e-12,1e-12)
# axis.set_xlim(150000,250000)
axis.grid(True, which='both',linewidth=3,linestyle='--')
axis.tick_params(width=3)
for spine in ['top','bottom','left','right']:
    axis.spines[spine].set_linewidth(3)

# Display plot
plt.subplots_adjust(hspace=0.4)
# mng = plt.get_current_fig_manager()
# mng.resize(*mng.window.maxsize())
# print("paso")
plt.show()