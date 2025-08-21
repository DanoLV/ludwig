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
vx = y1[:,1]
vy = y1[:,2]
vz = y1[:,3]

#------------------------------------------
# Plotting the Graph

# Vx
axis.plot(x1,vx,'g')
# Vy
axis.plot(x1,vy,'r')
# Vz
axis.plot(x1,vz,'b')

# Plot definitions
axis.set_ylabel('Vx/Vy/Vz')
axis.set_xlabel('Cycle')
# axis.set_ylim(-5.22e-4,-5.205e-4)
# axis.set_xlim(900000,)
axis.grid(True, which='both',linewidth=3,linestyle='--')
axis.tick_params(width=3)
for spine in ['top','bottom','left','right']:
    axis.spines[spine].set_linewidth(3)

# # Display plot
# plt.subplots_adjust(hspace=0.4)
# # mng = plt.get_current_fig_manager()
# # mng.resize(*mng.window.maxsize())
# # print("paso")
# plt.show()

plt.savefig("velocidades.png", dpi=600, bbox_inches='tight') 