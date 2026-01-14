#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import Rbf

# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", help="output file name")

# Parse command line arguments
try: 
    args = argParser.parse_args()
except:
    sys.exit("Could not read command line parameters")

fdata = args.i
if(fdata is None):
    sys.exit("Please specify an input data file name")
    
# Plot definitions
plt.rc('axes', labelsize=30) 
plt.rc('ytick', labelsize=24)
plt.rc('xtick', labelsize=24)

# subplots
figure, axis = plt.subplots(1) 

#------------------------------------------
# Read data file
y1 = np.genfromtxt(fdata,delimiter=';')

#------------------------------------------
# subplots data
print("size=%d\n",y1[:,0].size)
x1 = np.linspace(0,y1[:,0].size,y1[:,0].size,dtype=int)
psi = y1[:,1]
qrho = y1[:,5]

#------------------------------------------
# Plotting the Graph

# Vx
axis.plot(x1,psi,'g')
# Vy
# axis.plot(x1,qrho,'r')


# # Vx
# axis.plot(x1,vfx,'g',linestyle='dotted')
# # Vy
# axis.plot(x1,vfy,'r',linestyle='dotted')
# # Vz
# axis.plot(x1,vfz,'b',linestyle='dotted')

# # Vx
# axis.plot(x1,vxr,'g',linestyle='dashed')
# # Vy
# axis.plot(x1,vyr,'r',linestyle='dashed')
# # Vz
# axis.plot(x1,vzr,'b',linestyle='dashed')

# Plot definitions
# axis.set_ylabel('Psi/QRho')
axis.set_ylabel('Psi')
# axis.set_ylabel('QRho')
axis.set_xlabel('Cycle')
# axis.set_ylim(-5.22e-4,-5.205e-4)
axis.set_xlim(34000,35000)
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

plt.savefig("variacioncargacampo.png", dpi=600, bbox_inches='tight') 