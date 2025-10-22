#! /usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.spatial import ConvexHull
    
# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-nciclos", help="number of cycles")
argParser.add_argument("-npaso", help="cycles per step")
argParser.add_argument("-o", help="output file name")

# Parse command line arguments
try: 
    args = argParser.parse_args()
except:
    sys.exit("Could not read command line parameters")

cant = int(args.nciclos)
if(cant <= 0):
    sys.exit("Invalid number of cycles")

paso = int(args.npaso)
if(paso <= 0):
    sys.exit("Invalid step value")

fout = args.o
if(fout is None):
    sys.exit("Please specify an output file name")

# Open output file
fo = open(fout, "w")
separador = ' '

fo.write("cycle,vx,vy,vz,<vfx>,<vfy>,<vfz>,vxr,vyr,vyz\n")
fo.write("0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n")

for i in range(paso, cant + 1, paso):
    
    # Read data file
    archivo = "colloids-{ciclo:08d}.csv".format(ciclo=i)
    y = np.genfromtxt(archivo, dtype=np.float64, delimiter=',', skip_header=1)

    vx = y[4]
    vy = y[5]
    vz = y[6]
    
    # Read data file
    archivo = "vel-{ciclo:09d}.001-001".format(ciclo=i)
    y = np.genfromtxt(archivo, dtype=np.float64, delimiter=[23, 23, 23])

    vfx = np.mean(y[:,0])
    vfy = np.mean(y[:,1])
    vfz =np.mean(y[:,2])
    
    vxr=vx-vfx
    vyr=vy-vfy
    vzr=vz-vfz

    fo.write("{paso:d}".format(paso=i) + "," +
             "{x:e}".format(x=vx) + "," +
             "{y:e}".format(y=vy) + "," +
             "{z:e}".format(z=vz) + "," +
             "{x:e}".format(x=vfx) + "," +
             "{y:e}".format(y=vfy) + "," +
             "{z:e}".format(z=vfz) + "," + 
             "{x:e}".format(x=vxr) + "," +
             "{y:e}".format(y=vyr) + "," +
             "{z:e}".format(z=vzr) + 
             "\n")
  
# Close output file
fo.close()