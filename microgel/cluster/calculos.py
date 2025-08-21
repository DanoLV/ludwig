#! /usr/bin/env python

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

# Read command line parameters
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

fo.write("cycle,xcm,ycm,zcm,Ixcm,Iycm,Izcm,<lbondm>,density,volume\n")

for i in range(0, cant + 1, paso):
    # Read data file
    archivo = "colloids-{ciclo:08d}.csv".format(ciclo=i)
    y = np.genfromtxt(archivo, delimiter=',', skip_header=1)
    
    # Need at least 2 monomers
    if (len(y.shape) <= 1): 
        print("Simulation has less than 2 monomers")
        break

    nmon = y.shape[0]
    y = y[y[:, 0].argsort()]

    # Process data
    xcm = np.mean(y[:, 1])
    ycm = np.mean(y[:, 2])
    zcm = np.mean(y[:, 3])

    Ixcm = np.mean((y[:, 1] - xcm)**2)
    Iycm = np.mean((y[:, 2] - ycm)**2)
    Izcm = np.mean((y[:, 3] - zcm)**2)

    # coords: N x 3 array with particle positions
    hull = ConvexHull(y[:, 1:4])
    volume = hull.volume
    density = float(nmon / volume)

    # Compute average bond length
    total_bond_length = 0
    total_bonds = 0
    
    for j in y:
        nbonds = int(j[8])
        dist_sum = 0
        n = 9  # Column where bonded monomers start
        # Sum the bond distances
        for n in range(n, n + nbonds, 1):
            idbond = int(j[n]) - 1
            dist_sum += np.sqrt((j[1] - y[idbond, 1])**2 +
                                (j[2] - y[idbond, 2])**2 +
                                (j[3] - y[idbond, 3])**2)
            
        total_bond_length += dist_sum
        total_bonds += nbonds
        
    avg_bond_length = total_bond_length / total_bonds

    # Write output
    fo.write("{paso:d},".format(paso=i) +
             "{x:f},{y:f},{z:f},".format(x=xcm, y=ycm, z=zcm) +
             "{x:f},{y:f},{z:f},".format(x=Ixcm, y=Iycm, z=Izcm) +
             "{lbondm:f},".format(lbondm=avg_bond_length) +
             "{density:f},".format(density=density) +
             "{volume:f}\n".format(volume=volume))
  
# Close output file
fo.close()