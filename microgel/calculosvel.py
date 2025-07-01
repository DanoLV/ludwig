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

fo.write("cycle,vx\n")

for i in range(0, cant + 1, paso):
    # Read data file
    archivo = "colloids-{ciclo:08d}.csv".format(ciclo=i)
    y = np.genfromtxt(archivo, delimiter=',', skip_header=1)

    vx = y[4]

    fo.write("{paso:d}".format(paso=i) + "," +
             "{x:e}".format(x=vx) +
             "\n")
  
# Close output file
fo.close()