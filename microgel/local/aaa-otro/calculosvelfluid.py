#! /usr/bin/python3 
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.spatial import ConvexHull
import shutil

# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-nciclos", help="number of cycles")
argParser.add_argument("-ninicio", help="Number of cycle to start from")
argParser.add_argument("-npaso", help="cycles per step")
argParser.add_argument("-o", help="output file name")
argParser.add_argument("--idir", help="input directory", default='.')

# Parse command line arguments
try: 
    args = argParser.parse_args()
except:
    sys.exit("Could not read command line parameters")

cant = int(args.nciclos)
if(cant < 0):
    sys.exit("Invalid number of cycles")

inicio = int(args.ninicio)

print('cant:',cant)
print('inicio:',inicio)
if( inicio<0 or inicio>cant):
    sys.exit("Invalid number of cycle to start from")
    
paso = int(args.npaso)
if(paso <= 0):
    sys.exit("Invalid step value")

fout = args.o
if(fout is None):
    sys.exit("Please specify an output file name")

# Open output file in append mode
fo = open(fout, "a")
separador = ' '

# Write header and initial line only if starting from 0
if( inicio == 0):
    fo.write("cycle,c,y,z,vx,vy,vz,<vfx>,<vfy>,<vfz>,vxr,vyr,vyz\n")
    shutil.copy('vel-000000001.001-001', 'vel-000000000.001-001')

    # fo.write("0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n")

# Process data files from inicio+paso to cant (skip paso if inicio==0 since we already wrote initial values)
# This allows incremental processing: each call processes only new data
# start_iter = paso if inicio == 0 else inicio
start_iter = inicio
for i in range(start_iter, cant + 1, paso):
        
    # Read data file
    archivo = str(args.idir) + "/" + "colloids-{ciclo:08d}.csv".format(ciclo=i)
    print('idir:',str(args.idir))
    print('archivo:',archivo)
    y = np.genfromtxt(archivo, dtype=np.float64, delimiter=',', skip_header=1)

    xc = y[1]
    yc = y[2]
    zc = y[3]
    vx = y[4]
    vy = y[5]
    vz = y[6]
    
    # Read data file
    archivo = "./vel-{ciclo:09d}.001-001".format(ciclo=i)
    y = np.genfromtxt(archivo, dtype=np.float64, delimiter=[23, 23, 23])

    vfx = np.mean(y[:,0])
    vfy = np.mean(y[:,1])
    vfz =np.mean(y[:,2])
    
    vxr=vx-vfx
    vyr=vy-vfy
    vzr=vz-vfz

    fo.write("{paso:d}".format(paso=i) + "," +
             "{x:e}".format(x=xc) + "," +
             "{y:e}".format(y=yc) + "," +
             "{z:e}".format(z=zc) + "," +
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