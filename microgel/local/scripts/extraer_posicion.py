#! /usr/bin/python3
import numpy as np
import argparse
import sys
import os

# Set command line parameters
argParser = argparse.ArgumentParser(description='Extract particle position and velocity data from colloids-*.csv files')
argParser.add_argument("-nciclos", help="number of cycles", required=True)
argParser.add_argument("-ninicio", help="Number of cycle to start from", required=True)
argParser.add_argument("-npaso", help="cycles per step", required=True)
argParser.add_argument("-o", help="output file name", required=True)

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

# Write header and initial line only if starting from 0
if( inicio == 0):
    fo.write("cycle,x,y,z,vx,vy,vz\n")

    # Try to read initial file
    archivo = "colloids-{ciclo:08d}.csv".format(ciclo=0)
    if os.path.exists(archivo):
        y = np.genfromtxt(archivo, dtype=np.float64, delimiter=',', skip_header=1)
        x = y[1]
        y_pos = y[2]
        z = y[3]
        vx = y[4]
        vy = y[5]
        vz = y[6]

        fo.write("0," +
                 "{x:.15e}".format(x=x) + "," +
                 "{y:.15e}".format(y=y_pos) + "," +
                 "{z:.15e}".format(z=z) + "," +
                 "{vx:.15e}".format(vx=vx) + "," +
                 "{vy:.15e}".format(vy=vy) + "," +
                 "{vz:.15e}".format(vz=vz) +
                 "\n")
    else:
        fo.write("0,0.0,0.0,0.0,0.0,0.0,0.0\n")

# Process data files from inicio+paso to cant (skip paso if inicio==0 since we already wrote initial values)
# This allows incremental processing: each call processes only new data
start_iter = paso if inicio == 0 else inicio
for i in range(start_iter, cant + 1, paso):

    # Read data file
    archivo = "colloids-{ciclo:08d}.csv".format(ciclo=i)

    if not os.path.exists(archivo):
        print(f"Warning: File {archivo} not found, skipping...")
        continue

    y = np.genfromtxt(archivo, dtype=np.float64, delimiter=',', skip_header=1)

    # Extract position and velocity data (assuming single particle or first particle)
    x = y[1]
    y_pos = y[2]
    z = y[3]
    vx = y[4]
    vy = y[5]
    vz = y[6]

    fo.write("{paso:d}".format(paso=i) + "," +
             "{x:.15e}".format(x=x) + "," +
             "{y:.15e}".format(y=y_pos) + "," +
             "{z:.15e}".format(z=z) + "," +
             "{vx:.15e}".format(vx=vx) + "," +
             "{vy:.15e}".format(vy=vy) + "," +
             "{vz:.15e}".format(vz=vz) +
             "\n")

# Close output file
fo.close()

print(f"Data successfully extracted to {fout}")
