#! /usr/bin/python3
import numpy as np
import sys
import argparse

# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-nciclos", help="number of cycles", type=int, required=True)
argParser.add_argument("-ninicio", help="Number of cycle to start from", type=int, required=True)
argParser.add_argument("-npaso", help="cycles per step", type=int, required=True)
argParser.add_argument("-o", help="output file name", required=True)

args = argParser.parse_args()

cant = args.nciclos
inicio = args.ninicio
paso = args.npaso
fout = args.o

print(f'cant: {cant}')
print(f'inicio: {inicio}')

# Open output file in append mode
with open(fout, "a") as fo:
    # Write header only if starting from 0
    if inicio == 0:
        fo.write("cycle,vx,vy,vz,<vfx>,<vfy>,<vfz>,vxr,vyr,vyz\n")
        # Write initial values at cycle 0
        fo.write("0,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00,0.0e+00\n")

    # Process each time step (skip paso if inicio==0 since we already wrote initial values)
    start_iter = paso if inicio == 0 else inicio
    for i in range(start_iter, cant + 1, paso):
        # Read velocity field data
        archivo = f"vel-{i:09d}.001-001"
        try:
            y = np.genfromtxt(archivo, dtype=np.float64, delimiter=[23, 23, 23])

            # Calculate mean fluid velocity
            vfx = np.mean(y[:,0])
            vfy = np.mean(y[:,1])
            vfz = np.mean(y[:,2])

            # For fluid-only simulations, particle velocity is same as fluid
            vx = vfx
            vy = vfy
            vz = vfz

            # Relative velocity is zero
            vxr = 0.0
            vyr = 0.0
            vzr = 0.0

            # Write data
            fo.write(f"{i},{vx:e},{vy:e},{vz:e},{vfx:e},{vfy:e},{vfz:e},{vxr:e},{vyr:e},{vzr:e}\n")
        except Exception as e:
            print(f"Warning: Could not process step {i}: {e}")
            continue

print(f"Processed data saved to {fout}")
