#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import Rbf


def limites_y(axis,x_min,x_max):
    
    axis.set_xlim(x_min,x_max)
    
    # -------------------------
    # Ajustar Y automáticamente pero solo para el rango de X actual
    # -------------------------
    # x_min, x_max = axis.get_xlim()

    ys_in_range = []
    for line in axis.get_lines():
        xd = np.asarray(line.get_xdata())
        yd = np.asarray(line.get_ydata())
        # Seleccionar solo puntos con x dentro del rango y valores finitos en y
        mask = (xd >= x_min) & (xd <= x_max) & np.isfinite(yd)
        if np.any(mask):
            ys_in_range.append(yd[mask])

    if ys_in_range:
        all_y = np.hstack(ys_in_range)
        y_min, y_max = all_y.min(), all_y.max()
        if y_min == y_max:
            # Si todos los valores son iguales, expandimos un poco para visualizar
            pad = abs(y_min)*0.01 if y_min != 0 else 1e-6
        else:
            pad = 0.05 * (y_max - y_min)   # 5% padding
        axis.set_ylim(y_min - pad, y_max + pad)
    else:
        # Si no hay datos dentro del rango, mantener límites actuales o poner un fallback:
        # axis.set_ylim(-1,1)
        pass
    
    return 0

# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", help="input file name")
argParser.add_argument("--out_dir", help="Output directory", default=".")
argParser.add_argument("--fluid-only", action="store_true", help="plot only fluid average velocity")

# Parse command line arguments
try:
    args = argParser.parse_args()
except:
    sys.exit("Could not read command line parameters")

# input file   
fdata = args.i
if(fdata is None):
    sys.exit("Please specify an input data file name")

# Output file    
fout= args.out_dir + "/" + "velocidades.png"

# Plot definitions
plt.rc('axes', labelsize=22)
plt.rc('ytick', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('legend', fontsize=16)

#------------------------------------------
# Read data file
y1 = np.genfromtxt(fdata,delimiter=',',skip_header=1)

#------------------------------------------
# Extract data columns
x1 = y1[:,0]
vx = y1[:,4]
vy = y1[:,5]
vz = y1[:,6]
vfx = y1[:,7]
vfy = y1[:,8]
vfz = y1[:,9]
vxr = y1[:,10]
vyr = y1[:,11]
vzr = y1[:,12]

#-------------------------------------------------------------
# Plotting the Graph
#-------------------------------------------------------------

# x_min_set = x1.min() 
x_min_set = x1.min() - (x1.max() - x1.min())/50
# x_min_set = 5000 # -10 1000 5000 7500 10000 12000 15000 20000
x_max_set = x1.max() 
# x_max_set = 5000

if args.fluid_only:
    #-------------------------------------------------------------
    # VELOCIDAD DEL FLUIDO SOLAMENTE - UN SOLO PANEL
    #-------------------------------------------------------------
    figure, axis = plt.subplots(1, figsize=(12, 6))

    axis.plot(x1, vfx, 'g', linewidth=2.0, label='Vfx')
    axis.plot(x1, vfy, 'r', linewidth=2.0, label='Vfy')
    axis.plot(x1, vfz, 'b', linewidth=2.0, label='Vfz')
    
    limites_y(axis,x_min_set,x_max_set)

    axis.set_ylabel('Fluid Velocity', fontsize=22)
    axis.set_xlabel('Cycle', fontsize=22)
    axis.legend(loc='best')
    axis.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    axis.tick_params(width=2)
    for spine in ['top','bottom','left','right']:
        axis.spines[spine].set_linewidth(2)
    
else:
    #-------------------------------------------------------------
    # TRES PANELES SEPARADOS
    #-------------------------------------------------------------
    # figure, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    # figure, (ax3, ax2) = plt.subplots(2, 1, figsize=(12, 14), sharex=True)
    figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 14), sharex=True)
    #-------------------------------------------------------------
    # Panel 1: VELOCIDAD DE LA PARTICULA
    #-------------------------------------------------------------
    ax1.plot(x1, vx, 'g', linewidth=2.0, label='Vx (particle)')
    ax1.plot(x1, vy, 'r', linewidth=2.0, label='Vy (particle)')
    ax1.plot(x1, vz, 'b', linewidth=2.0, label='Vz (particle)')

    limites_y(ax1,x_min_set,x_max_set)
    
    ax1.set_ylabel('Particle Velocity', fontsize=22)
    ax1.set_title('Particle Velocity', fontsize=24, fontweight='bold', pad=15)
    ax1.set_xlabel('Step', fontsize=22)
    ax1.legend(loc='best')
    ax1.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    ax1.tick_params(width=2)
    ax1.tick_params(axis='x', labelbottom=True)
    for spine in ['top','bottom','left','right']:
        ax1.spines[spine].set_linewidth(2)

    #-------------------------------------------------------------
    # Panel 2: VELOCIDAD DEL FLUIDO
    #-------------------------------------------------------------
    ax2.plot(x1, vfx, 'g', linewidth=2.0, label='Vfx (fluid)')
    ax2.plot(x1, vfy, 'r', linewidth=2.0, label='Vfy (fluid)')
    ax2.plot(x1, vfz, 'b', linewidth=2.0, label='Vfz (fluid)')

    limites_y(ax2,x_min_set,x_max_set)

    ax2.set_ylabel('Fluid Velocity', fontsize=22)
    ax2.set_xlabel('Step', fontsize=22)
    ax2.set_title('Fluid Velocity (Average)', fontsize=24, fontweight='bold', pad=15)
    ax2.legend(loc='best')
    ax2.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    ax2.tick_params(width=2)
    for spine in ['top','bottom','left','right']:
        ax2.spines[spine].set_linewidth(2)

    #-------------------------------------------------------------
    # Panel 3: VELOCIDAD RELATIVA
    #-------------------------------------------------------------
    # ax3.plot(x1, vxr, 'g', linewidth=2.0, label='Vxr (relative)')
    # ax3.plot(x1, vyr, 'r', linewidth=2.0, label='Vyr (relative)')
    # ax3.plot(x1, vzr, 'b', linewidth=2.0, label='Vzr (relative)')
    
    # limites_y(ax3,x_min_set,x_max_set)

    # ax3.set_ylabel('Relative Velocity', fontsize=22)
    # ax3.set_xlabel('Cycle', fontsize=22)
    # ax3.set_title('Relative Velocity (Particle - Fluid)', fontsize=24, fontweight='bold', pad=15)
    # ax3.legend(loc='best')
    # ax3.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    # ax3.tick_params(width=2)
    # for spine in ['top','bottom','left','right']:
    #     ax3.spines[spine].set_linewidth(2)

    #-------------------------------------------------------------
    # Panel 3: VELOCIDAD RELATIVA Y DE PARTICULA SUPERPUESTAS
    #-------------------------------------------------------------
    # ax3.plot(x1, vxr, 'g', linewidth=2.0, linestyle='dashed', label='Vxr (relative)')
    # ax3.plot(x1, vyr, 'r', linewidth=2.0, linestyle='dashed', label='Vyr (relative)')
    # ax3.plot(x1, vzr, 'b', linewidth=2.0, linestyle='dashed', label='Vzr (relative)')
    # ax3.plot(x1, vx, 'g', linewidth=2.0, label='Vx (particle)')
    # ax3.plot(x1, vy, 'r', linewidth=2.0, label='Vy (particle)')
    # ax3.plot(x1, vz, 'b', linewidth=2.0, label='Vz (particle)')
    
    # limites_y(ax3,x_min_set,x_max_set)

    # ax3.set_ylabel('Particle and Relative Velocity', fontsize=22)
    # ax3.set_xlabel('Cycle', fontsize=22)
    # ax3.set_title('Particle and Relative Velocity (Particle - Fluid)', fontsize=24, fontweight='bold', pad=15)
    # ax3.legend(loc='best')
    # ax3.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    # ax3.tick_params(width=2)
    # for spine in ['top','bottom','left','right']:
    #     ax3.spines[spine].set_linewidth(2)
        
    # Adjust spacing between subplots
    plt.tight_layout()

plt.savefig(fout, dpi=300, bbox_inches='tight')
print(f"Plot saved to velocidades.png")
