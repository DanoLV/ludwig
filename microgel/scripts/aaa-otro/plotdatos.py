#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


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
argParser = argparse.ArgumentParser(description='Plot particle position or velocity data from CSV file')
argParser.add_argument("-i", help="input file name", required=True)
argParser.add_argument("-t", "--tipo",
                      choices=['posicion', 'velocidad'],
                      default='posicion',
                      help="type of data to plot: 'posicion' or 'velocidad' (default: posicion)")
argParser.add_argument("-o", "--output",
                      default=None,
                      help="output file name (default: auto-generated based on type)")

# Parse command line arguments
try:
    args = argParser.parse_args()
except:
    sys.exit("Could not read command line parameters")

fdata = args.i
if(fdata is None):
    sys.exit("Please specify an input data file name")

tipo = args.tipo

# Set output filename
if args.output:
    output_file = args.output
else:
    output_file = f"{tipo}.png"

# Plot definitions
plt.rc('axes', labelsize=22)
plt.rc('ytick', labelsize=18)
plt.rc('xtick', labelsize=18)
plt.rc('legend', fontsize=16)

#------------------------------------------
# Read data file
try:
    y1 = np.genfromtxt(fdata, delimiter=',', skip_header=1)
except:
    sys.exit(f"Error reading data file: {fdata}")

#------------------------------------------
# Extract data columns
ciclo = y1[:,0]
x = y1[:,1]
y = y1[:,2]
z = y1[:,3]
vx = y1[:,4]
vy = y1[:,5]
vz = y1[:,6]

#-------------------------------------------------------------
# Plotting the Graph
#-------------------------------------------------------------

x_min_set = ciclo.min()
x_max_set = ciclo.max()

if tipo == 'posicion':
    #-------------------------------------------------------------
    # GRAFICAR POSICIÓN
    #-------------------------------------------------------------
    figure, axis = plt.subplots(1, figsize=(12, 6))

    axis.plot(ciclo, x, 'g', linewidth=2.0, label='X')
    axis.plot(ciclo, y, 'r', linewidth=2.0, label='Y')
    axis.plot(ciclo, z, 'b', linewidth=2.0, label='Z')

    limites_y(axis, x_min_set, x_max_set)

    axis.set_ylabel('Position', fontsize=22)
    axis.set_xlabel('Cycle', fontsize=22)
    axis.set_title('Particle Position vs Cycle', fontsize=24, fontweight='bold', pad=15)
    axis.legend(loc='best')
    axis.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    axis.tick_params(width=2)
    for spine in ['top','bottom','left','right']:
        axis.spines[spine].set_linewidth(2)

elif tipo == 'velocidad':
    #-------------------------------------------------------------
    # GRAFICAR VELOCIDAD
    #-------------------------------------------------------------
    figure, axis = plt.subplots(1, figsize=(12, 6))

    axis.plot(ciclo, vx, 'g', linewidth=2.0, label='Vx')
    axis.plot(ciclo, vy, 'r', linewidth=2.0, label='Vy')
    axis.plot(ciclo, vz, 'b', linewidth=2.0, label='Vz')

    limites_y(axis, x_min_set, x_max_set)

    axis.set_ylabel('Velocity', fontsize=22)
    axis.set_xlabel('Cycle', fontsize=22)
    axis.set_title('Particle Velocity vs Cycle', fontsize=24, fontweight='bold', pad=15)
    axis.legend(loc='best')
    axis.grid(True, which='both', linewidth=2, linestyle='--', alpha=0.7)
    axis.tick_params(width=2)
    for spine in ['top','bottom','left','right']:
        axis.spines[spine].set_linewidth(2)

plt.tight_layout()
plt.savefig(output_file, dpi=600, bbox_inches='tight')
print(f"Plot saved to {output_file}")
