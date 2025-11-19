#!/usr/bin/env python3
"""
Script para comparar el campo eléctrico simulado con el campo teórico
de una carga puntual.

Para una carga puntual q en el vacío (o medio uniforme), el campo eléctrico es:
E = (q / (4πε)) * (r / |r|³)

donde r es el vector desde la carga al punto de observación.

En Ludwig, el sistema está discretizado en una malla y se usa el sistema de
unidades de lattice Boltzmann.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os

class PsiReader:
    """Lee y procesa archivos psi de Ludwig"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.psi = None

    def read(self):
        """Lee el archivo psi y lo organiza en una malla 3D"""
        data = np.loadtxt(self.filename)

        if len(data) != self.nx * self.ny * self.nz:
            raise ValueError(f"El archivo tiene {len(data)} valores, "
                           f"pero se esperaban {self.nx * self.ny * self.nz}")

        self.psi = data.reshape((self.nx, self.ny, self.nz))
        return self.psi

    def compute_electric_field(self):
        """Calcula E = -∇ψ"""
        if self.psi is None:
            raise ValueError("Primero debe leer el archivo con read()")

        grad_psi = np.gradient(self.psi)
        Ex = -grad_psi[0]
        Ey = -grad_psi[1]
        Ez = -grad_psi[2]

        return Ex, Ey, Ez


def theoretical_field(charge_pos, grid_size, q=1.0, epsilon=1.0, kt=1.0, scale_factor=None,
                     kappa=0.0, ionic_strength=None):
    """
    Calcula el campo eléctrico teórico de una carga puntual

    Args:
        charge_pos: Posición de la carga (x, y, z)
        grid_size: Tupla (nx, ny, nz)
        q: Carga (default: 1.0)
        epsilon: Permitividad (default: 1.0)
        kt: Energía térmica k_B*T (default: 1.0)
        scale_factor: Factor de escala manual. Si es None, usa q/(4πε*kt)
        kappa: Parámetro de Debye (inverso de longitud de Debye). Si > 0, usa potencial screened
        ionic_strength: Fuerza iónica (mol/L). Si se proporciona, calcula kappa automáticamente

    Potenciales:
        - kappa = 0 (vacío): φ = q/(4πε*kt·r)  →  E = q/(4πε*kt) · r/r³
        - kappa > 0 (electrolito): φ = q/(4πε*kt·r) · exp(-κr)  →  E = q/(4πε*kt) · exp(-κr) · (1/r² + κ/r) · r̂

    Returns:
        Ex_theory, Ey_theory, Ez_theory: Componentes del campo teórico
    """
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos

    # Calcular kappa desde ionic strength si se proporciona
    if ionic_strength is not None and ionic_strength > 0:
        # Longitud de Debye: λ_D = sqrt(ε·k_B·T / (2·N_A·e²·I))
        # κ = 1/λ_D ≈ 3.29 · sqrt(I)  [en unidades de nm^-1 para I en mol/L]
        # En unidades de lattice, necesitas conocer la escala de longitud
        kappa = 3.29 * np.sqrt(ionic_strength)  # Asume unidades en nm
        print(f"Calculado κ = {kappa:.4f} nm^-1 desde ionic_strength = {ionic_strength} mol/L")

    # Crear malla de coordenadas centradas en nodos (0.5, 1.5, 2.5, ...)
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Vector distancia desde la carga
    dx = X - xc
    dy = Y - yc
    dz = Z - zc

    # Distancia
    r_squared = dx**2 + dy**2 + dz**2
    r = np.sqrt(r_squared)

    # Evitar división por cero en la posición de la carga
    r = np.where(r < 0.01, 0.01, r)
    r_squared = r**2

    # Prefactor
    if scale_factor is None:
        prefactor = q / (4.0 * np.pi * epsilon * kt)
    else:
        prefactor = scale_factor

    # Campo eléctrico según el tipo de medio
    if kappa > 0:
        # Medio con iones: Potencial de Yukawa/Debye-Hückel
        # E = (q/(4πε)) · exp(-κr) · (κ/r + 1/r²) · r̂
        exp_factor = np.exp(-kappa * r)
        magnitude = prefactor * exp_factor * (kappa / r + 1.0 / r_squared)

        Ex_theory = magnitude * dx / r
        Ey_theory = magnitude * dy / r
        Ez_theory = magnitude * dz / r
    else:
        # Vacío: Coulomb simple
        # E = (q/(4πε)) · (r / r³)
        r_cubed = r_squared * r

        Ex_theory = prefactor * dx / r_cubed
        Ey_theory = prefactor * dy / r_cubed
        Ez_theory = prefactor * dz / r_cubed

    return Ex_theory, Ey_theory, Ez_theory


def fit_scale_factor(Ex_sim, Ey_sim, Ez_sim, charge_pos, grid_size, q=1.0, kappa=0.0, kt=1.0):
    """
    Encuentra el factor de escala óptimo comparando el campo simulado con el patrón teórico

    Args:
        Ex_sim, Ey_sim, Ez_sim: Componentes del campo simulado
        charge_pos: Posición de la carga
        grid_size: Tupla (nx, ny, nz)
        q: Carga
        kappa: Parámetro de Debye para screening iónico
        kt: Energía térmica k_B*T

    Returns:
        scale_factor: Factor de escala óptimo
    """
    # Calcular campo teórico con factor unitario
    Ex_unit, Ey_unit, Ez_unit = theoretical_field(charge_pos, grid_size, q=1.0, epsilon=1.0, kt=1.0,
                                                   scale_factor=1.0, kappa=kappa)

    # Verificar que las dimensiones coincidan
    if Ex_sim.shape != Ex_unit.shape:
        raise ValueError(f"Dimensiones incompatibles: Ex_sim={Ex_sim.shape}, Ex_unit={Ex_unit.shape}")

    # Crear máscara para excluir región cercana a la carga (donde la discretización afecta)
    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    r_squared = (X - xc)**2 + (Y - yc)**2 + (Z - zc)**2
    mask = r_squared > 4.0  # Excluir región dentro de 2 unidades de lattice

    # Calcular factor de escala mediante ajuste de mínimos cuadrados
    # scale_factor = sum(E_sim * E_unit) / sum(E_unit^2)
    numerator = (np.sum(Ex_sim[mask] * Ex_unit[mask]) +
                 np.sum(Ey_sim[mask] * Ey_unit[mask]) +
                 np.sum(Ez_sim[mask] * Ez_unit[mask]))

    denominator = (np.sum(Ex_unit[mask]**2) +
                   np.sum(Ey_unit[mask]**2) +
                   np.sum(Ez_unit[mask]**2))

    scale_factor = numerator / denominator

    return scale_factor


def compare_fields_line(Ex_sim, Ey_sim, Ez_sim, Ex_theory, Ey_theory, Ez_theory,
                        start, end, grid_size, num_points=100, output=None,
                        charge_pos=None, scale_factor=None, show_nodes=False, kappa=0.0, epsilon=1.0, kt=1.0):
    """
    Compara los campos simulado y teórico a lo largo de una línea

    Si charge_pos y scale_factor se proporcionan, el campo teórico se calcula
    directamente en los puntos de la línea sin interpolación para mejor resolución.

    Args:
        show_nodes: Si es True, muestra los nodos de Ludwig con marcadores y no interpola
        kappa: Parámetro de Debye para screening iónico (0 = vacío/Coulomb simple)
        epsilon: Permitividad del medio
        kt: Energía térmica k_B*T
    """
    from scipy.interpolate import RegularGridInterpolator

    nx, ny, nz = grid_size

    # Coordenadas de la malla
    x = np.arange(0.5, nx)
    y = np.arange(0.5, ny)
    z = np.arange(0.5, nz)

    # Crear interpoladores para el campo simulado
    interp_Ex_sim = RegularGridInterpolator((x, y, z), Ex_sim, bounds_error=False, fill_value=0)
    interp_Ey_sim = RegularGridInterpolator((x, y, z), Ey_sim, bounds_error=False, fill_value=0)
    interp_Ez_sim = RegularGridInterpolator((x, y, z), Ez_sim, bounds_error=False, fill_value=0)

    if show_nodes:
        # Modo show_nodes:
        # - Campo simulado: solo nodos de Ludwig con marcadores
        # - Campo teórico: alta resolución con línea continua

        # Encontrar nodos de Ludwig que están en la línea
        nodes_on_line = []
        distances_on_line = []

        #------------------------------------------------------------------------
        # for i in range(len(x)):
        #     node_pos = np.array([x[i], start[1], start[2]])
        #     # Verificar si el nodo está en la línea (mismo y, z)
        #     if np.abs(node_pos[1] - start[1]) < 0.01 and np.abs(node_pos[2] - start[2]) < 0.01:
        #         if start[0] <= node_pos[0] <= end[0]:
        #             nodes_on_line.append(i)
        #             distances_on_line.append(node_pos[0] - start[0])
        
        # Umbral de distancia para considerar que un nodo está "sobre la línea"
        tolerance = 0.01
        start = np.array(start)
        end = np.array(end)
        # Vector que define la dirección y longitud del segmento de línea
        line_vec = end - start
        line_len_sq = np.dot(line_vec, line_vec) # Longitud al cuadrado para evitar raíces cuadradas

        # Iterar a través de todos los nodos de la red (en las 3 dimensiones)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    # Posición del nodo actual
                    node_pos = np.array([x[i], y[j], z[k]])

                    # 1. Verificar si el nodo está dentro de los límites del segmento de línea
                    # Proyectamos el vector del inicio al nodo sobre el vector de la línea
                    vec_start_to_node = node_pos - start
                    # La proyección 't' nos dice qué tan lejos está el punto más cercano en la línea infinita
                    # t = 0 corresponde a 'start', t = 1 corresponde a 'end'
                    t = np.dot(vec_start_to_node, line_vec) / line_len_sq

                    # Si t está entre 0 y 1, la proyección del nodo cae dentro del segmento
                    if 0 <= t <= 1:
                        # 2. Calcular la distancia perpendicular del nodo al segmento
                        # El punto más cercano en la línea al nodo es start + t * line_vec
                        closest_point_on_line = start + t * line_vec
                        dist_sq = np.sum((node_pos - closest_point_on_line)**2) # Distancia al cuadrado

                        # 3. Si la distancia es menor que la tolerancia, el nodo está "en la línea"
                        if dist_sq < tolerance**2:
                            # Guardamos la posición del nodo
                            nodes_on_line.append(node_pos)
                            
                            # Calculamos la distancia a lo largo de la línea desde el punto 'start'
                            distance_along_line = np.sqrt(np.sum((closest_point_on_line - start)**2))
                            distances_on_line.append(distance_along_line)
                    
        #------------------------------------------------------------------------
        
        
        # # Extraer valores del campo simulado solo en los nodos
        # points_nodes = np.array([[x[i], start[1], start[2]] for i in nodes_on_line])
        # distance_nodes = np.array(distances_on_line)
        
        # VERIFICACIÓN DE SEGURIDAD: ¿Se encontraron nodos?
        if not nodes_on_line:
            print("ADVERTENCIA: No se encontraron nodos en la línea con la tolerancia actual.")
            print("Prueba a aumentar el valor de 'tolerance'. Saliendo de la función.")
            # Aquí puedes decidir qué hacer: salir, devolver un valor vacío, etc.
            # Por ejemplo, si esto está en una función, podrías hacer:
            return # Sale de la función compare_fields_line

        # Si la lista no está vacía
        # Convertir las listas de resultados a arrays de NumPy
        # 'nodes_on_line' ya contiene las coordenadas completas de los nodos encontrados.
        points_nodes = np.array(nodes_on_line)
        distance_nodes = np.array(distances_on_line)

        Ex_sim_line = interp_Ex_sim(points_nodes)
        Ey_sim_line = interp_Ey_sim(points_nodes)
        Ez_sim_line = interp_Ez_sim(points_nodes)

        # Para el campo teórico, usar alta resolución
        points_theory = np.array([np.linspace(start[i], end[i], num_points) for i in range(3)]).T
        distance_theory = np.sqrt(np.sum((points_theory - start)**2, axis=1))
    else:
        # Modo normal: ambos campos con interpolación
        # Crear puntos a lo largo de la línea
        points_nodes = np.array([np.linspace(start[i], end[i], num_points) for i in range(3)]).T
        distance_nodes = np.sqrt(np.sum((points_nodes - start)**2, axis=1))

        # Interpolar campo simulado
        Ex_sim_line = interp_Ex_sim(points_nodes)
        Ey_sim_line = interp_Ey_sim(points_nodes)
        Ez_sim_line = interp_Ez_sim(points_nodes)

        # Campo teórico usa los mismos puntos
        points_theory = points_nodes
        distance_theory = distance_nodes

    # Calcular campo teórico: directamente en los puntos o interpolado
    if charge_pos is not None and scale_factor is not None:
        # Calcular campo teórico directamente en cada punto de la línea (sin interpolación)
        xc, yc, zc = charge_pos

        # Vector distancia desde la carga a cada punto
        dx = points_theory[:, 0] - xc
        dy = points_theory[:, 1] - yc
        dz = points_theory[:, 2] - zc

        # Distancia
        r_squared = dx**2 + dy**2 + dz**2
        r = np.sqrt(r_squared)

        # Evitar división por cero en la posición de la carga
        r = np.where(r < 0.01, 0.01, r)
        r_squared = r**2

        # Campo eléctrico según el tipo de medio
        if kappa > 0:
            # Medio con iones: Potencial de Yukawa/Debye-Hückel
            exp_factor = np.exp(-kappa * r)
            magnitude = scale_factor * exp_factor * (kappa / r + 1.0 / r_squared)

            Ex_theory_line = magnitude * dx / r
            Ey_theory_line = magnitude * dy / r
            Ez_theory_line = magnitude * dz / r
        else:
            # Vacío: Coulomb simple
            r_cubed = r_squared * r

            Ex_theory_line = scale_factor * dx / r_cubed
            Ey_theory_line = scale_factor * dy / r_cubed
            Ez_theory_line = scale_factor * dz / r_cubed
    else:
        # Interpolar el campo teórico de la malla
        interp_Ex_theory = RegularGridInterpolator((x, y, z), Ex_theory, bounds_error=False, fill_value=0)
        interp_Ey_theory = RegularGridInterpolator((x, y, z), Ey_theory, bounds_error=False, fill_value=0)
        interp_Ez_theory = RegularGridInterpolator((x, y, z), Ez_theory, bounds_error=False, fill_value=0)

        Ex_theory_line = interp_Ex_theory(points_theory)
        Ey_theory_line = interp_Ey_theory(points_theory)
        Ez_theory_line = interp_Ez_theory(points_theory)

    # Magnitudes
    E_sim_mag = np.sqrt(Ex_sim_line**2 + Ey_sim_line**2 + Ez_sim_line**2)
    E_theory_mag = np.sqrt(Ex_theory_line**2 + Ey_theory_line**2 + Ez_theory_line**2)

    # Crear gráfico
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Estilo de plot según modo
    if show_nodes:
        # Modo nodos: solo marcadores, sin líneas
        sim_style = 'bo'  # puntos azules
        sim_label = 'Simulado (nodos)'
        sim_linewidth = 0
        sim_markersize = 6
    else:
        # Modo interpolado: líneas continuas con marcadores opcionales
        sim_style = 'b-'
        sim_label = 'Simulado'
        sim_linewidth = 2
        sim_markersize = 0

    # Componente X
    axes[0, 0].plot(distance_nodes, Ex_sim_line, sim_style, label=sim_label,
                    linewidth=sim_linewidth, markersize=sim_markersize)
    axes[0, 0].plot(distance_theory, Ex_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[0, 0].set_xlabel('Distancia', fontsize=12)
    axes[0, 0].set_ylabel('Ex', fontsize=12)
    axes[0, 0].set_title('Componente X del Campo Eléctrico', fontsize=14)
    axes[0, 0].legend(fontsize=10)
    axes[0, 0].grid(True, alpha=0.3)

    # Componente Y
    axes[0, 1].plot(distance_nodes, Ey_sim_line, sim_style, label=sim_label,
                    linewidth=sim_linewidth, markersize=sim_markersize)
    axes[0, 1].plot(distance_theory, Ey_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[0, 1].set_xlabel('Distancia', fontsize=12)
    axes[0, 1].set_ylabel('Ey', fontsize=12)
    axes[0, 1].set_title('Componente Y del Campo Eléctrico', fontsize=14)
    axes[0, 1].legend(fontsize=10)
    axes[0, 1].grid(True, alpha=0.3)

    # Componente Z
    axes[1, 0].plot(distance_nodes, Ez_sim_line, sim_style, label=sim_label,
                    linewidth=sim_linewidth, markersize=sim_markersize)
    axes[1, 0].plot(distance_theory, Ez_theory_line, 'r--', label='Teórico', linewidth=2)
    axes[1, 0].set_xlabel('Distancia', fontsize=12)
    axes[1, 0].set_ylabel('Ez', fontsize=12)
    axes[1, 0].set_title('Componente Z del Campo Eléctrico', fontsize=14)
    axes[1, 0].legend(fontsize=10)
    axes[1, 0].grid(True, alpha=0.3)

    # Magnitud
    axes[1, 1].plot(distance_nodes, E_sim_mag, sim_style, label=sim_label,
                    linewidth=sim_linewidth, markersize=sim_markersize)
    axes[1, 1].plot(distance_theory, E_theory_mag, 'r--', label='Teórico', linewidth=2)
    axes[1, 1].set_xlabel('Distancia', fontsize=12)
    axes[1, 1].set_ylabel('|E|', fontsize=12)
    axes[1, 1].set_title('Magnitud del Campo Eléctrico', fontsize=14)
    axes[1, 1].legend(fontsize=10)
    axes[1, 1].grid(True, alpha=0.3)

    # Ajustar límites del eje Y basados en los valores de simulación para cada subplot
    # Componente X
    y_min, y_max = np.min(Ex_sim_line), np.max(Ex_sim_line)
    y_range = y_max - y_min
    axes[0, 0].set_ylim(y_min - 0.05*y_range, y_max + 0.05*y_range)

    # Componente Y
    y_min, y_max = np.min(Ey_sim_line), np.max(Ey_sim_line)
    y_range = y_max - y_min
    axes[0, 1].set_ylim(y_min - 0.05*y_range, y_max + 0.05*y_range)

    # Componente Z
    y_min, y_max = np.min(Ez_sim_line), np.max(Ez_sim_line)
    y_range = y_max - y_min
    axes[1, 0].set_ylim(y_min - 0.05*y_range, y_max + 0.05*y_range)

    # Magnitud
    y_min, y_max = np.min(E_sim_mag), np.max(E_sim_mag)
    y_range = y_max - y_min
    axes[1, 1].set_ylim(y_min - 0.05*y_range, y_max + 0.05*y_range)

    nx, ny, nz = grid_size
    xc, yc, zc = charge_pos

    if (kappa==0.0):
        l_D="inf"
        plt.suptitle(f'Comparación Campo Simulado vs Teórico\n'
                     f'l_Debye={l_D}, K={kappa:.6f}, epsilon={epsilon:.6f}, kt={kt:.6f}, L=({nx:d}, {ny:d}, {nz:d})\n'
                     f'Carga puntual en: ({xc:.3f},{yc:.3f},{zc:.3f})'
                     f'  //  '
                     f'Línea desde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                     f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})',
                     fontsize=16, y=0.995)
    else:
        l_D=1/kappa
        plt.suptitle(f'Comparación Campo Simulado vs Teórico\n'
                    f'l_Debye={l_D:.6f}, K={kappa:.6f}, epsilon={epsilon:.6f}, kt={kt:.6f}, L=({nx:d}, {ny:d}, {nz:d})\n'
                    f'Carga puntual en: ({xc:.3f},{yc:.3f},{zc:.3f})'
                    f'  //  '
                    f'Línea desde ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) '
                    f'hasta ({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})',
                    fontsize=16, y=0.995)

    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en: {output}")
    else:
        plt.show()

    plt.close()

    # Calcular error relativo promedio
    # Si show_nodes está activo, necesitamos interpolar el campo teórico en los nodos para calcular el error
    if show_nodes and charge_pos is not None and scale_factor is not None:
        # Calcular campo teórico en los nodos para comparación
        xc, yc, zc = charge_pos
        dx_nodes = points_nodes[:, 0] - xc
        dy_nodes = points_nodes[:, 1] - yc
        dz_nodes = points_nodes[:, 2] - zc
        r_squared_nodes = dx_nodes**2 + dy_nodes**2 + dz_nodes**2
        r_nodes = np.sqrt(r_squared_nodes)
        r_nodes = np.where(r_nodes < 0.01, 0.01, r_nodes)
        r_squared_nodes = r_nodes**2

        if kappa > 0:
            # Medio con iones
            exp_factor_nodes = np.exp(-kappa * r_nodes)
            magnitude_nodes = scale_factor * exp_factor_nodes * (kappa / r_nodes + 1.0 / r_squared_nodes)
            Ex_theory_at_nodes = magnitude_nodes * dx_nodes / r_nodes
            Ey_theory_at_nodes = magnitude_nodes * dy_nodes / r_nodes
            Ez_theory_at_nodes = magnitude_nodes * dz_nodes / r_nodes
        else:
            # Vacío
            r_cubed_nodes = r_squared_nodes * r_nodes
            Ex_theory_at_nodes = scale_factor * dx_nodes / r_cubed_nodes
            Ey_theory_at_nodes = scale_factor * dy_nodes / r_cubed_nodes
            Ez_theory_at_nodes = scale_factor * dz_nodes / r_cubed_nodes

        E_theory_mag_at_nodes = np.sqrt(Ex_theory_at_nodes**2 + Ey_theory_at_nodes**2 + Ez_theory_at_nodes**2)

        error_x = np.mean(np.abs(Ex_sim_line - Ex_theory_at_nodes) / (np.abs(Ex_theory_at_nodes) + 1e-10))
        error_y = np.mean(np.abs(Ey_sim_line - Ey_theory_at_nodes) / (np.abs(Ey_theory_at_nodes) + 1e-10))
        error_z = np.mean(np.abs(Ez_sim_line - Ez_theory_at_nodes) / (np.abs(Ez_theory_at_nodes) + 1e-10))
        error_mag = np.mean(np.abs(E_sim_mag - E_theory_mag_at_nodes) / (E_theory_mag_at_nodes + 1e-10))
    else:
        # Modo normal: ambos en los mismos puntos
        error_x = np.mean(np.abs(Ex_sim_line - Ex_theory_line) / (np.abs(Ex_theory_line) + 1e-10))
        error_y = np.mean(np.abs(Ey_sim_line - Ey_theory_line) / (np.abs(Ey_theory_line) + 1e-10))
        error_z = np.mean(np.abs(Ez_sim_line - Ez_theory_line) / (np.abs(Ez_theory_line) + 1e-10))
        error_mag = np.mean(np.abs(E_sim_mag - E_theory_mag) / (E_theory_mag + 1e-10))

    print("\n" + "="*60)
    print("ANÁLISIS DE ERROR RELATIVO PROMEDIO")
    print("="*60)
    print(f"Error relativo Ex:  {error_x:.4f} ({error_x*100:.2f}%)")
    print(f"Error relativo Ey:  {error_y:.4f} ({error_y*100:.2f}%)")
    print(f"Error relativo Ez:  {error_z:.4f} ({error_z*100:.2f}%)")
    print(f"Error relativo |E|: {error_mag:.4f} ({error_mag*100:.2f}%)")
    print("="*60 + "\n")


def compare_fields_plane(Ex_sim, Ey_sim, Ez_sim, Ex_theory, Ey_theory, Ez_theory,
                         plane='xy', position=None, grid_size=(32, 32, 32), output=None):
    """
    Compara los campos simulado y teórico en un plano
    """
    nx, ny, nz = grid_size
    plane = plane.lower()

    if position is None:
        if plane == 'xy':
            position = nz // 2
        elif plane == 'xz':
            position = ny // 2
        elif plane == 'yz':
            position = nx // 2

    # Extraer slices
    if plane == 'xy':
        E_sim_mag = np.sqrt(Ex_sim[:, :, position]**2 +
                           Ey_sim[:, :, position]**2 +
                           Ez_sim[:, :, position]**2)
        E_theory_mag = np.sqrt(Ex_theory[:, :, position]**2 +
                              Ey_theory[:, :, position]**2 +
                              Ez_theory[:, :, position]**2)
        x = np.arange(0.5, nx)
        y = np.arange(0.5, ny)
        xlabel, ylabel = 'X', 'Y'
    elif plane == 'xz':
        E_sim_mag = np.sqrt(Ex_sim[:, position, :]**2 +
                           Ey_sim[:, position, :]**2 +
                           Ez_sim[:, position, :]**2)
        E_theory_mag = np.sqrt(Ex_theory[:, position, :]**2 +
                              Ey_theory[:, position, :]**2 +
                              Ez_theory[:, position, :]**2)
        x = np.arange(0.5, nx)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'X', 'Z'
    elif plane == 'yz':
        E_sim_mag = np.sqrt(Ex_sim[position, :, :]**2 +
                           Ey_sim[position, :, :]**2 +
                           Ez_sim[position, :, :]**2)
        E_theory_mag = np.sqrt(Ex_theory[position, :, :]**2 +
                              Ey_theory[position, :, :]**2 +
                              Ez_theory[position, :, :]**2)
        x = np.arange(0.5, ny)
        y = np.arange(0.5, nz)
        xlabel, ylabel = 'Y', 'Z'

    X, Y = np.meshgrid(x, y, indexing='ij')

    # Crear figura con 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Campo simulado
    im1 = axes[0].contourf(X, Y, E_sim_mag, levels=20, cmap='viridis')
    axes[0].set_xlabel(xlabel, fontsize=12)
    axes[0].set_ylabel(ylabel, fontsize=12)
    axes[0].set_title('Campo Simulado |E|', fontsize=14)
    axes[0].set_aspect('equal')
    plt.colorbar(im1, ax=axes[0])

    # Campo teórico
    im2 = axes[1].contourf(X, Y, E_theory_mag, levels=20, cmap='viridis')
    axes[1].set_xlabel(xlabel, fontsize=12)
    axes[1].set_ylabel(ylabel, fontsize=12)
    axes[1].set_title('Campo Teórico |E|', fontsize=14)
    axes[1].set_aspect('equal')
    plt.colorbar(im2, ax=axes[1])

    # Diferencia
    diff = np.abs(E_sim_mag - E_theory_mag)
    im3 = axes[2].contourf(X, Y, diff, levels=20, cmap='RdYlBu_r')
    axes[2].set_xlabel(xlabel, fontsize=12)
    axes[2].set_ylabel(ylabel, fontsize=12)
    axes[2].set_title('Diferencia Absoluta |E_sim - E_theory|', fontsize=14)
    axes[2].set_aspect('equal')
    plt.colorbar(im3, ax=axes[2])

    plt.suptitle(f'Comparación en Plano {plane.upper()} (posición {position})', fontsize=16)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado en: {output}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Comparar campo eléctrico simulado con campo teórico de carga puntual',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Comparar en vacío (Coulomb simple) a lo largo del eje X:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

2. Comparar con screening iónico (Debye-Hückel):
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --kappa 0.1 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

3. Comparar usando fuerza iónica:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --ionic-strength 0.001 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

4. Mostrar solo nodos de Ludwig con marcadores:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --kappa 0.1 \\
       -m line --start 0.5 16 16 --end 31.5 16 16 --show-nodes

5. Comparar en un plano XY:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --kappa 0.1 \\
       -m plane -p xy --position 16

6. Alta resolución con permitividad y guardar:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --epsilon 4.0 --kappa 0.05 \\
       -m line --start 0.5 16 16 --end 31.5 16 16 \\
       --num-points 500 -o comparison.png

7. Restar un campo externo constante antes de comparar:
   ./compare_field_theory.py -f psi-000030000.001-001 -s 32 32 32 \\
       --charge-pos 16 16 16 --charge 1.0 --kappa 0.1 \\
       --external-field 0.1 0.0 0.0 \\
       -m line --start 0.5 16 16 --end 31.5 16 16

Notas:
  - kappa = 0 (default): Coulomb simple (vacío, sin screening)
  - kappa > 0: Debye-Hückel (electrolito con screening iónico)
  - La longitud de Debye es λ_D = 1/κ
  - Para electrolitos: κ ≈ 0.329/λ_D(nm) en unidades SI
  - --external-field: Resta un campo externo del campo simulado antes de comparar
        """
    )

    parser.add_argument('-f', '--file', required=True,
                       help='Archivo psi a leer')
    parser.add_argument('-s', '--size', nargs=3, type=int, required=True,
                       metavar=('NX', 'NY', 'NZ'),
                       help='Tamaño de la malla')
    parser.add_argument('--charge-pos', nargs=3, type=float, required=True,
                       metavar=('X', 'Y', 'Z'),
                       help='Posición de la carga puntual')
    parser.add_argument('--charge', type=float, default=1.0,
                       help='Magnitud de la carga (default: 1.0)')
    parser.add_argument('--epsilon', type=float, default=1.0,
                       help='Permitividad del medio (default: 1.0)')
    parser.add_argument('--kt', type=float, default=1.0,
                       help='Energía térmica k_B*T (default: 1.0)')
    parser.add_argument('--kappa', type=float, default=0.0,
                       help='Parámetro de Debye κ (inverso de longitud de Debye) para screening iónico. '
                            'kappa=0 usa Coulomb simple (vacío). kappa>0 usa Debye-Hückel (electrolito). '
                            'En unidades de lattice. (default: 0.0)')
    parser.add_argument('--ionic-strength', type=float, default=None,
                       help='Fuerza iónica en mol/L. Si se proporciona, calcula kappa automáticamente. '
                            'Requiere conocer la escala de longitud del sistema.')

    parser.add_argument('-m', '--mode', choices=['line', 'plane'], required=True,
                       help='Modo de comparación: line o plane')

    # Para modo línea
    parser.add_argument('--start', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto inicial de la línea (obligatorio para mode=line)')
    parser.add_argument('--end', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Punto final de la línea (obligatorio para mode=line)')
    parser.add_argument('--num-points', type=int, default=100,
                       help='Número de puntos a lo largo de la línea (default: 100)')
    parser.add_argument('--show-nodes', action='store_true',
                       help='Mostrar solo los nodos de Ludwig con marcadores (sin interpolación)')

    # Para modo plano
    parser.add_argument('-p', '--plane', choices=['xy', 'xz', 'yz'],
                       help='Plano a comparar (obligatorio para mode=plane)')
    parser.add_argument('--position', type=int,
                       help='Posición del plano (índice). Default: centro')

    parser.add_argument('-o', '--output',
                       help='Archivo de salida (default: mostrar en pantalla)')

    parser.add_argument('--external-field', nargs=3, type=float, metavar=('EX', 'EY', 'EZ'),
                       help='Componentes del campo externo constante a restar (Ex, Ey, Ez). '
                            'Ejemplo: --external-field 0.1 0.0 0.0')

    args = parser.parse_args()

    # Validar argumentos
    if args.mode == 'line' and (args.start is None or args.end is None):
        parser.error("mode=line requiere --start y --end")
    if args.mode == 'plane' and args.plane is None:
        parser.error("mode=plane requiere --plane")

    if not os.path.exists(args.file):
        print(f"Error: No se encuentra el archivo {args.file}")
        sys.exit(1)

    # Leer archivo psi
    print(f"Leyendo archivo {args.file}...")
    reader = PsiReader(args.file, tuple(args.size))
    psi = reader.read()
    print(f"Archivo leído. Malla: {args.size[0]}x{args.size[1]}x{args.size[2]}")

    # Calcular campo simulado
    print("Calculando campo eléctrico simulado E = -∇ψ...")
    Ex_sim, Ey_sim, Ez_sim = reader.compute_electric_field()

    # Restar campo externo si se proporcionó
    if args.external_field is not None:
        Ex_ext, Ey_ext, Ez_ext = args.external_field

        print(f"\nRestando campo externo constante:")
        print(f"  Ex_ext = {Ex_ext}")
        print(f"  Ey_ext = {Ey_ext}")
        print(f"  Ez_ext = {Ez_ext}")

        # Restar campo externo constante
        Ex_sim = Ex_sim - Ex_ext
        Ey_sim = Ey_sim - Ey_ext
        Ez_sim = Ez_sim - Ez_ext
        print("Campo externo restado exitosamente.")

    # Determinar kappa (parámetro de screening)
    kappa = args.kappa
    if args.ionic_strength is not None:
        print(f"\nUsando fuerza iónica I = {args.ionic_strength} mol/L")
    if kappa > 0 or args.ionic_strength is not None:
        print(f"Modelo: Debye-Hückel (electrolito con screening iónico)")
    else:
        print(f"Modelo: Coulomb simple (vacío, sin screening)")

    # Ajustar factor de escala automáticamente
    print("\nAjustando factor de escala entre simulación y teoría...")
    scale_factor = fit_scale_factor(Ex_sim, Ey_sim, Ez_sim, args.charge_pos,
                                    tuple(args.size), q=args.charge, kappa=kappa, kt=args.kt)
    print(f"Factor de escala encontrado: {scale_factor:.6e}")
    print(f"Factor teórico q/(4πε*kt):   {args.charge / (4.0 * np.pi * args.epsilon * args.kt):.6e}")
    print(f"Ratio simulado/teórico:      {scale_factor / (args.charge / (4.0 * np.pi * args.epsilon * args.kt)):.4f}")

    # Calcular campo teórico con factor de escala ajustado
    print(f"\nCalculando campo teórico para carga q={args.charge} en "
          f"({args.charge_pos[0]}, {args.charge_pos[1]}, {args.charge_pos[2]}) "
          f"con factor de escala ajustado...")
    if kappa > 0:
        print(f"Parámetro de Debye κ = {kappa:.4f} (longitud de Debye λ_D = {1.0/kappa:.4f})")
    Ex_theory, Ey_theory, Ez_theory = theoretical_field(
        args.charge_pos, tuple(args.size), q=args.charge, epsilon=args.epsilon, kt=args.kt,
        scale_factor=scale_factor, kappa=kappa, ionic_strength=args.ionic_strength
    )

    # Comparar según el modo
    if args.mode == 'line':
        print(f"\nComparando a lo largo de línea desde {args.start} hasta {args.end}...")
        if args.show_nodes:
            print("Modo: Solo nodos de Ludwig (sin interpolación)")
        else:
            print(f"Calculando campo teórico en {args.num_points} puntos sin interpolación...")
        compare_fields_line(
            Ex_sim, Ey_sim, Ez_sim,
            Ex_theory, Ey_theory, Ez_theory,
            args.start, args.end, tuple(args.size),
            num_points=args.num_points,
            output=args.output,
            charge_pos=args.charge_pos,
            scale_factor=scale_factor,
            show_nodes=args.show_nodes,
            kappa=kappa,
            epsilon=args.epsilon,
            kt=args.kt
        )
    else:  # mode == 'plane'
        print(f"\nComparando en plano {args.plane.upper()}...")
        compare_fields_plane(
            Ex_sim, Ey_sim, Ez_sim,
            Ex_theory, Ey_theory, Ez_theory,
            plane=args.plane,
            position=args.position,
            grid_size=tuple(args.size),
            output=args.output
        )

    print("\n¡Comparación completada!")


if __name__ == '__main__':
    main()
