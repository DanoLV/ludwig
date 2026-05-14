#!/usr/bin/env python3
"""
Calcula el potencial químico por nodo y especie:
    μ_i(x) = kT * ln(ρ_i(x)) + z_i * e * ϕ(x)

donde:
  - ρ_i: densidad de la especie i (del archivo qsi, columnas 0 y 1)
  - ϕ:   potencial electrostático (del archivo psi, escalar por nodo)
  - z_i: carga de la especie i (+1 para cationes, -1 para aniones)
  - e=1, kT: parámetros de entrada

Salidas por paso:
  - Por pantalla: μ_i mínimo, máximo y delta por especie
  - mui-NNNNNNNNN.001-001 : mismo formato que qsi (N filas, 2 columnas)
  - muigrad-NNNNNNNNN.001-001 : gradiente de μ_i (N filas, 6 columnas)

Salidas globales (serie temporal):
  - mui_series.csv : estadísticas por paso (min, max, delta, grad_max por especie)
  - mui_delta_evolution.png : evolución temporal de Δμ_i por especie

Uso (un paso):
    python calc_mui.py -n 100 -L 32 -kt 0.00001 [-d directorio]

Uso (serie temporal):
    python calc_mui.py --n-start 100 --n-end 5600 --n-step 100 -L 32 -kt 0.00001

Formato de archivos de salida:
  mui: igual que qsi, N filas con 2 columnas (μ₊, μ₋)
  muigrad: N filas con 6 columnas (∇μ₊_x, ∇μ₊_y, ∇μ₊_z, ∇μ₋_x, ∇μ₋_y, ∇μ₋_z)

Cálculo opcional de theta (factor de Boltzmann ponderado por densidad):
    θ_k(x) = ρ_k(x) · exp(β·z_k·e·ϕ(x))
           = ρ_k(x) · exp(z_k · ϕ_Ludwig(x))   [ϕ_Ludwig ya en unidades e·ϕ/kT]

  - theta-NNNNNNNNN.001-001 : mismo formato que qsi (N filas, 2 columnas, θ₊ θ₋)

Activar con --theta.
"""

import numpy as np
import argparse
import sys
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Lectores de archivos
# ---------------------------------------------------------------------------

class QsiReader:
    """Lee archivos qsi de Ludwig: N filas, 2 columnas (ρ₊, ρ₋)"""

    def __init__(self, filename, grid_size, num_species=2):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size
        self.num_species = num_species

    def read(self):
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz
        if data.shape[0] != n_points:
            raise ValueError(
                f"El archivo qsi tiene {data.shape[0]} filas, "
                f"pero se esperaban {n_points}")
        if data.ndim == 1 or data.shape[1] != self.num_species:
            raise ValueError(
                f"Se esperaban {self.num_species} columnas en qsi, "
                f"encontradas: {data.shape[1] if data.ndim > 1 else 1}")
        # Organizar cada especie en malla 3D (orden x-major de Ludwig)
        species = []
        for i in range(self.num_species):
            species.append(data[:, i].reshape((self.nx, self.ny, self.nz)))
        return species  # lista de arrays 3D


class PsiReader:
    """Lee archivos psi de Ludwig: N valores escalares (un valor por nodo)"""

    def __init__(self, filename, grid_size):
        self.filename = filename
        self.nx, self.ny, self.nz = grid_size

    def read(self):
        data = np.loadtxt(self.filename)
        n_points = self.nx * self.ny * self.nz
        if data.ndim != 1 or len(data) != n_points:
            raise ValueError(
                f"El archivo psi tiene {data.size} valores, "
                f"pero se esperaban {n_points} (1 escalar por nodo)")
        return data.reshape((self.nx, self.ny, self.nz))


# ---------------------------------------------------------------------------
# Cálculo de μ_i
# ---------------------------------------------------------------------------

def calc_mui(rho_3d, psi_3d, z_i, kt, rho_min=1e-20):
    """
    μ_i = kT * ln(ρ_i) + kT * z_i * ϕ  =  kT * (ln(ρ_i) + z_i * ϕ)
    donde ϕ = e*ψ/kT (adimensional de Ludwig), así que z_i*e*ψ = kT * z_i * ϕ.
    Ambos términos en unidades de kT (energía).
    Los nodos con ρ_i <= 0 se marcan como NaN para evitar log(0).
    """
    rho_safe = np.where(rho_3d > rho_min, rho_3d, np.nan)
    mu = kt * (np.log(rho_safe) + z_i * psi_3d)
    return mu


# ---------------------------------------------------------------------------
# Cálculo de θ_k = ρ_k · exp(z_k · ϕ_Ludwig)
# ---------------------------------------------------------------------------

def calc_theta(rho_3d, psi_3d, z_k):
    """
    θ_k = ρ_k · exp(β·z_k·e·ϕ)
    ϕ_Ludwig ya está en unidades e·ϕ/kT, así que β·z_k·e·ϕ = z_k·ϕ_Ludwig.
    No hay riesgo de log(0): el exponencial siempre es finito.
    """
    return rho_3d * np.exp(z_k * psi_3d)


# ---------------------------------------------------------------------------
# Gradiente por diferencias finitas centradas (condiciones periódicas)
# ---------------------------------------------------------------------------

def gradient_periodic(field_3d):
    """
    Gradiente de un campo escalar 3D usando diferencias centradas con
    condiciones de contorno periódicas. Paso de malla = 1 (unidades de Ludwig).
    Devuelve (gx, gy, gz), cada uno array 3D de la misma forma que field_3d.
    """
    gx = (np.roll(field_3d, -1, axis=0) - np.roll(field_3d, 1, axis=0)) / 2.0
    gy = (np.roll(field_3d, -1, axis=1) - np.roll(field_3d, 1, axis=1)) / 2.0
    gz = (np.roll(field_3d, -1, axis=2) - np.roll(field_3d, 1, axis=2)) / 2.0
    return gx, gy, gz


# ---------------------------------------------------------------------------
# Escritura de archivos en formato qsi
# ---------------------------------------------------------------------------

def write_qsi_format(filename, arrays_3d):
    """
    Escribe en formato qsi: N filas, una columna por array.
    arrays_3d: lista de arrays 3D con la misma forma (nx, ny, nz).
    Orden x-major (igual que Ludwig).
    """
    nx, ny, nz = arrays_3d[0].shape
    n = nx * ny * nz
    cols = [a.reshape(n) for a in arrays_3d]
    data = np.column_stack(cols)
    np.savetxt(filename, data, fmt='%20.12e')
    print(f"  Guardado: {filename}")


# ---------------------------------------------------------------------------
# Procesado de un paso
# ---------------------------------------------------------------------------

def process_step(nstep, grid_size, kt, z0, z1, base_dir, rho_min, do_theta=False):
    """
    Lee qsi y psi para el paso dado, calcula μ_i y su gradiente,
    guarda los archivos de salida y devuelve un dict con estadísticas.
    Devuelve None si los archivos no existen.
    """
    nstep_str9 = f"{nstep:09d}"
    qsi_file = os.path.join(base_dir, f"qsi-{nstep_str9}.001-001")
    psi_file = os.path.join(base_dir, f"psi-{nstep_str9}.001-001")

    if not os.path.exists(qsi_file) or not os.path.exists(psi_file):
        print(f"  AVISO: paso {nstep} omitido (no se encontraron qsi o psi)")
        return None

    # Leer
    species = QsiReader(qsi_file, grid_size).read()
    rho0 = species[0]
    rho1 = species[1]
    psi_3d = PsiReader(psi_file, grid_size).read()

    # μ_i
    mu0 = calc_mui(rho0, psi_3d, z0, kt, rho_min=rho_min)
    mu1 = calc_mui(rho1, psi_3d, z1, kt, rho_min=rho_min)

    # Estadísticas de μ
    def mu_stats(arr):
        valid = arr[np.isfinite(arr)]
        if len(valid) == 0:
            return dict(min=np.nan, max=np.nan, delta=np.nan, mean=np.nan)
        mn, mx = float(valid.min()), float(valid.max())
        return dict(min=mn, max=mx, delta=mx - mn, mean=float(valid.mean()))

    s0 = mu_stats(mu0)
    s1 = mu_stats(mu1)

    # Guardar mui
    mui_file = os.path.join(base_dir, f"mui-{nstep_str9}.001-001")
    write_qsi_format(mui_file, [mu0, mu1])

    # Gradientes
    gx0, gy0, gz0 = gradient_periodic(mu0)
    gx1, gy1, gz1 = gradient_periodic(mu1)

    def grad_stats(gx, gy, gz):
        mag = np.sqrt(gx**2 + gy**2 + gz**2)
        valid = mag[np.isfinite(mag)]
        if len(valid) == 0:
            return dict(grad_min=np.nan, grad_max=np.nan, grad_mean=np.nan)
        return dict(grad_min=float(valid.min()),
                    grad_max=float(valid.max()),
                    grad_mean=float(valid.mean()))

    g0 = grad_stats(gx0, gy0, gz0)
    g1 = grad_stats(gx1, gy1, gz1)

    # Guardar muigrad
    muigrad_file = os.path.join(base_dir, f"muigrad-{nstep_str9}.001-001")
    write_qsi_format(muigrad_file, [gx0, gy0, gz0, gx1, gy1, gz1])

    # Imprimir resumen del paso
    print(f"  μ₊: min={s0['min']:+.6e}  max={s0['max']:+.6e}  Δ={s0['delta']:.6e}"
          f"  |∇μ₊|_max={g0['grad_max']:.6e}")
    print(f"  μ₋: min={s1['min']:+.6e}  max={s1['max']:+.6e}  Δ={s1['delta']:.6e}"
          f"  |∇μ₋|_max={g1['grad_max']:.6e}")

    row = dict(
        nstep=nstep,
        mu0_min=s0['min'],   mu0_max=s0['max'],   mu0_delta=s0['delta'],   mu0_mean=s0['mean'],
        mu0_grad_min=g0['grad_min'], mu0_grad_max=g0['grad_max'], mu0_grad_mean=g0['grad_mean'],
        mu1_min=s1['min'],   mu1_max=s1['max'],   mu1_delta=s1['delta'],   mu1_mean=s1['mean'],
        mu1_grad_min=g1['grad_min'], mu1_grad_max=g1['grad_max'], mu1_grad_mean=g1['grad_mean'],
    )

    # --- θ_k opcional ---
    if do_theta:
        th0 = calc_theta(rho0, psi_3d, z0)
        th1 = calc_theta(rho1, psi_3d, z1)

        def th_stats(arr, name):
            valid = arr[np.isfinite(arr)]
            if len(valid) == 0:
                mn = mx = delta = mean = np.nan
            else:
                mn, mx = float(valid.min()), float(valid.max())
                delta = mx - mn
                mean = float(valid.mean())
            print(f"  θ{name}: min={mn:.6e}  max={mx:.6e}  Δ={delta:.6e}  media={mean:.6e}")
            return dict(min=mn, max=mx, delta=delta, mean=mean)

        t0 = th_stats(th0, '₊')
        t1 = th_stats(th1, '₋')

        theta_file = os.path.join(base_dir, f"theta-{nstep_str9}.001-001")
        write_qsi_format(theta_file, [th0, th1])

        row.update(
            th0_min=t0['min'], th0_max=t0['max'], th0_delta=t0['delta'], th0_mean=t0['mean'],
            th1_min=t1['min'], th1_max=t1['max'], th1_delta=t1['delta'], th1_mean=t1['mean'],
        )

    return row


# ---------------------------------------------------------------------------
# CSV y gráfico de serie temporal
# ---------------------------------------------------------------------------

CSV_FIELDS_BASE = [
    'nstep',
    'mu0_min', 'mu0_max', 'mu0_delta', 'mu0_mean',
    'mu0_grad_min', 'mu0_grad_max', 'mu0_grad_mean',
    'mu1_min', 'mu1_max', 'mu1_delta', 'mu1_mean',
    'mu1_grad_min', 'mu1_grad_max', 'mu1_grad_mean',
]
CSV_FIELDS_THETA = [
    'th0_min', 'th0_max', 'th0_delta', 'th0_mean',
    'th1_min', 'th1_max', 'th1_delta', 'th1_mean',
]


def write_csv(rows, csv_path, do_theta=False):
    fields = CSV_FIELDS_BASE + (CSV_FIELDS_THETA if do_theta else [])
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: f"{row[k]:.10e}" if k != 'nstep' else row[k]
                             for k in fields})
    print(f"\nCSV guardado: {csv_path}")


def plot_delta_evolution(rows, plot_path, do_theta=False):
    steps = [r['nstep'] for r in rows]
    delta0 = [r['mu0_delta'] for r in rows]
    delta1 = [r['mu1_delta'] for r in rows]

    n_panels = 2 if do_theta else 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 5 * n_panels), squeeze=False)

    ax = axes[0, 0]
    ax.plot(steps, delta0, 'o-', color='tab:blue',  linewidth=1, markersize=1.5,
            label=r'$\Delta\mu_+$ (cationes, $z=+1$)')
    ax.plot(steps, delta1, 's-', color='tab:orange', linewidth=1, markersize=1.5,
            label=r'$\Delta\mu_-$ (aniones, $z=-1$)')
    ax.set_xlabel('Paso de tiempo', fontsize=13)
    ax.set_ylabel(r'$\Delta\mu_i = \mu_i^{max} - \mu_i^{min}$', fontsize=13)
    ax.set_title(r'Evolución temporal de $\Delta\mu_i$ por especie', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    if do_theta:
        dth0 = [r.get('th0_delta', np.nan) for r in rows]
        dth1 = [r.get('th1_delta', np.nan) for r in rows]
        ax2 = axes[1, 0]
        ax2.plot(steps, dth0, 'o-', color='tab:blue',  linewidth=1, markersize=1.5,
                 label=r'$\Delta\theta_+$ (cationes)')
        ax2.plot(steps, dth1, 's-', color='tab:orange', linewidth=1, markersize=1.5,
                 label=r'$\Delta\theta_-$ (aniones)')
        ax2.set_xlabel('Paso de tiempo', fontsize=13)
        ax2.set_ylabel(r'$\Delta\theta_k = \theta_k^{max} - \theta_k^{min}$', fontsize=13)
        ax2.set_title(r'Evolución temporal de $\Delta\theta_k$  '
                      r'[$\theta_k = \rho_k \cdot e^{z_k \phi}$]', fontsize=14)
        ax2.legend(fontsize=12)
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Gráfico guardado: {plot_path}")


# ---------------------------------------------------------------------------
# Gráfico --theta-only: θ₊ y θ₋ en ejes separados (serie temporal)
# ---------------------------------------------------------------------------

def plot_theta_only(rows, plot_path):
    """Dos paneles: θ₊(t) arriba, θ₋(t) abajo, con min/max/delta cada uno."""
    steps  = [r['nstep'] for r in rows]
    th0min = [r.get('th0_min',   np.nan) for r in rows]
    th0max = [r.get('th0_max',   np.nan) for r in rows]
    th0d   = [r.get('th0_delta', np.nan) for r in rows]
    th1min = [r.get('th1_min',   np.nan) for r in rows]
    th1max = [r.get('th1_max',   np.nan) for r in rows]
    th1d   = [r.get('th1_delta', np.nan) for r in rows]

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax0.plot(steps, th0max, '-', color='tab:blue',   linewidth=1, markersize=1.5,
             label=r'$\theta_+^{max}$')
    ax0.plot(steps, th0min, '-', color='tab:cyan',   linewidth=1, markersize=1.5,
             label=r'$\theta_+^{min}$')
    ax0.plot(steps, th0d,   '-', color='tab:green',  linewidth=1, markersize=1.5,
             label=r'$\Delta\theta_+$')
    ax0.set_ylabel(r'$\theta_+ = \rho_+ \cdot e^{z_+ \phi}$', fontsize=12)
    ax0.set_title(r'Evolución de $\theta_k$ por especie  '
                  r'[$\theta_k = \rho_k \cdot e^{z_k \phi}$]', fontsize=13)
    ax0.legend(fontsize=10)
    ax0.grid(True, alpha=0.3)

    ax1.plot(steps, th1max, '-', color='tab:orange', linewidth=1, markersize=1.5,
             label=r'$\theta_-^{max}$')
    ax1.plot(steps, th1min, '-', color='tab:red',    linewidth=1, markersize=1.5,
             label=r'$\theta_-^{min}$')
    ax1.plot(steps, th1d,   '-', color='tab:purple', linewidth=1, markersize=1.5,
             label=r'$\Delta\theta_-$')
    ax1.set_xlabel('Paso de tiempo', fontsize=12)
    ax1.set_ylabel(r'$\theta_- = \rho_- \cdot e^{z_- \phi}$', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Gráfico theta-only guardado: {plot_path}")


# ---------------------------------------------------------------------------
# Gráfico exp(βz_k eψ) vs ρ_k  (scatter por nodo, un paso)
# ---------------------------------------------------------------------------

def plot_boltzmann_vs_rho(rho0, rho1, psi_3d, z0, z1, nstep, base_dir, rho_min=1e-20):
    """
    Figura 3×2 con scatters por nodo para cada especie (izq=cationes, der=aniones):
      Fila 1: ρ_k  vs  exp(-β z_k e ϕ)        [1]
      Fila 2: ρ_k  vs  ϕ                       [2]
      Fila 3: ln(ρ_k) vs ϕ                     [3]
    Guarda 'boltzmann_vs_rho-NNNNNNNNN.png' en base_dir.
    """
    psi   = psi_3d.ravel()
    exp0  = np.exp(-z0 * psi)   # exp(-β z_k e ψ)
    exp1  = np.exp(-z1 * psi)
    r0    = rho0.ravel()
    r1    = rho1.ravel()

    mask0 = r0 > rho_min
    mask1 = r1 > rho_min

    KW = dict(s=3, alpha=0.3, linewidths=0)

    fig, axes = plt.subplots(3, 2, figsize=(12, 13))
    nstep_str9 = f"{nstep:09d}"
    fig.suptitle(f'Scatter por nodo — paso {nstep}', fontsize=14)

    # --- Fila 1: ρ_k vs exp(-z_k · ϕ) ---
    ax = axes[0, 0]
    ax.scatter(exp0[mask0], r0[mask0], color='tab:blue', **KW)
    ax.set_xlabel(r'$e^{-\beta z_+ e \phi}$', fontsize=12)
    ax.set_ylabel(r'$\rho_+$', fontsize=12)
    ax.set_title(rf'Cationes ($z_+={z0:+g}$)', fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.scatter(exp1[mask1], r1[mask1], color='tab:orange', **KW)
    ax.set_xlabel(r'$e^{-\beta z_- e \phi}$', fontsize=12)
    ax.set_ylabel(r'$\rho_-$', fontsize=12)
    ax.set_title(rf'Aniones ($z_-={z1:+g}$)', fontsize=11)
    ax.grid(True, alpha=0.3)

    # --- Fila 2: ρ_k vs ϕ ---
    ax = axes[1, 0]
    ax.scatter(psi[mask0], r0[mask0], color='tab:blue', **KW)
    ax.set_xlabel(r'$\phi$', fontsize=12)
    ax.set_ylabel(r'$\rho_+$', fontsize=12)
    ax.set_title(rf'Cationes ($z_+={z0:+g}$)', fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.scatter(psi[mask1], r1[mask1], color='tab:orange', **KW)
    ax.set_xlabel(r'$\phi$', fontsize=12)
    ax.set_ylabel(r'$\rho_-$', fontsize=12)
    ax.set_title(rf'Aniones ($z_-={z1:+g}$)', fontsize=11)
    ax.grid(True, alpha=0.3)

    # --- Fila 3: ln(ρ_k) vs ϕ  +  recta ln(θ) - β·z_k·e·ϕ ---
    def _ln_rho_panel(ax, psi_v, rho_v, mask, z_k, color, label_rho, label_z):
        ln_rho = np.log(rho_v[mask])
        phi_v  = psi_v[mask]

        ax.scatter(phi_v, ln_rho, color=color, **KW)

        # ln(θ) = Δln(ρ) / Δϕ  entre nodo de máx y mín densidad
        idx_max = np.argmax(rho_v[mask])
        idx_min = np.argmin(rho_v[mask])
        delta_ln_rho = ln_rho[idx_max] - ln_rho[idx_min]
        delta_phi    = phi_v[idx_max]  - phi_v[idx_min]
        if abs(delta_phi) > 1e-30:
            slope     = delta_ln_rho / delta_phi          # pendiente de la recta entre extremos
            intercept = ln_rho[idx_max] - slope * phi_v[idx_max]   # ln(θ): valor en ϕ=0
            ln_theta  = slope
        else:
            slope = intercept = ln_theta = np.nan

        phi_line    = np.linspace(phi_v.min(), phi_v.max(), 300)
        ln_rho_line = intercept + slope * phi_line

        ax.plot(phi_line, ln_rho_line, '-', color='gray', linewidth=0.5,
                label=rf'$\ln\rho = \ln\theta - z_k \phi$  ($\ln\theta={intercept:.4f}$)')

        ax.set_xlabel(r'$\phi$', fontsize=12)
        ax.set_ylabel(rf'$\ln(\rho_{{{label_rho}}})$', fontsize=12)
        ax.set_title(
            rf'{label_z}  —  $\ln\theta = {intercept:.4f}$  '
            rf'($\Delta\ln\rho/\Delta\phi = {ln_theta:.4f}$)',
            fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        return intercept, ln_theta

    ax = axes[2, 0]
    _ln_rho_panel(ax, psi, r0, mask0, z0, 'tab:blue',   '+', rf'Cationes ($z_+={z0:+g}$)')

    ax = axes[2, 1]
    _ln_rho_panel(ax, psi, r1, mask1, z1, 'tab:orange', '-', rf'Aniones ($z_-={z1:+g}$)')

    out = os.path.join(base_dir, f"boltzmann_vs_rho-{nstep_str9}.png")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfico Boltzmann vs ρ: {out}")


# ---------------------------------------------------------------------------
# Gráfico --ln-rho: ln(ρ_k) vs ϕ con recta de Boltzmann (figura independiente)
# ---------------------------------------------------------------------------

def plot_ln_rho(rho0, rho1, psi_3d, z0, z1, nstep, base_dir, rho_min=1e-20):
    """
    Figura 1×2: ln(ρ_k) vs ϕ para cada especie, con recta ln(θ) - z_k·ϕ superpuesta.
    Guarda 'ln_rho-NNNNNNNNN.png' en base_dir.
    """
    psi  = psi_3d.ravel()
    r0   = rho0.ravel()
    r1   = rho1.ravel()
    mask0 = r0 > rho_min
    mask1 = r1 > rho_min

    KW = dict(s=3, alpha=0.3, linewidths=0)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    nstep_str9 = f"{nstep:09d}"
    fig.suptitle(rf'$\ln(\rho_k)$ vs $\phi$ — paso {nstep}', fontsize=14)

    def _panel(ax, psi_v, rho_v, mask, z_k, color, sp_label, title_label):
        ln_rho = np.log(rho_v[mask])
        phi_v  = psi_v[mask]

        ax.scatter(phi_v, ln_rho, color=color, **KW, label=rf'$\ln(\rho_{{{sp_label}}})$')

        # ln(θ): intersección con ϕ=0 de la recta entre nodo de máx y mín densidad
        idx_max = np.argmax(rho_v[mask])
        idx_min = np.argmin(rho_v[mask])
        delta_phi = phi_v[idx_max] - phi_v[idx_min]
        if abs(delta_phi) > 1e-30:
            slope     = (ln_rho[idx_max] - ln_rho[idx_min]) / delta_phi
            intercept = ln_rho[idx_max] - slope * phi_v[idx_max]
        else:
            slope = intercept = np.nan

        phi_line = np.linspace(phi_v.min(), phi_v.max(), 300)
        ax.plot(phi_line, intercept + slope * phi_line, '-', color='gray', linewidth=0.5,
                label=rf'$\ln\rho = \ln\theta - z_k \phi$  ($\ln\theta={intercept:.4f}$)')

        ax.set_xlabel(r'$\phi$', fontsize=13)
        ax.set_ylabel(rf'$\ln(\rho_{{{sp_label}}})$', fontsize=13)
        ax.set_title(
            rf'{title_label}'
            rf'$\quad\ln\theta={intercept:.4f}$'
            rf'$\quad$pendiente$={slope:.4f}$',
            fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    _panel(axes[0], psi, r0, mask0, z0, 'tab:blue',   '+', rf'Cationes ($z_+={z0:+g}$) — ')
    _panel(axes[1], psi, r1, mask1, z1, 'tab:orange', '-', rf'Aniones ($z_-={z1:+g}$) — ')

    out = os.path.join(base_dir, f"ln_rho-{nstep_str9}.png")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfico ln(ρ) vs ϕ: {out}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Calcula potencial químico μ_i = kT*ln(ρ_i) + z_i*ϕ por nodo y especie',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  Un solo paso:
    python calc_mui.py -n 100 -L 32 -kt 0.00001

  Serie temporal:
    python calc_mui.py --n-start 100 --n-end 5600 --n-step 100 -L 32 -kt 0.00001
""")

    # Paso único (retrocompatible)
    parser.add_argument('-n', '--nstep', type=int, default=None,
                        help='Paso de tiempo único')
    # Serie temporal
    parser.add_argument('--n-start', type=int, default=None,
                        help='Paso inicial de la serie')
    parser.add_argument('--n-end', type=int, default=None,
                        help='Paso final de la serie')
    parser.add_argument('--n-step', type=int, default=None,
                        help='Incremento entre pasos de la serie')

    parser.add_argument('-L', '--grid-size', type=int, required=True,
                        help='Tamaño de la malla cúbica L×L×L')
    parser.add_argument('-kt', '--kt', type=float, required=True,
                        help='Energía térmica kT')
    parser.add_argument('-d', '--directory', default='.',
                        help='Directorio con los archivos de Ludwig (default: .)')
    parser.add_argument('--z0', type=float, default=1.0,
                        help='Carga de la especie 0 (cationes, default: +1)')
    parser.add_argument('--z1', type=float, default=-1.0,
                        help='Carga de la especie 1 (aniones, default: -1)')
    parser.add_argument('--rho-min', type=float, default=1e-20,
                        help='Densidad mínima para evitar log(0) (default: 1e-20)')
    parser.add_argument('--theta', action='store_true', default=False,
                        help='Calcular θ_k = ρ_k·exp(z_k·ϕ) y guardar theta-NNNNNNNNN.001-001')
    parser.add_argument('--theta-only', action='store_true', default=False,
                        help='Gráfico de θ_k por especie en ejes separados (implica --theta)')
    parser.add_argument('--boltzmann-plot', action='store_true', default=False,
                        help='Scatter de exp(βz_k eψ) vs ρ_k por nodo y especie')
    parser.add_argument('--ln-rho', action='store_true', default=False,
                        help='Scatter ln(ρ_k) vs ϕ con recta de Boltzmann (figura independiente)')
    parser.add_argument('--csv', default='mui_series.csv',
                        help='Nombre del archivo CSV de salida (default: mui_series.csv)')
    parser.add_argument('--plot', default='mui_delta_evolution.png',
                        help='Nombre del gráfico de salida (default: mui_delta_evolution.png)')

    args = parser.parse_args()

    # Determinar lista de pasos
    if args.nstep is not None:
        steps = [args.nstep]
    elif args.n_start is not None and args.n_end is not None and args.n_step is not None:
        steps = list(range(args.n_start, args.n_end + 1, args.n_step))
    else:
        print("ERROR: Indica -n PASO o bien --n-start, --n-end y --n-step")
        sys.exit(1)

    grid_size = (args.grid_size, args.grid_size, args.grid_size)
    base_dir = args.directory

    # --theta-only implica --theta
    if args.theta_only:
        args.theta = True

    print(f"Parámetros: L={args.grid_size}, kT={args.kt:.6e}, "
          f"z₊={args.z0}, z₋={args.z1}, rho_min={args.rho_min:.1e}"
          + ("  [--theta activado]" if args.theta else ""))
    print(f"Pasos a procesar: {steps[0]} → {steps[-1]} "
          f"({len(steps)} pasos)\n")

    rows = []
    for nstep in steps:
        print(f"Paso {nstep:9d}:")
        row = process_step(nstep, grid_size, args.kt, args.z0, args.z1,
                           base_dir, args.rho_min, do_theta=args.theta)
        if row is not None:
            rows.append(row)

        # Leer datos originales si algún gráfico de scatter los necesita
        if args.boltzmann_plot or args.ln_rho:
            nstep_str9 = f"{nstep:09d}"
            qsi_file = os.path.join(base_dir, f"qsi-{nstep_str9}.001-001")
            psi_file = os.path.join(base_dir, f"psi-{nstep_str9}.001-001")
            if os.path.exists(qsi_file) and os.path.exists(psi_file):
                sp  = QsiReader(qsi_file, grid_size).read()
                psi = PsiReader(psi_file, grid_size).read()
                if args.boltzmann_plot:
                    plot_boltzmann_vs_rho(sp[0], sp[1], psi,
                                          args.z0, args.z1, nstep,
                                          base_dir, rho_min=args.rho_min)
                if args.ln_rho:
                    plot_ln_rho(sp[0], sp[1], psi,
                                args.z0, args.z1, nstep,
                                base_dir, rho_min=args.rho_min)

    if not rows:
        print("No se procesó ningún paso.")
        sys.exit(1)

    # CSV y gráfico sólo si hay más de un paso o si se pidió explícitamente la serie
    if len(steps) > 1:
        csv_path  = os.path.join(base_dir, args.csv)
        plot_path = os.path.join(base_dir, args.plot)
        write_csv(rows, csv_path, do_theta=args.theta)
        if len(rows) > 1:
            plot_delta_evolution(rows, plot_path, do_theta=args.theta)
            if args.theta_only and args.theta:
                theta_plot = os.path.join(base_dir, 'theta_only_evolution.png')
                plot_theta_only(rows, theta_plot)
        else:
            print("Solo un paso procesado con éxito: se omite el gráfico.")
    else:
        # Un solo paso: mostrar estadísticas finales en pantalla
        r = rows[0]
        print(f"\nResumen paso {r['nstep']}:")
        print(f"  μ₊: min={r['mu0_min']:+.8e}  max={r['mu0_max']:+.8e}  "
              f"Δ={r['mu0_delta']:.8e}  media={r['mu0_mean']:+.8e}")
        print(f"  μ₋: min={r['mu1_min']:+.8e}  max={r['mu1_max']:+.8e}  "
              f"Δ={r['mu1_delta']:.8e}  media={r['mu1_mean']:+.8e}")
        print(f"  |∇μ₊|: min={r['mu0_grad_min']:.6e}  max={r['mu0_grad_max']:.6e}")
        print(f"  |∇μ₋|: min={r['mu1_grad_min']:.6e}  max={r['mu1_grad_max']:.6e}")
        if args.theta:
            print(f"  θ₊: min={r['th0_min']:.8e}  max={r['th0_max']:.8e}  Δ={r['th0_delta']:.8e}")
            print(f"  θ₋: min={r['th1_min']:.8e}  max={r['th1_max']:.8e}  Δ={r['th1_delta']:.8e}")

    print("\nListo.")


if __name__ == '__main__':
    main()
