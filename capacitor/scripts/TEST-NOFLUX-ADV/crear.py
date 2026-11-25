#!/usr/bin/env python3
"""
TEST-NOFLUX-ADV (2026-11-16, pedido del usuario): conservación de los iones con advección centred de NP,
LB térmico, paredes y electrolito, ANTES y DESPUÉS de agregar la condición de no-flujo (np_no_flux_boundary)
al driver de GPU de Nernst-Planck.

Motivo: en el driver de GPU, np_advective_fluxes_kernel no mira el mapa; un enlace fluido-pared lleva el
flujo W_p (u_i.c_p)/2 (n_i + n_pared)/2, que sale del nodo fluido y no entra en ninguno (la pared no se
actualiza). Se espera que el total de iones en el fluido NO se conserve antes y SÍ después.

Corridas (cada una con su Ludwig.exe copiado):
  <ver>/pared     caja 16x16x18 con paredes sólidas en z = 1 y z = 18 (porous_media_init wall_z), placas a psi = 0
  <ver>/sinpared  la misma sin paredes (periódica): control de regresión, antes y después deben ser idénticas
con <ver> = antes (exe 221f3229) o despues (exe nuevo). Sin partículas.

Uso:  python3 capacitor/scripts/TEST-NOFLUX-ADV/crear.py <ver> <exe> [caso ...]
Crea capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV/<ver>/{pared,sinpared}/ (no sobrescribe).
"""
import os
import shutil
import sys

ROOT = "/home/danolv/Sim/ludwig"
BASE = f"{ROOT}/capacitor/local/WALL-DIRICHLET/TEST-NOFLUX-ADV"

INPUT = """##############################################################################
#
#  TEST-NOFLUX-ADV/{ver}/{caso}  (crear.py, 2026-11-16)
#  Conservacion de iones con adveccion centred de NP + LB termico + electrolito.
#  {desc}
#  Exe: {exe_md5}. Comparar el total de iones en los nodos fluidos vs paso.
#
##############################################################################

N_start  0
N_cycles 2000            # corrida corta

##############################################################################
#  Sistema
##############################################################################

size {size}
grid 1_1_1
periodicity 1_1_1        # x, y periodicos; con paredes, z queda cerrado por el mapa

##############################################################################
#  Fluido (LB termico)
##############################################################################

fluid_rho0     0.8
viscosity      0.05
viscosity_bulk 0.05
lb_relaxation_scheme M10
lb_fluctuations on       # ruido termico: da velocidad normal a la pared en la primera capa de fluido
temperature    1.0e-5    # kT; u_rms ~ sqrt(kT/rho) = 3.5e-3
random_seed    8361235   # semilla del ruido del LB (igual antes y despues)

##############################################################################
#  Electrocinetica
##############################################################################

free_energy fe_electro
electric_e0 0.0_0.0_0.0
electrokinetics_init        uniform
electrokinetics_init_rho_el 1.11111111e-02   # electrolito 1:1 por especie (lambda_D = 3 con l_B = 0.398)
electrokinetics_z0 +1
electrokinetics_z1 -1
electrokinetics_eunit 1.0
electrokinetics_epsilon     2.0e4            # l_B = 1/(4 pi eps kT) = 0.398
# Difusividad: stencil 27 -> A0 = 6.137828; input d = 0.05 -> D_fisico = d*A0 = 0.3069
electrokinetics_d0 0.05
electrokinetics_d1 0.05
electrokinetics_rel_tol     1.0e-18
electrokinetics_abs_tol     1.0e-18
electrokinetics_multisteps  1
electrokinetics_solver_stencil 27
electrokinetics_solver_type petsc
electrokinetics_advection   centred          # el esquema bajo prueba (es el valor por defecto)
fe_force_method phi_gradmu_correction
{walls}
##############################################################################
#  Salida
##############################################################################

default_io_mode   mpiio
default_io_format ascii
freq_statistics 500
freq_psi_resid  500
freq_config     100      # qsi/psi/rho/vel cada 100 pasos: total de iones vs paso
vel_io_freq     100
config_at_end   yes
"""

WALLS = """
##############################################################################
#  Paredes solidas en z = 1 y z = 18 a psi = 0 (Asta et al. 2019)
##############################################################################

porous_media_init wall_z
electrokinetics_wall_axis z
electrokinetics_wall_potential 0.0_0.0
electrokinetics_electroneutral no   # exige paredes (sin paredes aborta)
"""

CASOS = {"pared": ("16_16_18", WALLS, "Paredes en z = 1 y z = 18 (fluido 2..17), placas a tierra."),
         "sinpared": ("16_16_16", "", "Sin paredes (periodica): control de regresion antes/despues.")}


def main():
    ver, exe = sys.argv[1], sys.argv[2]
    import hashlib
    md5 = hashlib.md5(open(exe, "rb").read()).hexdigest()[:8]
    casos = sys.argv[3:] or list(CASOS)
    for caso in casos:
        size, walls, desc = CASOS[caso]
        d = f"{BASE}/{ver}/{caso}"
        if os.path.exists(d):
            sys.exit(f"{d} ya existe: no se sobrescribe")
        txt = INPUT.format(ver=ver, caso=caso, desc=desc, exe_md5=md5, size=size, walls=walls)
        for l in txt.split("\n"):
            assert len(l.encode()) <= 120, l      # NKEY_LENGTH 128 (src/runtime.c)
        os.makedirs(f"{d}/logs")
        os.makedirs(f"{d}/proceced_data")
        open(f"{d}/input", "w").write(txt)
        shutil.copy2(exe, f"{d}/Ludwig.exe")
        shutil.copy2(f"{ROOT}/.petscrc", f"{d}/.petscrc")
        print("creada", d, "exe", md5)


if __name__ == "__main__":
    main()
