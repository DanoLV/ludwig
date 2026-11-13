#!/usr/bin/env python3
"""Baterias C-electrolito-hann{6,8}-LD{4,6,8}: una particula subgrid q=+1 fija en
x=y=16, a distancia d de la pared inferior, entre paredes a tierra, CON sal 1:1
(lambda_D = 4, 6, 8). Una corrida por posicion hasta el estacionario ionico.
Se parte de C-sin-fondo-hann08/ (misma caja 32x32x34, L=32, mismo bloque de particula).
Uso: python3 crear_baterias_electrolito.py
"""
import os, re, shutil
import numpy as np
SRC = "C-sin-fondo-hann08"
blk = open(os.path.join(SRC, "bloque_particula.txt")).readlines()
base = open(os.path.join(SRC, "input.base")).read()
LB4 = 10.0                                   # 4 pi l_B = beta/epsilon
for hann in (6, 8):
    for ld in (4, 6, 8):
        rhos = 1.0 / (2 * LB4 * 2 * np.pi * 0 + 8 * np.pi * (LB4 / (4 * np.pi)) * ld ** 2)  # (8 pi l_B lambda_D^2)^-1
        bat = "C-electrolito-hann%d-LD%d" % (hann, ld)
        os.makedirs(bat, exist_ok=True)
        dmin = hann / 2.0 - 0.5
        ds = sorted(set([dmin, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0]))
        ds = [x for x in ds if x >= dmin]
        for d in ds:
            z = d + 1.5
            dd = os.path.join(bat, "z_%05.2f" % z)
            for sub in ("logs", "proceced_data", "colloid_data"):
                os.makedirs(os.path.join(dd, sub), exist_ok=True)
            b = list(blk)
            b[0] = "%24d\n" % 1
            b[34] = "    %.15e    %.15e    %.15e\n" % (16.0, 16.0, z)
            b[45] = "    %.15e\n" % 1.0
            b[46] = "    %.15e\n" % 0.0
            with open(os.path.join(dd, "config.cds.init.001-001"), "w") as fh:
                fh.write("%22d\n" % 1); fh.writelines(b)
            t = base
            t = re.sub(r"^N_cycles .*$", "N_cycles 5000", t, flags=re.M)
            t = re.sub(r"^subgrid_hann_order .*$", "subgrid_hann_order %d      # orden Hann" % hann, t, flags=re.M)
            t = re.sub(r"^electrokinetics_init_rho_el .*$", "electrokinetics_init_rho_el %.8e" % rhos, t, flags=re.M)
            t = re.sub(r"^freq_psi_resid .*$", "freq_psi_resid 500", t, flags=re.M)
            t = re.sub(r"^freq_statistics .*$", "freq_statistics 500", t, flags=re.M)
            t = re.sub(r"^freq_config .*$", "freq_config 5000", t, flags=re.M)
            t = re.sub(r"^colloid_io_freq .*$", "colloid_io_freq 500", t, flags=re.M)
            t = re.sub(r"^pforce_io_freq .*$", "pforce_io_freq 500", t, flags=re.M)
            t += ("\n# --- Electrolito 1:1 (esta bateria) ---------------------------------------\n"
                  "# lambda_D = %d -> rho_s = 1/(8 pi l_B lambda_D^2) = %.6e por especie\n"
                  "# hidrodinamica queda ON: subgrid_force_from_particles exige hydro (assert). Con geometria plana u=0.\n"
                  "electrokinetics_z0 +1\nelectrokinetics_z1 -1\n"
                  "electrokinetics_eunit 1.0\n"
                  "# stencil 27 -> A0 = 6.137828 (memoria nernst_planck_A0_stencil_factor).\n"
                  "# d0/d1 = 0.05 NO corregido por A0: D_fisica efectiva = 0.05*6.137828 = 0.3069.\n"
                  "# Solo importa el ESTACIONARIO (no depende de D); D grande lo alcanza en ~5 veces menos pasos.\n"
                  "electrokinetics_d0 0.05\nelectrokinetics_d1 0.05\n" % (ld, rhos))
            open(os.path.join(dd, "input"), "w").write(t)
            shutil.copy("D-paper-asta/.petscrc", os.path.join(dd, ".petscrc"))
            lnk = os.path.join(dd, "Ludwig.exe")
            if not os.path.exists(lnk):
                os.symlink("/home/danolv/Sim/ludwig/src/Ludwig.exe", lnk)
        shutil.copy(os.path.join(SRC, "correr.sh"), os.path.join(bat, "correr.sh"))
        print(bat, "rho_s = %.4e" % rhos, "d =", ds)
