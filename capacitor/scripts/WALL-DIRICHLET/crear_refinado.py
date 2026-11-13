#!/usr/bin/env python3
"""
Refinamiento de la red para C-electrolito-hann8-LD6/z_09.50 (carga a d = 8 de la placa inferior, placas a tierra).

Mismo sistema fisico con paso de red dx = 1/r (r = n/8), de modo que Hann-n tenga el mismo ancho fisico (8) que Hann-8:
  - longitudes x r: caja 32 -> 32r (x, y), canal L = 32 -> 32r (size_z = 32r + 2), carga en (16r, 16r, 1.5 + 8r)
  - epsilon -> epsilon/r : q/(eps kT) en unidades de red = 10 r, asi psi(x_fisico) es el mismo (l_B -> r l_B)
  - rho_el -> rho_el/r^3 : densidad numerica por nodo; con eps/r da lambda_D -> r lambda_D
  - potenciales de pared, q, kT: sin cambios (psi = e phi / kT es adimensional)
  - D (por paso) igual -> el tiempo de difusion en pasos crece como r^2: N = 5000 r^2 redondeado hacia arriba
Solo se compara el estacionario, que no depende de D ni de la viscosidad (particula fija, u ~ 0).

Uso: python3 crear_refinado.py [n ...]   (por defecto n = 10 12; p. ej. "python3 crear_refinado.py 16" crea
     C-electrolito-hann16-LD12-ref2)
     python3 crear_refinado.py gauss [n ...]   kernel gaussiano en lugar de Hann, con r = n/8 (n = 8 incluido):
     sigma = 1.44612 r (misma varianza que Hann-8), corte a 5.5 sigma (soporte 11 sigma = 15.907 fisico; la mitad,
     7.95, queda dentro del canal: la carga esta a 8 de la placa). Crea C-electrolito-gauss-LD{6r}-ref{r}.
     Control de suavidad para el orden ~1 cerca del borde del soporte de Hann (Hann es solo C^1).
     python3 crear_refinado.py [gauss] q=0.25 [n ...]   carga q0 = 0.25 en lugar de 1 (sufijo -q0.25 en el nombre):
     prueba de linealidad. Si el piso C(z) del ajuste residuo = C + a dx^2 es la no linealidad de Poisson-Boltzmann
     (la analitica es Debye-Huckel linealizado), C/q escala como q^2 y con q = 0.25 baja 16 veces.
     python3 crear_refinado.py peskin6 [q=...] [n ...]   kernel de Peskin de 6 puntos (Bao et al. 2016, C^3) estirado por
     el entero s = 2r = n/4 (subgrid_peskin6_stretch): soporte 6s nodos = 12 fisico, sigma = 0.84503 s nodos = 1.69006
     fisico. n = 4, 8, 16 -> s = 1, 2, 4 (dx = 2, 1, 1/2). Crea C-electrolito-peskin6s{s}-LD{6r}-ref{r}.
     Control de suavidad con ancho fisico fijo (C^3, con mas derivadas continuas que Hann).
     python3 crear_refinado.py p=4 [q=...] [n ...]   Hann generalizado w ~ cos^p (subgrid_hann_power p, C^(p-1)), mismo
     soporte fisico 8 que Hann-8 (n = 8r nodos): para p = 4 sigma = 8 sqrt(1/12 - 5/(8 pi^2)) = 1.13155 fisico.
     Crea C-electrolito-hann{n}p{p}-LD{6r}-ref{r}. Control de suavidad (C^3) con el mismo soporte y las mismas redes.
"""
import math, os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "C-electrolito-hann8-LD6", "z_09.50")
EXE = os.path.join(HERE, "../../../src/Ludwig.exe")
EPS0, RHO0, L0, D0, SIG8 = 1.0e4, 1.38888889e-03, 32, 8.0, 1.44612
GCUT = 5.5                                            # corte del gaussiano en sigmas
SIGP6 = math.sqrt(59.0 / 60 - math.sqrt(29.0) / 20)  # sigma del Peskin-6 sin estirar (2do momento K, Bao et al. 2016)

HEADER = """\
##############################################################################
#
#  Refinamiento x{r:g} de C-electrolito-hann8-LD6/z_09.50 (crear_refinado.py)
#  Mismo sistema fisico con dx = 1/{r:g}: {kern},
#  L = {L} nodos (= 32 fisicos), lambda_D = {ld:g} nodos (= 6), d = {d:g} nodos (= 8).
#  Una carga q = +1 fija, placas a tierra, electrolito 1:1, sistema cerrado.
#
##############################################################################
"""

def sub(txt, pat, new):
    out, k = re.subn(pat, new, txt, flags=re.M)
    assert k == 1, pat
    return out

def make(n, gauss=False, q=1.0, peskin6=False, hp=2):
    r = n / 8.0; L = int(round(L0 * r)); s6 = int(round(2 * r))
    assert not peskin6 or abs(2 * r - s6) < 1e-12 and s6 >= 1, "peskin6: s = 2r debe ser entero"
    name = ("C-electrolito-gauss-LD%g-ref%g" % (6 * r, r) if gauss else
            "C-electrolito-peskin6s%d-LD%g-ref%g" % (s6, 6 * r, r) if peskin6 else
            "C-electrolito-hann%dp%d-LD%g-ref%g" % (n, hp, 6 * r, r) if hp != 2 else
            "C-electrolito-hann%d-LD%g-ref%g" % (n, 6 * r, r))
    if q != 1.0: name += "-q%g" % q
    zc = 1.5 + D0 * r; d = os.path.join(HERE, name, "z_%05.2f" % zc)
    nsteps = int(math.ceil(5000 * r * r / 1000.0) * 1000)
    for sd in ("", "logs", "colloid_data", "proceced_data"):
        os.makedirs(os.path.join(d, sd), exist_ok=True)
    t = open(os.path.join(BASE, "input")).read()
    kern = ("gaussiano sigma = %g nodos (= 1.44612 fisico), corte %g sigma" % (SIG8 * r, GCUT) if gauss else
            "Peskin-6 estirado s = %d (soporte %d nodos = 12 fisico, sigma = 1.69006 fisico)" % (s6, 6 * s6) if peskin6
            else "Hann-%d cos^%d (ancho fisico 8 = Hann-8, C^%d)" % (n, hp, hp - 1) if hp != 2
            else "Hann-%d (ancho fisico 8 = Hann-8)" % n)
    t = HEADER.format(r=r, kern=kern, L=L, ld=6 * r, d=D0 * r) + t[t.index("\nN_start"):]
    t = sub(t, r"^N_cycles .*$", "N_cycles %d             # 5000 r^2 (D por paso igual que Hann-8)" % nsteps)
    t = sub(t, r"^size .*$", "size %d_%d_%d          # 32r x 32r x (32r + 2 nodos de pared)" % (L, L, L + 2))
    t = sub(t, r"^electrokinetics_epsilon .*$",
            "electrokinetics_epsilon %.8e   # 1e4/r: q/(eps kT) = 10 r en unidades de red, l_B = %.6f" % (EPS0 / r, 0.7957747 * r))
    t = sub(t, r"^electrokinetics_init_rho_el .*$",
            "electrokinetics_init_rho_el %.8e   # 1.388889e-3/r^3 -> lambda_D = 6 r = %g" % (RHO0 / r ** 3, 6 * r))
    krange = math.ceil(GCUT * SIG8 * r) if gauss else 3 * s6 - 1 if peskin6 else math.ceil(n / 2)
    sigk = SIGP6 * s6 if peskin6 else 8 * r * math.sqrt(1 / 12 - 5 / (8 * math.pi ** 2)) if hp == 4 else SIG8 * r
    assert hp in (2, 4), "sigma del kernel solo calculada para p = 2, 4"
    t = sub(t, r"^pm_sr_sigma .*$", "pm_sr_sigma       %.5f      # = sigma_kernel (no se usa: pm_sr_correction no)" % sigk)
    t = sub(t, r"^pm_sr_rcut .*$", "pm_sr_rcut        %.1f      # 2 + krange (no se usa)" % (2 + krange))
    if gauss:
        t = sub(t, r"^subgrid_kernel .*$", "subgrid_kernel    gauss    # kernel de spread/gather: gaussiano truncado (C^inf salvo el corte)")
        t = sub(t, r"^subgrid_hann_order .*$",
                "subgrid_gauss_sigma   %.6f   # = 1.44612 r: misma varianza que Hann-8 en unidades fisicas\n"
                "subgrid_gauss_support %.6f   # = 2 x %g sigma (soporte completo); g - g(corte) = %.1e relativo"
                % (SIG8 * r, 2 * GCUT * SIG8 * r, GCUT, math.exp(-GCUT ** 2 / 2)))
    elif peskin6:
        t = sub(t, r"^subgrid_kernel .*$", "subgrid_kernel    peskin6  # kernel de spread/gather: Peskin 6 puntos (Bao et al. 2016), C^3")
        t = sub(t, r"^subgrid_hann_order .*$",
                "subgrid_peskin6_stretch %d   # s = 2r: w(u) = phi6(u/s)/s, soporte |u| < 3s = %d nodos (6 fisico),\n"
                "#                           sigma = 0.84503 s = %.5f nodos (1.69006 fisico); halo >= 3s - 1 = %d"
                % (s6, 3 * s6, SIGP6 * s6, krange))
    elif hp != 2:
        t = sub(t, r"^subgrid_hann_order .*$",
                "subgrid_hann_order %d      # soporte Hann: ancho fisico n/r = 8\n"
                "subgrid_hann_power %d       # w ~ cos^%d(pi u/n): C^%d, suma 1 exacta (n >= %d); sigma = %.5f nodos (%.5f fisico)"
                % (n, hp, hp, hp - 1, hp // 2 + 1, sigk, sigk / r))
    else:
        t = sub(t, r"^subgrid_hann_order .*$", "subgrid_hann_order %d      # orden Hann: ancho fisico n/r = 8" % n)
    t = sub(t, r"^vel_io_freq .*$", "vel_io_freq %d" % nsteps)
    t = sub(t, r"^freq_config .*$", "freq_config %d        # N/4: psi intermedios para ver la convergencia al estacionario" % (nsteps // 4))
    t = sub(t, r"^# lambda_D = 6 -> .*$",
            "# lambda_D = 6r = %g -> rho_s = 1/(8 pi l_B lambda_D^2) = %.6e por especie (l_B = %.6f)"
            % (6 * r, RHO0 / r ** 3, 0.7957747 * r))
    t = sub(t, r"^# Solo importa el ESTACIONARIO.*$",
            "# Solo importa el ESTACIONARIO (no depende de D). D por paso igual que Hann-8: en unidades fisicas el paso\n"
            "# de tiempo es 1/r^2 del de Hann-8, por eso N = 5000 r^2.")
    if q != 1.0:
        t += ("\n# --- Prueba de linealidad: q0 = %g (en config.cds.init) en lugar de 1. Si Debye-Huckel linealizado fuera\n"
              "# exacto, psi/q seria independiente de q; la parte no lineal de psi/q escala como q^2. ---\n" % q)
    open(os.path.join(d, "input"), "w").write(t)
    lines = open(os.path.join(BASE, "config.cds.init.001-001")).read().split("\n")
    assert lines[35].split() == ["1.600000000000000e+01", "1.600000000000000e+01", "9.500000000000000e+00"]
    lines[35] = "".join("{:>25.15e}".format(v) for v in (16 * r, 16 * r, zc))
    assert lines[46].split() == ["1.000000000000000e+00"]                  # q0
    lines[46] = "{:>25.15e}".format(q)
    open(os.path.join(d, "config.cds.init.001-001"), "w").write("\n".join(lines))
    shutil.copy(EXE, os.path.join(d, "Ludwig.exe")); shutil.copy(os.path.join(BASE, ".petscrc"), os.path.join(d, ".petscrc"))
    shutil.copy(os.path.join(HERE, "C-electrolito-hann8-LD6", "correr.sh"), os.path.join(HERE, name, "correr.sh"))
    print("%s: size %d_%d_%d, N = %d, carga q0 = %g en (%g, %g, %g)" % (name, L, L, L + 2, nsteps, q, 16 * r, 16 * r, zc))

if __name__ == "__main__":
    args = sys.argv[1:]; gauss = "gauss" in args; peskin6 = "peskin6" in args
    q = float(([a[2:] for a in args if a.startswith("q=")] or ["1"])[0])
    hp = int(([a[2:] for a in args if a.startswith("p=")] or ["2"])[0])
    args = [a for a in args if a not in ("gauss", "peskin6") and not a.startswith(("q=", "p="))]
    for n in ([int(a) for a in args] or ((8, 10, 12, 16) if gauss else (4, 8, 16) if peskin6 else (10, 12))):
        make(n, gauss, q, peskin6, hp)
