#!/usr/bin/env python3
"""
Perfiles de campo electrico (E), potencial (psi) y carga (rho/qsi) para un DIMERO.

- Perfil 1D: a lo largo de la recta que une las 2 particulas.
- Planos 2D: los dos planos que contienen esa recta (el que tiene la normal en
  cada una de las dos direcciones transversales).

Formato de datos Ludwig (ver plot_electric_field.py como base):
  - psi / rho / qsi : 1 columna, nx*ny*nz filas, orden z-major -> reshape (nx,ny,nz)
  - efield          : 3 columnas (Ex,Ey,Ez)
  - Nodo entero i (0-based) esta en la posicion fisica i + 0.5 (Lmin=0.5).
  - Posiciones de particulas: colloid_data/colloids-XXXXXXXX.csv (id,x,y,z,...)

Uso tipico:
  python3 plot_dimero_perfil.py -d <run_dir> -step 2000
  python3 plot_dimero_perfil.py -d <run_dir> -step 2000 -L 32 -o salida_
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")            # sin display (evita warnings Qt/wayland en WSL)
import matplotlib.pyplot as plt
import argparse, os, sys, glob


# ---------------------------------------------------------------- lectura ----
def load_scalar(path, L):
    d = np.loadtxt(path)
    if d.size == L**3:
        return d.reshape((L, L, L))  # z-major
    raise ValueError(f"{path}: {d.size} valores, esperaba {L**3}")


def load_charge(path, L, valency=(+1.0, -1.0)):
    """Carga neta desde qsi (2 columnas = especies rho0, rho1) o rho (1 col).

    qsi: carga_neta = z0*rho0 + z1*rho1  (por defecto z=+1,-1).
    """
    d = np.loadtxt(path)
    if d.ndim == 1 and d.size == L**3:
        return d.reshape((L, L, L))
    if d.ndim == 2 and d.shape[0] == L**3:
        nk = d.shape[1]
        z = np.asarray(valency[:nk], float)
        net = d @ z                     # combinacion lineal por especie
        return net.reshape((L, L, L))
    raise ValueError(f"{path}: forma {d.shape} inesperada para L={L}")


def load_vector(path, L):
    d = np.loadtxt(path)
    if d.shape[0] != L**3:
        raise ValueError(f"{path}: {d.shape[0]} filas, esperaba {L**3}")
    return d.reshape((L, L, L, 3))


def find_field(run_dir, prefix, step):
    """Busca <prefix>-<step con ceros>.001-001 de forma flexible."""
    pats = [os.path.join(run_dir, f"{prefix}-{step:09d}.001-001"),
            os.path.join(run_dir, f"{prefix}-{step}.001-001"),
            os.path.join(run_dir, f"{prefix}-*{step}*.001-001")]
    for p in pats:
        hits = sorted(glob.glob(p))
        hits = [h for h in hits if "metadata" not in h]
        if hits:
            return hits[0]
    # fallback: el ultimo disponible
    hits = sorted(h for h in glob.glob(os.path.join(run_dir, f"{prefix}-*.001-001"))
                  if "metadata" not in h)
    if hits:
        print(f"  aviso: no hallado {prefix} step {step}, uso {os.path.basename(hits[-1])}")
        return hits[-1]
    return None


def load_particles(run_dir, step):
    """Devuelve array (n,3) con posiciones de las particulas mas cercano a step."""
    cdir = os.path.join(run_dir, "colloid_data")
    cands = sorted(glob.glob(os.path.join(cdir, "colloids-*.csv")))
    if not cands:
        cands = sorted(glob.glob(os.path.join(run_dir, "colloids-*.csv")))
    if not cands:
        raise FileNotFoundError("No hay colloids-*.csv")

    def cyc(f):
        return int(os.path.basename(f).split("-")[1].split(".")[0])

    # elegir el csv con datos (>=1 particula) mas cercano al step por debajo
    best = None
    for f in cands:
        arr = np.genfromtxt(f, delimiter=",", skip_header=1)
        if arr.size == 0:
            continue
        arr = np.atleast_2d(arr)
        if arr.shape[1] < 4:
            continue
        if cyc(f) <= step or best is None:
            best = (f, arr)
    if best is None:
        raise ValueError("Ningun colloids csv con datos")
    f, arr = best
    print(f"  particulas de {os.path.basename(f)}: {arr.shape[0]}")
    return arr[:, 1:4]  # x,y,z


# --------------------------------------------------------- interpolacion 1D --
def sample_line(field, p0, p1, npts, comp=None):
    """
    Muestrea 'field' (LxLxL o LxLxLx3) a lo largo del segmento p0->p1 (coords
    fisicas). Convierte a indice de nodo restando 0.5 y usa interpolacion
    trilineal periodica. Devuelve (t, valores).

    comp: si field es vectorial, 'x'/'y'/'z'/'mag' o None (=mag).
    """
    L = field.shape[0]
    p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
    t = np.linspace(0.0, 1.0, npts)
    pts = p0[None, :] + t[:, None] * (p1 - p0)[None, :]
    idx = pts - 0.5  # posicion fisica i+0.5 <-> nodo i

    vals = np.array([trilinear(field, q) for q in idx])
    if field.ndim == 4:
        if comp in ("x", 0):   vals = vals[:, 0]
        elif comp in ("y", 1): vals = vals[:, 1]
        elif comp in ("z", 2): vals = vals[:, 2]
        else:                  vals = np.linalg.norm(vals, axis=1)
    dist = t * np.linalg.norm(p1 - p0)
    return dist, vals


def trilinear(field, q):
    """Interp. trilineal periodica en indice continuo q=(x,y,z)."""
    L = field.shape[0]
    i0 = np.floor(q).astype(int)
    f = q - i0
    out = 0.0 if field.ndim == 3 else np.zeros(3)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((f[0] if dx else 1 - f[0]) *
                     (f[1] if dy else 1 - f[1]) *
                     (f[2] if dz else 1 - f[2]))
                ii = (i0[0] + dx) % L
                jj = (i0[1] + dy) % L
                kk = (i0[2] + dz) % L
                out = out + w * field[ii, jj, kk]
    return out


# ------------------------------------------------------------ teoria D-H ----
def read_input_params(run_dir):
    """Lee epsilon, temperature, rho_el del fichero input de la corrida."""
    p = {"epsilon": 1.0e4, "beta": 1.0e5, "rho_el": 0.0, "e": 1.0}
    fin = os.path.join(run_dir, "input")
    if not os.path.exists(fin):
        return p
    with open(fin) as f:
        for ln in f:
            ln = ln.split("#")[0].strip()
            if not ln:
                continue
            key = ln.split()[0]
            try:
                val = float(ln.split()[-1])
            except (ValueError, IndexError):
                continue
            if key == "electrokinetics_epsilon":
                p["epsilon"] = val
            elif key == "temperature" and val > 0:
                p["beta"] = 1.0 / val
            elif key == "electrokinetics_init_rho_el":
                p["rho_el"] = val
            elif key == "electrokinetics_eunit":
                p["e"] = val
    return p


def read_particle_charge(run_dir, default=1.0):
    """Lee la carga q0 de la primera particula del config.cds.init.001-001.

    En el formato ascii, q0 es el 3er double despues de deltaphi; es mas robusto
    leer del nombre de carpeta (q_<valor>) si el parseo falla.
    """
    # intentar del nombre de carpeta: ..._q_1.0_...
    import re
    m = re.search(r"[-_]q[_-]([0-9.eE+-]+)", os.path.basename(run_dir.rstrip("/")))
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return default


def debye_kappa(params, valency=(+1.0, -1.0)):
    """kappa = sqrt(beta e^2 sum_i z_i^2 rho_i / epsilon).

    Con electrolito simetrico rho_el por especie: sum z^2 rho = (z0^2+z1^2)*rho_el.
    Devuelve (kappa, lambda_D). Si no hay sal, kappa=0 (Coulomb puro).
    """
    z = np.asarray(valency, float)
    s2 = float(np.sum(z**2)) * params["rho_el"]
    if s2 <= 0.0:
        return 0.0, np.inf
    kap2 = params["beta"] * params["e"]**2 * s2 / params["epsilon"]
    kap = np.sqrt(kap2)
    return kap, 1.0 / kap


def psi_dh_line(d_line, p0, pA, pB, u, params, q, kappa):
    """Potencial teorico Debye-Hueckel de 2 cargas puntuales q sobre la recta.

    UNIDADES DE LUDWIG: el solver resuelve  eps * lap(psi) = eunit*beta*rho
    (ver src/psi_sor.c:238), es decir  lap(psi) = -(beta e / eps) rho.  Por eso
    el potencial de Ludwig lleva un factor beta*e extra respecto al fisico.
    -> prefactor = beta * q e^2/(4 pi eps) = q * lB  (lB = long. de Bjerrum).
    Se evalua en los mismos puntos d_line del muestreo del dato.
    """
    pre = params["beta"] * q * params["e"]**2 / (4.0 * np.pi * params["epsilon"])
    # d_line=0 corresponde al punto p0 real de la recta muestreada
    pts = p0[None, :] + (d_line[:, None]) * u[None, :]
    out = np.zeros_like(d_line)
    for center in (pA, pB):
        r = np.linalg.norm(pts - center[None, :], axis=1)
        r = np.clip(r, 1e-6, None)          # evitar singularidad en el centro
        out += pre * np.exp(-kappa * r) / r
    return out


def Epar_dh_line(d_line, p0, pA, pB, u, params, q, kappa):
    """Componente de E paralela a la recta, teorica D-H (E = -grad psi).

    Para cada carga: E(r) = pre * (1+kappa r) exp(-kappa r)/r^2 * rhat.
    Proyectamos sobre u.

    UNIDADES: a diferencia de psi, el efield de Ludwig esta en unidades FISICAS
    (E_ludwig = -grad(psi_ludwig)/(beta*e), verificado numericamente). Por eso el
    prefactor de E NO lleva el factor beta -> pre = q e^2/(4 pi eps).
    """
    pre = q * params["e"]**2 / (4.0 * np.pi * params["epsilon"])
    pts = p0[None, :] + (d_line[:, None]) * u[None, :]
    Epar = np.zeros_like(d_line)
    for center in (pA, pB):
        dvec = pts - center[None, :]
        r = np.linalg.norm(dvec, axis=1)
        r = np.clip(r, 1e-6, None)
        mag = pre * (1.0 + kappa * r) * np.exp(-kappa * r) / r**2
        rhat = dvec / r[:, None]
        Epar += mag * (rhat @ u)            # proyeccion sobre la recta
    return Epar


# ------------------------------------------------ teoria carga Gaussiana ----
try:
    from scipy.special import erfc as _erfc
except ImportError:
    from math import erf as _merf
    _erfc = np.vectorize(lambda x: 1.0 - _merf(x))


def _gauss_dh_radial(r, kappa, sigma):
    """Factor radial del potencial D-H de una carga Gaussiana de ancho sigma.

    Para rho(r) = q * (2 pi sigma^2)^(-3/2) exp(-r^2/2 sigma^2), la solucion de
    (lap - kappa^2) psi = -rho/eps_efectivo es (sin el prefactor q/(4 pi eps)):

        f(r) = 1/(2r) [ e^{ k r} erfc( r/(s2) + k s/ (s2)*... ) ... ]

    Forma cerrada estandar (screened Gaussian, Yukawa suavizado):
        f(r) = 1/(2r) [ e^{-kappa r} erfc( (sigma^2 kappa - r)/(sqrt2 sigma) )
                       - e^{ kappa r} erfc( (sigma^2 kappa + r)/(sqrt2 sigma) ) ]
    Limite sigma->0 : f(r) -> e^{-kappa r}/r (puntual). Finito en r->0.
    """
    s = sigma
    r = np.clip(r, 1e-9, None)
    a = (s * s * kappa - r) / (np.sqrt(2.0) * s)
    b = (s * s * kappa + r) / (np.sqrt(2.0) * s)
    return (np.exp(-kappa * r) * _erfc(a)
            - np.exp(kappa * r) * _erfc(b)) / (2.0 * r)


def psi_gauss_line(d_line, p0, pA, pB, u, params, q, kappa, sigma):
    """Potencial D-H de 2 cargas Gaussianas (ancho sigma). Unidades de Ludwig
    (prefactor con beta, igual que psi_dh_line). Finito en los centros."""
    pre = params["beta"] * q * params["e"]**2 / (4.0 * np.pi * params["epsilon"])
    pts = p0[None, :] + (d_line[:, None]) * u[None, :]
    out = np.zeros_like(d_line)
    for center in (pA, pB):
        r = np.linalg.norm(pts - center[None, :], axis=1)
        out += pre * _gauss_dh_radial(r, kappa, sigma)
    return out


def Epar_gauss_line(d_line, p0, pA, pB, u, params, q, kappa, sigma):
    """E paralela de 2 cargas Gaussianas (E = -grad psi), unidades FISICAS
    (prefactor SIN beta, igual que Epar_dh_line). Derivada numerica de f(r)."""
    pre = q * params["e"]**2 / (4.0 * np.pi * params["epsilon"])
    pts = p0[None, :] + (d_line[:, None]) * u[None, :]
    Epar = np.zeros_like(d_line)
    h = 1e-4
    for center in (pA, pB):
        dvec = pts - center[None, :]
        r = np.linalg.norm(dvec, axis=1)
        # E_r = -d psi/dr = -pre * f'(r); derivada central de f
        fp = (_gauss_dh_radial(r + h, kappa, sigma)
              - _gauss_dh_radial(r - h, kappa, sigma)) / (2.0 * h)
        rsafe = np.clip(r, 1e-9, None)
        rhat = dvec / rsafe[:, None]
        Epar += (-pre * fp) * (rhat @ u)
    return Epar


def cloud_integral(L, centers, params, q, kappa, sigma=None):
    """Integra la nube de carga D-H teorica sobre TODA la caja LxLxL.

    rho_nube(r) = -eps*kappa^2 * psi_fisico,  psi = pre * f(r)  (pre con beta).
    Devuelve sum(rho_nube) sobre los nodos (dV=1). Sirve para fijar el fondo
    neutralizante correcto en caja periodica finita: si la nube integra a Qc,
    el fondo que falta para neutralidad es (-q_total - Qc)/V, NO -q_total/V.
    """
    pre = params["beta"] * q * params["e"]**2 / (4.0 * np.pi * params["epsilon"])
    xs = np.arange(L) + 0.5
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    rho = np.zeros((L, L, L))
    for c in centers:
        r = np.sqrt((X - c[0])**2 + (Y - c[1])**2 + (Z - c[2])**2)
        if sigma is None:
            r = np.clip(r, 1e-6, None)
            f = np.exp(-kappa * r) / r
        else:
            f = _gauss_dh_radial(r, kappa, sigma)
        psi_fis = pre * f / (params["beta"] * params["e"])
        rho += -params["epsilon"] * kappa**2 * psi_fis
    return rho.sum()


# --------------------------------------------------------------- graficos ----
def plot_line(run_dir, step, L, out_prefix, extend=4.0, npts=400,
              theory=True, q=None, full_box=False, gauss_sigma=None):
    pos = load_particles(run_dir, step)
    if pos.shape[0] < 2:
        print("  aviso: <2 particulas; uso la unica dos veces")
        pos = np.vstack([pos, pos])
    pA, pB = pos[0], pos[1]
    axis = pB - pA
    Lax = np.linalg.norm(axis)
    if Lax < 1e-9:
        raise ValueError("Las 2 particulas coinciden")
    u = axis / Lax
    if full_box:
        # recta centrada en el punto medio del dimero, de longitud = L (todo el
        # ancho de la caja a lo largo de la direccion del dimero). Periodica.
        mid = 0.5 * (pA + pB)
        p0 = mid - 0.5 * L * u
        p1 = mid + 0.5 * L * u
    else:
        p0 = pA - extend * u        # extender antes de A
        p1 = pB + extend * u        # y despues de B
    dA = np.dot(pA - p0, u)
    dB = np.dot(pB - p0, u)

    psi = load_scalar(find_field(run_dir, "psi", step), L)
    # carga: preferir 'qsi' (carga neta) si existe, si no 'rho'
    fq = find_field(run_dir, "qsi", step) or find_field(run_dir, "rho", step)
    rho = load_charge(fq, L)
    E = load_vector(find_field(run_dir, "efield", step), L)

    d, vpsi = sample_line(psi, p0, p1, npts)
    _, vrho = sample_line(rho, p0, p1, npts)
    _, vEmag = sample_line(E, p0, p1, npts, comp="mag")
    _, vEpar = sample_line(E, p0, p1, npts, comp=None) if False else (d, None)
    # componente de E a lo largo de la recta (proyeccion sobre u)
    Evec = np.array([trilinear(E, (p0 + t*(p1-p0) - 0.5))
                     for t in np.linspace(0, 1, npts)])
    vEpar = Evec @ u

    # ---- muestreo en los NODOS reales del grid sobre la recta ----
    # Las posiciones de nodo que caen sobre la recta: proyectamos los centros de
    # nodo (i+0.5) del eje dominante del dimero. dn = distancia sobre la recta.
    Ltot = np.linalg.norm(p1 - p0)
    dn = np.arange(np.ceil(-0.0), np.floor(Ltot) + 1e-9) + (0.5 - (np.dot(p0, u) % 1.0))
    dn = dn[(dn >= 0) & (dn <= Ltot)]
    pn = p0[None, :] + dn[:, None] * u[None, :]        # puntos sobre la recta
    idxn = pn - 0.5                                    # indice de nodo continuo
    npsi = np.array([trilinear(psi, q) for q in idxn])
    nrho = np.array([trilinear(rho, q) for q in idxn])
    nEvec = np.array([trilinear(E, q) for q in idxn])
    nEmag = np.linalg.norm(nEvec, axis=1)
    nEpar = nEvec @ u

    # ---- curvas teoricas Debye-Hueckel (carga puntual, 1/r apantallado) ----
    tstr = ""
    psi_th = rho_th = Epar_th = None
    use_gauss = gauss_sigma is not None
    if theory:
        params = read_input_params(run_dir)
        if q is None:
            q = read_particle_charge(run_dir)
        kappa, lamD = debye_kappa(params)
        if use_gauss:
            psi_th = psi_gauss_line(d, p0, pA, pB, u, params, q, kappa, gauss_sigma)
            Epar_th = Epar_gauss_line(d, p0, pA, pB, u, params, q, kappa, gauss_sigma)
            thlabel = f"D-H Gaussiana ($\\sigma$={gauss_sigma:g})"
        else:
            psi_th = psi_dh_line(d, p0, pA, pB, u, params, q, kappa)
            Epar_th = Epar_dh_line(d, p0, pA, pB, u, params, q, kappa)
            thlabel = "Debye-Hückel (puntual)"
        # carga inducida D-H: rho_fisica = -eps * kappa^2 * psi_fisico.
        # psi_th esta en unidades de Ludwig (factor beta*e); la carga rho del
        # output es fisica -> dividir por beta*e para volver a unidades de carga.
        psi_fis = psi_th / (params["beta"] * params["e"])
        rho_th = -params["epsilon"] * kappa**2 * psi_fis
        # FONDO NEUTRALIZANTE (caja periodica finita): el sistema es neutro, la
        # carga del fluido debe ser -q_total. La nube D-H de bano infinito
        # integrada SOBRE LA CAJA da Qc (< -q_total en magnitud si lambda_D ~ caja,
        # porque parte de la nube cae fuera). El deficit se reparte como fondo
        # uniforme:  rho_bg = (-q_total - Qc)/V.  (NO es -q_total/V: eso contaria
        # la neutralizacion dos veces, ya que la nube ya integra a Qc.)
        # Con lambda_D grande (poca sal) este fondo es visible en el campo lejano.
        q_total = q * pos.shape[0]
        Qc = cloud_integral(L, (pA, pB), params, q, kappa,
                            sigma=(gauss_sigma if use_gauss else None))
        rho_bg = (-q_total - Qc) / float(L**3)
        rho_th = rho_th + rho_bg
        # gauge: Ludwig fija <psi>=0 en la caja (fondo != 0). La teoria decae a 0.
        # Alinear en el CAMPO LEJANO (extremos de la recta, lejos de las cargas),
        # donde ambas curvas son planas y comparables. Se ignora la media global
        # (contaminada por la divergencia 1/r de la teoria puntual en los centros).
        nfar = max(3, len(d) // 20)
        far = np.r_[np.arange(nfar), np.arange(len(d) - nfar, len(d))]
        offset = np.mean(vpsi[far]) - np.mean(psi_th[far])
        psi_th = psi_th + offset
        tstr = (f"  |  {'Gauss ' if use_gauss else ''}q={q:g}, "
                + (f"$\\lambda_D$={lamD:.2f}, $\\kappa$={kappa:.3f}"
                   if np.isfinite(lamD) else "Coulomb (sin sal)")
                + (f", $\\sigma$={gauss_sigma:g}" if use_gauss else ""))

    def ylim_from(data_ref):
        """Rango Y con margen segun el DATO de Ludwig (no la teoria singular)."""
        lo, hi = np.min(data_ref), np.max(data_ref)
        pad = 0.15 * (hi - lo if hi > lo else abs(hi) + 1e-30)
        return lo - pad, hi + pad

    mk = dict(marker="o", ms=4, ls="none", mfc="none", mew=1.2, zorder=5)
    fig, ax = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    # -- psi --
    ax[0].plot(d, vpsi, "-", color="tab:blue", label="Ludwig")
    ax[0].plot(dn, npsi, color="tab:blue", label="nodos Ludwig", **mk)
    if theory:
        ax[0].plot(d, psi_th, "--", color="tab:orange", label=thlabel)
    ax[0].set_ylabel(r"$\psi$ (potencial)")
    ax[0].set_ylim(*ylim_from(vpsi))     # acotar a la escala del dato
    ax[0].legend(loc="upper right", fontsize=8)
    # -- rho --
    ax[1].plot(d, vrho, "-", color="tab:red", label="Ludwig")
    ax[1].plot(dn, nrho, color="tab:red", label="nodos Ludwig", **mk)
    if theory:
        ax[1].plot(d, rho_th, "--", color="tab:orange", label="D-H")
    ax[1].set_ylabel(r"$\rho$ (carga)")
    ax[1].set_ylim(*ylim_from(vrho))
    ax[1].legend(loc="upper right", fontsize=8)
    # -- E --
    ax[2].plot(d, vEmag, "-", color="k", label=r"$|E|$ Ludwig")
    ax[2].plot(dn, nEmag, color="k", **mk)
    ax[2].plot(d, vEpar, "--", color="tab:green", label=r"$E_\parallel$ Ludwig")
    ax[2].plot(dn, nEpar, color="tab:green", label="nodos Ludwig", **mk)
    if theory:
        ax[2].plot(d, Epar_th, ":", color="tab:orange", lw=2,
                   label=r"$E_\parallel$ D-H")
    ax[2].set_ylabel("Campo E")
    ax[2].set_xlabel("distancia a lo largo de la recta")
    ax[2].set_ylim(*ylim_from(np.concatenate([vEmag, vEpar])))
    ax[2].legend(loc="upper right", fontsize=8)
    for a in ax:
        a.axvline(dA, color="gray", ls=":", lw=1)
        a.axvline(dB, color="gray", ls=":", lw=1)
        a.grid(alpha=0.3)
        # marcador de posicion de las particulas SOBRE el eje X
        a.plot([dA, dB], [0, 0], marker="^", ms=11, ls="none",
               color="magenta", mec="k", mew=0.6,
               transform=a.get_xaxis_transform(), clip_on=False, zorder=10)
    # etiqueta de las particulas solo bajo el eje del panel inferior
    ax[2].annotate("P1", (dA, 0), xycoords=ax[2].get_xaxis_transform(),
                   xytext=(0, -22), textcoords="offset points",
                   ha="center", fontsize=8, color="magenta")
    ax[2].annotate("P2", (dB, 0), xycoords=ax[2].get_xaxis_transform(),
                   xytext=(0, -22), textcoords="offset points",
                   ha="center", fontsize=8, color="magenta")
    ax[0].set_title(f"Perfil sobre la recta del dimero (step {step})\n"
                    f"P1={np.round(pA,2)}  P2={np.round(pB,2)}" + tstr)
    fig.tight_layout()
    fn = f"{out_prefix}linea_step{step}.png"
    fig.savefig(fn, dpi=150)
    print("  guardado", fn)
    plt.close(fig)


def plot_planes(run_dir, step, L, out_prefix):
    """Dos planos que contienen la recta. Para un dimero alineado con X, esos
    planos son XY (z=cte) y XZ (y=cte) pasando por las particulas."""
    pos = load_particles(run_dir, step)
    pA, pB = pos[0], pos[1]
    mid = 0.5 * (pA + pB)
    axis = np.argmax(np.abs(pB - pA))   # eje dominante del dimero (0=x,1=y,2=z)

    psi = load_scalar(find_field(run_dir, "psi", step), L)
    fq = find_field(run_dir, "qsi", step) or find_field(run_dir, "rho", step)
    rho = load_charge(fq, L)
    E = load_vector(find_field(run_dir, "efield", step), L)
    Emag = np.linalg.norm(E, axis=3)

    # las dos direcciones transversales al eje del dimero
    trans = [a for a in (0, 1, 2) if a != axis]
    names = {0: "x", 1: "y", 2: "z"}

    for tnorm in trans:  # normal del plano = una de las transversales
        k = int(round(mid[tnorm] - 0.5)) % L   # indice de nodo del plano
        # ejes del plano: axis (largo del dimero) y la otra transversal
        other = [a for a in (0, 1, 2) if a not in (tnorm,)]
        # construir slices
        def take(vol):
            sl = [slice(None)] * 3
            sl[tnorm] = k
            return vol[tuple(sl)]  # queda 2D en (other[0], other[1])

        P = take(psi); R = take(rho); EM = take(Emag)
        a0, a1 = other  # ejes que quedan
        ext = [0.5, L + 0.5, 0.5, L + 0.5]  # coords fisicas aprox

        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        for a, dat, tit in zip(ax, [P, R, EM],
                               [r"$\psi$", r"$\rho$", r"$|E|$"]):
            im = a.imshow(dat.T, origin="lower", extent=ext, aspect="equal",
                          cmap="RdBu_r" if tit != r"$|E|$" else "viridis")
            a.set_title(f"{tit}  (plano {names[a0]}{names[a1]}, {names[tnorm]}={k+0.5:.1f})")
            a.set_xlabel(names[a0]); a.set_ylabel(names[a1])
            fig.colorbar(im, ax=a, fraction=0.046)
            # marcar particulas proyectadas
            a.plot([pA[a0], pB[a0]], [pA[a1], pB[a1]], "k+", ms=10, mew=2)
        fig.tight_layout()
        fn = f"{out_prefix}plano_{names[a0]}{names[a1]}_step{step}.png"
        fig.savefig(fn, dpi=150)
        print("  guardado", fn)
        plt.close(fig)


# ------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dir", required=True, help="directorio de la corrida")
    ap.add_argument("-step", type=int, default=2000, help="time-step a graficar")
    ap.add_argument("-L", type=int, default=32, help="lado del grid")
    ap.add_argument("-o", "--out", default="dimero_", help="prefijo de salida")
    ap.add_argument("--no-planes", action="store_true")
    ap.add_argument("--no-theory", action="store_true",
                    help="no superponer la curva teorica Debye-Hueckel")
    ap.add_argument("-q", "--charge", type=float, default=None,
                    help="carga de cada particula (default: del nombre de carpeta)")
    ap.add_argument("--full-box", action="store_true",
                    help="la recta cubre todo el ancho L de la caja (periodica), "
                         "no solo el tramo entre particulas +/- margen")
    ap.add_argument("--extend", type=float, default=4.0,
                    help="margen a cada lado de las particulas (si no --full-box)")
    ap.add_argument("--gaussian-sigma", type=float, default=None, metavar="SIGMA",
                    help="usar teoria D-H con cargas GAUSSIANAS de ancho SIGMA "
                         "(finita en los centros) en vez de cargas puntuales")
    a = ap.parse_args()

    run = a.dir
    outp = os.path.join(run, a.out)
    print("Perfil 1D...")
    plot_line(run, a.step, a.L, outp, extend=a.extend,
              theory=not a.no_theory, q=a.charge, full_box=a.full_box,
              gauss_sigma=a.gaussian_sigma)
    if not a.no_planes:
        print("Planos 2D...")
        plot_planes(run, a.step, a.L, outp)
    print("Listo.")


if __name__ == "__main__":
    main()
