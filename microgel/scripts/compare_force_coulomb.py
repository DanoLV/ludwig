#!/home/danolv/Sim/ludwig/.venv/bin/python3
"""
Compara la fuerza simulada entre dos colloides contra tres curvas teoricas:

  - Coulomb puro (sin apantallamiento):
        F(r) = q1 * q2 * e_unit^2 / (4 * pi * epsilon * r^2)

  - Debye-Huckel apantallado, carga puntual (nube ionica linealizada, con
    el kappa que ya calcula el propio codigo para esta corrida):
        f(d)    = q1 * q2 * e_unit^2 / (4 * pi * epsilon)
                  * exp(-kappa * d) * (kappa / d + 1 / d^2)

  - Debye-Huckel apantallado, carga gaussiana (--gaussian-sigma sigma):
    cada colloide modelado como una nube de carga gaussiana de ancho sigma
    (no puntual), resolviendo exactamente la ecuacion de Debye-Huckel
    linealizada para dos gaussianas de ancho combinado
    sigma_pair = sigma*sqrt(2). Potencial cerrado (obtenido via
    parametrizacion de Schwinger 1/(k^2+kappa^2) = integral de e^{-t(k^2+
    kappa^2)}dt en el espacio de Fourier, resuelto en forma cerrada):
        B = sigma_pair*kappa/sqrt(2), C = 1/(sigma_pair*sqrt(2))
        T1(d) = exp(-kappa*d) * erfc(B - C*d)
        T2(d) = exp(+kappa*d) * erfc(B + C*d)
        phi(d) = [q1*q2*e_unit^2/(8*pi*epsilon)] * exp(kappa^2*sigma_pair^2/2)
                 * (T1 - T2) / d
    El factor exp(kappa^2*sigma_pair^2/2) es imprescindible (viene de
    completar la integral) y se anula en los limites sigma->0 (recupera el
    potencial de carga puntual) y kappa->0 (recupera el potencial estandar
    de Coulomb con carga gaussiana, erf(d/(sigma_pair*sqrt(2)))/d) — por eso
    NO alcanza con chequear esos dos limites para validar la formula (un
    primer intento sin ese factor los pasaba igual). Se verifico numerica-
    mente contra la convolucion 3D directa antes de darla por buena. La
    fuerza f(d) = -dphi/dr se calcula por diferencias finitas centradas del
    potencial (una derivada analitica a mano resulto propensa a errores).

  Para Debye-Huckel (puntual y gaussiana) se suman ademas las imagenes
  periodicas de la particula fija a lo largo del eje x cada L (el dominio
  es periodico en las 3 direcciones, ver "periodicity 1_1_1" en el input):
  la particula movil en r "ve" copias de la carga fija en r - n*L para todo
  entero n, cada una repeliendola en la direccion que corresponda. La serie
  converge rapido porque el apantallamiento es exponencial (ver --n-images):
        F(r) = sum_{n=-N}^{N} sign(r - n*L) * f(|r - n*L|)
  Sin esto, F(r) no se anula en r=L/2 como exige la simetria de la caja
  periodica (la simulacion si lo hace).

Recorre todos los subdirectorios "pos_X_Y_Z" del directorio base, y para
cada uno toma, en un paso de tiempo dado, la fuerza sobre la particula
movil (componente en el eje que une ambas particulas) desde
proceced_data/particle_force.csv (escrito por Ludwig, ver
src/subgrid.c:subgrid_print_force). La particula fija es la de menor x
(la movil siempre se desplaza hacia +x desde su posicion inicial), leida
desde colloid_data/colloids-*.csv en el mismo paso.

Los parametros fisicos (carga q0, permitividad epsilon, carga unitaria,
kappa) se leen del log de cada simulacion (logs/output.txt) salvo que se
pasen por linea de comandos.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc


POS_RE = re.compile(r'pos_([\d.]+)_([\d.]+)_([\d.]+)')


def find_sim_dirs(base_dir):
    dirs = []
    for name in sorted(os.listdir(base_dir)):
        full = os.path.join(base_dir, name)
        if os.path.isdir(full) and POS_RE.search(name):
            dirs.append(full)
    return dirs


def parse_physics(log_path):
    """Extrae q0, epsilon, unit_charge y kappa del log de Ludwig (logs/output.txt)."""
    if not os.path.isfile(log_path):
        return {}

    text = open(log_path, errors='replace').read()

    def grab(pattern, cast=float):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else None

    epsilon = grab(r'Permittivity:\s+([0-9.eE+-]+)')
    unit_charge = grab(r'Unit charge:\s+([0-9.eE+-]+)')
    qtot = grab(r'qtot=([0-9.eE+-]+)')
    n_near = grab(r'n_near=(\d+)', cast=int)
    q0 = qtot / n_near if (qtot is not None and n_near) else unit_charge
    kappa = grab(r'Kappa calculado:\s*([0-9.eE+-]+)')
    box_l = grab(r'System size:\s+([0-9]+)\s+[0-9]+\s+[0-9]+', cast=float)

    return dict(epsilon=epsilon, unit_charge=unit_charge, q0=q0, kappa=kappa, box_l=box_l)


def read_force_csv(path, step):
    """Lee particle_force.csv (Step;Index;Fmod;Fx;Fy;Fz, decimal ',').

    Devuelve (step_usado, {index: (Fx, Fy, Fz)}) para el Step pedido, o el
    mayor Step disponible <= step si no hay una coincidencia exacta.
    """
    steps = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';')
            s = int(parts[0])
            idx = int(parts[1])
            fx, fy, fz = (float(p.replace(',', '.')) for p in parts[3:6])
            steps.setdefault(s, {})[idx] = (fx, fy, fz)

    if not steps:
        return None, {}
    if step in steps:
        return step, steps[step]

    below = [s for s in steps if s <= step]
    chosen = max(below) if below else min(steps)
    return chosen, steps[chosen]


def read_colloid_positions(sim_dir, step):
    """Lee colloid_data/colloids-<step>.csv (o el disponible mas cercano <= step).

    Devuelve (step_usado, {index: (x, y, z)}).
    """
    step_of = {}
    for fpath in glob.glob(os.path.join(sim_dir, 'colloid_data', 'colloids-*.csv')):
        m = re.search(r'colloids-(\d+)\.csv$', fpath)
        if m:
            step_of[int(m.group(1))] = fpath

    if not step_of:
        return None, {}

    if step in step_of:
        chosen = step
    else:
        below = [s for s in step_of if s <= step]
        chosen = max(below) if below else min(step_of)

    data = np.genfromtxt(step_of[chosen], delimiter=',', skip_header=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    positions = {int(row[0]): (row[1], row[2], row[3]) for row in data}
    return chosen, positions


def main():
    parser = argparse.ArgumentParser(
        description='Compara la fuerza simulada (particle_force.csv) contra '
                     'la fuerza teorica de Coulomb puro entre dos colloides.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplo (el grafico queda en el directorio base, junto a las simulaciones):
  ./compare_force_coulomb.py /home/danolv/Sim/ludwig/microgel/local/TEST-2_particulas_Ffile -n 3000
""")

    parser.add_argument('base_dir', nargs='?', default='.',
                         help='Directorio con los subdirectorios pos_X_Y_Z '
                              '(default: directorio actual)')
    parser.add_argument('-n', '--step', type=int, required=True,
                         help='Paso de tiempo a comparar (usa el ultimo '
                              'disponible <= step si no hay coincidencia exacta)')
    parser.add_argument('-q', '--charge', type=float, default=None,
                         help='Carga q0 de cada colloide (override; si no se '
                              'da, se lee del log de cada simulacion)')
    parser.add_argument('--epsilon', type=float, default=None,
                         help='Permitividad del medio (override; si no se da, '
                              'se lee del log de cada simulacion)')
    parser.add_argument('--kappa', type=float, default=None,
                         help='Parametro de Debye kappa=1/lambda_D (override; '
                              'si no se da, se toma la mediana de "Kappa '
                              'calculado" entre los logs de cada simulacion)')
    parser.add_argument('--no-dh', action='store_true',
                         help='No graficar la curva de Debye-Huckel apantallada')
    parser.add_argument('-L', '--box-size', type=float, default=None,
                         help='Tamano de la caja periodica en x (override; si '
                              'no se da, se lee "System size" del log)')
    parser.add_argument('--n-images', type=int, default=2,
                         help='Cantidad de imagenes periodicas a cada lado '
                              'sumadas en la curva de Debye-Huckel (default: 2; '
                              'la serie converge rapido por el apantallamiento '
                              'exponencial)')
    parser.add_argument('--no-images', action='store_true',
                         help='Graficar Debye-Huckel de una carga aislada, sin '
                              'sumar las imagenes periodicas de la caja')
    parser.add_argument('--gaussian-sigma', type=float, default=None,
                         help='Si se da, agrega la curva de Debye-Huckel para '
                              'dos cargas gaussianas de este ancho (sigma, por '
                              'particula, en unidades de red) en vez de '
                              'puntuales. Ancho combinado del par = '
                              'sigma*sqrt(2)')
    parser.add_argument('-o', '--output', default=None,
                         help='Archivo de salida del grafico (default: '
                              'force_vs_distance_step<N>.png dentro del '
                              'directorio base, junto a las simulaciones)')
    parser.add_argument('--kernel', default=None,
                        help='kernel label shown in the plot title, e.g. "Hann-8"')
    parser.add_argument('--relerr', action='store_true',
                        help='bottom panel: relative error (%%) against the Gaussian-Gaussian '
                             'force, instead of the sim/theory ratio')
    parser.add_argument('--linear', action='store_true',
                         help='Usar escala lineal en vez de log-log (default: log-log)')
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    output = args.output
    if output is None:
        output = os.path.join(base_dir, f'force_vs_distance_step{args.step}.png')

    sim_dirs = find_sim_dirs(base_dir)
    if not sim_dirs:
        sys.exit(f'No se encontraron subdirectorios pos_X_Y_Z en {base_dir}')

    rows = []
    for sim_dir in sim_dirs:
        name = os.path.basename(sim_dir)
        force_csv = os.path.join(sim_dir, 'proceced_data', 'particle_force.csv')
        log_path = os.path.join(sim_dir, 'logs', 'output.txt')

        if not os.path.isfile(force_csv):
            print(f'[SKIP] {name}: no existe particle_force.csv (simulacion vieja, '
                  f're-correr con el binario que lo escribe)')
            continue

        force_step, forces = read_force_csv(force_csv, args.step)
        if force_step is None:
            print(f'[SKIP] {name}: particle_force.csv vacio')
            continue

        pos_step, positions = read_colloid_positions(sim_dir, force_step)
        if len(positions) < 2:
            print(f'[SKIP] {name}: no se encontraron 2 colloides en colloid_data')
            continue

        # La particula fija siempre queda a menor x que la movil (esta se
        # desplaza hacia +x desde su posicion inicial). Usar min/max en vez
        # de "mas cercano a --fixed-x" evita ambiguedades cuando ambas
        # particulas quedan a la misma distancia de --fixed-x (p.ej. cuando
        # la separacion es tan chica que el offset de grilla las deja
        # simetricas alrededor de esa posicion nominal).
        idxs = list(positions.keys())
        fixed_idx = min(idxs, key=lambda i: positions[i][0])
        moving_idx = max(idxs, key=lambda i: positions[i][0])

        if moving_idx not in forces:
            print(f'[SKIP] {name}: sin fuerza para index {moving_idx} en step {force_step}')
            continue

        x_fixed = positions[fixed_idx][0]
        x_moving = positions[moving_idx][0]
        r = x_moving - x_fixed
        Fx_sim = forces[moving_idx][0]

        rows.append(dict(name=name, r=r, Fx_sim=Fx_sim,
                          force_step=force_step, pos_step=pos_step,
                          phys=parse_physics(log_path)))

        print(f'{name}: step usado={force_step} (colloids step={pos_step}) '
              f'r={r:.4f} Fx_sim={Fx_sim:.6e}')

    if not rows:
        sys.exit('No se pudo procesar ninguna simulacion.')

    rows.sort(key=lambda d: d['r'])

    r_arr = np.array([row['r'] for row in rows])
    Fsim_arr = np.array([row['Fx_sim'] for row in rows])

    # Parametros fisicos: override de linea de comandos, si no la mediana de
    # lo leido en los logs (deberian ser iguales en todas las simulaciones;
    # la mediana ademas ignora outliers puntuales, p.ej. una corrida donde
    # "Kappa calculado" se imprimio como 0 por algun problema puntual de esa
    # corrida particular).
    def phys_value(key, override, default=1.0):
        if override is not None:
            return override
        values = [row['phys'].get(key) for row in rows]
        values = [v for v in values if v is not None]
        return float(np.median(values)) if values else default

    q0 = phys_value('q0', args.charge, 1.0)
    unit_charge = phys_value('unit_charge', None, 1.0)
    epsilon = phys_value('epsilon', args.epsilon, 1.0)
    kappa = phys_value('kappa', args.kappa, 0.0)
    box_l = phys_value('box_l', args.box_size, 0.0)

    use_images = (not args.no_images) and box_l > 0
    if not args.no_images and box_l <= 0:
        print('AVISO: no se pudo determinar el tamano de la caja (L) para las '
              'imagenes periodicas; se usa Debye-Huckel de carga aislada.')

    q_eff = q0 * unit_charge

    def coulomb_force(r):
        return (q_eff ** 2) / (4.0 * np.pi * epsilon * r ** 2)

    def dh_force_isolated(d):
        return (q_eff ** 2) / (4.0 * np.pi * epsilon) * np.exp(-kappa * d) * (kappa / d + 1.0 / d ** 2)

    # Ancho combinado del par: si cada particula es una gaussiana de ancho
    # sigma (independiente), la convolucion de las dos equivale a una unica
    # gaussiana de ancho sigma*sqrt(2) (ver docstring del modulo).
    sigma_pair = args.gaussian_sigma * np.sqrt(2.0) if args.gaussian_sigma else 0.0

    def gaussian_dh_potential_isolated(d):
        """Potencial DH exacto (por unidad de q2) de una carga gaussiana de
        ancho combinado sigma_pair, verificado contra la convolucion 3D
        directa y contra los limites sigma->0 (carga puntual) y kappa->0
        (Coulomb con carga gaussiana: erf(d/(sigma_pair*sqrt(2)))/d).
        NOTA: el factor exp(kappa^2*sigma_pair^2/2) es imprescindible (sale
        de completar la integral de Schwinger-parametrizacion) — sin el, la
        formula da el limite correcto en sigma->0 y kappa->0 pero es
        incorrecta para cualquier r intermedio (fue verificado numericamente
        contra la convolucion 3D directa antes de incluirlo)."""
        B = sigma_pair * kappa / np.sqrt(2.0)
        C = 1.0 / (sigma_pair * np.sqrt(2.0))
        T1 = np.exp(-kappa * d) * erfc(B - C * d)
        T2 = np.exp(kappa * d) * erfc(B + C * d)
        pref = (q_eff ** 2) / (8.0 * np.pi * epsilon)
        return pref * np.exp(kappa ** 2 * sigma_pair ** 2 / 2.0) * (T1 - T2) / d

    def gaussian_dh_force_isolated(d, h=1e-4):
        """F = -dphi/dr, por diferencias finitas centradas del potencial
        exacto de arriba (evita re-derivar a mano la derivada analitica,
        que resulto propensa a errores de signo en pasos intermedios)."""
        return -(gaussian_dh_potential_isolated(d + h)
                 - gaussian_dh_potential_isolated(d - h)) / (2.0 * h)

    def periodic_sum(force_isolated_fn, r):
        """Suma sobre las imagenes periodicas de la particula fija (en 0,
        +-L, +-2L, ...) si use_images esta activo; si no, solo la carga
        real. Ver docstring del modulo."""
        r = np.asarray(r, dtype=float)
        if not use_images:
            return force_isolated_fn(r)
        total = np.zeros_like(r)
        for n in range(-args.n_images, args.n_images + 1):
            d = np.abs(r - n * box_l)
            sign = np.where(r > n * box_l, 1.0, -1.0)
            total = total + sign * force_isolated_fn(d)
        return total

    def dh_force(r):
        return periodic_sum(dh_force_isolated, r)

    def gaussian_dh_force(r):
        return periodic_sum(gaussian_dh_force_isolated, r)

    F_theory = coulomb_force(r_arr)
    F_dh = dh_force(r_arr)
    show_dh = kappa > 0 and not args.no_dh
    show_gauss = kappa > 0 and args.gaussian_sigma is not None
    F_gauss = gaussian_dh_force(r_arr) if show_gauss else None

    print()
    print(f'Parametros teoricos: q0={q0:.6g}, unit_charge={unit_charge:.6g}, '
          f'epsilon={epsilon:.6g}, kappa={kappa:.6g}'
          + (f' (lambda_D={1.0/kappa:.4g})' if kappa > 0 else '')
          + (f', L={box_l:.6g}, imagenes periodicas: +-{args.n_images}' if use_images else '')
          + (f', gaussian_sigma={args.gaussian_sigma:.4g} (pair={sigma_pair:.4g})' if show_gauss else ''))
    header = f'{"sim":45s} {"r":>8s} {"F_sim":>14s} {"F_coulomb":>14s} {"F_sim/F_coulomb":>16s}'
    if show_dh:
        header += f' {"F_DH":>14s} {"F_sim/F_DH":>12s}'
    if show_gauss:
        header += f' {"F_gauss_DH":>14s} {"F_sim/F_gauss":>14s}'
    print(header)
    F_gauss_iter = F_gauss if show_gauss else [None] * len(rows)
    for row, Fth, Fdh, Fg in zip(rows, F_theory, F_dh, F_gauss_iter):
        ratio = row['Fx_sim'] / Fth if Fth != 0 else float('nan')
        line = (f'{row["name"][:45]:45s} {row["r"]:8.3f} {row["Fx_sim"]:14.6e} '
                f'{Fth:14.6e} {ratio:16.4f}')
        if show_dh:
            ratio_dh = row['Fx_sim'] / Fdh if Fdh != 0 else float('nan')
            line += f' {Fdh:14.6e} {ratio_dh:12.4f}'
        if show_gauss:
            ratio_g = row['Fx_sim'] / Fg if Fg != 0 else float('nan')
            line += f' {Fg:14.6e} {ratio_g:14.4f}'
        print(line)

    # Grafico: fuerza vs distancia (arriba) + ratio sim/teoria (abajo)
    fig, (ax, ax_ratio) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]})

    ax.plot(r_arr, Fsim_arr, 'o', color='tab:blue', markersize=8,
            label='Simulation (particle_force.csv)', zorder=3)

    r_fine = np.linspace(r_arr.min(), r_arr.max(), 300)
    ax.plot(r_fine, coulomb_force(r_fine), 'r--', linewidth=2,
            label=r'Pure Coulomb: $F = q_1 q_2 / (4\pi\epsilon r^2)$')
    if show_dh:
        dh_label = r'Debye-Hückel ($\lambda_D$=' + f'{1.0/kappa:.3g}'
        dh_label += (f', {args.n_images} periodic images/side, L={box_l:.0f})'
                     if use_images else ')')
        ax.plot(r_fine, dh_force(r_fine), color='tab:green', linestyle='-.', linewidth=2,
                label=dh_label)
    if show_gauss:
        gauss_label = (r'Gaussian Debye-Hückel ($\sigma$=' + f'{args.gaussian_sigma:.3g}'
                       + r', $\sigma_{pair}$=' + f'{sigma_pair:.3g})')
        ax.plot(r_fine, gaussian_dh_force(r_fine), color='tab:purple', linestyle=':', linewidth=2.5,
                label=gauss_label)

    if not args.linear:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_ylabel(r'Force $F_x$ on the moving particle', fontsize=12)
    # Segunda linea del titulo: las condiciones que distinguen una corrida de
    # otra (apantallamiento y kernel), que es lo que se compara entre figuras.
    cond = []
    if show_dh and kappa > 0:
        cond.append(r'$\lambda_D$ = ' + f'{1.0/kappa:.4g}')
    elif show_dh:
        cond.append('no salt')
    if args.kernel:
        cond.append(f'{args.kernel} kernel')
    if show_gauss:
        cond.append(r'$\sigma$ = ' + f'{args.gaussian_sigma:.5g}')
    title = f'Simulated force vs theory (step={args.step})'
    if cond:
        title += '\n' + ',   '.join(cond)
    title += f'\nq={q0:.4g}, epsilon={epsilon:.4g}'
    ax.set_title(title, fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10)

    # Cerca de r=L/2 la fuerza teorica (con imagenes) y la simulada se
    # cancelan casi exactamente: el cociente queda dominado por ruido de
    # punto flotante y dispara a valores enormes. Se oculta ese punto del
    # panel de cociente (sigue impreso en la tabla de arriba).
    def masked_ratio(Fsim, Fth):
        ratio = Fsim / Fth
        noise_floor = 1e-6 * np.max(np.abs(Fth))
        ratio[np.abs(Fth) < noise_floor] = np.nan
        return ratio

    if args.relerr:
        # Bottom panel: relative error against the Gaussian-Gaussian force,
        # which is the correct physical reference for the pair of charge clouds
        # the kernel deposits on the mesh.
        if not show_gauss:
            raise SystemExit('--relerr requires the Gaussian reference '
                             '(pass --gaussian-sigma)')
        err = 100.0 * (Fsim_arr - F_gauss) / F_gauss
        noise_floor = 1e-6 * np.max(np.abs(F_gauss))
        err[np.abs(F_gauss) < noise_floor] = np.nan
        ax_ratio.plot(r_arr, err, 'o-', color='tab:purple', markersize=6,
                      label=r'$100\,(F_{sim}-F_{gauss})/F_{gauss}$')
        ax_ratio.axhline(0.0, color='k', linestyle='--', linewidth=1.0, alpha=0.6)
        for lim, sty in ((1.0, ':'), (-1.0, ':')):
            ax_ratio.axhline(lim, color='tab:gray', linestyle=sty, linewidth=0.9, alpha=0.7)
        ax_ratio.set_ylabel('relative error vs\nGaussian-Gaussian (%)', fontsize=10)
        finite = err[np.isfinite(err)]
        if finite.size:
            m = max(2.0, 1.15 * np.max(np.abs(finite)))
            ax_ratio.set_ylim(-m, m)
        ax_ratio.legend(fontsize=9, loc='best')
    else:
        ratio_arr = masked_ratio(Fsim_arr, F_theory)
        ax_ratio.plot(r_arr, ratio_arr, 'o-', color='tab:red', label='F_sim / F_pure_Coulomb')
        if show_dh:
            ratio_dh_arr = masked_ratio(Fsim_arr, F_dh)
            ax_ratio.plot(r_arr, ratio_dh_arr, 's-', color='tab:green', label='F_sim / F_DH')
        if show_gauss:
            ratio_g_arr = masked_ratio(Fsim_arr, F_gauss)
            ax_ratio.plot(r_arr, ratio_g_arr, '^-', color='tab:purple', label='F_sim / F_gaussian_DH')
        ax_ratio.axhline(1.0, color='k', linestyle='--', linewidth=1.0, alpha=0.6)
        ax_ratio.set_ylabel('F_sim / F_theory', fontsize=11)
        if show_dh or show_gauss:
            ax_ratio.legend(fontsize=9)
    if not args.linear:
        ax_ratio.set_xscale('log')
    ax_ratio.set_xlabel('Separation r between particles (lattice units)', fontsize=12)
    ax_ratio.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches='tight')
    print(f'\nGrafico guardado en: {output}')


if __name__ == '__main__':
    main()
