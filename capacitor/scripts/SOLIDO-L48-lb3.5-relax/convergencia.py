#!/usr/bin/env python3
"""Seguimiento de la relajacion: extrae de run.log, cada freq_statistics pasos,
la longitud minima/maxima de enlace, la energia armonica y la velocidad maxima
de los coloides. Uso: python3 convergencia.py [run.log]  (guarda convergencia.png/.csv)"""
import re, sys
import numpy as np

log = sys.argv[1] if len(sys.argv) > 1 else "run.log"
rows, cur = [], {}
num = r"([-+0-9.eE]+)"
for line in open(log, errors="replace"):
    m = re.search(r"Bond harmonic potential minimum r is:\s+" + num, line)
    if m: cur = {"rmin": float(m.group(1))}
    m = re.search(r"Bond harmonic potential maximum r is:\s+" + num, line)
    if m: cur["rmax"] = float(m.group(1))
    m = re.search(r"Bond harmonic potential energy is:\s+" + num, line)
    if m: cur["E"] = float(m.group(1))
    m = re.match(r"\[(minimum|maximum) \]\s+" + num + r"\s+" + num + r"\s+" + num, line)
    if m and cur:
        v = [abs(float(m.group(i))) for i in (2, 3, 4)]
        cur["vmax"] = max(cur.get("vmax", 0.0), max(v))
    m = re.match(r"\[fe\]\s+(\d+)\s", line)
    if m and "E" in cur:
        cur["step"] = int(m.group(1)); rows.append(cur); cur = {}

if not rows:
    sys.exit("Sin datos de estadisticas todavia en " + log)
a = np.array([[r["step"], r["rmin"], r["rmax"], r["E"], r.get("vmax", np.nan)] for r in rows])
np.savetxt("convergencia.csv", a, header="step rmin rmax E_bond vmax", fmt="%d %.7e %.7e %.7e %.7e")
print(f"{'paso':>7} {'rmin':>8} {'rmax':>8} {'E_enlace':>11} {'v_max':>10}")
for r in a[-15:]:
    print(f"{int(r[0]):7d} {r[1]:8.4f} {r[2]:8.4f} {r[3]:11.4e} {r[4]:10.3e}")
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(a[:, 0], a[:, 1], label="r min"); ax[0].plot(a[:, 0], a[:, 2], label="r max")
    ax[0].axhline(3.5, color="k", ls=":"); ax[0].set(xlabel="paso", ylabel="longitud de enlace"); ax[0].legend()
    ax[1].semilogy(a[:, 0], a[:, 3]); ax[1].set(xlabel="paso", ylabel="energia armonica")
    ax[2].semilogy(a[:, 0], a[:, 4]); ax[2].set(xlabel="paso", ylabel="|v| max coloides")
    fig.tight_layout(); fig.savefig("convergencia.png", dpi=120)
except Exception as e:
    print("(sin grafico:", e, ")")
