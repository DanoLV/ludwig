#!/usr/bin/env python3
"""
Continuation of the 4 x 7 viscosity runs of TEST-HIDRO-dx2 (2026-11-16, user's request: run 3000 more steps, check
for steady state, and if not reached run 3000 more; then redo the Delta mu curve with the steady state).

Variants (source of step N0):   old       ./hidro_si*              (exe 4f0c20e7, unweighted centred advection)
                                centred   ./advfix/hidro_si*       (exe 221f3229)
                                flla      ./advfix_flla/hidro_si*
                                flla_mass ./advfix_flla_mass/hidro_si*
Destination: cont<NEND>/<variant>/<run>/. Recipe of crear_continuacion_dxhalf.py: copies of input, config.cds.init,
the SAME Ludwig.exe as the source run, the *-metadata files, the files of step N0 only, and proceced_data/*.csv
truncated to step <= N0. Input: N_start N0, N_cycles NEND - N0 (nothing else changed). Ludwig restarts from
dist/psi/qsi/rho/vel of step N0 ("Re-starting simulation at step N0").
Note: flla/flla_mass store the LB link flux after each collision; right after a restart there is none yet, so the
first NP step of the continuation has no advection (np_mflux_valid_ = 0, src/nernst_planck.c), as at step 1.
The originals are not touched; existing destinations are skipped.

Usage: python3 scripts/crear_cont_visc.py N0 NEND [run ...]
   e.g. python3 scripts/crear_cont_visc.py 3000 6000            (all 28)
        python3 scripts/crear_cont_visc.py 6000 9000 flla/hidro_si_eta2.0 ...   (sources taken from cont6000/)
"""
import glob, os, re, shutil, sys

BASE = "/home/danolv/Sim/ludwig/microgel/local/TEST-HIDRO-dx2"   # scripts/ is a symlink: absolute paths
RUNS = ("hidro_si_eta0.02", "hidro_si", "hidro_si_eta0.1", "hidro_si_eta0.2", "hidro_si_eta0.5",
        "hidro_si_eta1.0", "hidro_si_eta2.0")
VARIANTS = {"old": BASE, "centred": BASE + "/advfix", "flla": BASE + "/advfix_flla",
            "flla_mass": BASE + "/advfix_flla_mass"}


def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2


n0, nend = int(sys.argv[1]), int(sys.argv[2])
sel = sys.argv[3:] or ["%s/%s" % (v, r) for v in VARIANTS for r in RUNS]
for vr in sel:
    v, r = vr.split("/")
    src = os.path.join(VARIANTS[v] if n0 == 3000 else os.path.join(BASE, "cont%d" % n0, v), r)
    dst = os.path.join(BASE, "cont%d" % nend, v, r)
    if os.path.exists(dst):
        print("%s: ya existe, no se modifica" % dst); continue
    files = glob.glob(os.path.join(src, "*-%09d.*" % n0)) + [os.path.join(src, "config.cds%08d.001-001" % n0)]
    for f in files:
        assert os.path.exists(f), f
    assert any("dist-" in f for f in files), src
    for sd in ("logs", "colloid_data", "proceced_data", "plots"):
        os.makedirs(os.path.join(dst, sd))
    for f in files + glob.glob(os.path.join(src, "*-metadata.*")) + \
            [os.path.join(src, x) for x in ("config.cds.init.001-001", "Ludwig.exe")]:
        shutil.copy2(f, dst)
    for c in glob.glob(os.path.join(src, "proceced_data", "*.csv")):
        out = []
        for l in open(c, errors="replace"):
            s = l.split(";")[0].strip()
            if l.startswith("#") or not s.isdigit() or int(s) <= n0:
                out.append(l)
        open(os.path.join(dst, "proceced_data", os.path.basename(c)), "w").writelines(out)
    t = open(os.path.join(src, "input")).read()
    t = sub(t, r"^N_start .*$", "N_start %d             # continuacion desde el paso %d (crear_cont_visc.py, 2026-11-16)"
            % (n0, n0))
    t = sub(t, r"^N_cycles .*$", "N_cycles %d            # pasos adicionales: termina en el paso %d" % (nend - n0, nend))
    open(os.path.join(dst, "input"), "w").write(t)
    print("%s -> %s (%d archivos del paso %d)" % (src, dst, len(files), n0))
