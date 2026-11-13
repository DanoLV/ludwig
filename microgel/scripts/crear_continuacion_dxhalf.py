#!/usr/bin/env python3
"""
Continuation of the dx = 1/2 two-charge runs SLIDES-LD8-hann16-corr_no-dxhalf/r_05, r_06, r_07 up to step 20000
(2026-11-10, user's request: continue the runs already done, no new runs). Recipe of
capacitor/local/WALL-DIRICHLET/crear_continuacion.py (Ludwig restarts from dist/psi/qsi/rho/vel of step N0 and
config.cds000N0, "Re-starting simulation at step N0").
  r_05: from step 12000 (end of its run), r_07: from 12000, r_06: from 10000 (its run died at 10775 by GPU out of memory;
  10000 is its last dump).
The originals are NOT touched: SLIDES-LD8-hann16-corr_no-dxhalf-cont/r_XX gets copies of input, config.cds.init,
Ludwig.exe (same exe, 4f0c20e7), the *-metadata files and the files of step N0 only (not the earlier dumps), and the
proceced_data/*.csv truncated to step <= N0 (Ludwig appends to existing csv files; r_06 had rows up to 10775).
Input: N_start N0, N_cycles 20000 - N0 (freq_config 2000 unchanged).
correr_cont.sh does not delete anything; log in run_cont.log; [OK] requires "finished normally" and the restart line.
Usage: python3 crear_continuacion_dxhalf.py   and then   cd SLIDES-LD8-hann16-corr_no-dxhalf-cont && ./correr_cont.sh 3 > bateria.log
"""
import os, re, glob, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SLIDES-LD8-hann16-corr_no-dxhalf")
BAT = SRC + "-cont"
NEND = 20000
N0 = {"r_05": 12000, "r_06": 10000, "r_07": 12000}

RUNNER = """\
#!/bin/bash
# Resumes every r_* from its step N0 (N in parallel). Deletes nothing; log in run_cont.log.
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
N=${1:-3}
uno() {
  cd "$1" || return 1
  if grep -aqs "finished normally" run_cont.log; then echo "[SKIP] $1"; return 0; fi
  n0=$(awk '/^N_start/ {print $2}' input)
  stdbuf -o0 ./Ludwig.exe > run_cont.log 2>&1
  grep -aq "finished normally" run_cont.log && grep -aq "Re-starting simulation at step $n0" run_cont.log \\
    && echo "[OK] $1" || echo "[FALLO] $1"
}
export -f uno
ls -d r_* | xargs -P "$N" -I{} bash -c 'uno {}'
echo BATERIA_LISTA
"""

def sub(t, pat, new):
    t2, n = re.subn(pat, new, t, flags=re.M)
    assert n == 1, pat
    return t2

if __name__ == "__main__":
    os.makedirs(BAT, exist_ok=True)
    for r, n0 in N0.items():
        src, dst = os.path.join(SRC, r), os.path.join(BAT, r)
        if os.path.exists(dst):
            print("%s: ya existe, no se modifica" % dst); continue
        files = glob.glob(os.path.join(src, "*-%09d.*" % n0)) + [os.path.join(src, "config.cds%08d.001-001" % n0)]
        for f in files:
            assert os.path.exists(f), f
        assert any("dist-" in f for f in files), src
        for sd in ("logs", "colloid_data", "proceced_data", "plots", "scripts"):
            os.makedirs(os.path.join(dst, sd), exist_ok=True)
        for f in files + glob.glob(os.path.join(src, "*-metadata.*")) + [os.path.join(src, x) for x in ("config.cds.init.001-001", "Ludwig.exe")]:
            shutil.copy2(f, dst)
        for c in glob.glob(os.path.join(src, "proceced_data", "*.csv")):
            out = []
            for l in open(c):
                v = l.split(";")[0].strip()
                if l.startswith("#") or not v.isdigit() or int(v) <= n0:
                    out.append(l)
            open(os.path.join(dst, "proceced_data", os.path.basename(c)), "w").writelines(out)
        t = open(os.path.join(src, "input")).read()
        t = sub(t, r"^N_start 0$", "N_start %d            # continuacion desde el paso %d (crear_continuacion_dxhalf.py, 2026-11-10)" % (n0, n0))
        t = sub(t, r"^N_cycles .*$", "N_cycles %d            # pasos adicionales: termina en el paso %d (t = %d en pasos de dx = 1)" % (NEND - n0, NEND, NEND // 4))
        open(os.path.join(dst, "input"), "w").write(t)
        print("%s -> %s (desde el paso %d, %d archivos del paso)" % (src, dst, n0, len(files)))
    open(os.path.join(BAT, "correr_cont.sh"), "w").write(RUNNER); os.chmod(os.path.join(BAT, "correr_cont.sh"), 0o755)
