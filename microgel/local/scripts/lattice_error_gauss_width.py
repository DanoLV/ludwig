import numpy as np
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'lattice_error_gauss_vs_hann.py')).read().split('print("Error de la red')[0])
print("Kernel gaussiano: error de la red contra SU PROPIO continuo, en r = sigma y r = 2 sigma")
print(f"{'sigma':>7} {'soporte (+-4s)':>15} {'nodos 3D':>9} {'err r=sigma':>12} {'err r=2sigma':>13}")
for s in (0.72, 1.084, 1.446, 2.0, 3.0, 4.0):
    sig = s
    T = int(np.ceil(4*s))
    kern = lambda x, T=T, s=s: np.where(np.abs(x) <= T, np.exp(-x*x/(2*s*s)), 0.0)
    ftg = lambda k, s=s: np.exp(-k*k*s*s/2)
    out = []
    for r in (s, 2*s):
        fg = lattice_force(kern, 24.0, 24.0 + r)
        cg = continuum_force(ftg, r)
        out.append(100*(fg/cg - 1))
    print(f"{s:7.3f} {2*T:15d} {(2*T)**3:9d} {out[0]:+11.2f}% {out[1]:+12.2f}%")
