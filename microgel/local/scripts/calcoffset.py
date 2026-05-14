#!/usr/bin/env python3
import numpy as np, math, time
import argparse
import sys
from functools import lru_cache

# Set command line parameters
argParser = argparse.ArgumentParser()
argParser.add_argument("-l", default='32', help="side length")

def phi_continuous(r):
    r = abs(r)
    if r < 1.0:
        return (1.0/8.0)*(3.0 - 2.0*r + math.sqrt(1.0 + 4.0*r - 4.0*r*r))
    elif r < 2.0:
        return (1.0/8.0)*(5.0 - 2.0*r - math.sqrt(-7.0 + 12.0*r - 4.0*r*r))
    else:
        return 0.0

@lru_cache(maxsize=None)
def phi_shifted_hat(k, s):
    # discrete sum r = -2,-1,0,1,2 of phi(r - s) * e^{-i k r}
    vals = [phi_continuous(r - s) * np.exp(-1j * k * r) for r in range(-2,3)]
    return sum(vals)

def compute_aL_tensorial(L, n_offsets=4):
    half = L//2
    # build list of k vectors (exclude k=0)
    k_list = []
    for nx in range(-half, half):
        for ny in range(-half, half):
            for nz in range(-half, half):
                if nx==0 and ny==0 and nz==0:
                    continue
                kx = 2.0*np.pi*nx / L
                ky = 2.0*np.pi*ny / L
                kz = 2.0*np.pi*nz / L
                k2 = kx*kx + ky*ky + kz*kz
                k_list.append((kx,ky,kz,k2))

    offsets = np.linspace(0.0, 1.0, n_offsets, endpoint=False)
    inv_sum = 0.0

    # sum over reciprocal vectors, average over offsets
    for kx,ky,kz,k2 in k_list:
        # trace of projector I - khat khat = 2 in 3D
        P_trace = 2.0
        sum_over_s = 0.0
        for sx in offsets:
            hx = phi_shifted_hat(kx, sx)
            for sy in offsets:
                hy = phi_shifted_hat(ky, sy)
                for sz in offsets:
                    hz = phi_shifted_hat(kz, sz)
                    dh = hx * hy * hz
                    sum_over_s += (abs(dh)**2) * P_trace
        avg = sum_over_s / (n_offsets**3)
        inv_sum += avg / k2

    inv_sum *= (6.0 * np.pi) / (L**3)
    aL = 1.0 / inv_sum
    return aL

if __name__ == "__main__":
    start = time.time()
    
    # Parse command line arguments
    try: 
        args = argParser.parse_args()
    except:
        sys.exit("Could not read command line parameters")

    L = int(args.l)
    if(L <= 0):
        sys.exit("Invalid length")

    val = compute_aL_tensorial(L, n_offsets=4)   # n_offsets=4..6 give good average
    elapsed = time.time() - start
    print(f"aL(L={L}) = {val:.15e}  (computed in {elapsed:.2f}s)")
