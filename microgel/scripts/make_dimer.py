#!/usr/bin/env python3
"""
Genera ficheros de configuracion de colonias (config.cds.init.001-001) para un
DIMERO: dos particulas subgrid cargadas, enlazadas por un bond harmonico.

Escribe el estado en el orden EXACTO de colloid_state_write_ascii (src/colloid.c):
un valor por linea, primero los enteros (32) y luego los doubles.

Uso:
    python3 make_dimer.py
Produce:
    config.cds.init.001-001.dimer_pp   (cargas +1 / +1, repulsivo)
    config.cds.init.001-001.dimer_pm   (cargas +1 / -1, atractivo)
"""

NBOND_MAX = 6       # ver src/colloid.h
NPAD_INT_TOTAL = 32 # el bloque de enteros se rellena hasta 32 ints
NPAD_DBL_TOTAL = 64 # el bloque de doubles se rellena hasta 512 bytes (64 doubles)

IF = "{:>24d}\n"        # isformat  (entero)
SF = "{:>24.15e}\n"     # sformat   (double escalar)
VF = "{:>24.15e}{:>24.15e}{:>24.15e}\n"  # vformat (vector 3)


def write_particle(fp, index, r, q0, q1, bond_to, nbonds, a0, ah, al,
                    epsilon, inter_type=0):
    ints = []
    ints.append(index)      # index
    ints.append(1)          # rebuild
    ints.append(nbonds)     # nbonds
    ints.append(0)          # nangles
    ints.append(0)          # isfixedr
    ints.append(0)          # isfixedv
    ints.append(0)          # isfixedw
    ints.append(1)          # isfixeds
    ints.append(0)          # type (no usado)
    # bond[NBOND_MAX]
    bonds = [0] * NBOND_MAX
    for i, b in enumerate(bond_to):
        bonds[i] = b
    ints.extend(bonds)
    ints.append(0)          # rng
    ints.extend([0, 0, 0])  # isfixedrxyz
    ints.extend([0, 0, 0])  # isfixedvxyz
    ints.append(inter_type) # inter_type
    ints.append(0)          # ioversion
    ints.append(2)          # bc  (2 = COLLOID_BC_SUBGRID; ver colloid.h)
    ints.append(0)          # shape
    ints.append(0)          # active
    ints.append(0)          # magnetic
    ints.append(0)          # attr
    # padding entero hasta 32 ints
    while len(ints) < NPAD_INT_TOTAL:
        ints.append(0)
    for v in ints:
        fp.write(IF.format(v))

    # doubles
    fp.write(SF.format(a0))
    fp.write(SF.format(ah))
    fp.write(VF.format(*r))              # r
    fp.write(VF.format(0.0, 0.0, 0.0))   # v
    fp.write(VF.format(0.0, 0.0, 0.0))   # w
    fp.write(VF.format(0.0, 0.0, 0.0))   # s
    fp.write(VF.format(1.0, 0.0, 0.0))   # m
    fp.write(SF.format(0.0))             # b1
    fp.write(SF.format(0.0))             # b2
    fp.write(SF.format(0.0))             # c
    fp.write(SF.format(0.0))             # h
    fp.write(VF.format(0.0, 0.0, 0.0))   # dr
    fp.write(SF.format(0.0))             # deltaphi
    fp.write(SF.format(q0))              # q0
    fp.write(SF.format(q1))              # q1
    fp.write(SF.format(epsilon))         # epsilon
    fp.write(SF.format(0.0))             # deltaq0
    fp.write(SF.format(0.0))             # deltaq1
    fp.write(SF.format(0.0))             # sa
    fp.write(SF.format(0.0))             # saf
    fp.write(SF.format(al))              # al
    fp.write(SF.format(0.0))             # elabc[0]
    fp.write(SF.format(0.0))             # elabc[1]
    fp.write(SF.format(0.0))             # elabc[2]
    fp.write(SF.format(1.0))             # quat[0]
    fp.write(SF.format(0.0))             # quat[1..3]
    fp.write(SF.format(0.0))
    fp.write(SF.format(0.0))
    fp.write(SF.format(1.0))             # quatold[0]
    fp.write(SF.format(0.0))
    fp.write(SF.format(0.0))
    fp.write(SF.format(0.0))
    # padding double
    written_dbl = 2 + 3*5 + 4 + 3 + 1 + 3 + 4 + 3 + 4 + 4
    # (a0,ah=2)(r,v,w,s,m=15)(b1,b2,c,h=4)(dr=3)(deltaphi=1)
    # (q0,q1,eps=3)(dq0,dq1,sa,saf=4)(al=1)(elabc=3)(quat=4)(quatold=4)
    # recompute exactly:
    written_dbl = 0
    written_dbl += 2          # a0, ah
    written_dbl += 3*5        # r v w s m
    written_dbl += 4          # b1 b2 c h
    written_dbl += 3          # dr
    written_dbl += 1          # deltaphi
    written_dbl += 3          # q0 q1 epsilon
    written_dbl += 4          # deltaq0 deltaq1 sa saf
    written_dbl += 1          # al
    written_dbl += 3          # elabc
    written_dbl += 4          # quat
    written_dbl += 4          # quatold
    for _ in range(NPAD_DBL_TOTAL - written_dbl):
        fp.write(SF.format(0.0))


def make_dimer(fname, q_a, q_b, r0_sep, center, a0, ah, al, epsilon):
    """Dos particulas separadas r0_sep en x, centradas en 'center'."""
    cx, cy, cz = center
    r1 = (cx - 0.5 * r0_sep, cy, cz)
    r2 = (cx + 0.5 * r0_sep, cy, cz)
    with open(fname, "w") as fp:
        # numero de colonias (una linea, formato isformat)
        fp.write(IF.format(2))
        write_particle(fp, 1, r1, q_a, 0.0, bond_to=[2], nbonds=1,
                       a0=a0, ah=ah, al=al, epsilon=epsilon)
        write_particle(fp, 2, r2, q_b, 0.0, bond_to=[1], nbonds=1,
                       a0=a0, ah=ah, al=al, epsilon=epsilon)
    print("Escrito", fname, "  particulas en", r1, "y", r2)


if __name__ == "__main__":
    # Parametros geometricos (caja 16^3)
    CENTER = (8.0, 8.0, 8.0)
    R0_SEP = 2.0     # separacion inicial = r0 del bond
    A0 = 0.05
    AH = 0.05
    AL = 0.5
    EPS = 1.0e4      # debe coincidir con electrokinetics_epsilon

    make_dimer("config.cds.init.001-001.dimer_pp",
               q_a=+1.0, q_b=+1.0, r0_sep=R0_SEP, center=CENTER,
               a0=A0, ah=AH, al=AL, epsilon=EPS)
    make_dimer("config.cds.init.001-001.dimer_pm",
               q_a=+1.0, q_b=-1.0, r0_sep=R0_SEP, center=CENTER,
               a0=A0, ah=AH, al=AL, epsilon=EPS)
