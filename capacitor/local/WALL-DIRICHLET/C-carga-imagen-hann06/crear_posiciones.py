#!/usr/bin/env python3
"""Crea z_XX.XX/ con input, config de UNA particula (q=+1, x=y=16) y enlaces."""
import os, numpy as np
blk = open('bloque_particula.txt').readlines()
for z in np.arange(2.0, 17.51, 0.5):
    d = 'z_%05.2f' % z
    os.makedirs(os.path.join(d, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(d, 'proceced_data'), exist_ok=True)
    os.makedirs(os.path.join(d, 'colloid_data'), exist_ok=True)
    b = list(blk)
    b[0] = '%24d\n' % 1                        # indice
    b[34] = '    %.15e    %.15e    %.15e\n' % (16.0, 16.0, z)
    b[45] = '    %.15e\n' % 1.0                 # q0
    b[46] = '    %.15e\n' % 0.0                 # q1
    with open(os.path.join(d, 'config.cds.init.001-001'), 'w') as fh:
        fh.write('%22d\n' % 1); fh.writelines(b)
    os.system('cp input.base %s/input' % d)
    os.system('cp /home/danolv/Sim/ludwig/microgel/local/.petscrc %s/' % d)
    if not os.path.exists(os.path.join(d, 'Ludwig.exe')):
        os.symlink('/home/danolv/Sim/ludwig/src/Ludwig.exe', os.path.join(d, 'Ludwig.exe'))
print("posiciones creadas")
