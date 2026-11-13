#!/usr/bin/env python3
"""Genera un directorio de corrida con N particulas en posiciones dadas."""
import sys, os, shutil
BASE='/home/danolv/Sim/ludwig/microgel/local/'
REF=BASE+'FFT-Force-hann08-corr_yes-sigma_0.40/FFT-Force-hann08-corr_yes-sigma_0.40-n_3000-s_50-L_48-pos_10.00_24.00_24.00-e_0.0_0.0_0.0-kT_0.00001-rho_0.8-eta_0.05-eps_1.0e4-rhoel_7.76e-04-stencil_27/'
POS={1:(20.30,23.70,24.10), 2:(25.10,24.30,23.90), 3:(22.70,28.10,25.30), 4:(24.10,24.90,28.30)}

def make(tag, ids, nsteps=3, q=1.0):
    d=os.path.join(os.path.dirname(os.path.abspath(__file__)), tag)
    os.makedirs(d+'/proceced_data', exist_ok=True); os.makedirs(d+'/logs', exist_ok=True)
    shutil.copy('/home/danolv/Sim/ludwig/src/Ludwig.exe', d+'/Ludwig.exe')
    cfg=open(REF+'config.cds.init.001-001').read().split('\n')
    block=cfg[1:69]                                   # plantilla de una particula
    out=['%22d' % len(ids)]
    for k,pid in enumerate(ids, start=1):
        b=list(block)
        b[0]='%24d' % k                               # indice de particula
        b[34]=''.join('%24.15e' % v for v in POS[pid])  # posicion
        b[45]='%24.15e' % q                              # carga q0
        out += b
    open(d+'/config.cds.init.001-001','w').write('\n'.join(out)+'\n')
    inp=open(REF+'input').read().split('\n')
    rep={'N_cycles':f'N_cycles {nsteps}', 'electrokinetics_init_rho_el':'electrokinetics_init_rho_el 0.0',
         'vel_io_freq':'vel_io_freq 1000', 'freq_config':'freq_config 1000', 'colloid_io_freq':'colloid_io_freq 1000',
         'config_at_end':'config_at_end no', 'freq_statistics':'freq_statistics 1000'}
    new=[]
    for l in inp:
        k=l.split()[0] if l.split() else ''
        new.append(rep.get(k, l))
    new.append('pm_sr_pair_ref measured')
    open(d+'/input','w').write('\n'.join(new))
    return d

if __name__=='__main__':
    tag=sys.argv[1]; ids=[int(c) for c in sys.argv[2]]
    q=float(sys.argv[3]) if len(sys.argv)>3 else 1.0
    print(make(tag, ids, q=q))
