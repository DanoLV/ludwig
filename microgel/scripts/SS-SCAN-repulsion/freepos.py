import sys
path,xm,xf = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
L=open(path).readlines()
pos=[]
for i,l in enumerate(L):
    s=l.split()
    if len(s)==3:
        try: v=[float(x) for x in s]
        except ValueError: continue
        if 1<v[0]<47 and 1<v[1]<47 and 1<v[2]<47: pos.append(i)
assert len(pos)==2, pos
# bloque: 32 enteros, luego a0, ah, luego la posicion -> isfixedr = 5to entero = pos-30
for b,x in zip(pos,(xm,xf)):
    L[b]=f"   {x:.15e}    2.400000000000000e+01    2.400000000000000e+01\n"
    L[b-30]=f"{0:24d}\n"
open(path,'w').writelines(L)
