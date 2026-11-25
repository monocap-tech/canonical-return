#!/usr/bin/env python3
"""Replay the frozen interval-matrix witnesses without any numerical optimizer.
This checks saved interval matrices, not the independent construction of entries.
"""
from pathlib import Path
import json
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
for parity in (1,2):
    tr=json.loads((ROOT/f'result_{parity}.json').read_text())
    s=int(tr['scale']); n=len(tr['corrected_energy']); assert n==44-parity
    A=tr['corrected_energy']; d=s//10**65; C=[];rad=[]
    for i,row in enumerate(A):
        R=[];r=0
        for j,(lo,hi) in enumerate(row):
            l=int(lo)//d;h=-((-int(hi))//d)
            if i==j:l-=10**65//4;h-=10**65//4
            m=(l+h)//2;R.append(m);r+=max(m-l,h-m)
        C.append(R);rad.append(r)
    for i in range(n):
        for j in range(i): assert C[i][j]==C[j][i]
        C[i][i]-=max(rad)
    prev=1
    for k in range(n):
        pivot=C[k][k];assert pivot>0
        for i in range(k+1,n):
            for j in range(i,n):
                x,r=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert r==0;C[i][j]=C[j][i]=x
        prev=pivot
    assert F(tr['energy_congruence_error'])<F(1,10**12)
    assert F(tr['band_residual_upper'])<F(1,10**96)
    print('finite trial margin 1/4 verified:',parity,n,'pivots')
    for prefix in ('extended','projected'):
        b=json.loads((ROOT/f'{prefix}_{parity}.json').read_text());assert int(b['scale'])==s
        z=[]
        for lo,hi in b['witness']['z']:
            assert lo==hi;z.append(int(lo))
        upper=0
        for i,row in enumerate(b['checked_band_matrix']):
            for j,(lo,hi) in enumerate(row):
                a=z[i]*z[j];upper+=a*(int(hi) if a>=0 else int(lo))
        u=F(upper,s**3);assert u < (-3 if parity==1 and prefix=='projected' else -2 if prefix=='projected' else -4)
        print(prefix,parity,'frozen rational witness upper:',float(u),'PASS')
print('PASS: saved interval matrices, no reconstruction of analytic inputs')
