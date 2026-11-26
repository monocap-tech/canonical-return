#!/usr/bin/env python3
"""CR-ALP-090: choose a capped Gamma trial from saved physical coefficients.
Floating computation is candidate selection only. The original response
evaluator measures the trial with gamma I+P, then 089 transfers exactly
to C=L+A, including original correction and infinite tail.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import numpy as np
from cr_alp_088_terminal_response import ids,read
from cr_alp_075_shifted_energy import floor
def main(k,source,old,work):
    I=ids(k);J=[j for j in I if j>85]
    r=read(old/'coeff.log');log=read(old/'response.log')
    assert set(r)=={'G',*(f'r{j}' for j in I)}
    def mid(v):return float(sum(v)/2)
    rhs=np.array([mid(r[f'r{j}']) for j in J])
    c=np.array([mid(log[f'pole{j}']) for j in J])
    st=mid(log['pole_st']);assert st>.003
    lines=source.read_text().splitlines();assert lines[0]==k and int(lines[3])==425
    t={}
    for row in lines[4:429]:
        j,l,h,a,b=row.split();t[int(j)]=mid((F(a),F(b)))
    assert sorted(t)==list(range(int(k),851,2))
    t=np.array([t.get(j,0.) for j in J])
    tau=(t@c)/st
    bands=json.loads((Path(__file__).parent/'CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    gamma=F(72,25)
    def dnew(b):return F(floor(F(b['gamma_lower'])*10**9),10**9)
    cap=next(dnew(b) for b in bands if b['start']==2049)
    diag=np.array([float(min(cap,dnew(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))) for j in J])
    A=diag-float(gamma);Acap=float(cap-gamma)
    z0=A*t-(A-Acap)*tau*c
    rhs=rhs-z0
    a=1.1;v=np.array(J)*np.pi/(2*a)
    sh=np.log([2,3,4,5,7,8,9]);wt=np.log([2,3,2,5,7,2,3])/np.sqrt([2,3,4,5,7,8,9])
    sn=np.sin(v[:,None]*sh);cs=np.cos(v[:,None]*sh);T=sn@wt
    num=v[None,:]*T[:,None]-v[:,None]*T[None,:];den=v[None,:]**2-v[:,None]**2
    np.fill_diagonal(den,1)
    P=-2/a*num/den;np.fill_diagonal(P,-((2*a-sh)/a*cs+sn/(a*v[:,None]))@wt)
    mat=P+np.diag(diag)
    z=np.linalg.solve(mat,np.column_stack([rhs,c]))
    x=z[:,0]-z[:,1]*(c@z[:,0])/(c@z[:,1])
    assert np.isfinite(x).all() and np.linalg.norm(x)<10
    def dec(n):return ('-' if n<0 else '')+f'{abs(n)//10**12}.{abs(n)%10**12:012d}'
    (work/'trial.txt').write_text(''.join(f'{j} {dec(t)} {dec(t)}\n' for j,t in zip(J,[int(round(y*10**12)) for y in x])))
    # The physical log is deliberately evaluated against L; 089 performs the exact C transfer.
    (work/'diag.txt').write_text(''.join(f'{j} 2.88 2.88\n' for j in J))
    print('NUMERICAL C TRIAL',k,'finite variational energy',float(rhs@x),'norm',float(np.linalg.norm(x)),flush=True)
if __name__=='__main__':main(sys.argv[1],*map(Path,sys.argv[2:]))
