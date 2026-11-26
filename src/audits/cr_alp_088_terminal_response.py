#!/usr/bin/env python3
"""CR-ALP-088: scalar exact response bounds for L=gamma I+P_Q.
Numerical finite trial selection only. Full physical residual norm, including
the infinite tail, is bounded by orthogonal projection and exact endpoints.
Comparison rejection is not a negative direction of the original Weil form.
"""
from fractions import Fraction as F
from pathlib import Path
import sys,json
from cr_alp_080_selected_second_moment import read,sqrtup
from cr_alp_085_obstruction_block import basis,mm,transpose
from cr_alp_086_terminal_moment import entrymap,input_for
from cr_alp_078_moment_receiver import sha
from cr_alp_074_balanced_selection import checked,HASHES,ELL
SHIFT=F(1,10**13)
def ids(k):return list(range(int(k),2049,2))
def trial(k,work):
    import numpy as np
    I=ids(k);J=[j for j in I if j>85]
    rows=read(work/'coeff.log');assert set(rows)=={'G',*(f'r{j}' for j in I)}
    r=np.array([float(sum(rows[f'r{j}'])/2) for j in J])
    a=1.1;v=np.array(J)*np.pi/(2*a)
    sh=np.log([2,3,4,5,7,8,9]);wt=np.log([2,3,2,5,7,2,3])/np.sqrt([2,3,4,5,7,8,9])
    sn=np.sin(v[:,None]*sh);cs=np.cos(v[:,None]*sh);T=sn@wt
    num=v[None,:]*T[:,None]-v[:,None]*T[None,:];den=v[None,:]**2-v[:,None]**2
    np.fill_diagonal(den,1)
    P=-2/a*num/den;np.fill_diagonal(P,-((2*a-sh)/a*cs+sn/(a*v[:,None]))@wt)
    A=P+2.88*np.eye(len(J));c=v/(v*v+.25)
    z=np.linalg.solve(A,np.column_stack([r,c]));x=z[:,0]-z[:,1]*(c@z[:,0])/(c@z[:,1])
    assert np.isfinite(x).all() and np.linalg.norm(x)<10
    def dec(n):return ('-' if n<0 else '')+f'{abs(n)//10**12}.{abs(n)%10**12:012d}'
    (work/'trial.txt').write_text(''.join(f'{j} {dec(t)} {dec(t)}\n' for j,t in zip(J,[int(round(y*10**12)) for y in x])))
    (work/'diag.txt').write_text(''.join(f'{j} 2.88 2.88\n' for j in J))
    print('NUMERICAL TRIAL',k,float(r@x),flush=True)
def sq(a):
    return (F(0) if a[0]<=0<=a[1] else min(x*x for x in a),max(x*x for x in a))
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    z=[x*y for x in a for y in b];return min(z),max(z)
def div(a,b):
    assert b[0]>0;return mul(a,(1/b[1],1/b[0]))
def receive(root,bank,prior,block,terminal,work,out):
    b=json.loads((block/'basis.json').read_text());sel=json.loads((terminal/'terminal_selection.json').read_text())
    E=entrymap(bank);report={}
    for k,J,H,D in basis(root,prior):
        p=sel[k];Z=[list(map(F,r)) for r in b[k]['basis']];q=list(map(F,p['q']))
        raw,eps=input_for(bank,E,k,J,Z,q);f=terminal/f'input_{k}_terminal.txt'
        assert raw==f.read_text() and sha(f)==p['input_sha256'] and eps==F(p['error_upper'])
        x=[sum(t*u for t,u in zip(row,q)) for row in Z]
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',HASHES[k]);scale=int(w['scale'])
        kl=kh=F(0)
        for i,xi in zip(J,x):
            for j,xj in zip(J,x):
                l,h=map(int,w['K'][i][j]);a,c=F(l,scale)*xi*xj,F(h,scale)*xi*xj
                kl+=min(a,c);kh+=max(a,c)
        r=read(work/k/'response.log');I=ids(k)
        assert set(r)=={'X2','G','E','Y2','Yb','pole_st',*(f'y{j}' for j in I),*(f'pole{j}' for j in I)}
        old=read(terminal/f'mode_{k}_terminal.log')['G']
        assert max(old[0],r['G'][0])<=min(old[1],r['G'][1])
        assert r['pole_st'][0]>F(3,1000) and 0<=r['X2'][1]<100
        Y,B=r['Y2'],r['Yb']
        for j in I:
            if j>85:break
            Y=sub(Y,sq(r[f'y{j}']));B=sub(B,mul(r[f'y{j}'],r[f'pole{j}']))
        Y=sub(Y,div(sq(B),r['pole_st']));assert Y[1]>=0
        # Trial has exact finite pole-null projection. Error in E due to
        # source proxy is 2||x||eps, and rational-shift action <= SHIFT.
        e=sqrtup(r['X2'][1])*(2*eps+SHIFT)
        lo=max(F(0),r['E'][0]-e)
        hi=r['E'][1]+e+(sqrtup(Y[1])+eps+SHIFT)**2/ELL
        assert lo<=hi
        report[k]={'corrected_energy':[str(kl),str(kh)],'inverse_interval':[str(lo),str(hi)],
          'margin_interval':[str(kl-hi),str(kh-lo)],'projected_response_error_upper':str(sqrtup(Y[1])+eps+SHIFT),
          'classification':'PASS_SCALAR_L_COMPARISON' if hi<kl else 'REJECT_SCALAR_L_COMPARISON' if lo>kh else 'UNDECIDED',
          'trial_sha256':sha(work/k/'trial.txt'),'response_log_sha256':sha(work/k/'response.log')}
    ans={'status':'CONDITIONAL_TERMINAL_RESPONSE_DIAGNOSTIC','retained_rank':83,'conditional_index_upper':2,'sectors':report,
         'premises':'Inherited physical source and numerical error contracts; exact finite pole-null trial and full infinite residual transport. Scalar tests alone cannot promote rank; rejection only concerns L comparator.'}
    out.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
if __name__=='__main__':
    if sys.argv[1]=='trial':trial(sys.argv[2],Path(sys.argv[3]))
    elif sys.argv[1]=='receive':receive(*map(Path,sys.argv[2:]))
    else:raise ValueError(sys.argv[1])
