#!/usr/bin/env python3
"""CR-ALP-089: transfer terminal physical trials to capped Gamma comparator.
No new quadrature. Original correction is enclosed directly by the SHA-pinned
CR-ALP-057 sine intervals and exact Q85 projection, not by assuming midpoint
attachment. Scalar diagnostics alone do not promote finite source rank.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from math import isqrt
from cr_alp_088_terminal_response import ids,read,sq,sub,mul,div,sqrtup,SHIFT
from cr_alp_085_obstruction_block import basis
from cr_alp_086_terminal_moment import entrymap,input_for
from cr_alp_074_balanced_selection import checked,HASHES
from cr_alp_078_moment_receiver import sha
from cr_alp_075_shifted_energy import floor,ceil
def P(x):return F(x),F(x)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def total(xs):
    lo=hi=F(0)
    for l,h in xs:lo+=l;hi+=h
    return lo,hi
def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def outlo(x):return str(F(floor(x*10**18),10**18))
def outhi(x):return str(F(ceil(x*10**18),10**18))
def interval(x):return [outlo(x[0]),outhi(x[1])]
def main(root,bank,prior,block,response,out):
    tbank=checked(root/'cr_alp_057/trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    bands=json.loads((Path(__file__).parent/'CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    gamma,s=F(72,25),F(21,8)
    def dnew(b):return F(floor(F(b['gamma_lower'])*10**9),10**9)
    cap=next(dnew(b) for b in bands if b['start']==2049)
    assert all(dnew(b)>=gamma for b in bands)
    assert all(dnew(b)>=cap for b in bands if b['start']>=2049)
    def d(j):return min(cap,dnew(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))
    dmin=min(min(cap,dnew(b))-s for b in bands);dtail=cap-s;Acap=cap-gamma
    assert 0<dmin<=dtail and Acap>0
    b=json.loads((block/'basis.json').read_text())
    E=entrymap(bank);report={}
    for k,J,H,D in basis(root,prior):
        Z=[list(map(F,r)) for r in b[k]['basis']]
        report[k]={}
        trials={}
        for tag in ('0','1','2','01','02','12'):
            W=response/k/tag
            q=[F(str(i) in tag) for i in range(3)]
            raw,eps=input_for(bank,E,k,J,Z,q);f=block/f'input_{k}_{tag}.txt'
            p=b[k]['inputs'][tag]
            assert raw==f.read_text() and sha(f)==p['input_sha256'] and eps==F(p['error_upper'])
            zsrc=[sum(t*u for t,u in zip(row,q)) for row in Z]
            w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',HASHES[k]);scale=int(w['scale'])
            K=total(mul(P(x*y),tuple(F(int(v),scale) for v in w['K'][i][j])) for i,x in zip(J,zsrc) for j,y in zip(J,zsrc))
            I=ids(k);hi=[j for j in I if j>85];low=[j for j in I if j<=85]
            r=read(W/'response.log');assert set(r)=={'X2','G','E','Y2','Yb','pole_st',*(f'y{j}' for j in I),*(f'pole{j}' for j in I)}
            assert overlap(r['G'],read(block/f'mode_{k}_{tag}.log')['G'])
            diag=[line.split() for line in (W/'diag.txt').read_text().splitlines()]
            assert len(diag)==len(hi) and all(int(j)==i and F(l)==F(h)==gamma for (j,l,h),i in zip(diag,hi))
            rows=[line.split() for line in (W/'trial.txt').read_text().splitlines()]
            assert len(rows)==len(hi) and all(int(j)==i and F(l)==F(h) for (j,l,h),i in zip(rows,hi))
            x0={int(j):P(l) for j,l,h in rows}
            trials[tag]={j:v[0] for j,v in x0.items()}
            if len(tag)==2:assert all(v==trials[tag[0]][j]+trials[tag[1]][j] for j,v in trials[tag].items())
            poles={j:r[f'pole{j}'] for j in I};y={j:r[f'y{j}'] for j in I}
            c=div(total(mul(x0[j],poles[j]) for j in hi),total(sq(poles[j]) for j in hi))
            x={j:sub(x0[j],mul(c,poles[j])) for j in hi}
            X2=total(sq(x[j]) for j in hi);assert overlap(X2,r['X2'])
            # Original full trial U=Q85 t with t finite through 850.
            # This identity is inherited from the 056/057 pole-null construction.
            tt=tbank['sectors'][k];ts=int(tbank['scale'])
            assert tt['ids']==list(range(int(k),851,2))
            assert len(tt['D'])==len(tt['ids'])
            def contract(row):
                assert len(row)==len(J)
                return total(mul(P(a),tuple(F(int(v),ts) for v in row[j])) for j,a in zip(J,zsrc))
            t={j:contract(row) for j,row in zip(tt['ids'],tt['D'])}
            for j in hi:t.setdefault(j,P(0))
            st=r['pole_st'];assert st[0]>F(3,1000)
            tau=div(total(mul(t[j],poles[j]) for j in hi),st)
            assert overlap(contract(tt['Db'][0]),sub(P(0),tau))
            U={j:sub(t[j],mul(tau,poles[j])) for j in hi}
            # Subtract constant-cap pole multiple, killed by Q85, so AU=Q85 z0.
            A={j:d(j)-gamma for j in hi}
            z0={j:sub(mul(P(A[j]),t[j]),mul(P(A[j]-Acap),mul(tau,poles[j]))) for j in hi}
            AU=total(mul(U[j],z0[j]) for j in hi);assert AU[1]>=0
            KC=add(K,(max(F(0),AU[0]),AU[1]))
            delta={j:mul(P(A[j]),x[j]) for j in hi}
            change={j:add(delta[j],z0[j]) for j in hi}
            # Reuse the existing unprojected y and trial energy.
            energy=sub(sub(r['E'],total(mul(x[j],delta[j]) for j in hi)),mul(P(2),total(mul(x[j],z0[j]) for j in hi)))
            Y2=add(sub(r['Y2'],mul(P(2),total(mul(y[j],change[j]) for j in hi))),total(sq(change[j]) for j in hi))
            Yb=sub(r['Yb'],total(mul(change[j],poles[j]) for j in hi))
            B=sub(Yb,total(mul(y[j],poles[j]) for j in low))
            lam=div(B,st)
            Q=sub(sub(Y2,total(sq(y[j]) for j in low)),div(sq(B),st));assert Q[1]>=0
            v={j:sub(sub(y[j],change[j]),mul(lam,poles[j])) for j in hi}
            vv={j:sq(v[j]) for j in hi}
            tail=sub(Q,total(vv.values()));assert tail[1]>=0
            # E plus weighted squared residual is one quadratic form in
            # source coordinates, because the common finite trial is linear.
            weights={j:1/(d(j)-s)-1/dtail for j in hi}
            assert min(weights.values())>=0
            penalty=add(mul(Q,P(1/dtail)),total(mul(vv[j],P(weights[j])) for j in hi))
            assert penalty[1]>=0
            errE=sqrtup(r['X2'][1])*(2*eps+SHIFT)
            err=(eps+SHIFT)*sqrtup(1/dmin)
            def sqrtlo(v):
                v=max(F(0),v)*10**60
                return F(isqrt(v.numerator//v.denominator),10**30)
            bounds=(energy[0]-errE+max(F(0),sqrtlo(penalty[0])-err)**2,
                    energy[1]+errE+(sqrtup(penalty[1])+err)**2)
            assert bounds[0]<=bounds[1]
            report[k][tag]={'majorant':interval(bounds),'energy':interval(energy),'weighted_proxy_residual':interval(penalty),
                            'source_sha256':sha(f),'response_sha256':sha(W/'response.log'),'trial_sha256':sha(W/'trial.txt')}
        print(k,'PASS common trial and six quadratic majorant intervals',flush=True)
    out.write_text(json.dumps({'status':'CONDITIONAL_CAPPED_BLOCK_QUADRATIC','sectors':report},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
