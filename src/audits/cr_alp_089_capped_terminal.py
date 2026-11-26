#!/usr/bin/env python3
"""CR-ALP-089: transfer terminal physical trials to capped Gamma comparator.
No new quadrature. Original correction is enclosed directly by the SHA-pinned
CR-ALP-057 sine intervals and exact Q85 projection, not by assuming midpoint
attachment. Scalar diagnostics alone do not promote finite source rank.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
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
def main(root,bank,prior,block,terminal,response,out):
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
    b=json.loads((block/'basis.json').read_text());sel=json.loads((terminal/'terminal_selection.json').read_text())
    prev=json.loads((response/'result.json').read_text());assert prev['retained_rank']==83
    E=entrymap(bank);report={}
    for k,J,H,D in basis(root,prior):
        p=sel[k];Z=[list(map(F,r)) for r in b[k]['basis']];q=list(map(F,p['q']))
        raw,eps=input_for(bank,E,k,J,Z,q);f=terminal/f'input_{k}_terminal.txt'
        assert raw==f.read_text() and sha(f)==p['input_sha256'] and eps==F(p['error_upper'])
        zsrc=[sum(t*u for t,u in zip(row,q)) for row in Z]
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',HASHES[k]);scale=int(w['scale'])
        K=total(mul(P(x*y),tuple(F(int(v),scale) for v in w['K'][i][j])) for i,x in zip(J,zsrc) for j,y in zip(J,zsrc))
        I=ids(k);hi=[j for j in I if j>85];low=[j for j in I if j<=85]
        r=read(response/k/'response.log');assert set(r)=={'X2','G','E','Y2','Yb','pole_st',*(f'y{j}' for j in I),*(f'pole{j}' for j in I)}
        assert overlap(r['G'],read(terminal/f'mode_{k}_terminal.log')['G'])
        assert sha(response/k/'response.log')==prev['sectors'][k]['response_log_sha256']
        assert sha(response/k/'trial.txt')==prev['sectors'][k]['trial_sha256']
        diag=[line.split() for line in (response/k/'diag.txt').read_text().splitlines()]
        assert len(diag)==len(hi) and all(int(j)==i and F(l)==F(h)==gamma for (j,l,h),i in zip(diag,hi))
        rows=[line.split() for line in (response/k/'trial.txt').read_text().splitlines()]
        assert len(rows)==len(hi) and all(int(j)==i and F(l)==F(h) for (j,l,h),i in zip(rows,hi))
        x0={int(j):P(l) for j,l,h in rows}
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
        # Equivalent weighted Parseval formula avoids subtracting interval
        # lower squares and adding upper squares at the same tail weight.
        penalty=Q[1]/dtail+sum(vv[j][1]*(1/(d(j)-s)-1/dtail) for j in hi)
        assert penalty>=0
        # True source and rational shift transport in the weighted inverse norm.
        errE=sqrtup(r['X2'][1])*(2*eps+SHIFT)
        lo=max(F(0),energy[0]-errE)
        up=energy[1]+errE+(sqrtup(penalty)+(eps+SHIFT)*sqrtup(1/dmin))**2
        assert lo<=up
        report[k]={'K_L':interval(K),'correction_energy':interval(AU),'K_C':interval(KC),
          'trial_energy':interval(energy),'Q_residual_norm2':interval(Q),'tail_norm2_upper':outhi(tail[1]),
          'weighted_residual_penalty_upper':outhi(penalty),'inverse_interval':[outlo(lo),outhi(up)],
          'margin_interval':[outlo(KC[0]-up),outhi(KC[1]-lo)],
          'classification':'PASS_SCALAR_C_COMPARISON' if up<KC[0] else 'REJECT_SCALAR_C_COMPARISON' if lo>KC[1] else 'UNDECIDED'}
        print(k,report[k],flush=True)
    ans={'status':'CONDITIONAL_CAPPED_TERMINAL_DIAGNOSTIC','retained_rank':83,'conditional_index_upper':2,
         'cap':str(cap),'dmin':str(dmin),'dtail':str(dtail),'sectors':report,
         'premises':'Inherited 056/057 exact correction attachment to Q85 finite sine precursor, 068 Gamma diagonal lower comparison, 088 physical and numerical enclosures; exact interval transfer; no full-sector rank claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
