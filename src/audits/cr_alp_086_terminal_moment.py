#!/usr/bin/env python3
"""CR-ALP-086: target one direction per remaining three-dimensional block.
Floating Schur/eigenvalue calculations propose directions and parameters only.
The final full-matrix lower comparison and inertia are exact.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import numpy as np
from cr_alp_085_obstruction_block import basis,combined,mm,transpose,inv,shift
from cr_alp_078_moment_receiver import sha
from cr_alp_075_shifted_energy import S,floor
from cr_alp_082_inertia_receiver import inertia
from cr_alp_080_selected_second_moment import receive as moments,selftest
def fmat(A):return np.array([[float(x) for x in row] for row in A])
def exact_spd(A):return inertia(A)==[len(A),0,0]
def data(root,prior,block):
    old=json.loads((block/'result.json').read_text());assert old['rank']==83
    saved=json.loads((block/'basis.json').read_text())
    for k,J,H,D in basis(root,prior):
        assert saved[k]['coordinates']==J
        Z=[list(map(F,r)) for r in saved[k]['basis']]
        A=[list(map(F,r)) for r in old['sectors'][k]['affine_gram_upper']]
        assert exact_spd(A)
        M=mm(transpose(Z),[[F(d,S)*x for x in r] for d,r in zip(D,Z)])
        assert exact_spd(M)
        W=mm(inv(M),transpose([[F(d,S)*x for x in r] for d,r in zip(D,Z)]))
        DP=mm(transpose(W),mm(M,W))
        yield k,J,H,D,Z,A,W,DP
def full(H,D,W,DP,B,eta):
    n=len(H);BP=mm(transpose(W),mm(B,W))
    C=[[floor(S*(F(H[i][j],S)-(1+eta)*(F(D[i],S)*(i==j)-DP[i][j])-(1+1/eta)*BP[i][j])) for j in range(n)] for i in range(n)]
    for i in range(n):C[i][i]-=n-1
    return C
def update(A,q,a,tau):
    v=[r[0] for r in mm(A,[[x] for x in q])]
    d=sum(x*y for x,y in zip(q,v));assert d>0 and a>=0
    return [[(1+tau)*(A[i][j]-v[i]*v[j]/d)+(1+1/tau)*a*v[i]*v[j]/d**2 for j in range(3)] for i in range(3)]
def threshold(H,D,Z,A):
    h=np.asarray(H,dtype=float)/S;d=np.asarray(D,dtype=float)/S;z=fmat(Z);a=fmat(A)
    Q,_=np.linalg.qr(np.sqrt(d[:,None])*z,mode='complete')
    C=Q[:,3:]/np.sqrt(d[:,None]);c=C.T@h@C
    cross=C.T@h@z;zz=z.T@h@z
    limit=np.linalg.eigvalsh(c)[0]-1;assert limit>0
    es=sorted({F(i,1000) for i in range(5,int(limit*1000),5)})
    ts=[F(i,100) for i in range(1,31)]+[F(i,20) for i in range(7,61)]
    L=np.linalg.cholesky(a);V=np.linalg.inv(L)
    effective=[]
    for e in es:
        ef=float(e);T=(zz-cross.T@np.linalg.solve(c-(1+ef)*np.eye(len(c)),cross))*ef/(1+ef)
        effective.append((e,T))
    proposals=[]
    for e,T in effective:
        ev,v=np.linalg.eigh(V@T@V.T);q=V.T@v[:,0];q=q/(2*np.linalg.norm(q))
        if ev[1]>1:proposals.append(q)
    # Include original coordinate directions as fallback.
    proposals.extend(np.eye(3)/2)
    best=None
    for qf in proposals:
        q=[F(int(round(x*10**12)),10**12) for x in qf];qf=np.array(list(map(float,q)))
        v=a@qf;d0=qf@v;assert d0>0
        AP=np.outer(v,v)/d0;w=v/d0
        rows=[]
        for eta,T in effective:
            for tau in ts:
                t=float(tau);M=T-(1+t)*(a-AP)
                if np.linalg.eigvalsh(M)[0]<=0:continue
                val=1/((1+1/t)*(w@np.linalg.solve(M,w)))
                if val>0:rows.append((val/d0,eta,tau,val))
        if rows:
            row=max(rows,key=lambda x:x[0])
            if best is None or row[0]>best[0]:best=(*row,q)
    assert best is not None,'No positive scalar receiving threshold found'
    return best
def entrymap(bank):
    E={(e['parity'],e['coordinate_zero_based']):e for e in json.loads((bank/'source_bank.json').read_text())['entries']}
    assert len(E)==85
    return E
def input_for(bank,E,k,J,Z,q):
    x=[sum(a*b for a,b in zip(r,q)) for r in Z]
    paths=[bank/E[k,j]['file'] for j in J]
    assert all(sha(p)==E[k,j]['sha256'] for p,j in zip(paths,J))
    eps=sum(abs(t)*(F(5,10**9)+F(E[k,j]['omitted_group_midpoint_norm_upper'])+shift) for t,j in zip(x,J))+shift
    return combined(paths,x),eps
def prepare(root,bank,prior,block,out):
    E=entrymap(bank);report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        ratio,eta,tau,val,q=threshold(H,D,Z,A)
        # Turn the floating threshold proposal into an exact positive comparison.
        # A tiny shrink is not assumed sufficient; back off until exact SPD.
        a=F(int(val*10**15),10**15)*F(99,100)
        for _ in range(20):
            assert a>0
            B=update(A,q,a,tau)
            if exact_spd(full(H,D,W,DP,B,eta)):break
            a*=F(9,10)
        else:raise AssertionError('No exact positive threshold found')
        raw,eps=input_for(bank,E,k,J,Z,q);f=out/f'input_{k}_terminal.txt';f.write_text(raw)
        base=sum(x*y for x,y in zip(q,[r[0] for r in mm(A,[[x] for x in q])]))
        report[k]={'q':list(map(str,q)),'eta':str(eta),'tau':str(tau),
           'exact_sufficient_inverse_ceiling':str(a),'prior_affine_ceiling':str(base),
           'threshold_fraction_of_affine':str(a/base),'input_sha256':sha(f),'error_upper':str(eps)}
        print('EXACT SUFFICIENT THRESHOLD',k,float(a), 'relative to affine',float(a/base),flush=True)
    (out/'terminal_selection.json').write_text(json.dumps(report,indent=2)+'\n')
def receive(root,bank,prior,block,out):
    E=entrymap(bank);sel=json.loads((out/'terminal_selection.json').read_text());report={}
    selftest()
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        p=sel[k];q=list(map(F,p['q']));eta,tau=F(p['eta']),F(p['tau'])
        bound=F(p['exact_sufficient_inverse_ceiling']);assert bound>0
        assert exact_spd(full(H,D,W,DP,update(A,q,bound,tau),eta))
        raw,eps=input_for(bank,E,k,J,Z,q);f=out/f'input_{k}_terminal.txt'
        assert raw==f.read_text() and sha(f)==p['input_sha256'] and eps==F(p['error_upper'])
        r=moments(out/f'mode_{k}_terminal.log',out/f'second_{k}_terminal.log',eps)
        base=sum(x*y for x,y in zip(q,[r[0] for r in mm(A,[[x] for x in q])]))
        assert base==F(p['prior_affine_ceiling'])
        a=min(F(r['inverse_upper']),base)
        n=len(H);best=n-1
        # Verify inherited rank 83 with a fixed saved eta in the unrefined block.
        assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=best
        C=full(H,D,W,DP,update(A,q,a,tau),eta)
        ii=inertia(C);best=max(best,ii[0])
        report[k]={'selection':p,'moments':r,'retained_inverse_upper':str(a),
           'receiving_inertia':ii,'rank':best,'pass_full_sector':best==n}
    rank=sum(v['rank'] for v in report.values())
    ans={'status':'CONDITIONAL_TERMINAL_MOMENT_RECEIVER','rank':rank,'conditional_index_upper':85-rank,
         'sectors':report,'premises':'Inherited source, continuum, physical numerical enclosure and documented shift transport premises; exact scalar threshold and final rational full matrix inertia. No unrestricted RH claim.'}
    (out/'result.json').write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
if __name__=='__main__':
    mode,*args=sys.argv[1:]
    (prepare if mode=='prepare' else receive)(*map(Path,args))
