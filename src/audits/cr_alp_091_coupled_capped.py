#!/usr/bin/env python3
"""CR-ALP-091: coupled capped-comparator receiver, no new quadrature.

B_C=(r_L-AU)^* C^-1 (r_L-AU) <= B_L + U^*AU by the
variational characterization / infimal convolution of L and A.
If B_C <= M and z* B_C z <= b, then
B_C <= ((1-w) M^-1 + w zz*/b)^-1 for 0<=w<1.
The Sherman-Morrison formula gives a rational upper comparison.
All floating work proposes eta,w only. Exact interval Gram enclosure
and congruence decide acceptance.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,sys
import numpy as np
from cr_alp_086_terminal_moment import data,full,exact_spd
from cr_alp_085_obstruction_block import mm,transpose
from cr_alp_089_capped_terminal import P,add,sub,mul,div,sq,total,overlap
from cr_alp_088_terminal_response import ids,read
from cr_alp_074_balanced_selection import checked
from cr_alp_075_shifted_energy import S,floor,ceil,positive
from cr_alp_082_inertia_receiver import inertia
Q=10**24
def sqrt_interval(a):
    assert 0<=a[0]<=a[1]
    def r(x):
        y=x*Q*Q
        return isqrt(y.numerator//y.denominator)
    return F(r(a[0]),Q),F(r(a[1])+1,Q)
def gram(rows,n):
    # Each row contains intervals for an actual weighted sine coefficient.
    # Exact integer midpoint Gram and row-sum error (including all cross terms).
    G=[[0]*n for _ in range(n)];E=[0]*n
    for row in rows:
        l=[floor(x[0]*Q) for x in row];h=[ceil(x[1]*Q) for x in row]
        m=[(a+b)//2 for a,b in zip(l,h)];r=[max(x-a,b-x) for x,a,b in zip(m,l,h)]
        R=sum(r);M=sum(map(abs,m))
        for i in range(n):
            E[i]+=abs(m[i])*R+r[i]*(M+R)
            for j in range(i,n):G[i][j]+=m[i]*m[j]
    L=[[0]*n for _ in range(n)];U=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i,n):
            v=G[i][j]*S//(Q*Q);L[i][j]=L[j][i]=v;U[i][j]=U[j][i]=v
        e=ceil(F(E[i]*S,Q*Q))+n
        L[i][i]-=e;U[i][i]+=e
    return L,U
def correction_gram(root,k,J,response):
    tbank=checked(root/'cr_alp_057/trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    tt=tbank['sectors'][k];scale=int(tbank['scale']);assert tt['ids']==list(range(int(k),851,2))
    n=len(J);I=[j for j in ids(k) if j>85];r=read(response/k/'response.log');st=r['pole_st'];assert st[0]>F(3,1000)
    p={j:r[f'pole{j}'] for j in I};zero=[P(0)]*n
    t={i:[tuple(F(int(v),scale) for v in row[j]) for j in J] for i,row in zip(tt['ids'],tt['D'])}
    for i in I:t.setdefault(i,zero)
    tau=[div(total(mul(t[i][j],p[i]) for i in I),st) for j in range(n)]
    assert all(overlap(tuple(F(int(v),scale) for v in tt['Db'][0][j]),sub(P(0),tau[k0])) for k0,j in enumerate(J))
    bands=json.loads((Path(__file__).parent/'CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    gamma=F(72,25)
    def dg(b):return F(floor(F(b['gamma_lower'])*10**9),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049)
    assert all(dg(b)>=gamma for b in bands) and all(dg(b)>=cap for b in bands if b['start']>=2049)
    Acap=cap-gamma;assert Acap>0
    def a(i):return min(cap,dg(next(b for b in bands if b['start']<=i and (b['end'] is None or i<=b['end']))))-gamma
    rows=[]
    for i in I:
        if a(i)==0:continue
        w=sqrt_interval(P(a(i)))
        rows.append([mul(w,sub(t[i][j],mul(tau[j],p[i]))) for j in range(n)])
    tail=sub(st,total(sq(p[i]) for i in I));assert tail[1]>=0
    w=sqrt_interval(mul(P(Acap),(max(F(0),tail[0]),tail[1])))
    rows.append([mul(w,x) for x in tau])
    return gram(rows,n),tail
def fmat(A):return np.asarray([[float(x) for x in r] for r in A])
def main(root,prior,block,terminal,response,capped,out):
    c=json.loads(capped.read_text());assert c['retained_rank']==83
    sel=json.loads((terminal/'terminal_selection.json').read_text());report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        n=len(J);assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=n-1
        (GL,GU),tail=correction_gram(root,k,J,response)
        K=[[F(H[i][j]+GL[i][j],S) for j in range(n)] for i in range(n)]
        assert exact_spd(K)
        q=list(map(F,sel[k]['q']));z=[sum(a*b for a,b in zip(r,q)) for r in Z]
        b=F(c['sectors'][k]['inverse_interval'][1]);assert b>0 and c['sectors'][k]['classification']=='PASS_SCALAR_C_COMPARISON'
        # Cross-check actual source correction enclosure against scalar 090.
        lo=sum(z[i]*F(GL[i][j],S)*z[j] for i in range(n) for j in range(n))
        hi=sum(z[i]*F(GU[i][j],S)*z[j] for i in range(n) for j in range(n))
        al,ah=map(F,c['sectors'][k]['correction_energy']);assert max(lo,al)<=min(hi,ah)
        BP=mm(transpose(W),mm(A,W))
        def major(eta):
            return [[(1+eta)*(F(D[i],S)*(i==j)-DP[i][j])+(1+1/eta)*BP[i][j]+F(GU[i][j],S) for j in range(n)] for i in range(n)]
        # Efficient float search over the two weights. eta controls the
        # inherited affine block; w controls the actual capped scalar.
        h=fmat(K);zf=np.array(list(map(float,z)));bf=float(b);proposals=[]
        for eta in [F(i,100) for i in range(1,61)]+[F(i,10) for i in range(7,21)]:
            M=major(eta);mf=fmat(M);v=mf@zf;d=float(zf@v)
            if np.linalg.eigvalsh(mf)[0]<=0:continue
            for w0 in range(0,991,5):
                w=w0/1000;a=1-w
                U=mf/a-w*np.outer(v,v)/(a*(a*bf+w*d))
                lam=float(np.linalg.eigvalsh(h-U)[0])
                proposals.append((lam,eta,F(w0,1000)))
        assert proposals
        proposals.sort(reverse=True)
        best=n-1;trials=[]
        for score,eta,w in proposals[:8]:
            M=major(eta);assert exact_spd(M)
            v=[r[0] for r in mm(M,[[x] for x in z])];d=sum(x*y for x,y in zip(z,v));assert d>0
            a=1-w;assert 0<a<=1
            U=[[M[i][j]/a-w*v[i]*v[j]/(a*(a*b+w*d)) for j in range(n)] for i in range(n)]
            C=[[floor((K[i][j]-U[i][j])*S) for j in range(n)] for i in range(n)]
            assert all(C[i][j]==C[j][i] for i in range(n) for j in range(n))
            for i in range(n):C[i][i]-=n-1
            ii=inertia(C);best=max(best,ii[0])
            trials.append({'eta':str(eta),'weight':str(w),'inertia':ii,'numeric_min_eigenvalue':score})
            if best==n:break
        report[k]={'rank':best,'pass_full_sector':best==n,'trials':trials,
                   'correction_gram_lower_grid':[[str(x) for x in r] for r in GL],
                   'correction_gram_upper_grid':[[str(x) for x in r] for r in GU],
                   'scale':str(S),'terminal_inverse_upper':str(b),
                   'tail_pole_norm2_interval':list(map(str,tail))}
        print('SECTOR',k,'RANK',best,'BEST NUMERICAL MARGIN',proposals[0][0],flush=True)
    rank=sum(v['rank'] for v in report.values())
    ans={'status':'CONDITIONAL_COUPLED_CAPPED_RECEIVER','rank':rank,'conditional_index_upper':85-rank,
        'sectors':report,'premises':'Original Q85 source and Gamma floor identities, inherited physical and numerical enclosures; 090 scalar bounds independently replayed in workflow; no unrestricted RH claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
