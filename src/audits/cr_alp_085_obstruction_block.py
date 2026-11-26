#!/usr/bin/env python3
"""CR-ALP-085: physical affine inverse Gram on the six comparison obstructions.
Only source-combination and trial-vector selection use floating point.
All endpoint transport and final inertia calculations use rational arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_075_shifted_energy import lower_matrix,S,floor,ceil
from cr_alp_074_norm_receiver import up
from cr_alp_078_moment_receiver import sha
from cr_alp_082_inertia_receiver import inertia
Q=10**80;U=F(1101,200);g=F(72,25);s=F(21,8);shift=F(1,10**13)
def dec(x):
    sign='-' if x<0 else '';x=abs(x)
    return sign+str(x//Q)+'.'+str(x%Q).zfill(80)
def enc(l,h):
    assert l<=h
    return dec(floor(l*Q))+' '+dec(ceil(h*Q))
def combined(inputs,z):
    rows=[p.read_text().splitlines() for p in inputs];n=len(rows[0]);assert all(len(r)==n for r in rows)
    out=[];i=0
    def same():
        nonlocal i
        vals=[r[i] for r in rows];assert len(set(vals))==1
        out.append(vals[0]);i+=1;return int(vals[0])
    def pair(vals):
        l=h=F(0)
        for a,x in zip(z,vals):
            u,v=map(F,x);assert u<=v
            l+=min(a*u,a*v);h+=max(a*u,a*v)
        return enc(l,h)
    def plain():
        nonlocal i
        out.append(pair([r[i].split() for r in rows]));i+=1
    assert same() in (1,2);plain();plain();assert same()==425
    for _ in range(425):
        vals=[r[i].split() for r in rows];assert len({x[0] for x in vals})==1 and all(len(x)==5 for x in vals)
        out.append(vals[0][0]+' '+pair([x[1:3] for x in vals])+' '+pair([x[3:5] for x in vals]));i+=1
    assert same()==128
    for _ in range(128):
        vals=[r[i].split() for r in rows];assert len({tuple(x[:2]) for x in vals})==1 and all(len(x)==4 for x in vals)
        out.append(' '.join(vals[0][:2])+' '+pair([x[2:] for x in vals]));i+=1
    assert same()==25
    for _ in range(25):
        vals=[r[i] for r in rows];assert len(set(vals))==1 and int(vals[0].split()[-1])==40
        out.append(vals[0]);i+=1
        for _ in range(40):plain()
    assert i==n
    return '\n'.join(out)+'\n'
def transpose(A):return list(map(list,zip(*A)))
def mm(A,B):
    T=transpose(B)
    return [[sum(x*y for x,y in zip(a,b)) for b in T] for a in A]
def inv(A):
    n=len(A);B=[list(map(F,r))+[F(i==j) for j in range(n)] for i,r in enumerate(A)]
    for i in range(n):
        j=next(j for j in range(i,n) if B[j][i]);B[i],B[j]=B[j],B[i]
        d=B[i][i];B[i]=[x/d for x in B[i]]
        for j in range(n):
            if i!=j:
                d=B[j][i];B[j]=[x-d*y for x,y in zip(B[j],B[i])]
    return [r[n:] for r in B]
def basis(root,prior):
    old=json.loads(prior.read_text());assert old['status']=='CONDITIONAL_TARGETED_INERTIA_RECEIVER' and old['rank']==79
    for k,h in HASHES.items():
        p=old['sectors'][k];J=p['coordinate_order'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',h);H=lower_matrix(w,J)
        t=p['full_direct_inertia'];D=list(map(int,t['diagonal_grid']));assert len(D)==n and min(D)>0
        a=[F(p['inverse_bounds'][str(j)]) for j in J]
        assert sum(x*S/F(d) for x,d in zip(a,D))<=1
        assert inertia([[x-(D[i] if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(H)])==[n-3,3,0]
        yield k,J,H,D
def prepare(root,bank,prior,out):
    import numpy as np
    E={(e['parity'],e['coordinate_zero_based']):e for e in json.loads((bank/'source_bank.json').read_text())['entries']}
    assert len(E)==85
    report={}
    for k,J,H,D in basis(root,prior):
        d=np.array(D,dtype=float)/S;C=np.asarray(H,dtype=float)/S/np.sqrt(d[:,None]*d[None,:])
        lam,V=np.linalg.eigh(C);assert lam[2]<1<lam[3]
        Z=V[:,:3]/np.sqrt(d[:,None]);Z/=4*np.linalg.norm(Z,axis=0)
        # Rounded rational directions; the later projection uses their exact Gram.
        Z=[[F(int(round(v*10**12)),10**12) for v in row] for row in Z]
        assert all(sum(x*x for x in col)<F(1,10) for col in transpose(Z))
        paths=[bank/E[k,j]['file'] for j in J]
        assert all(sha(p)==E[k,j]['sha256'] for p,j in zip(paths,J))
        vecs={str(i):[r[i] for r in Z] for i in range(3)}
        for i in range(3):
            for j in range(i+1,3):vecs[f'{i}{j}']=[r[i]+r[j] for r in Z]
        logs={}
        for tag,z in vecs.items():
            f=out/f'input_{k}_{tag}.txt';f.write_text(combined(paths,z))
            eps=sum(abs(x)*(F(5,10**9)+F(E[k,j]['omitted_group_midpoint_norm_upper'])+shift) for j,x in zip(J,z))+shift
            logs[tag]={'input_sha256':sha(f),'error_upper':str(eps)}
        report[k]={'coordinates':J,'basis':[[str(x) for x in r] for r in Z],
                   'comparison_generalized_eigenvalues':list(map(float,lam[:6])),'inputs':logs}
    (out/'basis.json').write_text(json.dumps(report,indent=2)+'\n')
def affine(p,e):
    rows={}
    for line in p.read_text().splitlines():
        n,l,h,d,q=line.split();assert d=='/' and q=='1e15' and n not in rows
        rows[n]=(F(int(l),10**15),F(int(h),10**15));assert rows[n][0]<=rows[n][1]
    assert set(rows)=={'G','M','penalty',*(f'corr{i}' for i in range(7))}
    l,h=rows['G'];assert h>=0
    err=e*(2*up(h)+e);merr=s*err+shift
    ml,mh=rows['M']
    return ((g*max(F(0),l-err)-mh-merr)/(ELL*U),(g*(h+err)-ml+merr)/(ELL*U))
def receive(root,bank,prior,out):
    E={(e['parity'],e['coordinate_zero_based']):e for e in json.loads((bank/'source_bank.json').read_text())['entries']}
    b=json.loads((out/'basis.json').read_text());report={}
    for k,J,H,D in basis(root,prior):
        sec=b[k];assert sec['coordinates']==J;n=len(J)
        Z=[list(map(F,r)) for r in sec['basis']]
        M=mm(transpose(Z),[[F(d,S)*x for x in r] for d,r in zip(D,Z)])
        W=mm(inv(M),transpose([[F(d,S)*x for x in r] for d,r in zip(D,Z)]))
        paths=[bank/E[k,j]['file'] for j in J];assert all(sha(p)==E[k,j]['sha256'] for p,j in zip(paths,J))
        v={str(i):[r[i] for r in Z] for i in range(3)}
        for i in range(3):
            for j in range(i+1,3):v[f'{i}{j}']=[r[i]+r[j] for r in Z]
        A={}
        for tag,z in v.items():
            f=out/f'input_{k}_{tag}.txt';assert sha(f)==sec['inputs'][tag]['input_sha256']
            assert f.read_text()==combined(paths,z)
            e=sum(abs(x)*(F(5,10**9)+F(E[k,j]['omitted_group_midpoint_norm_upper'])+shift) for j,x in zip(J,z))+shift
            assert e==F(sec['inputs'][tag]['error_upper'])
            A[tag]=affine(out/f'mode_{k}_{tag}.log',e)
        ints=[[None]*3 for _ in range(3)]
        for i in range(3):
            ints[i][i]=A[str(i)]
            for j in range(i+1,3):
                l,h=A[f'{i}{j}'];a,c=A[str(i)];d,e=A[str(j)]
                ints[i][j]=ints[j][i]=((l-c-e)/2,(h-a-d)/2)
        # Interval symmetric matrix upper enclosure via absolute row error.
        B=[[(l+h)/2 for l,h in row] for row in ints]
        for i in range(3):B[i][i]+=sum((h-l)/2 for l,h in ints[i])
        DP=mm(transpose(W),mm(M,W));BP=mm(transpose(W),mm(B,W))
        # D-orthogonal projection: ||R x||_inv^2 <=
        # (1+eta)||R(1-P)x||_inv^2 + (1+1/eta)||R P x||_inv^2.
        # The first term is bounded by D on the complement; second by B.
        best=n-3;trials=[]
        for eta in sorted({F(i,1000) for i in range(1,101)}|{F(i,20) for i in range(3,41)}|{F(3),F(5),F(10)}):
            # Outward lower rounding of full exact rational comparison.
            T=[[F(H[i][j],S)-(1+eta)*(F(D[i],S)*(i==j)-DP[i][j])-(1+1/eta)*BP[i][j] for j in range(n)] for i in range(n)]
            C=[[floor(x*S) for x in row] for row in T]
            for i in range(n):C[i][i]-=n-1
            # Fast numerical screening only; exact inertia for candidates.
            import numpy as np
            e=np.linalg.eigvalsh(np.asarray(C,dtype=float)/S)
            score=int((e>0).sum())
            if score>=best:
                ii=inertia(C)
                trials.append({'eta':str(eta),'inertia':ii})
                best=max(best,ii[0])
        report[k]={'rank':best,'affine_gram_upper':[[str(x) for x in r] for r in B],
                   'affine_intervals':{t:list(map(str,v)) for t,v in A.items()},'trials':trials,
                   'generalized_eigenvalues':sec['comparison_generalized_eigenvalues']}
    rank=sum(x['rank'] for x in report.values())
    ans={'status':'CONDITIONAL_OBSTRUCTION_BLOCK_DIAGNOSTIC','rank':rank,'conditional_index_upper':85-rank,'sectors':report,
         'premises':'Original fixed-collar continuum and numerical enclosures, source proxy transport and shift allowance for actual validated combinations; no unrestricted RH claim.'}
    (out/'result.json').write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
if __name__=='__main__':
    mode,*args=sys.argv[1:]
    (prepare if mode=='prepare' else receive)(*map(Path,args))
