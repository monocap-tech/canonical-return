#!/usr/bin/env python3
"""CR-ALP-078: complete coordinate moment bank and exact diagonal receiver.
Numerical optimization proposes rational diagonals only. Exact integer
pivots and outward rational sums decide acceptance. Inherited premises apply.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,hashlib
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_074_norm_receiver import up
from cr_alp_075_shifted_energy import lower_matrix,positive,S,ceil
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def moments(p,e):
    rows={}
    for line in p.read_text().splitlines():
        n,l,h,s,d=line.split();assert s=='/' and d=='1e15' and n not in rows
        l,h=F(int(l),10**15),F(int(h),10**15);assert l<=h
        rows[n]=(l,h)
    assert set(rows)=={'G','M','penalty',*(f'corr{i}' for i in range(7))}
    G=rows['G'][1];assert 0<=G<100
    eps=F(5,10**9)+F(e['omitted_group_midpoint_norm_upper'])+F(1,10**13)
    err=eps*(2*up(G)+eps)
    b=up(G)+eps
    Mlo=rows['M'][0]-F(21,8)*err-F(1,10**13)
    A=(F(72,25)*(G+err)-Mlo)/(ELL*F(1101,200));assert A>=0
    return {'norm_upper':str(b),'scalar_inverse_upper':str(b*b/ELL),
            'affine_inverse_upper':str(A),'inverse_upper':str(min(A,b*b/ELL)),
            'log_sha256':sha(p),'input_sha256':e['sha256']}
def candidate(H,ai):
    import numpy as np
    A=np.asarray(H,dtype=float)/S
    a=np.array([float(x) for x in ai]);n=len(a)
    p=np.full(n,max(1e-8,np.linalg.eigvalsh(A)[0]/4))
    def ch(q):
        if np.min(q)<=0:return None
        try:return np.linalg.cholesky(A-np.diag(q))
        except np.linalg.LinAlgError:return None
    def obj(q,L,mu):return np.sum(a/q)-2*mu*np.log(np.diag(L)).sum()
    best=None
    for mu in (.1,.03,.01,.003,.001,.0003,.0001,.00003,.00001,.000003,.000001):
        for _ in range(70):
            L=ch(p);assert L is not None
            inv=np.linalg.solve(L.T,np.linalg.solve(L,np.eye(n)))
            g=-a/p**2+mu*np.diag(inv)
            Hess=mu*inv**2+np.diag(2*a/p**3)
            step=np.linalg.solve(Hess,-g);descent=g@step
            if -descent<1e-13:break
            old=obj(p,L,mu);t=1.
            while t>1e-10:
                q=p+t*step;Lq=ch(q)
                if Lq is not None and obj(q,Lq,mu)<=old+.01*t*descent:
                    p=q;break
                t/=2
            if t<=1e-10:break
        P=[int(x*10**12)*10**6 for x in p]
        if min(P)<=0:continue
        if not positive([[v-(P[i] if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(H)]):continue
        # Outward sum on 10^-18 grid: no floating-point acceptance.
        theta=sum(ceil(a*S*S/F(d)) for a,d in zip(ai,P))
        row={'positive_pivots':True,'diagonal_grid':list(map(str,P)),
             'theta_upper_grid':str(theta),'scale':str(S),
             'pass_conditional':theta<S}
        if best is None or theta<int(best['theta_upper_grid']):best=row
        if theta<S:break
    return best or {'positive_pivots':False,'pass_conditional':False}
def main(p,bank,oldbank,oldresult,replayed,logs,out):
    old=json.loads(oldresult.read_text());assert old['rank']==69
    rank70=json.loads(replayed.read_text());assert rank70['status']=='PASS_CONDITIONAL_RANK70' and rank70['conditional_index_upper']==15
    assert F(rank70['theta_upper'])<1 and F(old['sectors']['1']['levels'][-1]['margin'])>0
    sel=json.loads((bank/'selection.json').read_text());assert sel['rank']==85
    entries=json.loads((bank/'source_bank.json').read_text())['entries'];assert len(entries)==85
    oldentries=json.loads((oldbank/'source_bank.json').read_text())['entries'];assert len(oldentries)==70
    E={(e['parity'],e['coordinate_zero_based']):e for e in entries}
    O={(e['parity'],e['coordinate_zero_based']):e for e in oldentries}
    assert len(E)==85 and len(O)==70 and O.keys()<=E.keys()
    for key,e in E.items():
        assert sha(bank/e['file'])==e['sha256'] and F(e['norm2'])==1
        if key in O:
            assert e['sha256']==O[key]['sha256'] and e['omitted_group_midpoint_norm_upper']==O[key]['omitted_group_midpoint_norm_upper']
    report={};total=0
    for k,h in HASHES.items():
        J=sel['sectors'][k]['indices_zero_based'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        cols=old['sectors'][k]['columns'];assert [c['coordinate'] for c in cols]==J[:35]
        a={}
        for c in cols:
            j=c['coordinate'];assert c['input_sha256']==E[k,j]['sha256']
            a[j]=F(c['norm_upper'])**2/ELL
        if k=='2':
            for j,x in rank70['refined_inverse_upper'].items():a[int(j)]=min(a[int(j)],F(x))
        rows={}
        for j in J[35:]:
            row=moments(logs/f'mode_{k}_{j}.log',E[k,j]);a[j]=F(row['inverse_upper']);rows[str(j)]=row
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        attempts={};best=35
        for r in range(n,35,-1):
            I=J[:r];v=candidate(lower_matrix(w,I),[a[j] for j in I])
            attempts[str(r)]=v
            if v['pass_conditional']:best=r;break
        total+=best
        report[k]={'coordinate_order':J,'new_moments':rows,'attempts':attempts,'best_conditional_rank':best}
    result={'status':'CONDITIONAL_FULL85_DIAGNOSTIC','rank':total,'conditional_index_upper':85-total,
            'sectors':report,'premises':'Authenticated original fixed-collar inputs, continuum identities, proxy and numerical enclosure premises; numerical optimization only proposes diagonals. Not a Lean or unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
