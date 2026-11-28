#!/usr/bin/env python3
"""CR-ALP-071 candidate: select and attach a two-column original even source bank.

No new physical integration or Schur positivity is asserted. The first column
is exactly the accepted 070 witness. The second is chosen using only the old
corrected energy K_L; this is a deterministic candidate-selection heuristic,
not an optimization of the true Schur form. Plus/minus banks are for correlated
polarization with a COMMON linear response trial, not separately optimized
response trials. See the matrix receiving contract emitted in bank.json.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib, json, sys

D=10**12
OUT=10**80
def read_checked(path,sha):
    raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==sha,path
    return json.loads(raw)
def encint(n):
    return ('-' if n<0 else '')+str(abs(n)//OUT)+'.'+str(abs(n)%OUT).zfill(80)
def enc(x):
    x=F(x)*OUT
    l=x.numerator//x.denominator
    h=-((-x.numerator)//x.denominator)
    return encint(l)+' '+encint(h)
def interval_dot(row,v,den):
    assert len(row)==len(v)
    l=h=m=0
    for (a,b),x in zip(row,v):
        a,b=int(a),int(b); assert a<=b
        p,q=a*x,b*x
        l+=min(p,q);h+=max(p,q);m+=((a+b)//2)*x
    assert l<=m<=h
    return F(l,den),F(h,den),F(m,den)
def maxabs(x):return max(abs(x[0]),abs(x[1]))
def ep(x):return [str(x[0]),str(x[1])]
def qdot(a,b):return F(sum(x*y for x,y in zip(a,b)),D**2)
def upper_root(x,scale=10**40):
    assert x>=0
    y=x*scale**2
    return F(isqrt(y.numerator//y.denominator)+1,scale)
def run(package,out):
    pre=package/'cr_alp_057'
    w=read_checked(package/'cr_alp_056/result_extended_2.json','4017c7a656e43fe9fa97da91c37f35a9eda5a78dd3efaa2863b40de71177ac36')
    a=read_checked(pre/'result.json','e0dbb8cd24322cd5f0a96a7d8c1b188c38b7f5d2f706d644b8268cbb6533a809')
    g=read_checked(pre/'grouped_proxy.json','12c7fbc679f244f0016cd71acd7256cf00830bb40d14fb2fedf452cab7706c0c')
    t=read_checked(pre/'trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    S=int(a['scale']);assert S==int(t['scale'])==int(g['scale'])==10**100
    a,g,t=(x['sectors']['2'] for x in (a,g,t))
    z=list(map(int,w['comparison_witness']['z']))
    K=w['K'];ks=int(w['scale']);n=len(z)
    assert n>=2 and len(K)==n and all(len(row)==n for row in K)
    assert 0<qdot(z,z)<F(1000001,1000000)
    assert t['ids']==list(range(2,851,2))
    def bilinear(u,v):
        l=h=0
        for i,row in enumerate(K):
            for j,(a,b) in enumerate(row):
                a,b=int(a),int(b);assert a<=b
                c=u[i]*v[j]
                l+=min(a*c,b*c);h+=max(a*c,b*c)
        return F(l,ks*D**2),F(h,ks*D**2)
    k11=bilinear(z,z)
    assert k11[0]>0
    choices=[]
    for j in range(n):
        e=[0]*n;e[j]=D
        if qdot(z,z)*qdot(e,e)-qdot(z,e)**2<=0:continue
        col=interval_dot([row[j] for row in K],z,ks*D)[:2]
        row=interval_dot(K[j],z,ks*D)[:2]
        k12=max(col[0],row[0]),min(col[1],row[1])
        assert k12[0]<=k12[1],('symmetry',j)
        k22=tuple(F(int(x),ks) for x in K[j][j])
        score=k22[0]-maxabs(k12)**2/k11[0]
        choices.append((score,j,k12,k22,e))
    assert choices
    score,j,k12,k22,e=max(choices,key=lambda x:(x[0],-x[1]))
    assert bilinear(z,e)[0]<=k12[0]<=k12[1]<=bilinear(z,e)[1]
    assert bilinear(e,e)==k22
    gramdet=qdot(z,z)*qdot(e,e)-qdot(z,e)**2
    assert gramdet>0
    out.mkdir(parents=True,exist_ok=True)
    groups=[x for x in g['groups'] if int(x['M'])<2**32]
    omitted=[x for x in g['groups'] if int(x['M'])>=2**32]
    banks={}
    for tag,v in [('first',z),('second',e),('plus',[x+y for x,y in zip(z,e)]),('minus',[x-y for x,y in zip(z,e)])]:
        f=interval_dot(a['pole_coefficient'],v,S*D)
        db=interval_dot(t['Db'][0],v,S*D)
        lows={i:interval_dot(row,v,S*D) for i,row in zip(a['source_ids'],a['source_sine_coefficients'])}
        trial={i:interval_dot(row,v,S*D) for i,row in zip(t['ids'],t['D'])}
        # Exact input format of archived 058 prepare.py; the midpoint is a proxy.
        lines=['2',enc(f[2]),enc(db[2]),str(len(trial))]
        for i in sorted(trial):
            lines.append(str(i)+' '+enc(lows[i][2] if i in lows else 0)+' '+enc(trial[i][2]))
        lines.append(str(len(g['low_amplitudes'])))
        for m,row in enumerate(g['low_amplitudes']):
            lines.append(enc(F(4*m+1,2))+' '+enc(interval_dot(row,v,S*D)[2]))
        lines.append(str(len(groups)))
        for b in groups:
            lines.append(enc(F(b['center']))+' '+str(len(b['amplitudes'])))
            lines.extend(enc(interval_dot(row,v,S*D)[2]) for row in b['amplitudes'])
        raw=('\n'.join(lines)+'\n').encode()
        if tag=='first':
            original=(package/'cr_alp_058/input_2.txt').read_bytes()
            assert hashlib.sha256(original).hexdigest()=='f70358f279994211c4231ac43cdce46a20ca619a1f6d985e95d440e58b34923b'
            assert raw==original,'frozen first proxy input mismatch'
        path=out/('input_2_'+tag+'.txt');path.write_bytes(raw)
        # A finite rational upper for the omitted proxy groups, with per-group
        # upward rounding. Same 2/sqrt(2c) bound as 058, no floating arithmetic.
        Q=10**40;tail=0
        for b in omitted:
            c=F(b['center']);x=2*c*Q**2
            root=isqrt(x.numerator//x.denominator);assert root>0
            term=F(2*Q,root)*sum(abs(interval_dot(row,v,S*D)[2]) for row in b['amplitudes'])
            q=term*Q;tail+=-((-q.numerator)//q.denominator)
        # Projection kills the full pole term and the low sine coefficients.
        # This tolerance concerns the exact midpoint, not a new quadrature.
        error2=sum(max(abs(h-m),abs(l-m))**2 for i,(l,h,m) in trial.items() if i>85)
        banks[tag]={'z_numerator':v,'norm2':str(qdot(v,v)),
                    'K_L_interval':ep(bilinear(v,v)),
                    'input_sha256':hashlib.sha256(raw).hexdigest(),
                    'correction_exact_intervals':{str(i):ep(q) for i,q in trial.items()},
                    'correction_pole_exact_interval':ep(db),
                    'correction_projected_midpoint_error_upper':str(upper_root(error2)),
                    'omitted_group_midpoint_norm_upper':str(F(tail,Q))}
    result={'status':'SOURCE_BANK_CANDIDATE_ONLY','dimension':2,'first_input_byte_match':True,
            'source_denominator':D,'second_coordinate_zero_based':j,
            'K_L_11':ep(k11),'K_L_12':ep(k12),'K_L_22':ep(k22),
            'K_L_schur_lower':str(score),'K_L_pair_positive_sufficient':score>0,
            'euclidean_gram_determinant':str(gramdet),'banks':banks,
            'receiving_contract':{
                'K_C':'K_L + U* A U, preserve full cross terms',
                'r_C':'R - C U = r_L - A U',
                'response':'common linear trial T, E=T*r_C+r_C*T-T*C*T, V=r_C-C*T',
                'inverse_identity':'r_C*C^-1*r_C = E + V*C^-1*V',
                'weighted_gram':'V*C^-1*V <= finite weighted Gram + full projected tail Gram/dtail, using C >= Q85 diag(d) Q85',
                'tail':'full projected Gram minus retained Gram; use common linear trial for plus/minus, not separate solves',
                'errors':'bound actual proxy/action/attachment errors columnwise and transport in matrix order; do not treat diagonal bounds as a cross Gram',
                'gate':'certify K_C minus full inverse upper matrix positive definite; no result from the source bank alone'}}
    (out/'bank.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='banks'},indent=2))
    print('PASS source-bank construction and byte-identical original first input; no Schur promotion')
if __name__=='__main__':
    assert len(sys.argv)==3
    run(Path(sys.argv[1]),Path(sys.argv[2]))
