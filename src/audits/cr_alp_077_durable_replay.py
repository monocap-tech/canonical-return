#!/usr/bin/env python3
"""Regenerate and receive CR-ALP-077 from durable Git archive inputs.

Conditional on the inherited continuum, uniform proxy and numerical enclosure
premises. Original local ZIP custody is independent of this new native replay.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_074_norm_receiver import parse,up
from cr_alp_075_shifted_energy import lower_matrix,positive,S
D=10**80
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fmt(x):
    z=x*D;assert z.denominator==1;n=z.numerator
    return ('-' if n<0 else '')+str(abs(n)//D)+'.'+str(abs(n)%D).zfill(80)
def interval(a,b,sgn):
    l,h=map(F,a);ll,hh=map(F,b);assert l<=h and ll<=hh
    return [fmt(l+ll if sgn==1 else l-hh),fmt(h+hh if sgn==1 else h-ll)]
def combine(a,b,sgn):
    a=a.splitlines();b=b.splitlines();assert len(a)==len(b)
    k=0;out=[]
    def same():
        nonlocal k
        assert a[k]==b[k];out.append(a[k]);k+=1;return int(a[k-1])
    def pair():
        nonlocal k
        out.append(' '.join(interval(a[k].split(),b[k].split(),sgn)));k+=1
    assert same()==2;pair();pair();n=same();assert n==425
    for _ in range(n):
        x,y=a[k].split(),b[k].split();assert x[0]==y[0] and len(x)==len(y)==5
        out.append(' '.join([x[0]]+interval(x[1:3],y[1:3],sgn)+interval(x[3:5],y[3:5],sgn)));k+=1
    n=same();assert n==128
    for _ in range(n):
        x,y=a[k].split(),b[k].split();assert x[:2]==y[:2] and len(x)==len(y)==4
        out.append(' '.join(x[:2]+interval(x[2:],y[2:],sgn)));k+=1
    n=same();assert n==25
    for _ in range(n):
        assert a[k]==b[k];out.append(a[k]);m=int(a[k].split()[-1]);assert m==40;k+=1
        for _ in range(m):pair()
    assert k==len(a)
    return '\n'.join(out)+'\n'
def sources(bank):
    entries=json.loads((bank/'source_bank.json').read_text())['entries']
    assert len(entries)==70
    E={e['coordinate_zero_based']:e for e in entries if e['parity']=='2'}
    assert len(E)==35
    for j,e in E.items():
        assert e['file']==f'input_2_{j}.txt' and F(e['norm2'])==1
        assert sha(bank/e['file'])==e['sha256']
    return E
def prepare(bank,out):
    E=sources(bank);out.mkdir(parents=True,exist_ok=True)
    for j in (14,9):
        for s,tag in [(1,'plus'),(-1,'minus')]:
            f=out/f'12_{j}_{tag}.txt'
            f.write_text(combine((bank/E[12]['file']).read_text(),(bank/E[j]['file']).read_text(),s))
    (out/'mode12.txt').write_bytes((bank/E[12]['file']).read_bytes())
    print('PASS four regenerated polarization inputs and mode12')
def receive(p,bank,logs,result,witness,out):
    E=sources(bank);data=json.loads(result.read_text());assert data['rank']==69
    W=json.loads(witness.read_text());J=W['coordinates'];assert len(J)==len(set(J))==35
    cols=data['sectors']['2']['columns'];assert J==[c['coordinate'] for c in cols]
    b={c['coordinate']:F(c['norm_upper']) for c in cols}
    for c in cols:assert c['input_sha256']==E[c['coordinate']]['sha256']
    # Sum of the two coordinate approximation errors is a valid transport
    # bound for the plus/minus proxy, including per-input displacement reserves.
    def eps(j):return F(5,10**9)+F(E[j]['omitted_group_midpoint_norm_upper'])+F(1,10**13)
    def normsq(tag,e):
        l,h=parse(logs/f'{tag}.log');q=up(h);err=e*(2*q+e)
        return max(F(0),l-err),h+err
    cross={}
    for j in (14,9):
        a,b1=normsq(f'12_{j}_plus',eps(12)+eps(j))
        c,d=normsq(f'12_{j}_minus',eps(12)+eps(j))
        cross[j]=((a-d)/4,(b1-c)/4)
    rows={}
    for line in (logs/'mode12.log').read_text().splitlines():
        n,l,h,s,d=line.split();assert s=='/' and d=='1e15' and n not in rows
        l,h=F(int(l),10**15),F(int(h),10**15);assert l<=h
        rows[n]=(l,h)
    assert set(rows)=={'G','M','penalty',*(f'corr{i}' for i in range(7))}
    G=rows['G'][1];assert G>=0
    err=eps(12)*(2*up(G)+eps(12))
    Gup=G+err;Mlo=rows['M'][0]-F(21,8)*err-F(1,10**13)
    A=(F(72,25)*Gup-Mlo)/(ELL*F(1101,200));assert A>=0
    refined={12:A};remainder={}
    for j,alpha in [(14,F(233661,250000)),(9,-F(143787,200000))]:
        l,h=cross[j]
        e2=b[j]**2-2*alpha*(l if alpha>=0 else h)+alpha**2*b[12]**2
        assert e2>=0
        remainder[j]=up(e2)
        refined[j]=(abs(alpha)*up(A)+up(e2/ELL))**2
    # Use authenticated original energy matrix, not the screen's stored lower matrix.
    w=checked(p/'cr_alp_056/result_extended_2.json',HASHES['2'])
    H=lower_matrix(w,J);Dg=list(map(int,W['results']['three_refined']['diagonal']))
    assert int(W['matrix_scale'])==S and len(Dg)==35 and all(d>0 for d in Dg)
    assert H==[list(map(int,r)) for r in W['energy_lower']]
    assert positive([[v-(Dg[i] if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(H)])
    ai=[min(b[j]**2/ELL,refined.get(j,b[j]**2/ELL)) for j in J]
    theta=sum(a/F(d,S) for a,d in zip(ai,Dg))
    odd=data['sectors']['1'];assert odd['best_conditional_rank']==35
    assert F(odd['levels'][-1]['margin'])>0
    ans={'status':'PASS_CONDITIONAL_RANK70' if theta<1 else 'UNDECIDED',
         'theta_upper':str(theta),'reserve_lower':str(1-theta),
         'conditional_index_upper':15 if theta<1 else 16,
         'cross_intervals':{str(j):list(map(str,x)) for j,x in cross.items()},
         'refined_inverse_upper':{str(j):str(a) for j,a in refined.items()},
         'remainder_norm_upper':{str(j):str(a) for j,a in remainder.items()},
         'premises':'Inherited original fixed-collar continuum identities, proxy and numerical enclosure premises; authenticated source and norm banks; no unrestricted RH claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
    assert theta<1
if __name__=='__main__':
    mode,*args=sys.argv[1:]
    if mode=='prepare':prepare(*map(Path,args))
    elif mode=='receive':receive(*map(Path,args))
    else:raise ValueError(mode)
