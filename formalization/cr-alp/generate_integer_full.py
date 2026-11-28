#!/usr/bin/env python3
"""LEAN-02B2: source-bound matrix literals, row shards and terminal native receivers."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, sys
sys.set_int_max_str_digits(0)

def load(path, sha):
    raw = Path(path).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == sha, path
    return json.loads(raw)

def ints(rows):
    return [[int(x) for x in row] for row in rows]

def product(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]

def transpose(A):
    return list(map(list,zip(*A)))

def main(out, key):
    w = load('custody/cr-alp-lean-02ar2-matrix-soundness/witnesses.json','d729ae089f673a2426cfe0e29934cf88ec2074f30fbbd8ca0614f37bd236ee57')['matrices']
    q = load('custody/cr-alp-123aq2-frozen-response/result.json','daafe6fb527d6e4cc0a9322bbc1ea8555b502192b826bacefda4a3839b918261')['sectors']
    a = load('custody/cr-alp-123al-metric-first/result.json','b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e')['sectors']
    count=0
    for k,s in q.items():
        n=s['dimension']; T=int(s['matrix_grid'])
        H=ints(s['corrected_energy_lower']); P=ints(s['prefix_inverse_gram_upper'])
        Q=int(a[k]['whitening_grid']); B=ints(a[k]['whitening_matrix'])
        p=F(s['prefix_relative_upper'])
        C=product(product(B,H),transpose(B))
        W=[[C[i][j]-T*Q*Q*int(i==j) for j in range(n)] for i in range(n)]
        G=[[p.numerator*H[i][j]-p.denominator*P[i][j] for j in range(n)] for i in range(n)]
        for name,M,D in [('H',H,T),('W',W,T*Q*Q),('P',G,T*p.denominator)]:
            v=w[k+name]
            assert v['dimension']==n and ints(v['matrix_numerators'])==M and int(v['matrix_denominator'])==D
            N=ints(v['factor_numerators']); S=int(v['factor_denominator']); d=F(v['delta']).denominator; p0=F(v['delta']).numerator
            assert all(t>0 for t in (D,S,d,p0))
            assert len(N)==n and all(len(row)==n for row in N)
            NN=product(N,transpose(N))
            for i in range(n):
                for j in range(n):
                    R=d*S*S*M[i][j]-d*D*NN[i][j]-p0*D*S*S*int(i==j)
                    assert 2*n*abs(R)<=p0*D*S*S
                    count+=1
    assert count==10839
    print('All six matrix literal identities and 10839 exact residual entries PASS')
    assert set(w)=={'1H','1W','1P','2H','2W','2P'}
    assert key in w
    v=w[key]; n=v['dimension']; assert n==(43 if key[0]=='1' else 42)
    M=ints(v['matrix_numerators']); N=ints(v['factor_numerators'])
    D=int(v['matrix_denominator']); S=int(v['factor_denominator']); delta=F(v['delta'])
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    tag='ALP'+key
    ns='CRALP.Cert'+key
    data=tag+'Data'
    files=[]
    def emit(name,lines):
        f=out/(name+'.lean')
        f.write_text('\n'.join(lines)+'\n')
        files.append({'module':name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
    def rowdefs(name,A):
        lines=[f'def {name} (i j : Fin {n}) : ℤ :=','  match i.val with']
        for i,row in enumerate(A):
            lines.append(f'  | {i} => (#[{", ".join(map(str,row))}] : Array ℤ)[j.val]!')
        lines.append('  | _ => 0')
        return lines
    options=['set_option maxRecDepth 20000','set_option maxHeartbeats 0']
    lines=['import IntegerWitness','import Mathlib.Tactic',f'namespace {ns}']+options
    lines+=rowdefs('M',M)+rowdefs('N',N)
    lines += [f'def D : ℤ := {D}',f'def S : ℤ := {S}',f'def d : ℤ := {delta.denominator}',f'def p : ℤ := {delta.numerator}',f'end {ns}']
    emit(data,lines)
    shards=[]
    for lo in range(0,n,5):
        hi=min(lo+5,n)
        name=f'{tag}Rows{lo:02d}'
        lines=[f'import {data}',f'namespace {ns}']+options
        for i in range(lo,hi):
            for j in range(n):
                lines += [f'theorem entry_{i}_{j} :',
                          f'    -(p * D * S^2) ≤ 2 * ({n} : ℤ) * intResidual M N D S d p {i} {j} ∧',
                          f'    2 * ({n} : ℤ) * intResidual M N D S d p {i} {j} ≤ p * D * S^2 := by decide']
            lines += [f'theorem row_{i} : rowOK M N D S d p {i} := by','  intro j','  fin_cases j']
            lines += [f'  · exact entry_{i}_{j}' for j in range(n)]
        lines += [f'end {ns}']
        lines += [f'#print axioms {ns}.row_{i}' for i in range(lo,hi)]
        emit(name,lines)
        shards.append({'module':name,'first':lo,'end_exclusive':hi})
    terminal=tag+'Terminal'
    lines=[f'import {s["module"]}' for s in shards]+[f'namespace {ns}','open scoped BigOperators']+options
    lines += ['theorem allRows : ∀ i, rowOK M N D S d p i := by','  intro i','  fin_cases i']
    lines += [f'  · exact row_{i}' for i in range(n)]
    lines += [f'theorem matrixLower (x : Fin {n} → ℝ) :',
              '    ((p : ℝ) / (2 * d)) * (∑ i, x i ^ 2) ≤',
              '      quad (fun i j => (M i j : ℝ) / D) x := by',
              '  exact intWitnessLower M N D S d p (by decide) (by decide) (by decide)',
              '    (by decide) (by decide) allRows x',
              f'end {ns}',f'#print axioms {ns}.allRows',f'#print axioms {ns}.matrixLower']
    emit(terminal,lines)
    # Check every generated entry and row index, not merely declaration counts.
    import re
    entries=[]; rows=[]
    for f in files:
        s=(out/(f['module']+'.lean')).read_text()
        assert not re.search(r'\b(?:sorry|admit|sorryAx|axiom|native_decide)\b',s)
        entries += [(int(i),int(j)) for i,j in re.findall(r'(?m)^theorem entry_(\d+)_(\d+) :$',s)]
        rows += [int(i) for i in re.findall(r'(?m)^theorem row_(\d+) :',s)]
    assert sorted(entries)==[(i,j) for i in range(n) for j in range(n)]
    assert sorted(rows)==list(range(n))
    manifest={'matrix':key,'dimension':n,'entries':n*n,'rows':n,'namespace':ns,
              'data':data,'shards':shards,'terminal':terminal,'files':files,
              'scales':{'D':str(D),'S':str(S),'d':str(delta.denominator),'p':str(delta.numerator)}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'{key}: {n*n} entry propositions, {n} rows, {len(shards)} shards, terminal coercivity candidate')
if __name__=='__main__':main(sys.argv[1],sys.argv[2])
