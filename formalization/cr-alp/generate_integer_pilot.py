#!/usr/bin/env python3
"""Exact input-to-Lean literal transfer. Pilot: final row of the actual cosh prefix matrix."""
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

def main(out):
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
    v=w['1P']; n=v['dimension']; assert n==43
    M=ints(v['matrix_numerators']); N=ints(v['factor_numerators'])
    D=int(v['matrix_denominator']); S=int(v['factor_denominator']); delta=F(v['delta'])
    def arr(A):
        return '#[\n'+',\n'.join('  #['+', '.join(map(str,row))+']' for row in A)+'\n]'
    lines=['import IntegerWitness','import Mathlib.Tactic','namespace CRALP.Pilot',
           'set_option maxRecDepth 20000','set_option maxHeartbeats 0',
           'def matrixRows : Array (Array ℤ) := '+arr(M),
           'def factorRows : Array (Array ℤ) := '+arr(N),
           'def M (i j : Fin 43) : ℤ := (matrixRows[i.val]!)[j.val]!',
           'def N (i j : Fin 43) : ℤ := (factorRows[i.val]!)[j.val]!',
           f'def D : ℤ := {D}',f'def S : ℤ := {S}',f'def d : ℤ := {delta.denominator}',f'def p : ℤ := {delta.numerator}']
    for j in range(n):
        lines += [f'theorem pilotEntry{j} :',
                  f'    -(p * D * S^2) ≤ 2 * (43 : ℤ) * intResidual M N D S d p 42 {j} ∧',
                  f'    2 * (43 : ℤ) * intResidual M N D S d p 42 {j} ≤ p * D * S^2 := by decide']
    lines += ['theorem pilotRow : rowOK M N D S d p 42 := by','  intro j','  fin_cases j']
    lines += [f'  · exact pilotEntry{j}' for j in range(n)]
    lines += ['end CRALP.Pilot','#print axioms CRALP.Pilot.pilotRow']
    Path(out).write_text('\n'.join(lines)+'\n')
    print('Generated 43 kernel entry propositions and row proof for cosh prefix Gram (row 42).')
if __name__=='__main__':main(sys.argv[1])
