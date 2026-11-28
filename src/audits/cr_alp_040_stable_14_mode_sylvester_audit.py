#!/usr/bin/env python3
"""CR-ALP-040: stable exact spectral enclosure of source-signed 14-mode H.

Inputs are interval matrices produced by the separate checked-int256
source-respecting Jensen-Gamma lower-matrix builder. Unlike recursive interval
LDL (which may explode on near-null Schur pivots), use a SINGLE exact common
integer scaling + Sylvester determinants of the CENTER matrix, with a direct
row-sum bound on ALL input interval uncertainties.

CLASSICAL SOURCE PREMISES REQUIRED: original Weil signed prime matrix,
physical digamma Gamma series, positivity and Jensen of Gamma kernels,
proper pole null basis. No Lean, 15000-mode Schur, full collar or RH.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json

ROOT=Path(__file__).resolve().parent
S=10**25
SIGMA=F(1,10**14)
THETA=F(1,10**16)
HFILE=ROOT/'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv'
INFILE=ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt'
assert HFILE.exists() and INFILE.exists()

def read_h(path):
    d={}
    for ln in path.read_text().splitlines():
        i,j,lo,hi=map(int,ln.split())
        assert lo<=hi and (i,j) not in d
        d[i,j]=(lo,hi)
    return d

def read_m(path):
    d={}
    for ln in path.read_text().splitlines():
        key,lo,hi=ln.split()
        if key.startswith('M_'):
            _,i,j=key.split('_'); d[int(i),int(j)]=(int(lo),int(hi))
    return d
H,M=read_h(HFILE),read_m(INFILE)
assert len(H)==98

def sym(d,i,j):
    a,b=d[i,j],d[j,i]
    lo=max(a[0],b[0]);hi=min(a[1],b[1])
    assert lo<=hi, ('wrong symmetric interval',i,j)
    return lo,hi

def det_bareiss(a):
    """Integer fraction-free determinant, checking divisions at each step."""
    n=len(a)
    t=[row[:] for row in a]
    den=1
    if n==1:return t[0][0]
    for k in range(n-1):
        p=t[k][k]
        assert p!=0, ('zero elimination pivot',k)
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=t[i][j]*p-t[i][k]*t[k][j]
                assert num%den==0,('nonexact Bareiss division',n,k)
                t[i][j]=num//den
        for i in range(k+1,n):t[i][k]=0
        den=p
    return t[-1][-1]

blocks=[]
for ids in ([3,5,7,9,11,13,1],[4,6,8,10,12,14,2]):
    n=len(ids)
    h=[[sym(H,i,j) for j in ids] for i in ids]
    m=[[sym(M,i,j) for j in ids] for i in ids]
    max_row_uncertainty=F(0)
    for r in range(n):
        # True (H-sigma M) minus center has absolute row sum at most this.
        row=sum((F(h[r][c][1]-h[r][c][0],2*S) +
                 SIGMA*F(m[r][c][1]-m[r][c][0],2*S)
                 for c in range(n)),F(0))
        if row>max_row_uncertainty:max_row_uncertainty=row
    assert max_row_uncertainty < THETA, (ids, max_row_uncertainty)
    # Let A0=C_H - 1e-14 C_M - 1e-16 I.
    # All entries represented exactly by integers in units 1/(2*10^39).
    scale=2*10**39
    mat=[[(h[r][c][0]+h[r][c][1])*10**14 -
          (m[r][c][0]+m[r][c][1]) -
          (2*10**23 if r==c else 0)
          for c in range(n)] for r in range(n)]
    assert all(mat[r][c]==mat[c][r] for r in range(n) for c in range(n))
    minors=[]
    for k in range(1,n+1):
        d=det_bareiss([row[:k] for row in mat[:k]])
        assert d>0,('negative exact centered Sylvester minor',ids,k,d)
        minors.append(d)
    # Positive exact A0 and symmetric interval norm ≤theta give
    # W_original|E14 >= H_midpoint >= sigma M by Loewner order.
    print('PASS',ids, 'all',n,'exact centered Sylvester determinants positive')
    print('  uncertainty row sum: ',max_row_uncertainty, ' < ', THETA)
    print('  min determinant signs:', ['+' if d>0 else '-' for d in minors])
    blocks.append({'indices':ids,'positive_centered_principal_minors':n,
                   'max_row_uncertainty':str(max_row_uncertainty),
                   'theta':str(THETA),'lower_physical_shift':str(SIGMA),
                   'determinants_sha256':sha256(('|'.join(map(str,minors))).encode()).hexdigest()})

out={
 'stage':'CR-ALP-040',
 'status':'EXACT INTEGER SYLVESTER + SPECTRAL INTERVAL CERTIFICATE PASS on genuine original signed E14',
 'theorem_conditional_on_classical_analytic_premises':
    'W_original(f,f) > (1/10^14) ||f||_L2^2 for all nonzero f in E14 (original physical pole null basis)',
 'gamma_lower':'First 30000 individually evaluated positive Gamma kernels + 11689 Jensen midpoint groups over [30000,10000000), remaining PSD kernels omitted on lower side',
 'matrix_interval_scale':'1/10^25',
 'method':'Exact midpoint Sylvester minors (fraction-free Bareiss) with full row-sum norm error budget, avoiding interval LDL dependency blowup',
 'input_sha256':sha256(INFILE.read_bytes()).hexdigest(),
 'gamma_matrix_sha256':sha256(HFILE.read_bytes()).hexdigest(),
 'blocks':blocks,
 'not_claimed':'Full 14-mode plus Q15000 mixed Schur, 14-to-15000 shell, whole fixed collar, Lean, RH'
}
(ROOT/'CR-ALP-040_certified_14_mode_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('CR-ALP-040 EXACT 14-MODE SOURCE-SIGNED FINITE CERTIFICATE PASS (not whole-collar positivity)')
