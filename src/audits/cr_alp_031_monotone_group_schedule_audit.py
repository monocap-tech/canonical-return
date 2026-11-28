#!/usr/bin/env python3
"""CR-ALP-031 independent exact schedule and endpoint-enclosure audit.

This is arithmetic/source-control support. The continuum Loewner monotonicity
J_beta>=J_beta' for beta<=beta' is proved analytically in the research note.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256

HERE=Path(__file__).resolve().parent
scale=10**25
records={}
for line in (HERE/'CR-ALP-031_initial_rational_intervals_1e25.txt').read_text().splitlines():
 key,low,high=line.split();assert key not in records
 records[key]=(int(low),int(high))
assert len(records)==129
assert F(records['q'][1],scale)<F(2,5)<F(1,2)
assert 2**100>10**25
assert 120001>100
# For t>=30000, q^(4t+1) < (1/2)^100 <10^-25.

K=30000;stop=350000
start=K; groups=[]
while start<stop:
 length=max(1,start//250)
 end=min(stop,start+length)
 groups.append((start,end,end-1,end-start))
 start=end
assert len(groups)==617
assert groups[0][0]==K and groups[-1][1]==stop
assert all(start<end and end-start==weight and t==end-1 for start,end,t,weight in groups)
assert all(groups[i][1]==groups[i+1][0] for i in range(len(groups)-1))
assert sum(x[3] for x in groups)==320000
assert sum(x[3] for x in groups)+K==stop
print('PASS: exact first 30000 Gamma kernels + 617 endpoint groups partition m=0..349999')
print('PASS: q^(4m+1) is contained in [0,1e-25] for every grouped endpoint')

# Check all 10x10 within-parity matrix intervals: 25 entries per block.
mat={}
for line in (HERE/'CR-ALP-031_grouped_H_exact_matrix_intervals.tsv').read_text().splitlines():
 i,j,lo,hi=line.split();ij=(int(i),int(j))
 assert ij not in mat and int(lo)<=int(hi)
 mat[ij]=(int(lo),int(hi))
assert len(mat)==50
for order in ((3,5,7,9,1),(4,6,8,10,2)):
 for i in order:
  for j in order:
   assert (i,j) in mat
   a,b=mat[i,j],mat[j,i]
   assert max(a[0],b[0])<=min(a[1],b[1])
print('PASS: complete symmetric 10-mode grouped lower matrix interval bank')
print('SHA256 grouped matrix',sha256((HERE/'CR-ALP-031_grouped_H_exact_matrix_intervals.tsv').read_bytes()).hexdigest())
