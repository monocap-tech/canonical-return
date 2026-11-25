#!/usr/bin/env python3
"""Exact outward integer evaluation of finite projected Chebyshev moments.
Affine upper bound valid for global inverse on the compressed residual; higher degrees only for compressed inverse.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json, hashlib, sys
import numpy as np
ROOT=Path(__file__).resolve().parent
S=10**24

def cdiv(n,d):return -((-n)//d)
def load(A,old):
 d=old//S;assert old==d*S
 l=np.array([[int(x[0])//d for x in row] for row in A],dtype=object)
 h=np.array([[cdiv(int(x[1]),d) for x in row] for row in A],dtype=object)
 m=(l+h)//2;return m,np.maximum(m-l,h-m)
def add(A,B):return A[0]+B[0],A[1]+B[1]
def neg(A):return -A[0],A[1]
def sub(A,B):return add(A,neg(B))
def trans(A):return A[0].T,A[1].T
def scale(A,f):
 f=F(f);num,den=f.numerator,f.denominator
 m0=A[0]*num;r0=A[1]*abs(num);m=m0//den
 return m,np.maximum(m-(m0-r0)//den,-((-(m0+r0))//den)-m)
def mul(A,B):
 a,r=A;b,s=B;m0=a@b;rad=abs(a)@s+r@abs(b)+r@s;m=m0//S
 return m,np.maximum(m-(m0-rad)//S,-((-(m0+rad))//S)-m)
def zero(shape):return np.zeros(shape,dtype=object),np.zeros(shape,dtype=object)
def identity(n):return np.eye(n,dtype=object)*S,np.zeros((n,n),dtype=object)
def symmetric(A):
 l=np.maximum(A[0]-A[1],(A[0]-A[1]).T);h=np.minimum(A[0]+A[1],(A[0]+A[1]).T)
 assert np.all(l<=h);m=(l+h)//2;return m,np.maximum(m-l,h-m)
def upper(A):return F(max(sum(row) for row in abs(A[0])+A[1]),S)
def cert(A,shift):
 n=len(A[0]);x=shift*S;assert x.denominator==1;C=A[0].tolist();err=max(sum(row) for row in A[1]);
 for i in range(n):C[i][i]-=int(x)+err
 prev=1
 for k in range(n):
  pivot=C[k][k]
  if pivot<=0:return k
  for i in range(k+1,n):
   for j in range(i,n):
    q,r=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert r==0;C[i][j]=C[j][i]=q
  prev=pivot
 return n
def bounds(A):return [[[str(m-r),str(m+r)] for m,r in zip(row,rad)] for row,rad in zip(*A)]
def form(A,z):return mul(trans(z),mul(A,z))


x=F(192,175);delta=F(21,8);gamma=F(72,25);ell=F(51,200);u=F(1101,200)
assert gamma-delta==ell and gamma+delta==u
for parity in map(int,sys.argv[1:] or ['1','2']):
 raw=(ROOT/f'bank_{parity}.json').read_bytes();b=json.loads(raw);old=int(b['scale']);P=load(b['P'],old);RF=load(b['R'],old);K=load(b['K'],old)
 D=scale(P,-1/delta);n=K[0].shape[0];I=identity(n);G=symmetric(mul(trans(RF),RF));M=symmetric(mul(trans(RF),mul(P,RF)))
 # Affine reciprocal bound on spectrum: s(L)=(gamma I-P)/(ell*u).
 A=symmetric(scale(sub(scale(G,gamma),M),1/(ell*u)))
 out={'parity':parity,'input_sha256':hashlib.sha256(raw).hexdigest(),'scale':str(S),'K':bounds(K),'G':bounds(G),'M':bounds(M),'tests':{}}
 z=load(b['witness']['z'],old)
 def test(label,J):
  U=symmetric(sub(K,J)); count=0
  for sh in map(F,['1/4','1/8','1/16','1/32','1/100','1/1000','0']):
   count=cert(U,sh)
   if count==n:break
  w=form(U,z);norm=form(I,z);assert norm[0][0,0]-norm[1][0,0]>0
  vals=[F(int(a),int(d)) for a in [w[0][0,0]-w[1][0,0],w[0][0,0]+w[1][0,0]] for d in [norm[0][0,0]-norm[1][0,0],norm[0][0,0]+norm[1][0,0]]]
  print('parity',parity,label,'pivots',count,'shift',str(sh),'frozen witness',float(min(vals)),float(max(vals)),'penalty row norm',float(upper(J)),flush=True)
  out['tests'][label]={'positive_pivots':count,'shift':str(sh),'witness_lower':str(min(vals)),'witness_upper':str(max(vals)),'penalty':bounds(J),'comparison':bounds(U),'penalty_row_upper':str(upper(J))}
 test('affine',A)
 # Degree-m residual polynomial (q has degree m-1); compressed only.
 Z0=zero(RF[0].shape); Z1=scale(RF,1/gamma);t0,t1=F(1),x
 for m in range(1,13):
  if m>1:
   t2=2*x*t1-t0
   Z2=sub(scale(add(mul(D,Z1),scale(RF,1/delta)),2*t1/t2),scale(Z0,t0/t2))
   Z0,Z1,t0,t1=Z1,Z2,t1,t2
  if m not in (4,8,12):continue
  J=symmetric(scale(symmetric(mul(trans(RF),Z1)),1/(1-1/t1)))
  test('chebyshev_'+str(m),J)
 (ROOT/f'moments_{parity}.json').write_text(json.dumps(out)+'\n')
