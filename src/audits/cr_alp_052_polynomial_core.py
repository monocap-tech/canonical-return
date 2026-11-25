#!/usr/bin/env python3
"""Exact outward integer evaluation of finite projected Chebyshev moments.
Does not evaluate the full infinite inverse. Analytic source intervals inherited.
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

