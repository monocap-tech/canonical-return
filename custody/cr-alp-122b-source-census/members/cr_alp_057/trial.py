#!/usr/bin/env python3
"""Reconstruct the CR-ALP-056 frozen trial as physical sine and complete pole coefficient intervals.
No sine-tail sweep, no new solve. Reuses predecessor projection algebra.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent
source=PARENT/'cr_alp_056/extended.py'
s=source.read_text().split('nu={j:')[0].replace('ROOT=Path(__file__).resolve().parent','ROOT=Path(__file__).resolve().parent.parent/"cr_alp_056"')
exec(s)
ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent
nu={j:j*NUF for j in range(1,851)}
def mat(A):return np.array(A,dtype=object)
def load(A):return mat([[I(int(l),int(h),True) for l,h in row] for row in A])
def save(A):return [[[str(v.l),str(v.h)] for v in row] for row in A]
def mul(A,B):
 def mr(A):
  M=np.array([[(x.l+x.h)//2 for x in row] for row in A],dtype=object)
  R=np.array([[max(m-x.l,x.h-m) for x,m in zip(row,ms)] for row,ms in zip(A,M)],dtype=object)
  return M,R
 ma,ra=mr(A);mb,rb=mr(B);mc=ma@mb;rc=abs(ma)@rb+ra@abs(mb)+ra@rb
 return mat([[I((x-r)//S,ceildiv(x+r,S),True) for x,r in zip(row,rad)] for row,rad in zip(mc,rc)])
result={}
for bl in bank['blocks']:
 k=bl['parity'];ids=bl['ids'];n=len(ids);b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4;BM=(E-q)/2+(A if k==1 else -A)
 def c(j):v=nu[j];return 2*b0*v/(rt*(v.sq()+qt))
 # Exact code from CR-ALP-056, with only the input path redirected to its archive.
 part=source.read_text().split(' low=list(range(k,86,2));')[1].split(' high=list(range(851')[0];part='low=list(range(k,86,2));'+part
 import textwrap
 part=part.replace('\n ','\n')
 part=part.replace("(ROOT/f'correction_{k}.json')", "(PARENT/'cr_alp_056'/f'correction_{k}.json')")
 exec(part)
 corraw=(PARENT/'cr_alp_056'/f'correction_{k}.json').read_bytes();assert hashlib.sha256(corraw).hexdigest()=={1:'7766934de00db1dfe5862c80b3ad9f8fe9e92caaa4aefef6590d4ef4f2c7b742',2:'75aff6016b1388eb2585f9bda3de0256565f16579858ee1d57b22db69a5d60b8'}[k]
 # Physical pole moment and low sine moments must enclose zero for each column.
 cm=mat([[c(j) for j in fullids]]);po=mul(cm,D)+Db*BM
 assert all(v.l<=0<=v.h for v in po[0]);
 for j in range(len(low)):
  for z in range(n):
   v=D[j,z]+c(low[j])*Db[0,z];assert v.l<=0<=v.h
 # Interval width is small, but exact pole-nullness follows from symbolic projections.
 print(k,'trial physical rows',len(fullids),'columns',n,'max pole interval magnitude',float(max(v.up() for v in po[0])),flush=True)
 result[str(k)]={'ids':fullids,'D':save(D),'Db':save(Db),'correction_sha256':hashlib.sha256(corraw).hexdigest(),'trial_predecessor_sha256':hashlib.sha256(pin).hexdigest()}
(ROOT/'trial.json').write_text(json.dumps({'scale':str(S),'sectors':result},separators=(',',':'))+'\n')
