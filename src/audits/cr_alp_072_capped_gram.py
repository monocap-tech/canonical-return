#!/usr/bin/env python3
"""CR-ALP-072: exact capped correction Gram and robust conditional receiver.
Usage: python script.py ARTIFACT_DIRECTORY GAMMA068_JSON
No new integration. Conditional on original correction construction, 068/070
analytical premises, and actual second old projected residual norm <= 1/5.
"""
import json, zipfile, hashlib, sys
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
S=10**30
N=2048
ids=range(86,N+1,2)
def P(x):return F(x),F(x)
Z=P(0)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    t=[x*y for x in a for y in b];return min(t),max(t)
def div(a,b):
    assert b[0]>0;return mul(a,(1/b[1],1/b[0]))
def rnd(a):return F(a[0]*S//1,S),-F(-a[1]*S//1,S)
def total(a):
    l=h=F(0)
    for x,y in a:l+=x;h+=y
    return rnd((l,h))
def up(x):
    assert x>=0;t=x*S*S;return F(isqrt(t.numerator//t.denominator)+1,S)
def readzip(p,h):
    assert hashlib.sha256(p.read_bytes()).hexdigest()==h
    z=zipfile.ZipFile(p);assert z.testzip() is None;return z
def run(p,gamma):
    b=readzip(p/'CR-ALP-071-common-trial-bank.zip','ebbc51e56e67e7438483928b3570a320b77df06b92a574ad2a0e998fdb4b0ead')
    raw=b.read('bank.json');assert hashlib.sha256(raw).hexdigest()=='b2d5f187293de20f05a9e510c1f54aeb1828fc2a20f55a0e4127546a7d99f71e';bank=json.loads(raw)
    r=readzip(p/'CR-ALP-071-physical-plus.zip','bd2266c3a2b3e177f48ba9abe2f650dc31a54292b54ca81ce5a80fbfbc23ac60')
    pole={};st=None
    for line in r.read('response_plus.log').decode().splitlines():
        v=line.split()
        if len(v)!=5 or v[3:]!=['/','1e15']:continue
        if v[0]=='pole_st':st=tuple(F(int(x),10**15) for x in v[1:3])
        elif v[0].startswith('pole'):pole[int(v[0][4:])]=tuple(F(int(x),10**15) for x in v[1:3])
    assert st and st[0]>0 and set(pole)==set(range(2,N+1,2))
    # Verify Git blob identity of original Gamma bank.
    raw=gamma.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()=='1a7fe8629daa8b8e6cec185a79b6746d6b28975d'
    bands=json.loads(raw)['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    floors=[(b['start'],b['end'],F(F(b['gamma_lower'])*10**9//1,10**9)) for b in bands]
    cap=next(d for l,h,d in floors if l==2049);gamma0=F(72,25);acap=cap-gamma0
    assert cap==F(555992307,125000000) and all(d>=gamma0 for _,_,d in floors)
    assert all(d>=cap for l,h,d in floors if l>=2049)
    assert min(d for _,_,d in floors)-F(21,8)==F(26535863,100000000)
    A={j:min(cap,next(d for l,h,d in floors if l<=j and (h is None or j<=h)))-gamma0 for j in ids}
    U={};z={}
    for tag in ['first','second']:
        t={int(j):rnd(tuple(map(F,lh))) for j,lh in bank['banks'][tag]['correction_exact_intervals'].items()}
        assert set(t)==set(range(2,851,2));t={j:t.get(j,Z) for j in ids}
        tau=rnd(div(total(mul(t[j],pole[j]) for j in ids),st))
        U[tag]={j:rnd(sub(t[j],mul(tau,pole[j]))) for j in ids}
        # Q85 kills the full constant-cap multiple of the pole.
        z[tag]={j:rnd(sub(mul(P(A[j]),t[j]),mul(P(A[j]-acap),mul(tau,pole[j])))) for j in ids}
    gram={}
    for a,b in [('first','first'),('first','second'),('second','second')]:
        q=total(mul(U[a][j],z[b][j]) for j in ids);r=total(mul(U[b][j],z[a][j]) for j in ids)
        q=max(q[0],r[0]),min(q[1],r[1]);assert q[0]<=q[1];gram[a,b]=q
    k22=add(tuple(map(F,bank['K_L_22'])),gram['second','second'])
    k12=add(tuple(map(F,bank['K_L_12'])),gram['first','second'])
    a2=gram['second','second'][1];assert a2>=0
    m1=F(61214171350249519,10**18);B1=F(16059528681514303,62500000000000000);dmin=F(26535863,100000000)
    q=F(1,5);B2=(q*up(1/dmin)+up(a2))**2
    m2=k22[0]-B2;c=max(map(abs,k12))+up(B1*B2);det=m1*m2-c*c
    assert m2>0 and det>0
    return {'status':'PASS_CONDITIONAL','a11':list(map(str,gram['first','first'])),'a12':list(map(str,gram['first','second'])),'a22':list(map(str,gram['second','second'])),'KC12':list(map(str,k12)),'KC22':list(map(str,k22)),'norm_threshold':str(q),'B2_upper':str(rnd(P(B2))[1]),'S22_lower':str(rnd(P(m2))[0]),'cross_upper':str(rnd(P(c))[1]),'det_lower':str(rnd(P(det))[0])}
if __name__=='__main__':print(json.dumps(run(Path(sys.argv[1]),Path(sys.argv[2])),indent=2))
