#!/usr/bin/env python3
"""Exact diagnostic on preserved 071 trial and polarization artifacts.
Not a proof of applicability of the original forward-error enclosure."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,zipfile,json,sys
def readzip(path,sha):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
    z=zipfile.ZipFile(path);assert z.testzip() is None
    return z
def parse(z,tag):
    t=iter(z.read('input_2_'+tag+'.txt').decode().split())
    assert next(t)=='2'
    def rd():
        l,h=F(next(t)),F(next(t));assert l<=h
        return l,h
    f,db=rd(),rd();n=int(next(t));assert n==425
    ids=[];s=[];x=[]
    for _ in range(n):
        ids.append(int(next(t)));s.append(rd());x.append(rd())
    assert ids==list(range(2,851,2))
    n=int(next(t));assert n==128
    lc=[];low=[]
    for _ in range(n):lc.append(rd());low.append(rd())
    ng=int(next(t));assert ng==25
    centers=[];g=[]
    for _ in range(ng):
        centers.append(rd());n=int(next(t));g.append([rd() for _ in range(n)])
    assert next(t,None) is None
    return {'pole':[f],'correction_pole':[db],'source_sine':s,'correction_sine':x,'low_amplitudes':low,'group_amplitudes':[p for gg in g for p in gg]},(ids,lc,centers,list(map(len,g)))
def mag(v):return max(map(abs,v))
def getG(raw):
    q=[l.split() for l in raw.decode().splitlines() if l.startswith('G ')]
    assert len(q)==1 and q[0][3:]==['/','1e15']
    return tuple(F(int(t),10**15) for t in q[0][1:3])
def main(p):
    bank=readzip(p/'CR-ALP-071-common-trial-bank.zip','ebbc51e56e67e7438483928b3570a320b77df06b92a574ad2a0e998fdb4b0ead')
    a,ga=parse(bank,'first');b,gb=parse(bank,'second');assert ga==gb
    ratios={};exceptions={}
    for k in a:
        A=list(map(mag,a[k]));B=list(map(mag,b[k]))
        assert sum(B)<=sum(A)
        ratios[k]=str(sum(B)/sum(A)) if sum(A) else '0'
        exceptions[k]=sum(y>x for x,y in zip(A,B))
    plus=readzip(p/'CR-ALP-071-physical-plus.zip','bd2266c3a2b3e177f48ba9abe2f650dc31a54292b54ca81ce5a80fbfbc23ac60')
    minus=readzip(p/'CR-ALP-071-physical-minus.zip','e4e56a4eef044eefd021cd651bdc4dae98b5d82b19856035fb42d39fada8e458')
    P=getG(plus.read('response_plus.log'));M=getG(minus.read('response_minus.log'))
    A=F(136960221478274,10**15),F(136991883280279,10**15)
    B=getG(bank.read('second_coeff.log'))
    # This is a consistency diagnostic only until midpoint projection
    # linearity/error transport and all direction-specific envelopes are checked.
    q=(P[0]+M[0])/2-A[1],(P[1]+M[1])/2-A[0]
    c=(P[0]-M[1])/4,(P[1]-M[0])/4
    assert max(q[0],B[0])<=min(q[1],B[1])
    return {'status':'PASS_DIAGNOSTIC_ONLY','identical_input_geometry':True,'aggregate_budget_ratios':ratios,'entrywise_dominance_exceptions':exceptions,'polarized_G2_interval':list(map(str,q)),'direct_G2_interval':list(map(str,B)),'polarized_cross_interval':list(map(str,c))}
if __name__=='__main__':print(json.dumps(main(Path(sys.argv[1])),indent=2))
