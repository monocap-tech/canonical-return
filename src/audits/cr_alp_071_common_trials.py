#!/usr/bin/env python3
"""Freeze common linear response trials for the CR-ALP-071 pair.

The second numerical solve is candidate selection only. Its acceptance will
require new physical response logs and original-source error transport.
The first trial and diagonal are the original authenticated 067 files.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, sys

IDS=list(range(86,2049,2))
S=10**12
def load(path):
    rows={}
    for line in path.read_text().splitlines():
        j,l,h=line.split();j=int(j)
        assert j not in rows and F(l)==F(h)
        rows[j]=F(l)
    assert list(rows)==IDS
    return rows
def encode(x):
    n=x*S;assert n.denominator==1
    n=n.numerator
    return ('-' if n<0 else '')+f'{abs(n)//S}.{abs(n)%S:012d}'
def write(path,rows):
    path.write_text(''.join(f'{j} {encode(rows[j])} {encode(rows[j])}\n' for j in IDS))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def main(bank,first,second):
    a=load(first/'weighted_trial.txt')
    b=load(second/'trial.txt')
    assert sha(first/'weighted_trial.txt')=='27670ee547ec2f5fc30fd3adb155378267e82df6f49bbaf90701c31f1fe4c15e'
    assert sha(first/'weighted_diag.txt')=='38024951a1daf15f1a19c0298d793443458e9033e5b5620af3178a8378728126'
    assert load(first/'weighted_diag.txt')==load(second/'diag.txt')
    # The same exact finite pole projection is linear; no separate solve for +/−.
    for tag,v in [('first',a),('second',b),('plus',{j:a[j]+b[j] for j in IDS}),('minus',{j:a[j]-b[j] for j in IDS})]:
        write(bank/f'trial_{tag}.txt',v)
    assert load(bank/'trial_first.txt')==a
    p=load(bank/'trial_plus.txt');m=load(bank/'trial_minus.txt')
    assert all(p[j]==a[j]+b[j] and m[j]==a[j]-b[j] for j in IDS)
    (bank/'diag.txt').write_bytes((first/'weighted_diag.txt').read_bytes())
    print('PASS exact common response trial construction; no physical or Schur claim')
if __name__=='__main__':
    assert len(sys.argv)==4
    main(*(Path(p) for p in sys.argv[1:]))
