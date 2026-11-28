#!/usr/bin/env python3
"""Verify the stored conservative shifted terminal comparators, NOT source/action generation."""
from fractions import Fraction as F
from pathlib import Path
import json
d = json.loads(Path(__file__).with_name("CA-CB-terminal-matrices.json").read_text())
assert d["scale"] == 10**12
assert [(r["name"], r["margin"], r["n"]) for r in d["records"]] == [
    ("CA even 1","9229/10000",1), ("CA even 2","2749/5000",2),
    ("CA even 3","1427/5000",3), ("CB odd 12","2253/50000",12),
    ("CA odd 11","153/1000",11)]
for r in d["records"]:
    n = r["n"]
    t = r["upper_triangle"]
    assert len(t) == n*(n+1)//2 and all(type(x) is int for x in t)
    A = [[F(0) for _ in range(n)] for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i,n):
            A[i][j] = A[j][i] = F(t[k], d["scale"])
            k += 1
    for i in range(n):
        p = A[i][i]
        assert p > 0, (r["name"],i,p)
        for j in range(i+1,n):
            for k in range(j,n):
                A[k][j] -= A[k][i]*A[j][i]/p
                A[j][k] = A[k][j]
    print(r["name"],r["margin"],"PASS")
print("TERMINAL COMPARATOR PASS ONLY; no source/action replay")
