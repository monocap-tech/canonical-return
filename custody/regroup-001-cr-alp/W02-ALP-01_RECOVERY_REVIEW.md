# W02-ALP-01 — Recovery and acceptance review
Owner explicitly authorized recovery and acceptance review only, stop before executing 123AE. Governing entrypoint e96ebfce615154a5207aafd7b5b6f8e4f99e646e; R05 27213dfdc2ec6ba78ee25045c4ceec6529da3060; predecessor 7c0a1982eb25732c5f20d041e11e29ffc2d170ec. 
Disposition: **BYTE/STRUCTURAL ACCEPTANCE; CONDITIONAL RECEIVER READY FOR SEPARATE RELEASE.** No new mathematical execution, selftest, response solve, optimization, CI dispatch, proof promotion, or source mutation. Mathematical head 123AE staged; latest executed 123AC inconclusive; strongest conditional theorem 122P/R. Both pole signs open. Delivery unchanged: monocap-tech/canonical-return, research/cr-alp-122-fixed-collar-transfer-staged. This is administrative evidence only.

## 1. Independent recovery
Direct fetch_file of 30 MB result returned empty content with correct blob metadata; not used as byte proof. Original Actions artifact 15170167880, run 45933049682, was downloaded as actual ZIP and materialized locally; Python hashlib and zipfile verified full bytes and CRC. Frozen source/custody pin 185bc2a86ac577df7afd0d32dd6a89a97125e223; original manifest blob a7dcdea2322de3d6e7023a197d73d39139904eb3. ZIP 8506138 bytes, SHA256 32c2f2e01f228d56f54b8072f5b47515969dd9bfba43dd9d5cdbbc55d583e231. Both members match frozen Git object IDs recomputed from the actual bytes (Git blob header plus SHA1), as well as SHA256SUMS:
- result.json 30378190 bytes, SHA256 a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473, Git blob b4174a176e0318c1c7b9a562fe89889729e8eedc.
- receiver.log 8659 bytes, SHA256 aaf7c8ce388a258fc0d6760006dcef1b30df54730a3631227cbd7ba3fe678db1, Git blob d618a2af557070a2addd7f1fdbb5ca8e0f92b0ce.
The artifact is independent readback of the same runner output, not an independent numerical replication. Existing Git blobs are the durable material; no new upload of duplicate 30 MB output necessary.

Structural validation via Python JSON and exact Fraction comparisons (not numerical receiver execution):
| sector | columns | E850 bank rows | prefix rows | prefix interval pairs | bank coordinates | tail receipts | recorded overlap checks |
|---|---:|---:|---:|---:|---:|---:|---:|
| cosh | 43 | 382 | 981 | 42183 | 16426 | 43 | 129 |
| sinh | 42 | 383 | 982 | 41244 | 16086 | 42 | 126 |
Verified status, exact sector keys, coordinate permutations, square Dbar/LDL shapes, positive rational pivots, bank and prefix parity index sequences, integer coordinates, grids 10^40 and 10^60, nonnegative interval radii, full tail lengths and nine signed coefficient intervals, three 32-order moment arrays with ordered endpoints, nonnegative remainder charges. Recorded 255 direct overlaps are not reexecuted. Both classifications remain SUFFICIENT_BOUND_INCONCLUSIVE.

## 2. Dependencies and runtime
Independent direct text SHA256 readbacks (ASCII, SHA256 implementation checked against abc test vector):
- 121 off-branch 7c9dc0301df765ab31f9332e9e8c76015b59894c, custody/cr-alp-074-077/replayed-121/result.json, 319828 bytes, blob 2b55baf0e81ac9a296276f396f02633553e58b28, SHA256 92b5f6bc64c74b7910540c1e39c9171345771eb8dff96e370e2e708c73d34603 matches its own manifest. Rank 85; both terminal grids recorded full_pass. Exact integer matrix canonical serialization (convert decimal strings to integer tokens, not JSON strings or JS Number): cosh ef6a2ce2b37654f5785c37d4d793030e843eeb075a1d0380bf6242cc283aa402; sinh 42c495d86072e538250e00069268a726555ad3fb17c803d5db7651904da32c3d. Both match 122P. Initial string-serialization hashes differed; corrected serialization matches; not a source discrepancy.
- 082 result SHA256 d96f10b1f83d0243a47b877001e8cd6049e12b6993abae761d1ce680e9edcfc7, blob bffbeba9ad6a3268898cbb0c19c1b80cdac9d3b7, matches its manifest.
- 091 result SHA256 c0a399bc0128e7935791978c818790c0bf43b8fb3c004e05b4c4e87e7d6ec4d7, blob c64264282b63b49ed0ec64a34bf7fa69e02267a8, matches its manifest.
- 122P result SHA256 f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de, blob 7fb998716bdb5a8f0c178585f1f8ad00e7d3fe37, matches pin.
- 123Z joint-repaired SHA256 7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897, blob 420f19e94763870bc982b0972632dd21be2bc44a, matches pin.
121 off-branch reachability confirmed, not relocated. Historical archive reassembler source 1757db4ce01a5cbbf88036552f468dc05c368734 read; full archive not reassembled here. Source 046 blob 554e7ab54da4545987d89aa660a4abe6c951cb63, gamma bands 068 blob 1a7fe8629daa8b8e6cec185a79b6746d6b28975d, 123AC blob 44620577d22ad4307142462db89053c730ca9136, 123Y a5812874b3924e1b76574211e13ab3ab2f59604d. 123AC run job 140289849262 logs directly read: ubuntu-26.04 image 20261116.202.1, NumPy 2.2.6 compiled cp314 wheel SHA256 f5764bcb46deaf79d8e1bf5247a8c7b3ff72f7f063e67459bdbffa3888e25ac2; OPENBLAS_NUM_THREADS=2 from pinned workflow e0a015db02dcb373c65d9d13668919d3e0aa6f74. Input manifest checks and final output manifest checks recorded successful. Exact Python patch version, all compiler/library versions and full prior environment not independently reconstructed. NumPy proposals untrusted; interval arithmetic acceptance separate. Historical 058 VALIDATION.md blob f45e7cc9c523a41525842728626be7fcd896fb14 specifies binary128 RN, MPFR 192-bit and directed conversions, no fast-math, g++ -O2 -std=gnu++17 ... -lquadmath -l:libmpfr.so.6 -lgmp; MPFR not bundled. Independent physical-contract validation remains open.
122O blob 64f9a1c79e45b199c2555e390dc70ef8a4f61649 read, original IV103 Weil-Lab staging 018f7e59f4129ce521cd3ae86d1afd59d3b64b5f, source blob 51a58fff565bb27453ab19fe276d49621b951563, successful run 45723601649. Source identity only, not positivity. E8 normalization and source distinction retained.

## 3. Static 123AE review
Read staged executable blob a5913d6b3462cad622f579d31e1ba0f7cec29aea, derivation blob a7dff7a3d41f39fba6ee525c06570ac317d36d5c, 123AD proof blob 010a4bfc6075f05bcbacdfbc0e045f906d82f618, 123F 55e4d78b1ef8e54f34bd1f60e51919fdb7145808, 123AB 374773999ff7f5140536c398eace37a4d8180a49. Imports read: 075 positive d8958b1f457f711b1ac14562fb6335752cc6284f; 123Z namespace df402965b51e309ba9045f170ec8ca47e0fc8fe5; transitive 046, 123Y, 074 balanced selection cd452949470a7cc049c5d23a9bae7040867cd878 and 075 nested selection 31bf1774e7f8fe98b0a9bd35e7ffc27a3af43741. Review only; no import or execution.
- Prefix uses stored ZL midpoint/radii and sqrt(Delta/d_j) with outward boxes. Gram is exact integer; diagonal reserve covers entrywise spectral rounding; strict fraction-free positive shifts supply upper bound.
- Nine-function tail matches signed coefficient order, common frequency scale and capped denominator. Existing 123AD variation/Abel reserves and row sums give Loewner upper Gram, preserving cross-column signs; coefficient box and full remainder norm charged. No new tail estimate asserted.
- Joint prefix/tail midpoint bound and stacked perturbation reserve are consistent with the 123AE derivation. Complete Q85 pole-null pairing cancels constant-cap tail in weighted pair; full pole mass (not truncated) used for denominator upper bound. Midpoint rank-one subtraction is weakened by outward perturbation reserve. Critical inherited exact Q85 projection and comparator hypotheses must hold.
- Minimum retains old trace, split spectral, joint and projected bounds; finite charge U,E unchanged. Positive lower bound if achieved would be conditional only. No result claimed.
- Static interface weakness: tail loop zip could silently truncate malformed receipt; authenticated input has exactly n entries, separately verified. Do not generalize to another input. Existing script takes expected SHA as argument; future command must use frozen digest above, not self-selected digest. No script edit under this grant.
No new blocking discrepancy identified by scoped static review. This is not formal proof certification, synthetic preflight PASS, fresh physical enclosure audit or guarantee of numerical improvement. All 122R A-D and original 123Y/AC signed-tail and comparator premises retained.

## 4. Reproduction of byte custody and structural acceptance
Download existing artifact 15170167880 (no workflow) as 123AC-frozen.zip. Use the following standalone inspection (no imports from research sources, no receiver execution).
```python
import hashlib, zipfile, json, sys
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
p=Path(sys.argv[1])
assert p.stat().st_size==8506138
assert hashlib.sha256(p.read_bytes()).hexdigest()=="32c2f2e01f228d56f54b8072f5b47515969dd9bfba43dd9d5cdbbc55d583e231"
with zipfile.ZipFile(p) as z:
    assert set(z.namelist())=={"receiver.log","result.json"} and z.testzip() is None
    for name,size,sha,blob in [
        ("result.json",30378190,"a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473","b4174a176e0318c1c7b9a562fe89889729e8eedc"),
        ("receiver.log",8659,"aaf7c8ce388a258fc0d6760006dcef1b30df54730a3631227cbd7ba3fe678db1","d618a2af557070a2addd7f1fdbb5ca8e0f92b0ce")]:
        b=z.read(name)
        assert len(b)==size and hashlib.sha256(b).hexdigest()==sha
        assert hashlib.sha1(f"blob {size}\0".encode()+b).hexdigest()==blob
    d=json.loads(z.read("result.json"))
assert d["status"]=="CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER" and set(d["sectors"])=={"1","2"}
for k,s in d["sectors"].items():
    k=int(k); n=44-k; start=87 if k==1 else 86
    bank=list(range(start,851,2)); prefix=list(range(start,2049,2))
    assert sorted(s["coordinate_order"])==list(range(n))
    assert s["bank_indices"]==bank
    assert len(s["response_transformed_coordinates"])==len(bank)
    assert all(len(row)==n and all(type(x)==int for x in row) for row in s["response_transformed_coordinates"])
    assert s["response_coordinate_grid"]==10**40 and s["prefix_grid"]==10**60
    assert [r[0] for r in s["prefix_mid_radius"]]==prefix
    assert all(len(row)==n and all(len(v)==2 and all(type(x)==int for x in v) and v[1]>=0 for v in row) for _,row in s["prefix_mid_radius"])
    assert len(s["LDL_pivots"])==n and all(F(v)>0 for v in s["LDL_pivots"])
    for key in ["Dbar","LDL_lower"]:
        assert len(s[key])==n and all(len(row)==n for row in s[key])
    assert len(s["column_prefix_inverse_uppers"])==n and len(s["tail_receipts"])==n and s["direct_prefix_overlap_checks"]==3*n
    for t in s["tail_receipts"]:
        assert len(t["pole"]["cosine_coefficients"])==7
        for lo,hi in [t["pole"]["constant"],t["pole"]["log_coefficient"],*t["pole"]["cosine_coefficients"],t["source_pole"],t["combined_pole"],t["alpha"],t["leading_combined_moment"]]:
            assert F(lo)<=F(hi)
        assert set(t["moments"])=={"centered_combined","gamma_sine","prime_sine"}
        assert all(len(m)==32 and all(F(lo)<=F(hi) for lo,hi in m) for m in t["moments"].values())
        assert all(F(t[q])>=0 for q in ["sine_series_norm_upper","tail_norm_upper"])
        assert F(t["pole"]["pole_remainder_norm_upper"])>=0
    print(k,n,len(bank),len(prefix),len(s["tail_receipts"]),"structural PASS")
```
## 5. Four-axis disposition and stop
Mathematical: conditional 122P/R unchanged; both pole signs open; original ALP phase method STOP, arbitrary support out of scope. Computational: 123AC authenticated, inconclusive; 123AE staged only. Custodial: principal 30 MB result gap resolved; 121 and named downstream source receipts independently read; full historical clean replay and physical premise acceptance still provisional. Operational: stop before 123AE, await separate explicit owner release for bounded existing post-receiver only; no expansion on inconclusive result. No cross-repository writes or promotion. This report is written to existing branch as administrative append-only evidence [skip ci]; no workflow file edited. Branch-wide workflow trigger semantics not exhaustively reaudited; no dispatch. 
