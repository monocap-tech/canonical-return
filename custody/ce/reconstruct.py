import csv, hashlib, re, shutil, zipfile
from pathlib import Path, PurePosixPath
repo=Path(__file__).resolve().parents[2]
out=repo/'ce-output'; out.mkdir(exist_ok=True)
archives=out/'archives'; archives.mkdir(exist_ok=True)
work=out/'pristine'; work.mkdir(exist_ok=True)
parts=repo/'custody/T1D-archives/parts'
def digest(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
known={}
for line in (parts/'README.md').read_text().splitlines():
    m=re.fullmatch(r'\| (CR-INT-T1D-[^ ]+\.zip) \| (\d+) \| (\d+) \| .([0-9a-f]{64}). \|',line)
    if m: known[m[1]]=(int(m[2]),int(m[3]),m[4])
assert len(known)==8,known
with (repo/'custody/T1D-archive-inventory.tsv').open() as f:
    inv={r['archive']:int(r['bytes']) for r in csv.DictReader(f,delimiter='\t') if not r['archive'].endswith('custody-recovery.zip')}
assert len(inv)==56
assert {p.name for p in parts.iterdir() if p.is_dir()}==set(inv)
af=[]
with (out/'archive-verification.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['name','bytes','parts','sha256','original_hash_status','crc'])
    for name,size in sorted(inv.items()):
        pp=sorted((parts/name).glob('part-*.bin'))
        assert len(pp)==(size+779999)//780000,name
        assert all(p.name==f'part-{i:04d}.bin' and p.stat().st_size==min(780000,size-i*780000) for i,p in enumerate(pp)),name
        target=archives/name
        with target.open('wb') as g:
            for p in pp:
                with p.open('rb') as h:shutil.copyfileobj(h,g)
        assert target.stat().st_size==size
        sha=digest(target)
        if name in known: assert known[name]==(size,len(pp),sha),name
        with zipfile.ZipFile(target) as z:
            assert z.testzip() is None,name
            af.extend((name,x) for x in z.namelist() if PurePosixPath(x).name=='canonical_return_t1d_af_neumann_ritz.py')
        w.writerow([name,size,len(pp),sha,'MATCH' if name in known else 'NEW_DIGEST_NOT_INDEPENDENTLY_PINNED','PASS']);f.flush()
        print('ZIP PASS',name,size,sha,flush=True)
(out/'af-locations.txt').write_text('\n'.join(f'{n}\t{x}' for n,x in af)+'\n')
print('AF source locations',af,flush=True)
for stage in ['BP','BS','BT','BU','CA','CB','BK']:
    name=f'CR-INT-T1D-{stage}-reproducibility.zip'
    dest=work if stage in ['CA','CB'] else work/f'cr_int_t1d_{stage.lower()}'
    dest.mkdir(exist_ok=True)
    with zipfile.ZipFile(archives/name) as z:
        seen=set()
        for info in z.infolist():
            rel=PurePosixPath(info.filename)
            assert not rel.is_absolute() and '..' not in rel.parts and '\\' not in info.filename,info.filename
            assert info.filename not in seen,info.filename
            seen.add(info.filename)
            z.extract(info,dest)
for stage,manifest in [('bp','bp_sha256.txt'),('bs','bs_sha256.txt'),('bt','bt_sha256.txt'),('bu','bu_sha256.txt'),('ca','ca_sha256.txt'),('cb','cb_sha256.txt')]:
    folder=work/f'cr_int_t1d_{stage}'
    for line in (folder/manifest).read_text().splitlines():
        h,n=line.split(maxsplit=1);assert digest(folder/n.strip())==h,(stage,n)
    print('ORIGINAL MANIFEST PASS',stage,flush=True)
ca=work/'cr_int_t1d_ca'
expected={n:h for h,n in (line.split(maxsplit=1) for line in (ca/'verified_input_sha256.txt').read_text().splitlines())}
for n in ['even/physical_dictionary.npz','even/exceptional.npz','even/pole.npz']:
    rel='cr_int_t1d_bk/'+n
    assert digest(work/rel)==expected[rel],rel
    print('BK ORIGINAL HASH PASS',rel,flush=True)
