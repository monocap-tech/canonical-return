"""Verify and reconstruct exact original ZIP; no mathematical execution."""
import hashlib,json,sys,zipfile
from pathlib import Path
p=Path(__file__).resolve().parent
d=json.loads((p/'recovery.json').read_text())
out=Path(sys.argv[1])
with out.open('wb') as f:
    for part in d['parts']:
        b=(p/part['path']).read_bytes()
        assert len(b)==part['bytes'] and hashlib.sha256(b).hexdigest()==part['sha256']
        f.write(b)
with out.open('rb') as f: assert hashlib.file_digest(f,'sha256').hexdigest()==d['original_archive_sha256']
assert out.stat().st_size==d['original_archive_bytes']
with zipfile.ZipFile(out) as z:
    assert z.testzip() is None and set(z.namelist())==set(d['members'])
    for n,m in d['members'].items():
        with z.open(n) as f: assert hashlib.file_digest(f,'sha256').hexdigest()==m['sha256']
        assert z.getinfo(n).file_size==m['bytes']
print('Exact original ZIP and member readback PASS')
