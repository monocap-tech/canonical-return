#!/usr/bin/env python3
"""Preserve authenticated original Actions ZIPs, never promote proof status."""
import hashlib,json,os,subprocess,sys,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
inventory=root/'research/master-return-099c/custody_inventory.json'
assert subprocess.check_output(['git','hash-object',str(inventory)],text=True).strip()=='7e544f48f44ac42581af77a12c39aa80f7a7e80d'
d=json.loads(inventory.read_text());assert d['repository']==os.environ['GITHUB_REPOSITORY']=='monocap-tech/canonical-return'
assert len(d['artifacts'])==20 and d['total_bytes']==128209266
target=root/'archive/master-return-088-099'
assert not target.exists();target.mkdir()
work=Path(sys.argv[1]);work.mkdir(parents=True,exist_ok=True)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
records=[]
for a in d['artifacts']:
 label=f"{a['id']}-{a['name']}";assert all(c.isalnum() or c in '-_' for c in label)
 f=work/(label+'.zip')
 with f.open('wb') as fd:
  subprocess.run(['gh','api',f"repos/{d['repository']}/actions/artifacts/{a['id']}/zip"],stdout=fd,check=True)
 assert f.stat().st_size==a['bytes'] and sha(f)==a['sha256'],label
 with zipfile.ZipFile(f) as z:assert z.testzip() is None,label
 folder=target/label;folder.mkdir()
 parts=[];whole=hashlib.sha256();total=0
 with f.open('rb') as fd:
  while chunk:=fd.read(16*1024*1024):
   p=folder/f'archive.zip.part{len(parts):03d}';p.write_bytes(chunk)
   parts.append({'path':str(p.relative_to(target)),'bytes':len(chunk),'sha256':sha(p)})
 for p in parts:
  data=(target/p['path']).read_bytes();assert len(data)==p['bytes'] and hashlib.sha256(data).hexdigest()==p['sha256']
  whole.update(data);total+=len(data)
 assert whole.hexdigest()==a['sha256'] and total==a['bytes']
 records.append({**a,'parts':parts});print('custody source authenticated',label,total,len(parts),flush=True)
manifest={'schema':'return100-custody-v1','scope':'byte-identical preservation of original CI ZIPs; not mathematical proof promotion','source_inventory_git_blob':'7e544f48f44ac42581af77a12c39aa80f7a7e80d','archives':records,'total_bytes':sum(a['bytes'] for a in records)}
(target/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(target/'README.md').write_text('# RETURN088–099 evidence custody\n\nOriginal Actions ZIPs split into ordered parts of at most 16 MiB. Each source ZIP was checked against its pinned size, SHA-256 and CRC before splitting; the concatenation was checked again. See manifest.json. Reconstruct by concatenating archive.zip.part* in each directory, then verify the complete ZIP SHA-256. No proof promotion or merge.\n')
print('RETURN100 ALL 20 ARCHIVES AUTHENTICATED AND REASSEMBLED',flush=True)
