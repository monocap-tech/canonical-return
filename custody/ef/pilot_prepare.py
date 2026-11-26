from pathlib import Path
import json,shutil,hashlib,zipfile
import numpy as np
root=Path(__file__).resolve().parents[2];w=root/'ce-output/pristine';out=root/'custody/ef/generated';p=w/'cr_int_t1d_ef_pilot';p.mkdir(exist_ok=False)
src=root/'custody/ee/generated/sources/cr_int_t1d_ca_odd/bq_action.py'
assert src.read_bytes()==(w/'cr_int_t1d_ca_odd/bq_action.py').read_bytes()
shutil.copyfile(src,p/'bq_action.py')
shutil.copyfile(w/'cr_int_t1d_ca_odd/bn_energy.py',p/'bn_energy.py')
shutil.copyfile(root/'custody/dv/generated/merged/odd_merged_corrector.npz',p/'odd_corrector.npz')
a=np.load(p/'odd_corrector.npz');assert a['coefficient'].shape==(491,352)
# Recover original head construction scripts from authenticated archives.
# Preserve all matches with their provenance rather than selecting an unverified variant.
dest=out/'head-sources';dest.mkdir(exist_ok=True)
records=[]
for stage in ['BO','BP']:
    archive=root/'ce-output/archives'/f'CR-INT-T1D-{stage}-reproducibility.zip'
    with zipfile.ZipFile(archive) as z:
        for member in z.namelist():
            name=Path(member).name
            if name not in {'bp_head.py','bo_head.py','bo_head_extend.py','bo_normalize.py','bo_tail.py','bs_symbol.py'}:continue
            data=z.read(member);target=dest/stage/member
            assert not Path(member).is_absolute() and '..' not in Path(member).parts
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
            records.append(dict(archive=archive.name,member=member,sha256=hashlib.sha256(data).hexdigest()))
# Symbol generator may live in a later archive. Preserve literal source if available.
for f in sorted(w.rglob('*symbol*.py')):
    target=dest/'recovered'/f.relative_to(w);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,target)
    records.append(dict(recovered=str(f.relative_to(w)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(out/'source-inventory.json').write_text(json.dumps(records,indent=2)+'\n')
print('EF SOURCE INVENTORY',records,flush=True)
