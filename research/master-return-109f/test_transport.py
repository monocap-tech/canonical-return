import json, pathlib, tempfile, subprocess, sys
from unittest.mock import patch
import transport as T

def case(fault):
    with tempfile.TemporaryDirectory() as tmp:
        root=pathlib.Path(tmp)
        pins={}
        for i in range(3):
            files={f'p{j}':{'blob':T.digest(str(j).encode()),'content':str(j)} for j in range(i*7,min((i+1)*7,20))}
            if i==2:
                if fault=='duplicate': files['p0']={'blob':T.digest(b'0'),'content':'0'}
                if fault=='missing': files.pop('p19')
                if fault=='entry': files['p19']['content']='wrong'
            b=json.dumps(dict(canonical='wrong' if fault=='revision' and i==2 else T.CANONICAL,files=files)).encode()
            name=f'bundle-{i}.json'; (root/name).write_bytes(b); pins[name]=T.digest(b)
        if fault=='bundle': pins['bundle-0.json']='wrong'
        scripts={'prepare_f_assembly.py':T.digest(b'fixture')}; (root/'prepare_f_assembly.py').write_bytes(b'wrong' if fault=='script' else b'fixture')
        with patch.multiple(T,BUNDLES=pins,SCRIPTS=scripts):
            try: files=T.verified_files(root)
            except AssertionError: assert fault
            else: assert not fault and len(files)==20
for f in [None,'duplicate','missing','entry','revision','bundle','script']:
    case(f); print('PASS',f)
p=subprocess.run([sys.executable,'-O',str(pathlib.Path(T.__file__))],capture_output=True)
assert p.returncode and b'verification requires assertions enabled' in p.stderr
print('PASS optimized')
print('8/8 PASS')
