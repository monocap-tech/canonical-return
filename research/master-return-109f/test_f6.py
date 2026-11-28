import contextlib, importlib.util, io, json, os, pathlib, subprocess, sys, tempfile, unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT = pathlib.Path(__file__).parent

def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

P, A = load('prepare_f_assembly'), load('audit_f_stage')

class TestGenerator(unittest.TestCase):
    def run_case(self, fault=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            base, out = root/'base.py', root/'assembled.py'
            added = {f'supplier{i}': {'blob': f'blob{i}'} for i in range(31)}
            root_entries = [{'path': f'supplier{i}', 'blob': f'blob{i}', 'inE2': False} for i in range(4)]
            initial = {f'base{i}': {} for i in range(283)}
            entries = {k: {'blob': v['blob'], 'required_by': ['IV']} for k,v in added.items()}
            if fault == 'delta': entries['unexpected'] = {'blob':'x','required_by':['IV']}
            if fault == 'supplier_blob': entries['supplier0']['blob'] = 'wrong'
            if fault == 'lineage': entries['supplier5']['required_by'] = ['RECON']
            if fault == 'collision': collision = "pathlib.Path('WeilDefect/MasterReturn/M0.lean').touch()\n"
            else: collision = ''
            base.write_text('import pathlib,json,hashlib\n' +
                'def digest(b):\n    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\\0"+b).hexdigest()\n' +
                f'manifest={initial!r}\nentries={entries!r}\n' +
                'def visit(p, side):\n    manifest.update(entries)\n' +
                '# No source files are written until the entire selected closure passes.\n' +
                "pathlib.Path('WeilDefect/MasterReturn').mkdir(parents=True)\n"+collision)
            files, mods = {}, []
            for i in range(16):
                names = [f'M{i}.t{j}' for j in range(7 if i < 15 else 6)]
                p = f'research/master-return-109f/M{i}.lean'
                b = ('import Mathlib\n'+'\n'.join('#print axioms '+n for n in names)+'\n').encode()
                files[p] = b
                mods.append({'path':p,'blob':P.digest(b),'names':names,'imports':['Mathlib']})
            if fault == 'names': mods[0]['names'][0] = 'wrong'
            if fault == 'source_blob': mods[0]['blob'] = 'wrong'
            if fault == 'imports': mods[0]['imports'] = ['Wrong']
            data = {'modules':mods,'total':111,'unique':111}
            roots = {'roots':root_entries,'base_manifest_count':283,'IV':'i','RECON':'r'}
            disc = {'pending':[],'added':added,'base_count':283,'IV':'i','RECON':'r'}
            for n,obj in [('F_NAMED_AUDIT_MANIFEST',data),('F_SUPPLIER_ROOTS',roots),('F_CLOSURE_DISCOVERY',disc)]:
                files[str(P.DIR/(n+'.json'))] = json.dumps(obj).encode()
            def fake_git(args): return files[args[-1].split(':',1)[1]]
            consts = {'BASE_BLOB':P.digest(base.read_bytes()) if fault!='base' else 'wrong',
                      'MANIFEST_BLOB':P.digest(files[str(P.DIR/'F_NAMED_AUDIT_MANIFEST.json')]),
                      'ROOTS_BLOB':P.digest(files[str(P.DIR/'F_SUPPLIER_ROOTS.json')]),
                      'DISCOVERY_BLOB':P.digest(files[str(P.DIR/'F_CLOSURE_DISCOVERY.json')])}
            with patch.multiple(P, **consts), patch.object(P.subprocess,'check_output',side_effect=fake_git), patch.object(sys,'argv',['gen',str(base),str(root),str(out)]), contextlib.redirect_stdout(io.StringIO()):
                if fault in ('base','names','source_blob','imports'):
                    with self.assertRaises(AssertionError): P.main()
                    self.assertFalse(out.exists()); return
                P.main()
            p = subprocess.run([sys.executable,str(out)],cwd=root,capture_output=True,text=True)
            if fault:
                self.assertNotEqual(p.returncode,0,p.stdout)
                self.assertFalse((root/'return109f-closure.json').exists()); return
            self.assertEqual(p.returncode,0,p.stderr)
            rec=json.loads((root/'return109f-closure.json').read_text())
            self.assertEqual(len(rec['staged']),16); self.assertEqual(rec['named_targets'],111)
            for m in mods:
                self.assertEqual(P.digest((root/'WeilDefect/MasterReturn'/pathlib.Path(m['path']).name).read_bytes()),m['blob'])
    def test_valid(self): self.run_case()
    def test_base_drift(self): self.run_case('base')
    def test_name_drift(self): self.run_case('names')
    def test_source_drift(self): self.run_case('source_blob')
    def test_import_drift(self): self.run_case('imports')
    def test_unexpected_delta(self): self.run_case('delta')
    def test_supplier_blob(self): self.run_case('supplier_blob')
    def test_lineage(self): self.run_case('lineage')
    def test_destination_collision(self): self.run_case('collision')

class TestAudit(unittest.TestCase):
    def run_case(self, fault=None):
        with tempfile.TemporaryDirectory() as tmp:
            old=os.getcwd(); os.chdir(tmp)
            try:
                source=pathlib.Path('WeilDefect/MasterReturn/M.lean'); source.parent.mkdir(parents=True)
                names=['M.good','M.other']; source.write_text('\n'.join('#print axioms '+x for x in names)+'\n')
                blob=A.digest(source.read_bytes())
                plan={'allowed_axioms':['propext','Classical.choice','Quot.sound'],'stages':[{'stage':'test','count':2,'modules':[{'module':'WeilDefect.MasterReturn.M','blob':blob,'targets':names}]}]}
                pathlib.Path('plan.json').write_text(json.dumps(plan)); pathlib.Path('return109f-closure.json').write_text(json.dumps({'named_targets':111,'staged':{str(source):{'blob':blob}}}))
                if fault=='source': source.write_text('changed')
                if fault=='closure': pathlib.Path('return109f-closure.json').write_text(json.dumps({'named_targets':111,'staged':{str(source):{'blob':'wrong'}}}))
                def fake_run(cmd,stdout,stderr):
                    if cmd[1]=='build': return SimpleNamespace(returncode=1 if fault=='build' else 0)
                    log="'M.good' depends on axioms: [propext,\n Classical.choice, Quot.sound]\n'M.other' does not depend on any axioms\n"
                    if fault=='forbidden': log=log.replace('Quot.sound','sorryAx')
                    if fault=='missing': log=log.split("'M.other'")[0]
                    if fault=='duplicate': log+="'M.other' does not depend on any axioms\n"
                    stdout.write(log.encode()); return SimpleNamespace(returncode=1 if fault=='audit' else 0)
                with patch.object(A,'PLAN_BLOB',A.digest(pathlib.Path('plan.json').read_bytes()) if fault!='plan' else 'wrong'), patch.object(A.subprocess,'run',side_effect=fake_run), patch.object(sys,'argv',['audit','plan.json','test','evidence']),contextlib.redirect_stdout(io.StringIO()):
                    if fault:
                        with self.assertRaises(AssertionError): A.main()
                    else: A.main()
                rec=pathlib.Path('evidence/receipt.json')
                if fault=='plan': self.assertFalse(rec.exists()); return
                self.assertEqual(json.loads(rec.read_text())['complete'],fault is None)
            finally: os.chdir(old)
    def test_valid(self): self.run_case()
    def test_plan_drift(self): self.run_case('plan')
    def test_source_drift(self): self.run_case('source')
    def test_closure_drift(self): self.run_case('closure')
    def test_build_error(self): self.run_case('build')
    def test_audit_error(self): self.run_case('audit')
    def test_forbidden_axiom(self): self.run_case('forbidden')
    def test_missing_audit(self): self.run_case('missing')
    def test_duplicate_audit(self): self.run_case('duplicate')

if __name__=='__main__': unittest.main(verbosity=2)
