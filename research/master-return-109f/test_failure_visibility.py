import contextlib, io, json, os, pathlib, sys, tempfile
from types import SimpleNamespace
from unittest.mock import patch
import audit_f_stage as A

for fail_stage in ['build', 'audit']:
    with tempfile.TemporaryDirectory() as tmp:
        old = os.getcwd()
        os.chdir(tmp)
        try:
            p = pathlib.Path('WeilDefect/MasterReturn/M.lean')
            p.parent.mkdir(parents=True)
            p.write_text('#print axioms M.good\n')
            e = p.parent / 'LogCarrierBoundedAnalysis.lean'
            e.write_bytes(b'fixture E')
            blob, eb = A.digest(p.read_bytes()), A.digest(e.read_bytes())
            plan = dict(allowed_axioms=['propext','Classical.choice','Quot.sound'], stages=[dict(stage='test',count=1,modules=[dict(module='WeilDefect.MasterReturn.M',blob=blob,targets=['M.good'])])])
            pathlib.Path('plan.json').write_text(json.dumps(plan))
            pathlib.Path('return109f-closure.json').write_text(json.dumps(dict(named_targets=111,admitted_E_blob=eb,base_assembler_blob=A.BASE_BLOB,supplier_count=314,staged={str(p):dict(blob=blob)})))
            def run(cmd,stdout,stderr):
                phase = 'build' if cmd[1]=='build' else 'audit'
                if phase==fail_stage:
                    stdout.write(b'error: fixture root diagnostic\n'+b'x'*20000+b'\nfixture tail\n')
                    return SimpleNamespace(returncode=1)
                return SimpleNamespace(returncode=0)
            output = io.StringIO()
            with patch.multiple(A,PLAN_BLOB=A.digest(pathlib.Path('plan.json').read_bytes()),E_BLOB=eb),patch.object(A.subprocess,'run',side_effect=run),patch.object(sys,'argv',['audit','plan.json','test','evidence']),contextlib.redirect_stdout(output):
                try: A.main()
                except AssertionError: pass
                else: raise AssertionError('failure was accepted')
            visible = output.getvalue()
            assert 'error: fixture root diagnostic' in visible
            assert 'fixture tail' in visible and len(visible)<28500
            assert not json.loads(pathlib.Path('evidence/receipt.json').read_text())['complete']
            assert pathlib.Path('evidence/WeilDefect.MasterReturn.M.'+fail_stage+'.log').stat().st_size>20000
            print('PASS bounded diagnostic visibility and fail-closed receipt:',fail_stage)
        finally: os.chdir(old)
print('2/2 PASS')
