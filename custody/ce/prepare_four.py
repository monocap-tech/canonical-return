from pathlib import Path
import shutil
repo=Path(__file__).resolve().parents[2];out=repo/'ce-output';src=out/'pristine';dst=out/'trial4'
for name in ['cr_int_t1d_bp','cr_int_t1d_bs','cr_int_t1d_bt','cr_int_t1d_bu','cr_int_t1d_ca','cr_int_t1d_cb','cr_int_t1d_bk']:
    shutil.copytree(src/name,dst/name)
ca=dst/'cr_int_t1d_ca';p=ca/'bx_prepare.py';s=p.read_text()
assert s.count('[:,:3]')==5
p.write_text(s.replace('[:,:3]','[:,:4]'))
for n in ['bu_trial.npz','even_joint_proposal.npz','raw_part0.npz','raw_part1.npz','raw_part2.npz','joint_circle.npz','even_joint_Uaction.npz']:
    (ca/n).unlink()
shutil.copyfile(repo/'custody/ce/ce_four.py',dst/'ce_four.py')
print('FOUR COLUMN SCRATCH PREPARED',flush=True)
