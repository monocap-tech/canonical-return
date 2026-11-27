# 123AP initial run — failed energy overlap gate

Run 46030975840, job 140614374142, activation 8ec6026e950ea1987c803c991637f5f66f6c1592, conclusion failure. Authentication passed. Independent coefficient generation reached j=2048, cosh explicit prefix completed. Failure at the joint energy overlap/positivity assertion (the combined assertion does not distinguish which failed). No result accepted; no repository result custody from this run. Artifact 15214233394, 1016 bytes, digest sha256:79a85b4b34d99c57456b7fa123ab07f9b562ec1666f8372f5f2a232e69f71318. Preserve original failed execution; do not interpret this as a negative physical witness. Diagnostic successor must record actual independently reconstructed energy interval, saved energy interval, signed differences and component comparisons before any acceptance.

```text
2026-11-27T16:11:03.1071055Z digamma complete 2020
2026-11-27T16:11:04.3705689Z digamma complete 2040
2026-11-27T16:11:12.1181913Z direct prefix 1 201
2026-11-27T16:11:15.0506448Z direct prefix 1 401
2026-11-27T16:11:17.9829433Z direct prefix 1 601
2026-11-27T16:11:20.9101567Z direct prefix 1 801
2026-11-27T16:11:23.8095838Z direct prefix 1 1001
2026-11-27T16:11:26.7090395Z direct prefix 1 1201
2026-11-27T16:11:29.5992441Z direct prefix 1 1401
2026-11-27T16:11:32.4836150Z direct prefix 1 1601
2026-11-27T16:11:35.3597912Z direct prefix 1 1801
2026-11-27T16:11:38.2281030Z direct prefix 1 2001
2026-11-27T16:11:39.9601141Z Traceback (most recent call last):
2026-11-27T16:11:39.9611029Z   File "/home/runner/work/canonical-return/canonical-return/src/audits/cr_alp_123ap_frozen_pole_trial_audit.py", line 197, in <module>
2026-11-27T16:11:39.9612257Z     if __name__=='__main__':main(*sys.argv[1:])
2026-11-27T16:11:39.9612726Z                             ^^^^^^^^^^^^^^^^^^^
2026-11-27T16:11:39.9616764Z   File "/home/runner/work/canonical-return/canonical-return/src/audits/cr_alp_123ap_frozen_pole_trial_audit.py", line 138, in main
2026-11-27T16:11:39.9618185Z     assert overlap(energy,z['joint_energy_interval']) and energy.l>0,(k,'energy')
2026-11-27T16:11:39.9618977Z            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2026-11-27T16:11:39.9619500Z AssertionError: (1, 'energy')
2026-11-27T16:11:40.0204715Z ##[error]Process completed with exit code 1.
2026-11-27T16:11:40.0289676Z ##[group]Run actions/upload-artifact@v4
2026-11-27T16:11:40.0290007Z with:
2026-11-27T16:11:40.0290245Z   name: CR-ALP-123AP-frozen-pole-trial
2026-11-27T16:11:40.0290549Z   path: /tmp/alp123ap-evidence/
```
