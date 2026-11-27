# CR-ALP-123AQ — Failed implementation run

Run 46067602745, job 140736593893, activation 3d13e7b6249f17961f09856240affe7d2d84833a. Input authentication succeeded; all cosh prefix rows and 43 columns reached their printed completion, implying 981*43=42183 prefix overlaps, 43*96=4128 moment overlaps, 129 sparse direct-prefix overlaps, and 172 pole-coordinate overlaps before the matrix stage. These are execution-log implications, not accepted output matrices. The next line at source line 186 calls mat(...), but the earlier loop 'for mat,key in [(H,"H"),(Mass,"M")]' binds mat locally to the Mass list. The error is a Python name collision, not a failed matrix-order inequality. Sinh sector not reached. No positivity promotion. The runner did not preserve result.json to the branch.

Artifact 15231389846, 2132 bytes, API SHA256 6bb5f1f9ed1fce067495f9f4b3b50fb0c13be5b948fd9fa0224e7e7fb56e90ae. The compact diagnostic is preserved here; raw runner artifact remains external and expiring.

```text
2026-11-27T21:09:52.3663449Z column 1 35 tail 0.20363965294346045 prefix 0.012258365582425197
2026-11-27T21:09:53.1681287Z column 1 36 tail 0.209885933218047 prefix 0.004815700712146425
2026-11-27T21:09:53.9144164Z column 1 37 tail 0.1580470237227367 prefix 0.006956314770254028
2026-11-27T21:09:54.6550255Z column 1 38 tail 0.2832901520407571 prefix 0.009509846014262757
2026-11-27T21:09:55.3949241Z column 1 39 tail 0.14304325991696848 prefix 0.001771517626184685
2026-11-27T21:09:56.1428725Z column 1 40 tail 0.22249943329253666 prefix 0.004776206238522403
2026-11-27T21:09:56.9079916Z column 1 41 tail 0.16943997732843807 prefix 0.002019388101945511
2026-11-27T21:09:57.6408393Z column 1 42 tail 0.1343811163586072 prefix 0.005337794296990371
2026-11-27T21:09:57.7048594Z Traceback (most recent call last):
2026-11-27T21:09:57.7057131Z   File "/home/runner/work/canonical-return/canonical-return/src/audits/cr_alp_123aq_frozen_response_audit.py", line 253, in <module>
2026-11-27T21:09:57.7058569Z     if __name__=='__main__':main(*sys.argv[1:])
2026-11-27T21:09:57.7059248Z                             ^^^^^^^^^^^^^^^^^^^
2026-11-27T21:09:57.7060306Z   File "/home/runner/work/canonical-return/canonical-return/src/audits/cr_alp_123aq_frozen_response_audit.py", line 186, in main
2026-11-27T21:09:57.7061483Z     Vmat=mat([[ib(v) for v in row] for row in Vt])
2026-11-27T21:09:57.7062038Z          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2026-11-27T21:09:57.7062589Z TypeError: 'list' object is not callable
2026-11-27T21:09:57.8015914Z ##[error]Process completed with exit code 1.
2026-11-27T21:09:57.8084549Z ##[group]Run actions/upload-artifact@v4
2026-11-27T21:09:57.8085050Z with:
2026-11-27T21:09:57.8085279Z   name: CR-ALP-123AQ-frozen-response
```

Bounded repair: rename the locally bound finite-form matrix, retain all original assertions, no new response solve, and preflight imported callable names with a static AST binding census before execution. Do not interpret the prior run as an accepted corrected-energy/Gram certificate.
