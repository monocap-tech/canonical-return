# Recovery scan 001 — 59 byte-identical inherited files restored

The pinned CA/CB archive members were scanned by basename and SHA-256 to restore files under their recorded original paths, including inherited scripts copied forward into later gates. **59 of 578** now verified locally; **519** remain absent, zero mismatches. The five immediate BP/BT/BU/BS blockers are still missing. Re-executed local independent verifier: 166/166 local manifest entries PASS, 59/578 inherited PASS, five terminal matrices and seven matched arrays PASS. No source/action replay or CC advancement.

- recovery_scan_001.tsv preserves all 578 path/hash/status rows.
- recover_t1d_dependencies.py (Python 3.11+) scans a folder of actual ZIP bytes, verifies ZIP integrity and archive SHA-256, restores only exact matching expected bytes, refuses overwrite, emits archive and per-path reports. It does not check unrecorded or recursive dependencies, or extract nested ZIPs. Tested with synthetic matched/missing/mismatched files and the actual CA/CB inputs.
- The recovered 59 bytes remain recoverable from the two original CA/CB archives; this commit contains their inventory, not the original ZIPs.

Next obtain original BP, BS, BT and BU ZIPs, then the full ancestor inventory; use GitHub release assets or LFS for durable raw-byte storage and verify downloads. No promotion until actual clean-environment replay. Conditional bound 225 unchanged.
