# CR-ALP-055–059 repository recovery

This directory preserves the **original CR-ALP-058 cumulative ZIP bytes**, split into 90 Git blobs because the repository writer could not upload the large binary as one file. The ZIP contains the complete 055–058 replay tree (93 entries), including the predecessor banks and inputs absent from the earlier partial 058 text recovery. It is a draft research package, not a Lean or RH acceptance claim. The independent 059 threshold script is at `src/audits/cr_alp_059_moment_threshold_audit.py` and does not evaluate the second moment.

From a fresh checkout of the CR-ALP draft branch:

```sh
python3 src/audits/cr_alp_058_archive/reassemble.py /tmp/CR-ALP-058_global_physical_witness_moments.zip
mkdir -p /tmp/cr-alp-replay
cd /tmp/cr-alp-replay
python3 -m zipfile -e /tmp/CR-ALP-058_global_physical_witness_moments.zip .
python3 cr_alp_055/check_saved.py
python3 cr_alp_056/check_saved.py
python3 cr_alp_056/check_saved.py _extended
python3 cr_alp_057/verify_group.py
python3 cr_alp_057/assemble.py
python3 cr_alp_058/validate_geometry.py
python3 cr_alp_058/check_saved.py
```

The reassembler checks size, SHA-256 and every ZIP CRC before installing the output. Recovered original ZIP size: **62,872,042 bytes**; SHA-256: `f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2`. Each blob is 700,000 bytes except the last. The ZIP itself contains stage SHA256SUMS files. The original Library item is `libfile_cd3e46516dbc8191a0734cb871c2faf6`, version 2; the distinct 055, 056 and 057 archives passed CRC checks with respective SHA-256 `42f00009952266dbf323016a9c3e1aff466d57fc3878d0f1b5dd3bb9eb48051b`, `7428cc9e624503a1dad962652028b814bc9dcc6bbb321dfa607f03bec7fc17b4`, and `9c21a614f100d754ae712a9657c3f24778eb6e574d1c0bbe1867b2d182aca21f`.

Recovered locally on 2026-11-24 Pacific time: all four original ZIP CRC checks passed; embedded 056–058 SHA256SUMS checked without missing or mismatched entries; the saved 055, 056 (including extended), 057 group and assembly, and 058 geometry and endpoint checks passed. These checks replay saved arithmetic and documented error budgets; they do not independently prove every inherited analytic identity or compute 059's second moment. Keep this PR draft and nonpromoting.
