#!/usr/bin/env python3
"""Verify custody ZIP identity and recorded ledger bytes; never promotes source replay."""
import hashlib, sys, zipfile
from pathlib import Path
p = Path(sys.argv[1])
with p.open("rb") as f:
    assert hashlib.file_digest(f,"sha256").hexdigest() == "2cd5e9004dcb9508d075b2959311ef9be69aa04805ea83aebb3d29219019056f"
assert p.stat().st_size == 67253180
expected = {
    "originals/CR-INT-T1D-CA-reproducibility(1).zip":"08510d58bcdc3dac6923e60d7a8439d24f124de71460297490528c0a21181209",
    "originals/CR-INT-T1D-CB-reproducibility(1).zip":"7225e8a6870aaaf7b87f996cc74e7200ed10f79d87c8746e65d3d165ba917911",
    "custody/inherited_sha256.txt":"c76225ab58cab65cecddb0b01a4ad05fa9a67fc1fd9a00f106b537e9c2eb2397",
    "custody/missing_inherited.tsv":"6d8e153678078ef716aebafc28c4e4b64ad5c2ad3873e9bfaa5521aaf4645844",
    "custody/verify_custody.py":"f537a42d53ea7b31c6d6dd2ee01e436c875a860b1598e1977af2a79e325b10d3",
    "custody/replay.log":"fb629eb33c7660beb7b79dc3bc1f099175992795d52213aca75e59e7c69f580d",
}
with zipfile.ZipFile(p) as z:
    assert z.testzip() is None
    for name,h in expected.items():
        with z.open(name) as f:
            assert hashlib.file_digest(f,"sha256").hexdigest() == h, name
        print("PASS",name)
print("RECOVERY PACK HASH PASS; inherited bytes NOT verified")
