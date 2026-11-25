# T1D archive parts: eight recovered archives

Each directory below contains consecutive `part-0000.bin`, `part-0001.bin`, ... chunks of the original ZIP, at most 780000 bytes each. Concatenate them in lexical order to reconstruct the named ZIP and check SHA-256 and exact size against this table. The chunks are raw bytes; none is a standalone ZIP.

| Archive | Bytes | Parts | SHA-256 |
| --- | ---: | ---: | --- |
| CR-INT-T1D-BI-even-projected-checkpoint.zip | 111791850 | 144 | `649b2ff7d909fa85bb3319489a25d236fb9285fb9c08f42e0f8240ff0e0291e9` |
| CR-INT-T1D-BI-odd-projected-checkpoint.zip | 115503553 | 149 | `184ea2bf757fd6f85958ac5e8df296c207cf73115e65bb91d6de9686852a8d7c` |
| CR-INT-T1D-BP-reproducibility.zip | 8887439 | 12 | `e82efaa2bd8360221d574e43d0fa5eb55f7de281f180d83eca6f1adc951624be` |
| CR-INT-T1D-BS-reproducibility.zip | 8112156 | 11 | `892a7c143a79cf87c608dc3bf3eb89f30b8f3f7f33c2864a83d25ce24559ff8c` |
| CR-INT-T1D-BT-reproducibility.zip | 4884132 | 7 | `a67ec45b190f9de3a69aca2abcaa933229efaa2371f1672a4f9f40848c7f2aea` |
| CR-INT-T1D-BU-reproducibility.zip | 35864643 | 46 | `55ab94bd18a384b019fd66d0de0989076b0a1a726556469db1418fab93800bb2` |
| CR-INT-T1D-CA-reproducibility.zip | 45359073 | 59 | `08510d58bcdc3dac6923e60d7a8439d24f124de71460297490528c0a21181209` |
| CR-INT-T1D-CB-reproducibility.zip | 21832298 | 28 | `7225e8a6870aaaf7b87f996cc74e7200ed10f79d87c8746e65d3d165ba917911` |

Example from the repository root:

```bash
name=CR-INT-T1D-BP-reproducibility.zip
cat "custody/T1D-archives/parts/$name"/part-*.bin > "$name"
sha256sum "$name"
```

These hashes attest to archive identity after reconstruction; mathematical or source-level verification remains a separate task.
