#!/usr/bin/env python3
"""Reassemble and verify the archived CR-ALP-055–058 package."""
from __future__ import annotations

import hashlib
import pathlib
import sys
import zipfile

HERE = pathlib.Path(__file__).resolve().parent
PARTS = HERE / "parts"
NAME = "CR-ALP-058_global_physical_witness_moments.zip"
SIZE = 62_872_042
SHA256 = "f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2"
COUNT = 90
PART_SIZE = 700_000


def main() -> None:
    output = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / NAME
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_name(output.name + ".partial")
    digest = hashlib.sha256()
    total = 0
    try:
        with temp.open("wb") as dest:
            for i in range(COUNT):
                part = PARTS / f"part-{i:03d}.bin"
                data = part.read_bytes()
                expected_size = PART_SIZE if i < COUNT - 1 else SIZE - PART_SIZE * (COUNT - 1)
                if len(data) != expected_size:
                    raise ValueError(f"{part.name}: expected {expected_size} bytes, got {len(data)}")
                dest.write(data)
                digest.update(data)
                total += len(data)
        if total != SIZE or digest.hexdigest() != SHA256:
            raise ValueError(f"archive digest/size mismatch: {total} bytes, {digest.hexdigest()}")
        with zipfile.ZipFile(temp) as archive:
            bad = archive.testzip()
            if bad is not None:
                raise ValueError(f"ZIP CRC failed: {bad}")
            print(f"Verified {NAME}: {total} bytes, SHA256 {SHA256}, {len(archive.namelist())} entries")
        temp.replace(output)
        print(output)
    finally:
        temp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
