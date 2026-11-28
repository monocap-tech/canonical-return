#!/usr/bin/env python3
"""Restore the pinned CA/CB inherited files from ZIP archives by exact SHA-256.
No overwrite of mismatched files. This checks the recorded ledger only, not full
transitive closure. Original CA/CB archives should be extracted separately.
"""
import argparse
import csv
import hashlib
import shutil
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path, PurePosixPath


def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives', type=Path, help='folder containing original ZIPs')
    p.add_argument('workspace', type=Path, help='extracted CA/CB sibling workspace')
    p.add_argument('--ledger', type=Path, required=True)
    p.add_argument('--report', type=Path, default=Path('recovery_report.tsv'))
    args = p.parse_args()
    if not args.archives.is_dir():
        p.error('archive folder does not exist')
    args.workspace.mkdir(parents=True, exist_ok=True)
    root = args.workspace.resolve()
    rows = list(csv.DictReader(args.ledger.open(encoding='utf-8', newline=''), delimiter='\t'))
    expected = {}
    for row in rows:
        name, sha = row['relative_path'], row['sha256']
        rel = PurePosixPath(name)
        if rel.is_absolute() or '..' in rel.parts or '\\' in name or ':' in name or len(sha) != 64:
            raise ValueError(f'Invalid ledger entry: {name}')
        if name in expected and expected[name] != sha:
            raise ValueError(f'Conflicting ledger entry: {name}')
        expected[name] = sha
    status, origins, by_base = {}, {}, defaultdict(list)
    for name, sha in expected.items():
        dest = root / name
        if not dest.resolve().is_relative_to(root):
            raise ValueError(f'Unsafe destination: {dest}')
        status[name] = ('PASS' if digest(dest) == sha else 'MISMATCH') if dest.is_file() else ('MISMATCH' if dest.exists() else 'MISSING')
        if status[name] == 'MISSING':
            by_base[PurePosixPath(name).name].append(name)
    archives = sorted(args.archives.rglob('*.zip'))
    inventory = args.report.with_name(args.report.stem + '_archives.tsv')
    with inventory.open('w', encoding='utf-8', newline='') as log:
        w = csv.writer(log, delimiter='\t')
        w.writerow(['archive', 'bytes', 'sha256', 'crc'])
        for archive in archives:
            sha = digest(archive)
            try:
                with zipfile.ZipFile(archive) as z:
                    bad = z.testzip()
                    w.writerow([str(archive), archive.stat().st_size, sha, 'PASS' if bad is None else f'FAIL:{bad}'])
                    if bad is not None:
                        continue
                    for member in z.infolist():
                        if member.is_dir():
                            continue
                        names = [n for n in by_base.get(PurePosixPath(member.filename).name, []) if status[n] == 'MISSING']
                        if not names:
                            continue
                        with tempfile.TemporaryDirectory() as temp:
                            tmp = Path(temp) / 'payload'
                            with z.open(member) as src, tmp.open('wb') as dst:
                                shutil.copyfileobj(src, dst)
                            h = digest(tmp)
                            for name in names:
                                if expected[name] != h:
                                    continue
                                dest = root / name
                                dest.parent.mkdir(parents=True, exist_ok=True)
                                if not dest.resolve().is_relative_to(root):
                                    raise ValueError(f'Unsafe destination: {dest}')
                                with tmp.open('rb') as src, dest.open('xb') as dst:
                                    shutil.copyfileobj(src, dst)
                                if digest(dest) != h:
                                    raise IOError(f'Copy verification failed: {dest}')
                                status[name] = 'PASS'
                                origins[name] = f'{archive}!{member.filename}'
            except (OSError, zipfile.BadZipFile, RuntimeError) as e:
                w.writerow([str(archive), archive.stat().st_size, sha, f'ERROR:{e}'])
            log.flush()
    with args.report.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['relative_path', 'sha256', 'status', 'recovered_from'])
        for name, sha in sorted(expected.items()):
            w.writerow([name, sha, status[name], origins.get(name, '')])
    for s in ('PASS', 'MISSING', 'MISMATCH'):
        print(f'{s}: {sum(v == s for v in status.values())}')
    print('Reports:', args.report, inventory)
    print('Recorded ledger only; not a full source replay or proof of dependency closure.')
    return 0 if all(v == 'PASS' for v in status.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
