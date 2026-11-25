#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$root"
sha256sum -c SHA256SUMS.parts
out="${1:-reconstructed}"
mkdir -p "$out"
cat RETURN083/RETURN083_reproduction.zip.part* > "$out/RETURN083_reproduction.zip"
cat RETURN086/RETURN086_audit.zip.part* > "$out/RETURN086_audit.zip"
cp RETURN084_audit.zip RETURN087_audit.zip "$out/"
(cd "$out" && sha256sum -c "$root/SHA256SUMS.archives")
