#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sha256sum -c SHA256SUMS.parts
out="${1:-reconstructed}"
mkdir -p "$out"
cat RETURN083/RETURN083_reproduction.zip.part* > "$out/RETURN083_reproduction.zip"
cat RETURN086/RETURN086_audit.zip.part* > "$out/RETURN086_audit.zip"
cp RETURN084_audit.zip RETURN087_audit.zip "$out/"
(cd "$out" && sha256sum -c ../SHA256SUMS.archives)
