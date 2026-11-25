# MASTER-RETURN 083–087 archive custody

This branch preserves four retained ZIPs byte for byte. `RETURN084_audit.zip` and
`RETURN087_audit.zip` are single files. `RETURN083_reproduction.zip` and
`RETURN086_audit.zip` are split into numbered consecutive parts solely to fit
the repository upload limit. The part files are raw byte slices, with no
recompression or content changes.

From this directory run `bash reconstruct.sh` to verify every part, reconstruct
the two ZIPs under `reconstructed/`, and verify all four complete ZIP hashes.
An optional output directory can be passed as the first argument.

`SHA256SUMS.parts` records each raw part. `SHA256SUMS.archives` records the
complete original ZIPs, including the two single-file archives here.

This is source custody only. It does not promote a numerical or mathematical
claim, establish global Weil positivity, or prove RH.
