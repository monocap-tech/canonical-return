# RETURN088–099 evidence custody

Original Actions ZIPs split into ordered parts of at most 16 MiB. Each source ZIP was checked against its pinned size, SHA-256 and CRC before splitting; the concatenation was checked again. See manifest.json. Reconstruct by concatenating archive.zip.part* in each directory, then verify the complete ZIP SHA-256. No proof promotion or merge.
