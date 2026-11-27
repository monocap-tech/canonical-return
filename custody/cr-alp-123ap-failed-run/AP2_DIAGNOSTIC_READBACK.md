# CR-ALP-123AP2 — Recovered diagnostic and precision failure

Run 46060236238, job 140712443174, activation 3aa344dcbd78fa0478a054ef703b84ba2732d5b1. Both sectors completed the direct and polarization calculations; final JSON validation failed on a literal backslash-n appended to the otherwise complete JSON document. No accepted result was pushed by the runner.

Artifact 15228472238 downloaded independently, 11493 bytes, SHA256 aec5b17d6cda14d62dd880a5e49af673f3b8cc3df36789b26903cd101171763b (matches API digest), ZIP CRC PASS. Raw members:
- preflight.log 240 bytes SHA256 f917153f5ce279dc0a1c4f13d78bbf6f1f596e549e62605968e569b9b783885e
- receiver.log 3089 bytes SHA256 55e32fbfd5215e3b77168a3c1ac64a95935260318c14996b904604cc1488d477
- result.json 37276 bytes SHA256 40214e0a0c2d397b1d227f3792a312d8e36645e266c7ade1a28d2cbb113216fd

Used JSONDecoder.raw_decode and verified the entire unparsed suffix is exactly two literal characters backslash and n. This is diagnostic recovery, not acceptance of the failed runner. Full raw archive currently retained locally and as the expiring Actions artifact; this note is permanent custody of the compact diagnostic, not of the full raw payload.

Both original and final energy intervals overlap their historical intervals, all 85 conditioned residual intervals overlap, and all 6 pole-coordinate overlap checks pass. But the direct final intervals are enormously wide: cosh approximately [-544546.5634,544546.5634], sinh [-24359314.6218,24359314.6218]. Both direct positivity flags false. Historical final energies approximately 1.885132235e-38 and 3.104056888e-33. Overlap is not meaningful independent confirmation at these widths; no physical negative witness or independently established positive margin.

A concrete precision defect in the new 10^180 grid choice is exposed by the actual 046 digamma code: p is repeatedly rounded to a fixed absolute grid, then multiplied by Bernoulli coefficients in the Euler-Maclaurin loop. Exact rational recomputation of B_180 by the source recurrence gives 10^185 < |B_180| < 10^186, approximately 8.54328935788337e185. A single grid-unit error multiplied by this coefficient is approximately 854328.94, on the scale of the observed catastrophic overestimation. This explains why a lower grid is not an appropriate independent precision test for this implementation. This is a specific mechanism for enclosure inflation, not a proof that all other discrepancies are absent.

Next bounded repair: retain the independent recurrence shift 768, restore grid 10^340 (or use 10^400 with separately checked Taylor truncation and explicit outward reserves); fix the literal JSON suffix; record interval widths before accepting overlaps. Preserve failed receipts. Do not initiate a dependent response-column audit or promote the fixed-collar theorem until a narrow direct replay passes. No additional workflow activated in this check-in.
