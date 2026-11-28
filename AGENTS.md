# Agent entrypoint — Canonical Return

This file is a provider-neutral restart guide. The governing claim and change rules remain [governance/CLAIMS.md](governance/CLAIMS.md) and [governance/WORKFLOW.md](governance/WORKFLOW.md). Read them before editing. This repository owns Canonical Return mathematics and its native evidence. Do not edit another repository's acceptance ledger through this one.

## Start from live evidence

1. Fetch all refs and inspect `main`, active `research/` branches, PRs, Actions runs, run attempts, jobs, artifacts, and custody commits. The README describes the initial import; it is not a complete live status board.
2. Read the native claim and workflow rules, then select a named CR obligation with an exact branch and commit. Preserve historical snapshots; append corrections with source links.
3. For a claimed PASS, check the exact source commit, compiler output, named axiom audit, original input identities, and remote evidence readback. A job, artifact, or PR label alone does not admit a theorem.
4. If a CR result is proposed for Weil-Lab Phase IV, prepare the full source and form-preservation packet required by [Weil-Lab issue #935](https://github.com/monocap-tech/weil-lab/issues/935). The receiving IV owner decides ADOPT, INVESTIGATE, RETURN, or HOLD. A CR result does not change an IV or formalization gate on its own.

A fresh agent may use `git`, GitHub's REST API or `gh`, Python and the pinned runner/Lean setup present on the selected branch. No ChatGPT session, Library file, or proprietary connector is required to read the repository and replay repository-hosted checks. Check external artifact availability before relying on it; preserve hashes and remote commit paths. Never commit credentials.

At a handoff, record the source commit, run/attempt/job/failed step, artifact digest, native receipt, claim level, unresolved work, and next permitted CR operation in durable repository custody. For project-wide indexing and the full restart packet format, see [provenance's continuation protocol](https://github.com/monocap-tech/provenance/blob/main/protocols/PROVIDER_NEUTRAL_CONTINUATION_v0_1.md) if access is available. Provenance is private; this local guide and the native contracts remain sufficient to work within CR when access is unavailable. Cross-repository transfer stays blocked until the receiving repository verifies it.

**Status guard:** Neither a finite positive Gram calculation nor a scoped Lean pilot is global Weil positivity or RH. Keep conditional and open premises explicit.
