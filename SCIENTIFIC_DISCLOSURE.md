# Scientific Disclosure Charter

**Adaptive Investigative Governance in AI-Assisted Mathematical Research**  
*A longitudinal case study of coordinated reasoning, formal verification and research continuity through an investigation of Weil positivity and the Riemann hypothesis.*

**Document:** SD-CHARTER · Revision 0.1 · 2026-11-27  
**Standing:** Planning and editorial governance; manuscripts not released; no new mathematical claim admitted.  
**Custody:** `monocap-tech/canonical-return` is the disclosure entry point. The research repositories and their existing lineages remain authoritative for their own mathematical claims.

## 1. Purpose and jurisdiction

Produce a coherent, independently examinable scientific disclosure of the mathematical research program and its AI-assisted investigative methodology. The package must be useful to mathematicians, formal-methods researchers, AI research engineers, and institutional scientific leadership independently of any prospective collaboration or RH outcome.

The disclosure initiative is an editorial and evidence-reconciliation function, **not an eighth mathematical investigation**. It consumes source-pinned evidence from the seven existing lineages (MASTER-RETURN, MASTER-RECON, IV106, IV-AM, CR-ALP, CR-INT, NF618). It does not reassign their mathematical responsibilities, authorize computational work, activate runners, alter admission ledgers, or confer new theorem standing. Closed comparative work, including CROSS-PW/CROSS-SD, is treated as source material rather than an active research front.

The project does not claim a proof of RH, global Weil positivity, completed structural integration, or established priority for the composite criterion. Institutional interest, affiliation, endorsement, peer review and external reproducibility must never be implied without evidence.

Existing [claim classifications](governance/CLAIMS.md) and [research workflow](governance/WORKFLOW.md) remain in force. Historical archives are preserved; corrections enter as explicit revisions, not retroactive substitutions.

## 2. Document series and deliberation mandates

| ID | Deliverable | Principal audience | Target scope (indicative) |
|---|---|---|---|
| SD-00 | Executive Scientific Brief | Scientific leadership | 2–3 pages |
| SD-01 | Mathematical Research Program and Current Results | Mathematicians | 15–25 pages |
| SD-02 | Adaptive Investigative Governance | AI researchers | 15–25 pages |
| SD-03 | Multi-Agent Research Architecture | Research engineers | 10–20 pages |
| SD-04 | Formal Verification and Research Provenance | Formal methods specialists | 15–30 pages |
| SD-05 | Longitudinal Operational Findings | AI infrastructure teams | 8–15 pages |
| SD-06 | Mathematical Originality and Prior-Art Assessment | Scientific reviewers | 10–15 pages |
| SD-07 | Reproducibility and Evidence Compendium | Independent auditors | As required |

All eight documents are **planned, not drafted or released**. Page ranges are editorial targets, not quotas.

### SD-00 — Executive Scientific Brief
Decision-oriented account of the research question, substantive results and limitations, why the investigative method merits examination, and what institutional dialogue could address. It must be understandable without specialized Weil-theory knowledge. Gate: every factual claim is traceable to a reviewed entry in SD-01 through SD-07; no unsupported claims about capability, efficiency or proof proximity.

### SD-01 — Mathematical Research Program and Current Results
Reconstruct the signed explicit-formula foundation, defect morphology, support filtration, native carrier reconciliation, complementary estimates, pole restoration, exceptional-sector control and the proposed localized globalization criterion. Include definitions, explicit quantifiers, hypotheses, theorem/proof statements, counterexamples to insufficient methods, and a dependency graph. Gate: reconcile normalization, signs, domains, conditional premises, formal scope and remaining obligations against pinned sources. Keep finite or fixed-collar positivity distinct from arbitrary support. A derived RH-equivalent receiver is not a proof of its required estimate.

### SD-02 — Adaptive Investigative Governance
Define swarmflow as distributed investigative activity and distinguish it from coherent, adaptively governed mathematical progression. Develop jurisdiction, delegated continuation, evidence-driven regrouping, method stops, conflict precedence, failure recovery and preservation of epistemic state. Use the W02 governance interruption as an explicit case study, with contradictory instructions and recovery evidence where available. Gate: distinguish observed operating practice from proposed doctrine; do not assert causal superiority without a defensible comparison.

### SD-03 — Multi-Agent Research Architecture
Describe implemented thread organization, inter-lineage handoffs, repository topology, claim acceptance, execution environments, context recovery and the human decision boundary. Separate actual implementation from proposed API or infrastructure architecture. Gate: reproduce the control/data flow and identify its failure boundaries without exposing credentials, private conversations or sensitive operational details.

### SD-04 — Formal Verification and Research Provenance
Explain analytical proof, exact arithmetic, directed enclosure, Lean kernel acceptance, axiom audits, independent replay and immutable custody. Include representative successful and failed checks with actual scope and dependency pins. Gate: every formal claim identifies the declaration, source revision, toolchain, accepted axioms and audit evidence; every computational claim states what was checked and what was not.

### SD-05 — Longitudinal Operational Findings
Document available evidence on context exhaustion, recovery loops, custody failures, interrupted work, usage warnings and governance repair. Use measured durations or counts only with documented sampling and denominators. Distinguish observed UI warnings from unverified platform-side causes; do not claim to know why a safety or abuse-detection system intervened. Gate: disclose missing data, confounders, selection bias and limits of generalization; redact personal and security-sensitive material.

### SD-06 — Mathematical Originality and Prior-Art Assessment
Compare exact theorem statements, hypothesis strength, conclusions and proof ingredients against dated primary sources. Treat the CROSS-PW composite criterion as a candidate contribution, not established novelty; separately assess fixed detector construction, moment/pole constraints, translated signed energy and localized Landau reception. Include closest precedents, negative search results with search coverage, and open attribution questions. Gate: verify authorship, dates, bibliographic identifiers and theorem correspondence directly; no priority assertion from absence of a search hit.

### SD-07 — Reproducibility and Evidence Compendium
Supply source pins, declaration identifiers, dependency manifests, environment versions, runnable commands, expected outputs and exact result hashes. Record failed and unavailable replays. Gate: an independent reader can retrieve and rerun each claimed reproducible artifact within its declared scope; otherwise the limitation is explicit. An ephemeral CI artifact alone is not durable custody.

## 3. Publication and representation contract

Every released scientific manuscript has three synchronized representations:
- Complete repository-readable Markdown, with meaningful headings, mathematical equations, tables and source links; not merely an abstract.
- Editable LaTeX source and bibliography, in conventional AMS mathematical style for mathematical papers and a readable technical-report style for governance and systems papers.
- Publication-quality PDF with embedded fonts, selectable text, numbered equations/theorems where appropriate, references, accessible navigation and clear version/standing on its title page.

Use a single editorial source of truth per manuscript and a reproducible conversion/build recipe. Compare equation numbering, claims, references, captions and limitations across representations. Record and review any unavoidable format-specific divergence. Do not pretend that a manuscript is a published journal article or carries an external artifact-evaluation badge.

Proposed repository structure (files not yet created unless explicitly listed in the delivery record):

```text
README.md
SCIENTIFIC_DISCLOSURE.md
EXECUTIVE_BRIEF.md                 # SD-00 once released
CITATION.cff                       # after author/metadata approval
disclosure/
  DISCLOSURE_INDEX.md
  CLAIMS_REGISTER.md
  EVIDENCE_REGISTER.md
  REVIEW_STATUS.md
  CHANGELOG.md
  SD-00/ ... SD-07/
    MANUSCRIPT.md
    manuscript.tex
    references.bib
    manuscript.pdf
    CLAIMS.md
    EVIDENCE.json
    SHA256SUMS
```

The canonical-return README links this charter prominently now; the final disclosure index and executive brief become the primary reader entry points upon release. A reciprocal link from weil-lab is a separate editorial change subject to that repository's owner and claim boundaries. Do not relocate or duplicate authoritative research payloads merely to make the disclosure easier to browse.

## 4. Claim and evidence control

A claim record must carry: stable claim ID; exact assertion and scope; assumptions and unclosed premises; standing; research owner; immutable repository/commit/path and, where appropriate, blob and SHA-256 identities; declaration or computation identifiers; reproduction instructions; verification result; review status; and the disclosure revision using it.

Maintain separate axes rather than collapsing them into one PASS label:

| Axis | Examples |
|---|---|
| Mathematical standing | Imported result, derived, conditional, open, disproven method |
| Computational evidence | Exploratory numerical diagnostic, exact rational check, directed enclosure |
| Formal standing | Not attempted, staged, failed, kernel accepted with specified axioms |
| Custody | Local only, ephemeral artifact, remote source/result verified, complete replay package |
| Review | Author/internal review, independent replay, external mathematical review, peer review |

A successful workflow is not automatically a successful theorem. A conditional positive margin is not unconditional positivity. GitHub custody is not independent reproduction. A failed sufficient bound is not a counterexample to the original theorem. Preserve these distinctions in the executive brief as well as technical appendices.

## 5. Release process

| Gate | Required evidence | Initial state |
|---|---|---|
| D0 — Charter | Scope, jurisdictions, document plan, claim rules | Charter recorded; owner review pending |
| D1 — Evidence freeze | Immutable per-lineage snapshot and completeness inventory | Open |
| D2 — Deliberation | Document-specific thesis, outline, audience, claim list, exclusions | Open |
| D3 — Drafting | Complete manuscript and references | Open |
| D4 — Adversarial review | Claim-by-claim challenge, prior-art audit, limitation reconciliation | Open |
| D5 — Reproduction | Source/readback audit and executable evidence checks | Open |
| D6 — Typesetting | Markdown/TeX/PDF parity, clean build and visual review | Open |
| D7 — Author approval | Authorship, AI-use disclosure, licensing, privacy and release approval | Open |
| D8 — Release | Immutable version, checksums, links and remote readback | Open |

D2–D6 may iterate; no failed gate may be hidden by later formatting. Review each manuscript individually and reconcile shared claims across the entire package before release. Only the owner may approve authorship, public release metadata and institutional representation. No external communications are authorized by this charter.

Suggested publication sequence: v0.1 initial methodological and mathematical disclosure; v0.2 consolidated verified snapshot; v0.3 post-regrouping integration; v1.0 mature independently reviewable package. Versions describe disclosure maturity, not RH completion. Freeze a date and source commit for each version; never silently update a released PDF.

## 6. Immediate work package

D1 is the next editorial task. Recover the actual admitted and investigative heads of all seven lineages, governing amendments, terminal comparative records, exact claim ledgers and executable evidence. Reconcile contradictory or stale descriptions and mark inaccessible material. Freeze the snapshot with a source manifest before composing technical conclusions.

D2 then produces a separate deliberation record for each of SD-00 through SD-07, with document architecture, source map, intended claims, counterclaims, open questions and reviewer checklist. Draft SD-01 through SD-07 from this evidence; derive the executive brief from reviewed material rather than letting it dictate technical claims.

The package should eventually support institutional outreach, but no recipient, timing, endorsement or collaboration is assumed. Scientific disclosure stands on its own merits.

## 7. Delivery boundary

This charter records the plan and release discipline only. It makes no new mathematical deduction, activates no workflow, changes no research acceptance ledger, selects no legal author identity or license, and claims no completed manuscripts or typeset PDFs. Its repository visibility is not approval of future document content.
