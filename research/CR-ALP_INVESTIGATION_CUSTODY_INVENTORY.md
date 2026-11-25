# CR-ALP investigation-wide custody inventory — nonpromoting

Snapshot of draft PR #3 at commit 6de57fd76045bd15c3824b6ac3afeaf5d185410f. 346 changed paths. This is a path inventory and recovery record, not independent mathematical verification. All stages 001–058 have at least one research note; 059 has its threshold audit script but no measured second moment. No promotion to Weil-Lab, Phase IV or RH.

## Stage inventory

| Stage | Research paths | Other named paths | Custody status |
|---|---:|---:|---|
| 001 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 002 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 003 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 004 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 005 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 006 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 007 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 008 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 009 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 010 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 011 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 012 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 013 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 014 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 015 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 016 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 017 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 018 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 019 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 020 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 021 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 022 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 023 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 024 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 025 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 026 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 027 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 028 | 1 | 4 | Named Git paths present; independent complete replay not checked |
| 029 | 1 | 2 | Named Git paths present; independent complete replay not checked |
| 030 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 031 | 1 | 10 | Named Git paths present; independent complete replay not checked |
| 032 | 1 | 2 | Named Git paths present; independent complete replay not checked |
| 033 | 1 | 11 | Named Git paths present; independent complete replay not checked |
| 034 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 035 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 036 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 037 | 1 | 1 | Named Git paths present; independent complete replay not checked |
| 038 | 2 | 7 | Named Git paths present; independent complete replay not checked |
| 039 | 1 | 10 | Named Git paths present; independent complete replay not checked |
| 040 | 1 | 12 | Named Git paths present; independent complete replay not checked |
| 041 | 1 | 3 | Named Git paths present; independent complete replay not checked |
| 042 | 1 | 3 | Named Git paths present; independent complete replay not checked |
| 043 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 044 | 1 | 8 | Named Git paths present; independent complete replay not checked |
| 045 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 046 | 1 | 9 | Named Git paths present; independent complete replay not checked |
| 047 | 1 | 4 | Named Git paths present; independent complete replay not checked |
| 048 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 049 | 1 | 6 | Named Git paths present; independent complete replay not checked |
| 050 | 1 | 5 | Named Git paths present; independent complete replay not checked |
| 051 | 1 | 7 | Named Git paths present; independent complete replay not checked |
| 052 | 1 | 6 | Named Git paths present; independent complete replay not checked |
| 053 | 1 | 5 | Named Git paths present; independent complete replay not checked |
| 054 | 1 | 9 | Named Git paths present; independent complete replay not checked |
| 055 | 1 | 0 | Cumulative 058 archive; 93 entries confirmed in central directory |
| 056 | 1 | 0 | Cumulative 058 archive; 93 entries confirmed in central directory |
| 057 | 1 | 0 | Cumulative 058 archive; 93 entries confirmed in central directory |
| 058 | 1 | 98 | Cumulative 058 archive; 93 entries confirmed in central directory |
| 059 | 0 | 1 | Threshold audit only; N not evaluated |

## Recovered cumulative archive

See `src/audits/cr_alp_058_archive/README.md` and `reassemble.py`. The 90 part paths are present. Independently fetched part-089.bin from pinned commit, decoded its 572042 bytes and parsed the central directory: all 93 entries present in the directory, including 055–058 inputs, selected predecessor banks from 045,046,049,053,054, and 058 input_1.txt, input_2.txt, mesh.txt, gauss.txt, prepare.py, validated.cpp, logs and SHA256SUMS. This is independent directory inspection, not a fresh full ZIP checksum or execution of every entry. Prior recovery README reports SHA256 f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2 and successful saved replays. The earlier partial recovery README is historical and superseded for cumulative 055–058 archive custody.

## Standalone archive catalog (first page; not exhaustive for earlier stages)

| Library filename | Bytes |
|---|---:|
| CR-ALP-058_global_physical_witness_moments.zip | 62872042 |
| CR-ALP-057_global_physical_residual_proxy.zip | 62618105 |
| CR-ALP-056_corrected_Q850_residual_certificate.zip | 48340364 |
| CR-ALP-055_complete_E850_bridge_certificate.zip | 46666228 |
| CR-ALP-054_two_moment_structured_848_certificate.zip | 25107539 |
| CR-ALP-053_global_bridge_augmented_certificate.zip | 3555110 |
| CR-ALP-052_outside_band_bridge_falsifier.zip | 13866553 |
| CR-ALP-051_polynomial_moments_598_mode_certificate.zip | 13625276 |
| CR-ALP-050_second_structured_correction_certificate.zip | 3445742 |
| CR-ALP-049_trial_inverse_residual_falsifier.zip | 2611510 |
| CR-ALP-048_structured_prime_tail_certificate.zip | 1570886 |
| CR-ALP-047_scalar_falsifier_and_103_mode_certificate.zip | 1573537 |
| CR-ALP-046_complete_E85_finite_certificate.zip | 1567762 |
| CR-ALP-045_layered_gamma_and_mixed_certificate.zip | 64059 |
| CR-ALP-044_fitted_matrix_172_mixed_certificate.zip | 32775 |
| CR-ALP-043_reduced_cutoff_mixed_certificate.zip | 33069 |
| CR-ALP-042_support_aware_172_mode_tail_certificate.zip | 18627 |
| CR-ALP-041_energy_conditioned_14_mixed_certificate.zip | 13479 |
| CR-ALP-040_first_hierarchical_shell_certificate.zip | 41361 |
| CR-ALP-039_complete_twelve_plus_infinite_tail_certificate.zip | 33093 |
| CR-ALP-038_complete_signed_interior_infinite_tail_certificate.zip | 35558 |
| CR-ALP-037_uniform_original_signed_tail_package.zip | 10618 |
| CR-ALP-036_source_signed_prime_gamma_tail_package.zip | 18089 |
| CR-ALP-035_rank_two_tail_projection_package.zip | 7740 |
| CR-ALP-034_twelve_mode_certificate.zip | 26460 |
| CR-ALP-033_partial_eleven_mode_certificate.zip | 26593 |
| CR-ALP-032_full_reproducibility_package.zip | 21374 |
| CR-ALP-031_complete_ten_mode_certificate.zip | 29804 |
| CR-ALP-030_complete_eight_mode_certificate.zip | 29115 |
| CR-ALP-029_certified_six_mode_package.zip | 12539 |
| CR-ALP-028_research_and_audits.zip | 14455 |

Raw-byte materialization of the sampled standalone 029, 052 and 054 archives was denied in this session. Earlier Git sources and selected banks do not establish that all original standalone archives are mirrored byte-for-byte. Outstanding: reconcile each standalone package against Git paths and hashes, transfer any missing immutable evidence, and execute independent fresh-checkout replays. No claims of missing mathematics or lost files based on access denial alone.
