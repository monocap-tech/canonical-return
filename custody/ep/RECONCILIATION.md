# EP custody reconciliation — 2026-11-26 Pacific

**Status: recovered source snapshot; mathematical progression suspended.** This commit carries the exact 8,020,408-byte local EP archive in 31 ordered binary parts and a complete member digest inventory. It does not activate a workflow, endorse the numerical proposals, or change the conditional index bound **≤205**. Global positivity and RH are unresolved.

## Remote baseline

The staged EP branch `research/cr-int-t1d-ep-targeted-corrector-staged` was at EO evidence commit `600107472ade3316addb19eb2b7fd9d7bcc84799` (tree `2d78dbb5e0370fb00b5c4e900bbaa8c499a1f432`), parent `307e753bc23e81ccb2138eff97f2f7804498ee64`. EO workflow run [45907071606](https://github.com/monocap-tech/canonical-return/actions/runs/45907071606) succeeded on that parent; its evidence commit contains `custody/eo/screen.py`, `generated/screen.json`, `candidate32.npz`, `candidate32_positive_t.npz`, and `.github/workflows/t1d-eo-screen.yml`. This is a numerical screen with omitted exact residual/enclosure terms. EN's evidence attached note at the baseline explicitly retains ≤205 and says EN has not executed. No EP source, EP workflow, or EP branch run was found at the baseline.

## Dependency graph and investigative increments

- Hosted CE → CK → CQ → CR → DA → DG → DK → DM → DV → DZ → EG → EN → EO are the parent workflow artifact lineage. The thirteen parent artifact identifiers, sizes and SHA-256 values are in `PARENT-ARTIFACTS.json` and `manifest.json`. Their 810 uncompressed entries are matched to Git blobs or restored under `parent-ce-recovered`; see the verification receipts.
- EO stored candidates → EP authenticated candidate, finite physical action and archived inverse diagnostics → exceptional dictionary reconstruction and sampled multiplier diagnostics → recovered selection, receiver transport and enlarged receiver → exact conditional receiver norm and kernel audit → kernel energy and receiving → coupled kernel and joint reoptimization → row bank, scale and integer joint replay (the last bank is preserved separately in the Weil Library `joint_replay_bundle.zip`; this EP ZIP preserves the EP scripts/results).
- The local sequence includes failed or nonpromoting trials: 32-direction EO proxy screen, numerical least-squares collective action, artificial lambda pair tradeoff, exceptional dictionary identifiability and sampled multiplier replay, receiver-only and coupled kernel trials, and joint kernel reoptimization. Their saved diagnostics cannot be promoted to a full directed action/tail/physical inverse certificate. In particular the two joint comparison minima replay as negative (approximately -3.24897e-7 and -5.97256e-7). The integer joint bound is conditional on its bank and inherited analytic input, not Lean verification or new positivity.
- Every one of the 64 locally preserved items (source scripts, numerical arrays, reports, logs, checksums and parent references) is listed below. The ZIP passes CRC, and all 61 entries named in `SHA256SUMS` match; README and parent reference JSON are listed in the manifest. The archive SHA-256 is `7b5c0f48b0b63e719be161e4197c603be40d74aee051538fa086791624bada21`.

## Restore and verify

From repository root at this immutable commit, run `cat custody/ep/archive/EP-local-custody.zip.part??? > EP-local-custody.zip`, compare `sha256sum EP-local-custody.zip` with `7b5c0f48b0b63e719be161e4197c603be40d74aee051538fa086791624bada21`, then `unzip -t EP-local-custody.zip`. After extraction compare every member against `manifest.json` SHA-256. Never run the mathematical scripts as a substitute for this custody check.

## Open custody and verification

All thirteen parent workflow ZIPs were downloaded from Actions and matched their recorded SHA-256 and ZIP CRC. A member-by-member Git blob comparison matched 800 of 810 uncompressed entries to the EO baseline tree; ten missing CE files were added to this custody branch. See `PARENT-VERIFY.json` and `PARENT-TREE-COMPARISON.json`. The original ZIP container bytes remain in time-limited Actions storage; all their 810 member payloads now have repository blob custody. The EP archive has no EP runner execution, independent directed physical inverse or tail check, exact exceptional dictionary, complete mathematical acceptance or Lean proof. The provenance of the previously reported rejected EP write is not preserved in the archive; its exact rejection response is unverified. A direct `git ls-remote` to github.com from this environment returns HTTP 403, while the authenticated repository write interface accepts branch and blob creation; this recovery uses that interface. The baseline EP branch itself is intentionally left at EO; this custody branch stages the recoverable EP evidence.

## Exhaustive archive inventory

| Member | Bytes | SHA-256 |
| --- | ---: | --- |
| EP-EO-candidate-audit.json | 1397 | `93c49e04da39b1b27bf79884b71ca08c39ece963d1f5f2e37b272712b0ae7399` |
| EP-action-increment-audit.json | 51358 | `de11e3c592fc8af686a55d57c780473b9ff74c00c83be84352d7b955a05bd559` |
| EP-archived-inverse-replay.json | 2925 | `e27df500bb67bb2fab2297bd7690fb45695b9ab1159edb2e1efcfc1a80331584` |
| EP-collective-action-audit.json | 8402 | `ca4725efb312af787127e9fc18669142f1cf517520208f4e3defb3271d6f1d25` |
| EP-energy-increment-audit.json | 44917 | `bbdd3238d010fd69ab7e57791974747b72a3019a8c1c39a5479cf6e1e7f09faf` |
| EP-enlarged-receiver.json | 3317 | `5838cf91e883f5fc969cc81abe7d60782a7bb5bb011974bfaf272096597b3895` |
| EP-enlarged-receiver.npz | 2434432 | `99b901da0815ea1fe97f421af5f150d47d7287c3e25aa1343c4ee153b671dc29` |
| EP-exceptional-expanded-identifiability.json | 24602 | `e1cb643192ef6331ac8129201bbcbb813cff02a08fd09fdbb9931f5faef6f425` |
| EP-exceptional-gram-audit.json | 3387 | `a769dd217bd0cb0c2821b24cdf4cbebc5c368d9d368387d062f21933354cba45` |
| EP-exceptional-identifiability.json | 20018 | `80258475fac3195fea3e0d5e1a617bc76cbce7d39ff3982cc67a2d249969a96d` |
| EP-expanded-exceptional-replay.json | 7958 | `077a3be0f0c4503a151fa82fe35e5468f3f332fedf237512459598d2a8546acc` |
| EP-extended-metric-attachment.json | 6970 | `b84418469fec2680390c0b6feb27a16f8c598c372440c836c10e726ec8416803` |
| EP-kernel-coupled-block-0.npz | 50517 | `feef5afa171934780ec555554091369df31204e0bb6ceac738c4d39e5586ebc4` |
| EP-kernel-coupled-block-1.npz | 50495 | `9f55404e1175acb7ef07ae07c51b1c365cf7e55bc51c62e097aae9670bcddcf0` |
| EP-kernel-coupled.json | 5687 | `7a827dbf87f7a9247211953b3143784c4aae32143fba3334bfdb61e6e6339b6e` |
| EP-kernel-energy.json | 961 | `3897a7e574e3f207317164014d91785d8e2f09745606bd95f6ba7163351c34a2` |
| EP-kernel-energy.npz | 62038 | `bfc3438d0aed6b57e069aa9d5ba05cda6b0c60a98d1366c3ff67e30530a4b196` |
| EP-kernel-joint-reopt-block-0.npz | 212557 | `fd1bed41d050f2e0d795ba83e75ba07f3417ad51cca0faba3e79b13f5df64376` |
| EP-kernel-joint-reopt-block-1.npz | 212577 | `55ee47bec34fe66059a5068d78877d61aa6f7c63f2919f453503f4d124e0fcea` |
| EP-kernel-joint-reopt.json | 7438 | `c437ef15f3129226533df1e400a157357fa2051b0d0a77b87d4a0b0d797439e4` |
| EP-kernel-receiving-block-0.npz | 2044 | `9968bd0b5142fde6aa83df7c7b2a8af47cb3c3b4a1e799810121b262cf016e68` |
| EP-kernel-receiving-block-1.npz | 2045 | `e5b943b99f50ace998c5049128cb69b91623ece58645bbd310d415f847b494ef` |
| EP-kernel-receiving.json | 4057 | `8500adc9cc57dee672f4909ca747aaaefeee97321591164d008ec9a225581cfe` |
| EP-pair-tradeoff-audit.json | 735432 | `d3432deeb481bc02a1a8a8d9b33d25586ce9f9b2fa907a486bb4e727a482f041` |
| EP-physical-action-audit.json | 2559 | `6926ab8a82ec20872e06410ccfda37b6789e4fa2a442293147032b3a3179cd26` |
| EP-receiver-kernel-audit.json | 3309 | `8cc7e27f9862ce6b6eb0f5911332e57a62e25feaddd8e5931c734e3a76ad637b` |
| EP-receiver-kernel-proposal.npz | 4178445 | `eadd3e6d04796b741608ab9c09fd5c9f7b41254fce55047750fda144ee0764ea` |
| EP-receiver-norm-audit.json | 4393 | `e1f0afec5cccc09724ef2e4677599adfca1a6bb3d74a834314055de86ae43e08` |
| EP-receiver-norm-congruences.npz | 279718 | `0e221fb831de3bc2c11976bc1c16231b2e89c61a76dafbc2d0caf7f627bbdde4` |
| EP-receiver-transport.json | 11414 | `7a1e05c50af540d471959b9b50b532cd7159aed2d85aca707366df75cfd015dd` |
| EP-recovered-selection.json | 11590 | `9b47c57729506ee8ade437316fbc2cf889628794174ce8d0c7cffaf1083b6398` |
| EP-recovered-selection.npz | 363698 | `755dc57a4e406bbe132028d6b03e9d54f087e0247e33ea9b7eb27581509ff0ea` |
| EP-sampled-multiplier-replay.json | 3822 | `847efe33305a58cdcb2daf9d83648dd7f1cb1daad406bab53e9341c3ec81c7b5` |
| EP-trial-subspace-audit.json | 9817 | `54ad26fb922fdb4b7b88d5138ae242c585520cafb95aaf279d8bbe88d0a7003a` |
| ep_action_increment_audit.py | 3697 | `7bee3dd3a6484cbed498ac614d917dd096a960564e2dfed3ddd30ed66a641884` |
| ep_archived_inverse_replay.py | 4108 | `013bb53ae9af2a964da1d766a31d0047d256b78b3a8ccd0a21bd9dfecdda45df` |
| ep_collective_action_audit.log | 8402 | `ca4725efb312af787127e9fc18669142f1cf517520208f4e3defb3271d6f1d25` |
| ep_collective_action_audit.py | 4120 | `019deb60c9cf8df490a9c63ef872bba4d16573e7ab16a89f30c95f6611142fe6` |
| ep_energy_increment_audit.py | 3412 | `042f20bcd60ddbc51bae4ed9521d3b766d0795554cb31ed6b3cb719e7a558d18` |
| ep_enlarged_receiver.py | 4051 | `eb6553e38111e0e1c67d91f0c35abca37f55b41adff97165ba7762ecab10cf5f` |
| ep_exceptional_expanded_identifiability.py | 4176 | `9a55433340a4b7adf4c076957c38c6ff7ebd336ab2d3f0d7f71634b10bcf33ff` |
| ep_exceptional_gram_audit.py | 1795 | `6ea83e82ca79f44e974c312eb4afa66d72dddbce833646c5b1a42721de44b5df` |
| ep_exceptional_identifiability.py | 4071 | `fe56f6bfcd19e4f7a908e40bc926822fa0e22d7aca1817a0ab1d5c569f0af2c8` |
| ep_expanded_exceptional_replay.log | 8190 | `5f586f6d7a4af9810247ed1d336160a23451d7d9ae3c755fe930d12d7ee8b6cf` |
| ep_expanded_exceptional_replay.py | 9022 | `ac5110780b57a10e8379b21a4e86d5781a3474600d11886408f93e4aa7b07bdc` |
| ep_extended_metric_attachment.py | 3913 | `38b328033d6f83b8b4e620b8b26b99741811a2f8c25ebc8474b4fd3a1f635214` |
| ep_kernel_coupled.log | 9468 | `305caf1d42fe0d50e13c1100c2db52f1b90b4ee561a0dba9a80a2b4c6550ea7a` |
| ep_kernel_coupled.py | 3545 | `d5ef7c10e3511fc504735e40d58f25a0fec11cc1ee9911892695cb85600a830c` |
| ep_kernel_joint_reopt.log | 12087 | `32d11cfee94acd6417d10da1fc8efa797d78cb84056f03b806c86039e29e577a` |
| ep_kernel_joint_reopt.py | 2796 | `ad5901d329a575d4dc8609e7f11659f3dce0d2e1d9437ba7d2407b8b6491f1ca` |
| ep_kernel_receiving.py | 7625 | `75bc39950768c91a2a99121398f6d3b5978b20ddc21dc16028cc5df736e9121d` |
| ep_pair_tradeoff_audit.py | 4296 | `8ee584068b4c5ca8273536de17a9b1285e7c7d0fd8304f4dabbbabd34ac39f40` |
| ep_physical_action_audit.py | 4452 | `51671a8d812336bd0fd717a83e28973b27243b592a53205b9b362f8fd40add67` |
| ep_physical_preflight.py | 3742 | `cc11327c29db314dcbbde14c7d1e4b066401cbf875f02efecaec0f682a95beee` |
| ep_receiver_kernel_audit.py | 6831 | `087ba74257f9265cf7159600f71786811a54d5e04493d29092d56520b73b0d6d` |
| ep_receiver_norm_audit.py | 4461 | `36a04b89ee1b79a935b41952f2247ecfd86881d1e644d95f52209f6f05a60e36` |
| ep_receiver_transport.py | 4145 | `b467bc1a1661a018089059cf1de2a53a8c653367cdf68a4df677f62c2e209645` |
| ep_recovered_selection.py | 3067 | `8b5cc89ac655c828d37eb56d90b1656e7492aa0348199bad0b8b44fe44f54efb` |
| ep_sampled_multiplier_replay.py | 8327 | `0870110e19b0e1308f8265599ab2314fc761af8df8a655655843cbf867429ec1` |
| ep_trial_subspace_audit.log | 9817 | `54ad26fb922fdb4b7b88d5138ae242c585520cafb95aaf279d8bbe88d0a7003a` |
| ep_trial_subspace_audit.py | 2635 | `d9714d8c61986c8db6049cadc8d28b3b3f237364bd036d416f4743ddc8e9c780` |
| README.md | 1198 | `2fff05237702d3eca2a842612bd28ddc6557493caae40db3e394c491729db22f` |
| SHA256SUMS | 5818 | `a9d3a7595ee39ee9d3a5a5163b9dc9170d2532fed0e20ec94729d4194f95ff77` |
| PARENT-ARTIFACTS.json | 2498 | `fddd04d5d9f663389a9923a00773ae2872b277159aadf58ffa867ffbf106f873` |
