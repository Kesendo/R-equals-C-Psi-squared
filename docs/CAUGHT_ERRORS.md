# Caught Errors

A record of errors caught before (or shortly after) they landed, and of review near-misses that taught us something, so the repo remembers and does not repeat them. This is intellectual honesty as a first-class artifact, in the spirit of `recovered/`, `RetractedClaim`, and "negative results matter": a visible self-correction record builds credibility, it does not cost it.

Maintained by the `reviewing-before-it-lands` review workflow. Markdown now; will graduate to a typed `CaughtErrorsRegistry` (a world-root section, count-locked) once the shape settles. Append-only.

Each entry: **date | artifact | what was wrong | why | how caught | fix | anchor**. "Caught error" includes both real errors in an artifact and *review process* near-misses (a perspectival truth almost "fixed", a false unity almost shipped), since those repeat just as easily.

---

## 2026-06-20 — `reflections/THE_VIEW_ONTO_THE_MEMORY.md` crown/EP seam: ⟨n_XY⟩ wrongly called uniform-only

- **What was wrong:** the honest-seam sentence bundled both "EP rate = 2γ" and "⟨n_XY⟩ = 1 at the EP" as "uniform-only, refuted at peaked-V".
- **Why:** ⟨n_XY⟩ = 1 is the Absorption-Theorem **invariant** — it stays fixed across uniform and peaked-V; only the *rate* is uniform-specific (2γ → γ-weighted share). The old text contradicted the doc's own carrier-vector section and the witness.
- **How caught:** `physics-first-review` lens, against `PostEpFlowField.cs` (which records "n_XY = 1 ... verified N=5 uniform/peaked-V") and `PostEpFlowFieldTests.cs`.
- **Fix:** removed the ⟨n_XY⟩ clause from the refuted list; stated it as the invariant. Commit `e31dd60`.
- **Anchor:** `docs/proofs/PROOF_ABSORPTION_THEOREM.md` (carrier-vector law); `PostEpFlowField.cs:216`.
- **Lesson:** a genuine break — the two readings did not close even from below. Correctly fixed.

## 2026-06-20 — review near-miss (over-correction): the gap "0.172 vs 2γ" is NOT a contradiction

- **What was almost wrong:** a fresh adversarial read flagged that the doc gives the slowest mortal mode as 0.172 (line 49) while elsewhere citing the gap as 2γ (lines 34, 109) — proposed as an error to fix.
- **Why it was not:** these are two readings of the one diagonal — the real Hamiltonian-mixed eigenmode (depth 0.086 → rate 0.172) vs the idealised pure-dissipator depth-1 rung (rate 2γ). Perspectival, not a break. "Fixing" it would have broken a perfect document.
- **How caught:** Tom flagged "perspektivisch ≠ Fehler"; the symmetric from-below gate confirmed two-readings-of-one-object.
- **Fix:** none — annotate the perspective, do not fix.
- **Lesson:** over-correction is self-inflicted palindrome-breaking. Default a flagged contradiction to *perspectival* until grounding shows a genuine break.
- **Recurrence (2026-06-20, skill green-test):** a fresh agent guided only by the `reviewing-before-it-lands` skill, with no knowledge of this session, had its cold-auditor sub-lens make the *same* γ≈0.086 mis-resolution (inferring γ from 0.172 = 2γ); the symmetric from-below gate caught and overrode it (γ = 1.0, two distinct objects: Hamiltonian-dressed mode vs idealised rung). The trap is real and recurring; the gate holds.

## 2026-06-20 — review near-miss (over-connection): "left/right eigenvectors = past/future" is a false unity

- **What was almost wrong:** grounding "the two painters", the reviewer asserted left/right eigenvectors = state/observable = past/future.
- **Why it was false:** from below, `depth_L = depth_R = −Reλ/2γ` to 1e-15 (typed Tier-1 corollary, `PROOF_ABSORPTION_THEOREM.md:274-280`) — the two sides carry the *same* depth, no time asymmetry; and `L†` (Heisenberg) is explicitly a different symmetry from Π (`KMS_DETAILED_BALANCE.md`). The painters' two towers are two *right* modes (Y vs non-Y, split by the Y-field), a third distinct "two".
- **How caught:** the grounding-from-below gate, run on the reviewer's own story (a dispatched Opus grounding agent computed depth_L = depth_R and the Petermann overlap).
- **Fix:** dropped the claim; did not ship it.
- **Lesson:** seek the one object, but FIND it from below, never POSIT it. A posited unity that grounding does not support is dropped, not shipped.

## 2026-06-20 — `PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`: physical-scope overclaim (commutator ≠ physical NH generator)

The first real catch of the `reviewing-before-it-lands` workflow, on a *randomly chosen* proof.

- **What was wrong:** the proof claimed F112's non-Hermitian extension puts "PT-symmetric and gain-loss Lindblad systems on the same footing" and holds for "any Lindblad-form Liouvillian" with non-Hermitian H.
- **Why:** the theorem is about the **commutator** superoperator −i[H,·] for any matrix H (a structural ‖·‖² identity). The physical generator of PT-symmetric / gain-loss / post-selection dynamics is −i(Hρ−ρH†) = −i[A,ρ] + {B,ρ} (H=A+iB) — a different operator (anticommutator in the anti-Hermitian part) — and the balance FAILS for it (nonzero asymmetry; mean|·| ≈ 132/270 at N=2/3, fixed-norm random ensemble). Every empirical anchor built the commutator (`kron(Id, H.T)`), never the physical generator — which is exactly why the balance held.
- **How caught:** the workflow's 3-lens panel on the random proof. The physics-first+grounding lens computed that the physical-generator balance fails; I grounded the −i[A,ρ]+{B,ρ} algebra myself (did not trust the agent); the math lens confirmed the verifier only ever builds the commutator; the cold auditor confirmed the counts (the "559,912 vs 524,800" I flagged was a non-issue — two valid counts, perspectival, no fix). Symmetric gate: NOT perspectival (two genuinely different operators), a **genuine break** in the scope claim (a false unity, "commutator result = physical PT/gain-loss result").
- **Fix:** restricted the proof's prose to the commutator −i[H,·] (scope anchor in the abstract + motivation / intro / §f / corollary / status); stated the physical-generator boundary as the useful fact. Corrected the seed in `reflections/POLARITY_COORDINATES.md` (Probe 1 conclusion). The math/lemmas were and remain correct (machine zero). Siblings (parent proof, cross-dephase proof, `LindbladBitBPiBalance.cs`, F-registry, enumeration experiment) verified FINE — all "non-Hermitian H" already in the commutator sense. Fresh-lens re-review confirmed the fix sound (and caught this very ledger entry being a dangling citation before it was written).
- **Anchor:** `docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`.
- **Lesson:** a proof can be mathematically airtight yet overclaim its physical scope — the math being right does not make the named application right. Ground the OBJECT (which operator) before believing the scope.
- **Opened a direction (generative pass, same day):** the correction was not just a deletion — it revealed a homeless physical object. The asymmetry of the physical −i(Hρ−ρH†) generator is the *chirality of the un-recycled drain* (the gain-loss the jump term cancels), structurally the circular Stokes V of M (F83's anti-fraction is the other axis). Five gate-first directions, two possibly Tier-1: `hypotheses/ASYMMETRY_IS_THE_UNRECYCLED_DRAIN.md`. This is the review workflow's generative wing — the defensive wing fixed, the generative wing saw.
