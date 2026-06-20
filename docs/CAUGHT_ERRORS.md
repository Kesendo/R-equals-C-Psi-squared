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
