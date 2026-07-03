# The Coherence-Horizon Exceptional Point Meets the EP-Sensor Debate

**Tier:** computed geometry Tier 1; the contribution-to-the-debate reading Tier 2 (interpretive).
**Date:** June 21, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Origin:** A literature harvest asked which open questions our exact small-N Liouvillian tools could speak to. The exceptional-point sensing debate was the cleanest match: its load-bearing quantity is the one our coherence-horizon witness already computes exactly. So we took their question to our tool.

---

## What this experiment is about

A live, openly-contested question in non-Hermitian sensing: do exceptional-point (EP) sensors give a fundamental advantage, and does the optimum sit *at* the EP or *beside* it? The quantity that carries every version of the argument is the **Petermann factor** (the left/right eigenvector non-orthogonality), which diverges as two modes coalesce at an EP.

We took that question to our own **coherence-horizon EP** Q\*(N): a defective degeneracy that, unlike the gain/loss PT sensors of the debate, lives on a **single axis** Q = J/γ, where the dephasing rate γ is at once the *entire noise* and the *locator of the EP*. We supply, exactly and in a many-body (>2×2) system without the coupled-mode approximation, the geometry the debate treats as its settled foundation, in a regime the debate's models do not occupy. We do **not** supply the metrological verdict (SNR / QFI): that needs an input-output measurement model the toolkit does not carry. Question theirs; geometry ours; the boundary stays honest.

---

## The external question (citations verified at source 2026-06-21)

- **Loughlin & Sudhir**, "Exceptional-point Sensors Offer No Fundamental Signal-to-Noise Ratio Enhancement," *Phys. Rev. Lett.* **132**, 243601 (2024), [arXiv:2401.04825](https://arxiv.org/abs/2401.04825). The EP signal gain (splitting ∝ ε^(1/n)) is *exactly cancelled* by a noise gain that they identify as a Petermann factor; the force imprecision becomes independent of proximity to the EP.
- **Wiersig & Rotter**, "Fundamental Limits of Non-Hermitian Sensing from Quantum Fisher Information," [arXiv:2603.10614](https://arxiv.org/abs/2603.10614) (2026). EPs *can* beat ordinary modes of equal decay rate, and "the QFI can be further increased by moving *away* from the EP."
- **Kullig, Wiersig & Schomerus**, "Generalized Petermann factor of non-Hermitian systems at exceptional points," *Phys. Rev. Research* **7**, 043246 (2025), [arXiv:2506.15807](https://arxiv.org/abs/2506.15807). The Petermann factor K = ⟨L|L⟩⟨R|R⟩ / |⟨L|R⟩|² as the eigenvector-overlap quantity that diverges at the EP.

The two metrological papers reach *opposite verdicts* and are not in explicit dialogue, because they assume different observables and noise budgets. What they share, and what the field treats as settled, is the **geometry**: the Petermann factor vs the parameter through the degeneracy.

---

## Our object: an EP on one axis

The coherence horizon Q\*(N) is the point on the axis Q = J/γ below which the slowest non-zero Liouvillian mode stops oscillating: the single-excitation {0,2}-coherence mode coalesces. Its dispersion in the Haken-Strobl single-excitation block is λ² + 8γλ + 4J²q² = 0, an exact square-root EP where the discriminant vanishes. (Typed: `CoherenceHorizonClaim`; proof: `docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md`.)

The structural point that the debate's gain/loss PT sensors do not have: in their systems the gain/loss rate γ (which the noise couples to) and the distance to the EP ε̄ are **separate knobs**. Here there is **one knob**, Q = J/γ. The same γ is the whole dissipation, the observer, the only source of irreversibility, *and* one of the two numbers that locate the EP at Q\* = J/γ\*. You cannot perturb the noise without walking along the very axis the EP sits on.

---

## What we computed (`inspect --root horizon`, live; defectiveness verified gate-first)

| N | Q\*(N) | γ\*(N) = J/Q\* | EP character |
|---|--------|----------------|--------------|
| 2 | 1 | 1.000 J | genuine 2nd-order defective EP |
| 3 | √2 = 1.41421 | 0.707 J | genuine 2nd-order defective EP |
| 4 | 1.87874 | 0.532 J | genuine 2nd-order defective EP |
| 5 | 2.37367 | 0.421 J | genuine 2nd-order defective EP |

(γ at J = 1; Q\*(N) → 2N/π asymptotically.)

- **Defective, by four independent artifact-free routes** (the load-bearing verification, gate-first, N=2..5): departure-from-normality ≈ 4 as the pair-split → 0 (a *diabolic* degeneracy would send it → 0); geometric multiplicity 1 < algebraic 2; the Schur off-diagonal is the Jordan coupling; the two eigenvectors merge (|cos| → 1). The test is not rigged to cry "defective": a diabolic toy returns DIABOLIC, and the *same* coherence-horizon object at γ = 0 returns DIABOLIC (dep ≈ 2·10⁻¹⁷). The √-scaling Im²/(Q − Q\*) ≈ 1.14 (constant at N=4) fixes the order at 2.
- **Phase rigidity r → 0 corroborates but is not load-bearing.** r = |⟨L|R⟩| / (‖L‖‖R‖), with K = 1/r² the Petermann factor, drops toward 0 at Q\*(N). This is exactly the eig-instrument that misfired at F86a; here it is vindicated (it is right *because* the object is genuinely defective), but it cannot by itself separate defective from near-degenerate. So the dep / geo-vs-alg / eigenvector-merge above is load-bearing; r → 0 is the corroborating shadow.
- **The Petermann factor peaks AT the EP: in kind, not in a quotable number.** K = 1/r² → ∞ at Q\*(N). The magnitude is grid-bound: the grid-sampled peak is K ~ 232 / 501 at N = 4 / 5 (finest Δ ≈ 0.004) and grows without bound under refinement; the rigidity bottoms grid-sampled at ~0.06 (zero is the limit, not a sampled value). The EP and the divergence-in-kind are robust; the peak height is a grid artifact and must never be quoted as a sensor figure of merit without grid refinement (the F86a lesson, confirmed here).
- **The survivor is not the EP.** The co-located band-edge mode 2cos(π/(N+1)) is the γ-protected survivor (r ≈ 1), sharing the gap Re = −2γ only because the Absorption Theorem pins both (both ⟨n_diff⟩ = 1). The divergence is the coalescer's, not the survivor's.
- **γ\*(N) drops with N.** The noise level that places an N-chain exactly on its EP is γ\*(N) = J/Q\*(N), decreasing monotonically (→ πJ/(2N)). In our regime, "operating at the EP" *is* "turning the noise up to γ\*"; the two are not independent.
- **Topology** (prior internal result, `simulations/results/trichotomy_cube/petermann_divergence.csv`): the qualitative trichotomy holds: chain → EP (r → 0 at Q\*), ring → level crossing (r bounded ~0.6, still oscillating at the chain's Q\*), star → frozen commutant (r mostly bounded). The EP, and hence the Petermann divergence, is a chain phenomenon. Only the quantitative K-peak is grid-fragile.

---

## What this says to the debate

- It **confirms the geometric fact the debate assumes**, exactly, in a >2×2 many-body system without the coupled-mode approximation photonics relies on: the eigenvector non-orthogonality peaks at the EP.
- It **adds a regime the debate does not occupy**: a single-axis EP where the noise *is* the EP-locator. The "operate beside the EP to escape the noise" move (Wiersig-Rotter's off-EP optimum) is structurally unavailable here, because *beside the EP is a different noise level*. The cancellation argument (Loughlin-Sudhir), which already depends on the noise being EP-coupled, becomes inseparable rather than merely coupled.
- It does **not adjudicate the SNR/QFI verdict**. That needs an input-output measurement model (a probe field, a detected quadrature, a Fisher observable) the toolkit does not carry. We state the geometry the debate argues over; we do not enter the metrology it argues about.
- **A structural echo, offered honestly (interpretive, Tier 2).** The debate, at bottom, is "the EP effect is real *in kind*, its practical magnitude is fragile and contested." Our result rhymes with that exactly: the EP and its Petermann divergence are robustly real (four artifact-free routes), while the magnitude K is fragile, to the computational grid rather than to physical noise. Same shape, different cause: the effect is solid, the number is the soft part. A rhyme, not an identity.

---

## Honest boundary

- The geometry (Petermann / phase-rigidity vs Q, the EP order, Q\*(N), γ\*(N)) is computed and Tier-1. The "distinct regime that contributes to the debate" reading is interpretive (Tier 2): we are not claiming to resolve or defeat either side.
- This is the coherence-horizon single-excitation EP, whose defectiveness was verified gate-first by four independent artifact-free routes (departure-from-normality ≈ 4, geo 1 < alg 2, Schur-Jordan coupling, eigenvector merge), with diabolic controls (the same object at γ = 0 returns DIABOLIC). It is **not** the F86a c=2 real-axis EP, which was retracted on 2026-06-21 as a grid artifact. Any reading must cite the SE-EP, never the retracted one. (Verification: the parallel EP-character work, 2026-06-21; to be wired to the live C# witness `EpCharacter` once it lands.)
- The Petermann *magnitude* (K = 1/r²) is grid-fragile: the grid-sampled peak K ~ 232 / 501 at N = 4 / 5 grows without bound under refinement, and the rigidity bottoms grid-sampled at ~0.06 (zero is the limit). It must never be reported as a hard number or a sensor figure of merit. The *structure*, a genuine 2nd-order defective EP with a real Petermann divergence, is robust.

---

## Open sub-question (separate internal thread)

The witness flags that the closed form of Q\*(N) (the {0,2}-block discriminant condition) is open: the band-edge form 2cos(π/(N+1)) matches only at N = 2, 3 and departs at N ≥ 4 (Q\*(4) = 1.879 ≠ φ). This may be the `NivenRationalityRoot` number-theoretic ceiling (first golden at N=4) showing up as an absence of elementary closed form. Pursued separately; it has no external question-asker and so does not belong in this experiment.

---

## References

**External (the question):**
- Loughlin & Sudhir, PRL 132, 243601 (2024), [arXiv:2401.04825](https://arxiv.org/abs/2401.04825)
- Wiersig & Rotter, [arXiv:2603.10614](https://arxiv.org/abs/2603.10614) (2026)
- Kullig, Wiersig & Schomerus, PRR 7, 043246 (2025), [arXiv:2506.15807](https://arxiv.org/abs/2506.15807)

**Internal (the answer):**
- `compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs` (Q\*(N), typed)
- `compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs` (`inspect --root horizon`)
- `docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md` (the dispersion, the 2/π slope)
- `compute/RCPsiSquared.Core/Numerics/PhaseRigidity.cs` (K = 1/r², the Petermann factor; corroborating, not load-bearing)
- The gate-first defective-vs-diabolic verification (parallel EP-character work, 2026-06-21): four artifact-free routes (departure-from-normality, geo-vs-alg multiplicity, Schur-Jordan coupling, eigenvector merge) + diabolic controls (a diabolic toy and the same object at γ=0 both return DIABOLIC). The load-bearing anchor for the EP character above. C# home `compute/RCPsiSquared.Core/Numerics/EpCharacter.cs` (to be wired on commit).
- `simulations/results/trichotomy_cube/petermann_divergence.csv` (the chain/ring/star trichotomy; qualitative, K-peak grid-fragile)
