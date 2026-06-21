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

## What we computed (`inspect --root horizon`, live)

| N | Q\*(N) | γ\*(N) = J/Q\* | EP verdict (phase rigidity r) |
|---|--------|----------------|-------------------------------|
| 2 | 1 | 1.000 J | r → 0.001, genuine 2nd-order EP |
| 3 | √2 = 1.41422 | 0.707 J | r → 0.002, genuine 2nd-order EP |
| 4 | 1.87855 | 0.532 J | r → 0.012, genuine 2nd-order EP |
| 5 | 2.37216 | 0.422 J | r → 0.000, genuine 2nd-order EP |

(γ at J = 1; Q\*(N) → 2N/π asymptotically.)

- **The Petermann factor peaks AT the EP.** Our phase rigidity r = |⟨L|R⟩| / (‖L‖‖R‖) is exactly the inverse-root of the Petermann factor: K = 1/r². At Q\*(N) the coalescing {0,2}-coherence mode has r → 0, so K → ∞. The √-scaling Im²/(Q − Q\*) = 1.14 (constant at N=4) certifies a clean 2nd-order EP.
- **The survivor is not the EP.** The co-located band-edge mode 2cos(π/(N+1)) is the γ-protected survivor (r ≈ 1), sharing the gap Re = −2γ only because the Absorption Theorem pins both (both ⟨n_diff⟩ = 1). The divergence is the coalescer's, not the survivor's.
- **γ\*(N) drops with N.** The noise level that places an N-chain exactly on its EP is γ\*(N) = J/Q\*(N), decreasing monotonically (→ πJ/(2N)). In our regime, "operating at the EP" *is* "turning the noise up to γ\*"; the two are not independent.
- **Topology** (prior internal result, `simulations/results/trichotomy_cube/petermann_divergence.csv`): chain → EP (K peaks at Q\*), ring → level crossing (r bounded, no coalescence), star → frozen commutant (r ≈ 1 everywhere). The EP, and hence the Petermann divergence, is a chain phenomenon.

---

## What this says to the debate

- It **confirms the geometric fact the debate assumes**, exactly, in a >2×2 many-body system without the coupled-mode approximation photonics relies on: the eigenvector non-orthogonality peaks at the EP.
- It **adds a regime the debate does not occupy**: a single-axis EP where the noise *is* the EP-locator. The "operate beside the EP to escape the noise" move (Wiersig-Rotter's off-EP optimum) is structurally unavailable here, because *beside the EP is a different noise level*. The cancellation argument (Loughlin-Sudhir), which already depends on the noise being EP-coupled, becomes inseparable rather than merely coupled.
- It does **not adjudicate the SNR/QFI verdict**. That needs an input-output measurement model (a probe field, a detected quadrature, a Fisher observable) the toolkit does not carry. We state the geometry the debate argues over; we do not enter the metrology it argues about.

---

## Honest boundary

- The geometry (Petermann / phase-rigidity vs Q, the EP order, Q\*(N), γ\*(N)) is computed and Tier-1. The "distinct regime that contributes to the debate" reading is interpretive (Tier 2): we are not claiming to resolve or defeat either side.
- This is the coherence-horizon single-excitation EP, which **survived** an artifact-free defectiveness verification (r → 0 with √-scaling, departure-from-normality finite). It is **not** the F86a c=2 real-axis EP, which was retracted on 2026-06-21 as a grid artifact. Any reading must cite the SE-EP, never the retracted one.
- The Petermann *magnitude* (K = 1/r², with r ~ 10⁻³ giving K ~ 10⁶ at the grid resolution) is grid-fragile and must not be reported as a hard number. The *structure* (r → 0, clean 2nd-order √-scaling) is robust.

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
- `compute/RCPsiSquared.Core/Numerics/PhaseRigidity.cs` (K = 1/r², the Petermann factor)
- `simulations/results/trichotomy_cube/petermann_divergence.csv` (the chain/ring/star trichotomy)
