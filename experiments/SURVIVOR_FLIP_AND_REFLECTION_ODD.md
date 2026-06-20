# The Survivor is Spin-Flip-Odd and Reflection-Odd: Darkness is Antisymmetry

**Status:** Empirically verified, gate-first, at N = 4 and N = 6 for the open chain and the periodic ring, across every Q tested (`simulations/_survivor_particle_hole_mirror.py`, dense per-sector eigendecomposition). The four facts (X^⊗N-odd, R-odd, not staggered-PH-fixed, dark + real) are one fact about the half-filling density standing wave, read through three operators. This re-reads, in the survivor's own frame, the spectrum-wide R-parity block structure of [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md).
**Date:** 2026-06-19
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Depends on:**
- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): `Re λ = −2γ·⟨n_XY⟩`, the rate-equals-light-content law that makes "dark" a measurable number.
- [PROOF_DIFFUSION_RAYLEIGH_CLOSURE](../docs/proofs/PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md) / [`SurvivorDiffusionGradientClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs): the survivor is predominantly a real diagonal density mode `n(j)` (a standing wave), its decay carried by a subdominant Hamming-2 coherence admixture. The density profile is the lowest antisymmetric Neumann harmonic.
- [`SurvivalIncompletenessMirrorClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/SurvivalIncompletenessMirrorClaim.cs): the half-filling survivor (the slowest `(p,p)` mode) and its `⟨n_XY⟩ ~ Q²/N²` darkness, `inspect --root survivor`.
- [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md): the whole-spectrum R-parity decomposition, and (Steps 4, 6) the prior analysis of the Bogoliubov particle-hole operator on the XY chain.
- F1² = Π² = X^⊗N (the global spin-flip is the square of the palindrome conjugation; F1 entry in [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md)); F61 (the distinct n_XY-parity Z₂).
- Verifier: [`simulations/_survivor_particle_hole_mirror.py`](../simulations/_survivor_particle_hole_mirror.py)

---

## What this is about

Watch a Heisenberg spin chain and the watching wears its rhythms away. One coherence outlives all the others: the survivor, the chain's longest-lived memory, the note still sounding after the rest have gone quiet. We had already learned its shape. In the regime where the watching is strong, the survivor is not a fast-beating wave at all but a single slow swell of *density*, smooth and even across the chain, half the spins up and half down, spread so thin that the watching can barely read it. That is why it lives long: a pure population is invisible to Z-dephasing, and this survivor is mostly population.

This note answers a narrower question: what is the *shape* of that swell, exactly, and what does its shape tell us. The answer is that the swell is **antisymmetric**. It rises smoothly from low at one end of the chain to high at the other, passing through zero in the middle, and it does this in the gentlest possible way (one half-wavelength across the whole chain). Three things that looked like three separate facts about the survivor are really one fact about this antisymmetric shape, seen from three sides. Flip every spin (swap up and down everywhere): an antisymmetric, half-filled swell turns into its own negative. Reflect the chain end-for-end: an antisymmetric swell turns into its own negative again. Both operations send the survivor to minus-itself, so it is "odd" under both. And the reason a closely related operation from the condensed-matter toolbox (the staggered particle-hole symmetry, the one a physicist would reach for first) *fails* to fix the survivor is also about this shape, sharpened by which couplings the chain has. The survivor is dark because it is smooth and even; it is odd because that even smoothness is, more precisely, antisymmetric.

## Abstract

For the half-filling **survivor** (the slowest non-stationary Liouvillian mode of the Heisenberg `XX+YY+ZZ` chain or ring under uniform Z-dephasing, in the strong-dephasing regime below the coherence horizon), four robust facts hold for chain and ring at every Q tested:

1. The survivor is an eigenvector of the **global spin-flip** `X^⊗N = ⊗_l X_l` with eigenvalue **−1** (X-odd).
2. It is an eigenvector of **spatial reflection** `R: l → N−1−l` with eigenvalue **−1** (R-odd).
3. The **staggered (bipartite) particle-hole** operator `U = X^⊗N · ∏_{l odd} Z_l` does **not** fix it: the survivor is not an eigenvector of `U`, and `U` does not even commute with the Heisenberg Liouvillian.
4. It is **dark**: `⟨n_XY⟩ ≪ 1`, falling as `~Q²/N²` and `→ 0` as `Q = J/γ → 0`, with a real eigenvalue (`|Im| = 0`, overdamped, non-oscillating).

These are not four coincidences. The survivor's dominant part is a real diagonal density operator `M ≈ diag(n_a)` whose single-site occupation profile is the lowest antisymmetric Neumann standing wave `n(j) ∝ cos(π(j−½)/N)` (numerically `[−0.74, −0.46, +0.46, +0.74]` at N=4, antisymmetric about the chain centre). Facts 1, 2, 4 are three readings of that antisymmetry; fact 3 is the one place where the Heisenberg `ZZ` term matters. Gate-verified: under both `X^⊗N` and `R` the survivor returns eigenvalue exactly −1, and the diagonal density vector itself satisfies `n · n_complement / ‖n‖² = −1` and `n · n_reflected / ‖n‖² = −1` (machine precision, N=4).

## 1. The setup and the measured result

Open chain or periodic ring, `H = (J/2) Σ_{(i,j)∈bonds} (X_iX_j + Y_iY_j + Z_iZ_j)` (Heisenberg), uniform Z-dephasing at rate `γ`, Liouvillian `L = −i[H,·] + γ Σ_l (Z_l · Z_l − ·)`. The Liouvillian is block-diagonal in the joint-excitation-number pair `(p_c, p_r) = (popcount ket, popcount bra)` because both `H` (`U(1)`-conserving) and the dissipator (diagonal) preserve it. The survivor is the slowest strictly-decaying eigenmode over all blocks.

Running the gate-first verifier (`_survivor_particle_hole_mirror.py`) on the N=6 Heisenberg chain:

| Q = J/γ | survivor sector | Re λ | \|Im\| | ⟨n_XY⟩ | X^⊗N | R | U (staggered PH) |
|---|---|---|---|---|---|---|---|
| 1 | (3,3) half-filling | −0.405 | 0.000 | 0.203 | **−1** (fid 1.000) | **−1** (fid 1.000) | not an eigenvector (fid 0.82) |
| 2 | (3,3) half-filling | −0.636 | 0.000 | 0.636 | **−1** (fid 1.000) | **−1** (fid 1.000) | not an eigenvector (fid 0.55) |

The intertwiner gate (`G1`) confirms `X^⊗N` and `R` both **commute** with `L` exactly (so the survivor *can* have a definite parity under each), while `U` is "neither" (mixed: it neither commutes nor anticommutes with the Heisenberg `L`). The same signs hold on the ring. Control runs pin the regime: in pure XY (no `ZZ`) or on the star, the survivor leaves half-filling, so the half-filling premise is specifically a Heisenberg-on-dispersive-matter, low-Q statement; and at Q=20 the survivor has moved to the `(0,1)` band edge (`⟨n_XY⟩ = 1`, oscillating), the regime above the coherence horizon where this whole picture is replaced (cf. [PROOF_CHAIN_GAP_DOMINANCE](../docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md)).

## 2. Why dark: the survivor is a density standing wave

By the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), `Re λ = −2γ⟨n_XY⟩`: the decay rate is `2γ` per unit of "light" content `⟨n_XY⟩` (the average number of `X`/`Y` Pauli factors), and pure `{I,Z}` operators are dark (rate 0). The survivor is *almost* such an operator. As shown in [PROOF_DIFFUSION_RAYLEIGH_CLOSURE](../docs/proofs/PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md), it is **predominantly a real diagonal density mode**: a dominant `{I,Z}`-diagonal carrying the single-site occupation profile `n(j)`, dressed by a subdominant Hamming-distance-2 coherence admixture. The diagonal is dark; the decay runs entirely through the small admixture the density stirs up. As `Q → 0` the admixture vanishes and `⟨n_XY⟩ → 0`: the strict diffusion limit is the dark, density-only mode. So "dark + real" is just "this is a slowly-relaxing population, not a precessing coherence"; the eigenvalue is real because the dominant content does not rotate.

The shape of `n(j)` is fixed by the diffusion: under strong dephasing the slow density obeys a classical diffusion on the chain with **no-flux (Neumann) ends**, whose slowest non-trivial mode is the lowest Neumann harmonic

    n(j) ∝ cos(π (j − ½) / N).

This is the swell: low at one end, high at the other, zero in the middle, **antisymmetric about the chain centre**. Numerically (the dominant diagonal of the survivor at N=4, mean-removed): `n = [−0.74, −0.46, +0.46, +0.74]`. The two oddnesses are now immediate.

## 3. Why X-odd and R-odd: two readings of the same antisymmetry

**Reflection `R: l → N−1−l`.** A Neumann harmonic with one half-wavelength is antisymmetric: `n(N−1−j) = −n(j)`. As an operator the survivor's diagonal reverses sign under site-reversal, so `R` returns eigenvalue **−1**. (Verified directly: the diagonal density vector satisfies `n · n_reflected / ‖n‖² = −1.0` at N=4.) This is the same chain reflection that commutes with `L` ([`ChainMirror`](../compute/RCPsiSquared.Core/Symmetry/ChainMirror.cs) / [PROOF_C1_MIRROR_SYMMETRY](../docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md)) and organizes the entire dephased-chain spectrum into R-even and R-odd blocks ([SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md)); the survivor sits squarely in the R-**odd** block, and this oddness coincides with its `X^⊗N`-oddness because the F71 reflection eigenspaces and the half-filling complement read the same antisymmetric profile.

**Global spin-flip `X^⊗N`.** `X_l` flips spin `l` up↔down; `X^⊗N` flips all of them, which on a computational basis state is the bit-complement `a → ã`. At **half-filling**, the complement `ã` is exactly the **particle-hole conjugate**: every occupied site becomes empty and vice-versa, so the occupation profile flips sign, `n(j) → −n(j)`. The survivor's diagonal is therefore complement-odd, and `X^⊗N` returns eigenvalue **−1**. (Verified: `n · n_complement / ‖n‖² = −1.0` at N=4.) This works *because* the survivor lives at half-filling, which is where the complement closes the profile onto its own negative; off half-filling the complement maps to a different sector entirely.

A structural footnote worth keeping straight. `X^⊗N` is the project's **F1² = Π²** (the square of the palindrome conjugation Π; F1 in [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md)), the global charge-conjugation that pairs sectors `(p_c, p_r) ↔ (N−p_c, N−p_r)`. The survivor's `X^⊗N`-oddness is **not** the same thing as `n_XY`-parity (the F61 selection rule, which is conjugation by the *other* global string `Z^⊗N`). The survivor has **even** `n_XY` (it mixes `n_XY = 0` and `2`, the `{0,2}` family) and is **odd** under `X^⊗N`. Two different Z₂'s, two different signs; only the spin-flip one reads the density profile's antisymmetry.

## 4. Why the staggered particle-hole fails: the ZZ term

The natural condensed-matter move at half-filling is the **staggered (bipartite) particle-hole** transformation `U = X^⊗N · ∏_{l odd} Z_l`, i.e. complement the occupations and dress with a sublattice `(−1)^l` sign (the chiral `K = diag((−1)^sublattice)` of AZ class BDI, [`ChiralKClaim`](../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) / [BIPARTITE_CHIRALITY_DIAGONAL_CELL](BIPARTITE_CHIRALITY_DIAGONAL_CELL.md), whose defining property is `KHK = −H` for the *tight-binding hopping*). One might expect *this* to be the survivor's defining self-mirror. It is not, and the reason is exactly which couplings the chain carries.

`U` implements the free-fermion `E → −E` symmetry: it sends the hopping `X_iX_{i+1}+Y_iY_{i+1} → −(X_iX_{i+1}+Y_iY_{i+1})`. That is a clean symmetry of the **pure XY** chain (it flips the single-particle dispersion `ε_k → −ε_k`). But the Heisenberg chain also has the diagonal `Z_iZ_{i+1}` term, which the staggered `(−1)^l` dressing does **not** send to its negative. So `U` is not a symmetry of the Heisenberg `L` at all: the verifier measures `[U, L] ≠ 0` ("neither", fidelity-to-eigenvector only 0.55–0.82, falling further as Q grows and the coherence dressing thickens). The survivor is `X^⊗N`-odd and `R`-odd, but it is simply *not* a `U`-eigenstate. The plain global flip `X^⊗N` (no sublattice phase) is the operator that closes on the survivor, because `X^⊗N` reads the *density profile's* antisymmetry, which the `ZZ` term respects, rather than the *hopping's* sign, which it does not.

This also reconciles cleanly with [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md) Step 6, which studied the *Bogoliubov*-mode particle-hole `U_PH` on the XY chain and found its conjugation lift to operator space *preserves* R-parity (the Hilbert-space anticommutation `{U_PH, R} = 0` is "squared away" by the bra-ket product). That is a different operator (a momentum-space free-fermion construction on pure XY) and a different question (the whole spectrum's R-block isospectrality). Here, on Heisenberg matter, the real-space staggered `U` does not even survive to operator space as a symmetry, because the `ZZ` term breaks it before any lifting question arises.

## 5. The {0,2} connection, in one line

The survivor lives entirely on **diagonal** popcount sectors `(p,p)` (equal particle number in bra and ket). Within such a sector the only coherences are at Hamming distance 0 (the dark density diagonal) and Hamming distance 2 (the rate-bearing admixture, since a single `XX`/`YY` hop changes two bits). So the survivor is built from exactly the **`{0,2}` disagreement classes**: HD-0 dark population + HD-2 light coherence. The off-diagonal `(p, p±1)` sectors, which would carry HD-1 single-magnon coherences, are pinned at the `2γ` floor by F50 and cannot host a slow mode. "Lives in `{0,2}`", "dark", and "real" are the same statement: the survivor is the slow relaxation of a half-filling **density** wave, whose only handle for the dissipator is the faint two-site coherence it stirs up.

## 6. Scope and what is and isn't new

- **Regime.** Heisenberg (`XX+YY+ZZ`), dispersive topology (chain or ring), strong dephasing below the coherence horizon `Q*(N)`, where the survivor is the half-filling `(N/2, N/2)` density mode. Above the horizon the survivor hands over to the rigid, bright, oscillating `(0,1)` band edge (`⟨n_XY⟩ = 1`), and facts 1–4 no longer apply. On the **star** (no dispersion) the survivor sits at the popcount boundary, not half-filling, and is outside this picture (`SurvivalIncompletenessMirrorClaim`'s counterexample).
- **Odd N.** There is no exact `(N/2, N/2)` sector, so `X^⊗N` has no fixed half-filling sector to act within; the clean fixed-point test is an even-N statement (the verifier skips odd N gracefully). The reflection-oddness, being purely about the spatial profile, is not so restricted.
- **Already known / typed.** The survivor's identity and darkness (`SurvivalIncompletenessMirrorClaim`, Tier 1 candidate); the antisymmetric density profile and the diffusion mechanism (`SurvivorDiffusionGradientClaim` / F123, Tier 1 candidate); the spectrum-wide R-parity block structure and the Bogoliubov-`U_PH` analysis ([SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md)); `X^⊗N = F1² = Π²` and its distinctness from `n_XY`-parity (F1, F61).
- **What this note adds.** The *unified* reading: that X^⊗N-odd, R-odd, not-staggered-PH-fixed, and dark+real are one fact about the half-filling antisymmetric density standing wave; the clean separation of the two relevant Z₂'s (spin-flip vs `n_XY`-parity); and the sharp reason the textbook staggered particle-hole *fails* here (the `ZZ` term it cannot negate), distinguishing the real-space staggered `U` from the XY Bogoliubov `U_PH`. This is a reading-level synthesis over typed Tier-1-candidate parents, not a new derivation.
- **The seam to stitch.** The survivor's *sector* (half-filling = the `X^⊗N` self-paired `(N/2, N/2)`, [`XGlobalChargeConjugationPairing`](../compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs)) is typed, and the *operators* `X^⊗N`, `R`, `K` are each typed; but the *assignment of the survivor's eigenvalue* under them (`−1`, `−1`, none) is not yet a typed Claim or live witness. It is the natural typing candidate here, a witness in the spirit of `QuditPartialPalindromeWitness` with parents `XGlobalChargeConjugationPairing` + `ChainMirror`/`F71` + `SurvivalIncompletenessMirrorClaim`. The deepest survivor arc, [`felt_time_dimensions`](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs) (open), already frames the same darkness as "the survivor lives in the disagreement kernel `k=0`, the diagonal's blind spot" and notes "the survivor also has `ω=0`"; this note is the symmetry-axis reading of that same blind spot.

---

## Cross-references

- Verifier: [`simulations/_survivor_particle_hole_mirror.py`](../simulations/_survivor_particle_hole_mirror.py) (gate-first; `X → −1`, `R → −1`, `U_PH` not-fixed, with the math-lens guard that `X` and `U` genuinely complement diagonals).
- Survivor identity + darkness: [`SurvivalIncompletenessMirrorClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/SurvivalIncompletenessMirrorClaim.cs) (`inspect --root survivor`), [`simulations/carbon/incompleteness_survivor.py`](../simulations/carbon/incompleteness_survivor.py), [`reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md`](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md).
- The density standing wave: [PROOF_DIFFUSION_RAYLEIGH_CLOSURE](../docs/proofs/PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md) / [`SurvivorDiffusionGradientClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs) (`inspect --root gradient`), F123 in [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md).
- Rate = light content: [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md).
- Spectrum-wide R-parity + the Bogoliubov particle-hole: [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md), [MAJORANA_AXIS_MODES](MAJORANA_AXIS_MODES.md).
- Band-edge handover (the regime where this picture ends): [PROOF_CHAIN_GAP_DOMINANCE](../docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md), [PROOF_RING_GAP_DOMINANCE](../docs/proofs/PROOF_RING_GAP_DOMINANCE.md), [`CoherenceHorizonClaim`](../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs).
- The two global Z₂ strings: F1 (`X^⊗N = Π²`) and F61 (`n_XY` parity, `Z^⊗N`) in [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md); the staggered/bipartite chiral `K`: [BIPARTITE_CHIRALITY_DIAGONAL_CELL](BIPARTITE_CHIRALITY_DIAGONAL_CELL.md).
- The `{0,2}` disagreement reading: [`reflections/ON_THE_ONE_DIAGONAL.md`](../reflections/ON_THE_ONE_DIAGONAL.md), [CHAIN_GAP_SECTOR_DIAGNOSTIC](CHAIN_GAP_SECTOR_DIAGNOSTIC.md).
