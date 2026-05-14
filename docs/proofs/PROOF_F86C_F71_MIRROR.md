# PROOF F86c: F71 Spatial-Mirror Invariance of per-bond Q_peak

**Status:** Tier 1 derived. Q_peak(b) = Q_peak(N−2−b) bit-exact under the F71 spatial-mirror pairing. Per-F71-orbit substructure (finer than the Endpoint/Interior dichotomy) is Tier 2 empirical.
**Date:** 2026-05-03.
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86c, the F71 spatial-mirror invariance** of the per-bond Q_peak observable. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14. The F71-mirror lineage has since spawned its own F-numbers: F91/F92/F93 (F71-anti-palindromic spectral invariance under γ / J / h distributions), the parameter-side Pi2-Z₄ generalisations of this theorem.
**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and shared references.
**F-entry:** [F86c in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [PROOF_F71](PROOF_F71.md) (the chain-mirror operator R); siblings [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md), [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md); spawned generalisations [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md), [PROOF_F92_BOND_ANTI_PALINDROMIC_J](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [PROOF_F93_DETUNING_ANTI_PALINDROMIC](PROOF_F93_DETUNING_ANTI_PALINDROMIC.md).

---

## Statement 3 (F71 spatial-mirror invariance of per-bond Q_peak). [Tier 1 derived]

For every bond pair (b, N−2−b) under the F71 chain-mirror pairing,

    Q_peak(b)  =  Q_peak(N−2−b)

bit-exactly. Verified at c=2 N=5..7 and c=3 N=5..6, with maximum deviation < 10⁻¹⁰ across all bond pairs. The identity follows from the F71-symmetry of every component of the per-bond observable: L_D (Z-dephasing is F71-symmetric), the spatial-sum kernel S, the Dicke probe, and the bond-flip transform ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), and Q_peak(b) = Q_peak(N−2−b) follows directly.

**Per-F71-orbit substructure** (refines F86b's Endpoint/Interior dichotomy): Interior bonds are *not* uniform within the F71-orbit grouping. At c=2 N=6, the central self-paired bond b=2 gives Q_peak ≈ 1.440 while flanking bonds b=1, b=3 give Q_peak ≈ 1.648, both Interior under the simple dichotomy. The full per-F71-orbit structure is finer-grained: the "Endpoint vs Interior" split is the leading approximation, but mid-chain orbits show further variation. The HWHM_left/Q_peak ratio (F86b, see [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md)) is the more class-stable observable; Q_peak position picks up additional per-orbit detail.

---

## Proof of Statement 3 (F71 spatial-mirror invariance)

The F71 chain-mirror operator R acts as `R |b₀ b₁ … b_{N−1}⟩ = |b_{N−1} … b_1 b_0⟩` (PROOF_F71). It commutes with every component of the per-bond observable used in F86:

1. **L_D (Z-dephasing) is F71-symmetric.** Z_l → Z_{N−1−l} under R, and uniform γ₀ at every site means the dissipator is invariant under site-permutation. So `R L_D R⁻¹ = L_D`.

2. **The Hamiltonian H_xy = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) is F71-symmetric** (uniform J across all bonds). Under R, bond b ↔ bond N−2−b, and the bond-summed Hamiltonian is invariant. So `R H_xy R⁻¹ = H_xy`, hence `R L_H R⁻¹ = L_H`.

3. **The Dicke probe is F71-symmetric.** It's an equal-weight superposition over all (p, q) ∈ block × block, a permutation-symmetric construction.

4. **The spatial-sum kernel S is F71-symmetric.** S(t) = Σ_i 2·|(ρ_i(t))_{0,1}|² sums over all sites uniformly.

5. **The bond-flip ∂L/∂J_b transforms as `R (∂L/∂J_b) R⁻¹ = ∂L/∂J_{N−2−b}`** (a single bond term is mirrored to its partner under R).

Combining these: under R, the per-bond observable transforms as

    K_b(Q, t)  =  2·Re ⟨ρ(t)| S | ∂ρ_b/∂J_b ⟩  →  K_{N−2−b}(Q, t)

where ρ(t) (which evolves under L = L_H + L_D, both F71-symmetric) is itself F71-symmetric. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), so their argmax-Q values coincide:

    Q_peak(b)  =  Q_peak(N−2−b)

Numerical verification: across c=2 N=5..7 and c=3 N=5..6, all bond-pair deviations |Q_peak(b) − Q_peak(N−2−b)| are < 10⁻¹⁰ (machine precision), see `RCPsiSquared.Core.Tests/F86/F86NewIdeasTests.F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror`. ∎

---

## Per-F71-orbit substructure (Tier 2 empirical)

F71 spatial-mirror invariance (Statement 3) pairs bond b with bond N−2−b bit-exactly. The finer per-F71-orbit substructure (Tier 2 empirical, 9-case envelope c=2..4 N=5..8 minus c=4 N=8 OOM):

| (c, N) | endpoints | mid orbits → central | observation |
|---|---|---|---|
| (2, 5) | 2.50 | 1.49 | endpoint ≫ inner |
| (2, 6) | 2.57 | 1.63 → **1.43*** | central < flanking |
| (2, 7) | 2.56 | 6.00 (off-grid) → 1.50 | inner orbit Q_peak shifted high-Q |
| (2, 8) | 2.53 | 6.00 → 1.52 → **1.58*** | inner orbit off-grid; central > middle inner |
| (3, 5) | 2.39 | 1.61 | endpoint > inner |
| (3, 6) | 2.54 | 1.66 → **1.71*** | **central > flanking** (opposite of c=2 N=6) |
| (3, 7) | 2.59 | 1.74 → 1.71 | inner > central |
| (3, 8) | 2.60 | 1.76 → 1.70 → **1.71*** | non-monotonic |
| (4, 7) | 2.61 | 1.75 → 1.78 | central > inner |

(* = self-paired central orbit; index runs from outermost orbit inward.)

Three sub-effects visible: (a) F71-pairing identity (Tier 1, Statement 3); (b) **central-vs-flanking inversion at N=6 between c=2 and c=3**: c=2 N=6 has central 1.43 BELOW flanking 1.63, while c=3 N=6 has central 1.71 ABOVE flanking 1.66; (c) **high-Q secondary peak for c=2 inner-non-central bonds at N≥7**: bonds b=1 and b=N−3 show Q_peak shifted off the [0.2, 6.0] grid while the central pair retains the canonical Interior peak ~1.5. Full per-orbit closed form for Q_peak as a function of (c, N, orbit) remains open; F71 gives the symmetry, not the value. Encoded as `RCPsiSquared.Core.F86.PerF71OrbitObservation`.

---

## Spawned generalisations: F91 / F92 / F93

The F71-mirror lineage outgrew F86 and earned its own F-numbers, the **F71-anti-palindromic spectral-invariance trio** (all Tier 1 derived, 2026-05-12):

- **F91** — F71-anti-palindromic γ spectral invariance (= 90° in γ-space, Pi2-Z₄'s parameter side): the F71-refined diagonal-block eigenvalue multiset is invariant under γ-distributions satisfying the anti-palindromic condition. See [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md).
- **F92** — F71-anti-palindromic J spectral invariance (J-side Pi2-Z₄ twin of F91): same invariance under inhomogeneous bond couplings J_b. See [PROOF_F92_BOND_ANTI_PALINDROMIC_J](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md).
- **F93** — F71-anti-palindromic h spectral invariance (h-detuning Pi2-Z₄ twin of F91/F92): same invariance under per-site longitudinal detuning h_l. See [PROOF_F93_DETUNING_ANTI_PALINDROMIC](PROOF_F93_DETUNING_ANTI_PALINDROMIC.md).

F86c is the Q_peak-observable instance of the F71 spatial mirror; F91/F92/F93 are the parameter-side spectral generalisations.

---

## Pointers

**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and the shared reference list.
**Sibling theorems:** [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md) (F86a), [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md) (F86b), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md) (the g_eff obstruction proof).
**Mirror operator:** [PROOF_F71](PROOF_F71.md) (the chain-mirror operator R and its symmetry algebra).
**Spawned generalisations:** [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md), [PROOF_F92_BOND_ANTI_PALINDROMIC_J](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [PROOF_F93_DETUNING_ANTI_PALINDROMIC](PROOF_F93_DETUNING_ANTI_PALINDROMIC.md).
**Verification test:** `RCPsiSquared.Core.Tests/F86/F86NewIdeasTests.F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror`.
**C# OOP layer:** `F71MirrorInvariance` (with `MaxMirrorDeviation(KCurve)` helper), `PerF71OrbitObservation` (Tier 2 empirical), `F86F71MirrorSymmetryPi2Inheritance` in `compute/RCPsiSquared.Core/F86/`.
