# The N=4 Twin-Scalar Crossing and the Residual Additivity Hypothesis

**Status:** Exact N=4 twin-scalar restriction and Tier-2 residual interpretation. The N=4 path-3 character is Tier-1-derived (`F89Path3OcticEpClaim`); the sampled Delta response is finite-N evidence, not proof of integrability causality or an all-N protection theorem.
**Date:** 2026-06-22
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Origin:** the generative ("why") pass after the F89-octic diabolic *character* correction landed (Plan A + Plan B, master). The typed layer proves *that* it is diabolic; this asks *why* it is diabolic rather than the generic defective. See `docs/CAUGHT_ERRORS.md` (the EP-character trilogy) and `compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs`.

## The question

The F89 path-3 octic degeneracy (the (single-excitation SE, double-excitation DE) coherence block of the 4-site XY chain under Z-dephasing, at q_EP = √((−1+√13)/6) ≈ 0.659, λ_EP = −4γ + 2iJ) is **diabolic**: the two eigenvalues coalesce but the eigenvectors stay independent (geometric multiplicity g1 = algebraic multiplicity g2 = 2, departure-from-normality dep = 0, the 2×2 restriction L|₂D = λ·I, no Jordan block). This is settled and confirmed artifact-free (`F89Path3OcticEpClaim`, `inspect --root f89octic`).

But a one-parameter coalescence of a **non-normal** operator is *generically defective* (the eigenvectors also coalesce, a Jordan block / square-root EP, codimension 1); a diabolic point normally needs codimension ≥ 3 (the von Neumann-Wigner count). So a single knob `q` producing a diabolic point is non-generic; it needs a reason. The exact double discriminant factor locates the degeneracy and, after pair isolation, supports the analytic-crossing reading; it does not decide semisimplicity. The typed record now gives the twin-scalar restriction and live character diagnostics as the load-bearing route. This note asks why that scalar restriction occurs in an *irreducible* octic.

## The exact local restriction and the conditional interpretation

The N=4 twin-scalar restriction supplies semisimplicity at the stated point. Free-fermion additivity supplies its Hamiltonian multiplet, but the dephasing restriction is a separate condition. The sampled XXZ anisotropies below change the character. Since uniform XXZ remains Bethe-integrable, this control probes departure from free-fermion additivity, not loss of Hamiltonian integrability.

### The mechanism (from below)

`L` on the (SE,DE) block is built from two pieces: the hopping `H_eff` and the Z-dephasing `D`. At q_EP, **both restrict to scalars on the 2D coalescing eigenspace**, so their sum is scalar, no off-diagonal Jordan coupling:

- `H_eff|₂D = 2iJ·I` (scalar). **Why:** the XY model is free-fermion, so the two-excitation (DE) energies are exact **sums** of one-excitation (SE) energies, E_(a,b) = ε_a + ε_b (verified: M_DE spectrum = the 6 pairwise sums of M_SE = 4J·cos(kπ/5), to 1.3e-15). The coalescing pair descends from one **4-fold-degenerate free-fermion multiplet**, on which the hopping is a multiple of the identity.
- `D|₂D = −4γ·I` (scalar). **Why:** the dephasing rate is −2γ·n_diff (n_diff = the number of sites where the bra and ket basis labels differ) = −6γ + 4γ·p (p = overlap fraction), so p = ½ ⟺ rate −4γ. q_EP sits exactly at this **overlap-balanced rate midpoint**, the AT-spectral midpoint (AT = Absorption Theorem; the midpoint of the absorption rungs −2γ and −6γ) between the rate-2γ (overlap) F_a modes and the rate-6γ (no-overlap) F_b modes (`F89PathKAtLockMechanismClaim`; the F_a/F_b pair is the overlap↔no-overlap mirror, fixed point at 4γ). At p = ½ the dephasing-projector is exactly ½·I on the pair, i.e. a uniform shift, not a coupling.

So **L|₂D = (−4γ + 2iJ)·I = λ_EP·I**, semisimple. The discriminant double-zero is the *algebraic shadow* of this twin scalarity, not an independent cause.

### The decisive gate (XXZ-Δ breaks free-fermion additivity)

Probe `simulations/f89_zz_break_gate.py` (gate-first; Stage 0 reproduces the Δ=0 diabolic point and confirms the Pauli build equals the committed reference + is a genuine sub-block of the full 256² Liouvillian, both to 0.0e+00). Turning on Δ (the ZZ term makes the two magnons interact, so E_(a,b) ≠ ε_a + ε_b):

| Δ | q* | min pair-dist | g1 | g2 | dep | character |
|---|---|---|---|---|---|---|
| 0.00 | 0.658983 | 6.8e-15 | 2 | 2 | 0.000 | **DIABOLIC** |
| 0.02 | 0.660249 | 3.4e-05 | 1 | 2 | 0.022 | DEFECTIVE (Jordan) |
| 0.05 | 0.662459 | 7.5e-05 | 1 | 2 | 0.056 | DEFECTIVE |
| 0.10 | 0.667060 | 8.6e-05 | 1 | 2 | 0.112 | DEFECTIVE |
| 0.20 | 0.644962 | 1.0e-04 | 1 | 2 | 0.165 | DEFECTIVE |
| 0.50 | 0.639578 | 4.1e-04 | 1 | 2 | 0.422 | DEFECTIVE |

At the sampled nonzero Δ values, g1 changes 2 → 1 (g2 stays 2), the eigenvector-merge |cos| changes 0.60 → 1.0000, and dep grows approximately linearly over the tested range. A complex-q locator reduces the split to ~1e-8 and finds the defective partner on the real-q axis, consistent with the surviving pseudo-Hermiticity Σ L Σ = L†. The restriction probe finds `H_eff|₂D = 2iJ·I` at Δ=0 and non-scalar restrictions at the tested nonzero Δ values (diagonal deviation 0 → 0.30), while the dephasing half stays scalar at the −4γ midpoint. These finite-N measurements support the conditional residual twin-scalar interpretation; they neither isolate integrability as a cause nor establish the response at unsampled Δ or N.

## What it is NOT (candidates the gate refuted)

The following controls distinguish alternative explanations at the N=4 point; excluding them does not prove integrability causality:

- **Not a commuting-symmetry separation.** No symmetry puts the two modes in different sectors so they cannot couple: the site-reflection R commutes but both modes are R = +1; the overlap↔no-overlap involution does not even commute ([S,L] ≈ 19.6); the chiral Σ gives Σ L Σ = L† and an antilinear PT, but relates λ_EP to its *conjugate* −4γ−2iJ, not to its degenerate partner.
- **Not a bare free-fermion crossing "diabolic by construction."** The EP frequency Im(λ_EP)/J = 2 is *absent* from the bare γ=0 difference set; the pair is born inside a 4-fold free-fermion multiplet (at Im/J = −1+√5) and is split by dephasing, *re-coalescing* only at the tuned q_EP. Free-fermion supplies the degenerate origin, not a ready-made crossing.
- **Not a factorisation/resultant.** The octic is irreducible over Q(i), so it is not two factors whose roots cross (`F89Path3OcticGaloisClaim`; refuted in `docs/CAUGHT_ERRORS.md`).
- **Not Π / palindrome / Kramers.** Π pairs eigen*values* (the spectral palindrome), not eigen*vectors* at a crossing; no anti-unitary T² = −1 Kramers structure is present. The mirror fixes the line and pairing; the N=4 overlap/twin-scalar restriction supplies the on-line position. The mirror alone supplies neither that position nor semisimplicity; see "The line vs the silence" below and `experiments/F89_BRANCH_LOCUS_PALINDROME.md`.

## Where it sits in the framework

The N=4 Delta control probes departure from free-fermion additivity: uniform XXZ remains Bethe-integrable. The measured character flip and loss of the twin-scalar restriction support the Tier-2 hypothesis, but this control does not isolate integrability causality. Absence of level repulsion likewise does not prove integrability (see `experiments/RANDOM_MATRIX_THEORY.md`). The local observable is semisimplicity of the coalescing pair, not a spectral-statistics classification of the Hamiltonian.

The **EP-character trilogy** (`docs/CAUGHT_ERRORS.md`) distinguishes F86a (near-EP, no coalescence), the coherence-horizon defective √-EP, and the F89-octic diabolic. These character readings do not establish a shared integrability cause; the N=4 diabolic has its own twin-scalar restriction, including the overlap p=½ midpoint.

**The line vs the silence (the branch-locus palindrome).** The −4γ AT-midpoint above is exactly the palindrome's mirror centre Re λ = −σ = −4. The palindrome proof is independent of this Tier-2 hypothesis: the exact block identity proves the branch-locus pairing (`experiments/F89_BRANCH_LOCUS_PALINDROME.md`, typed `F89BranchLocusPalindromeClaim`, live at `inspect --root branchpalindrome`). The sampled N=4 Delta control stays on the line yet becomes defective; it is a counterexample separating mirror placement from semisimple character, not a premise of the palindrome proof or proof of general integrability protection. The N=4 character fact has its own twin-scalar restriction. Plain-words sibling: `reflections/ON_WHO_WATCHES_WHOM.md`.

## Open / next

- **Inventory boundary.** The N=6 inventory is complete: the parity-labelled direct-t Route B inventory reconciles 266 Diabolic; 0 Defective loci (133 per parity). Character provenance is 118 exact Hermitian-axis readings, 148 stable numerical EpCharacter readings, and 0 executed exact-rank certificates; exact algebraic boxes do not make numerical character exact. N=7 remains the inventory boundary. See [the current Route B result](../docs/THE_DOUBLE_ROOT.md).
- **Residual mechanism.** [The codimension-1 additivity proof](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §8 separates automatic AT-locked crossings from the conditional twin-scalar residual regime. The D-scalar half is an extra condition, supplied at the N=4 point and in the stated N=5 tests (`TwinScalarDHalfTests`); imaginary-q Hermiticity is a separate structural route to semisimplicity. Neither the completed N=6 character inventory nor the sampled Delta controls prove an all-N residual mechanism. Restriction tests require the HS-orthonormal coherence basis, not an unnormalized orbit compression.

## Anchors
- Typed *what*: `compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs` (Tier 1 derived, the Correction block); live `inspect --root f89octic` (`F89OcticCharacterWitness`); `compute/RCPsiSquared.Core/F89PathK/F89Path3OcticBlock.cs`.
- The AT-lock / overlap-no-overlap mirror: `F89PathKAtLockMechanismClaim`, `F89Path3SeDeFactorisationClaim`; `experiments/F89_PATH_K_GALOIS.md` (§ Path-3 octic diabolic-degeneracy location).
- The *why* probes (gate-first, this hypothesis's evidence): `simulations/f89_zz_break_gate.py` (the decisive Δ-break), `simulations/f89_why_diabolic_probe.py` (the twin-scalar / free-fermion-additivity mechanism).
- The *line vs silence* downstream (uses the Δ-gate): `experiments/F89_BRANCH_LOCUS_PALINDROME.md` (Tier 1), `F89BranchLocusPalindromeClaim`, `inspect --root branchpalindrome`, `reflections/ON_WHO_WATCHES_WHOM.md`. Order/sheet corroboration only: `inspect --root galoismonodromy` (an identity eigenvalue loop says the pair is single-valued; it does not independently decide Jordan character).
