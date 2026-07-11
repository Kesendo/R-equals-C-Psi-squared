# The Octic Crossing Is Diabolic Because the Magnons Don't Interact

**Status:** Tier 2: the diabolic→defective flip under an XXZ anisotropy is a **decisive gate-first computation** (N=4 path-3, machine precision); the "free-fermion integrability is the protection" reading is the (well-grounded) interpretation. The diabolic *fact* it explains is itself Tier-1-derived (`F89Path3OcticEpClaim`).
**Date:** 2026-06-22
**Updated:** 2026-06-25 (full re-read: every formula and number re-verified bit-exact from the committed machinery; term glosses for a cold reader; the branch-locus palindrome connection added, the now-closed markdown-sweep item retired).
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Origin:** the generative ("why") pass after the F89-octic diabolic *character* correction landed (Plan A + Plan B, master). The typed layer proves *that* it is diabolic; this asks *why* it is diabolic rather than the generic defective. See `docs/CAUGHT_ERRORS.md` (the EP-character trilogy) and `compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs`.

## The question

The F89 path-3 octic degeneracy (the (single-excitation SE, double-excitation DE) coherence block of the 4-site XY chain under Z-dephasing, at q_EP = √((−1+√13)/6) ≈ 0.659, λ_EP = −4γ + 2iJ) is **diabolic**: the two eigenvalues coalesce but the eigenvectors stay independent (geometric multiplicity g1 = algebraic multiplicity g2 = 2, departure-from-normality dep = 0, the 2×2 restriction L|₂D = λ·I, no Jordan block). This is settled and confirmed artifact-free (`F89Path3OcticEpClaim`, `inspect --root f89octic`).

But a one-parameter coalescence of a **non-normal** operator is *generically defective* (the eigenvectors also coalesce, a Jordan block / square-root EP, codimension 1); a diabolic point normally needs codimension ≥ 3 (the von Neumann-Wigner count). So a single knob `q` producing a diabolic point is non-generic; it needs a reason. The typed record says only "intrinsic to the double-zero of the discriminant," which *restates* semisimplicity (an even-multiplicity discriminant factor `(3q⁴+q²−1)²`) rather than explaining why an *irreducible* octic has one.

## The answer: free-fermion integrability

**The XY chain is free-fermion integrable, and that is the protection.** The moment the magnons are made to interact (an XXZ ZZ-anisotropy Δ ≠ 0), the degeneracy turns **defective**, at every Δ tested, instantly.

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

The instant Δ ≠ 0: g1 collapses 2 → 1 (g2 stays 2 = a Jordan block), the eigenvector-merge |cos| jumps 0.60 → 1.0000, and dep turns on and grows ~linearly with Δ. A complex-q locator confirms a **genuine** defective EP (split driven to ~1e-8) that stays **on the real-q axis** (pinned by the surviving pseudo-Hermiticity Σ L Σ = L†) but flips character diabolic → defective. The mechanism probe shows exactly the predicted cause: `H_eff|₂D` is scalar (2iJ·I) at Δ=0 and **non-scalar at every Δ > 0** (diagonal deviation grows 0 → 0.30), while the dephasing half stays scalar at the −4γ midpoint. Break free-fermion additivity → the H_eff-scalar half dies → L|₂D acquires the Jordan coupling.

**Generic XXZ gives the ordinary defective EP; the diabolic point is the integrable special case.**

## What it is NOT (candidates the gate refuted)

The honest record: each was proposed and killed from below, which is why "integrability" survives by elimination:

- **Not a commuting-symmetry separation.** No symmetry puts the two modes in different sectors so they cannot couple: the site-reflection R commutes but both modes are R = +1; the overlap↔no-overlap involution does not even commute ([S,L] ≈ 19.6); the chiral Σ gives Σ L Σ = L† and an antilinear PT, but relates λ_EP to its *conjugate* −4γ−2iJ, not to its degenerate partner.
- **Not a bare free-fermion crossing "diabolic by construction."** The EP frequency Im(λ_EP)/J = 2 is *absent* from the bare γ=0 difference set; the pair is born inside a 4-fold free-fermion multiplet (at Im/J = −1+√5) and is split by dephasing, *re-coalescing* only at the tuned q_EP. Free-fermion supplies the degenerate origin, not a ready-made crossing.
- **Not a factorisation/resultant.** The octic is irreducible over Q(i), so it is not two factors whose roots cross (`F89Path3OcticGaloisClaim`; refuted in `docs/CAUGHT_ERRORS.md`).
- **Not Π / palindrome / Kramers.** Π pairs eigen*values* (the spectral palindrome), not eigen*vectors* at a crossing; no anti-unitary T² = −1 Kramers structure is present. (The palindrome does fix the diabolic's *position*, on the mirror line Re λ = −4, the AT-midpoint being the palindrome centre, but not its silence; see "The line vs the silence" below and `experiments/F89_BRANCH_LOCUS_PALINDROME.md`.)

## Where it sits in the framework

This is a concrete, from-below-verified instance of the repo's recurring **integrability → no level repulsion → crossings allowed** theme (cf. `experiments/RANDOM_MATRIX_THEORY.md`; the EP-vs-level-crossing dichotomy in `docs/THE_TRICHOTOMY_SEEN.md` and `experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md`). The novelty here is *making the diabolic↔defective boundary itself the observable*: tuning the single integrability-breaking knob Δ flips the EP character, and the from-below mechanism (twin scalar restriction, free-fermion multiplet + p=½ midpoint) says exactly why.

It also sharpens the **EP-character trilogy** (`docs/CAUGHT_ERRORS.md`): F86a (near-EP, no coalescence) / coherence-horizon (defective √-EP) / F89-octic (diabolic). The trilogy's three different outcomes now have a unifying axis: the F89 octic is diabolic precisely because it is the *integrable* member; integrability-breaking would carry it to the coherence-horizon-style defective EP.

**The line vs the silence (the branch-locus palindrome).** The −4γ AT-midpoint above is exactly the palindrome's mirror centre Re λ = −σ = −4. The branch-locus result (`experiments/F89_BRANCH_LOCUS_PALINDROME.md`, the typed `F89BranchLocusPalindromeClaim`, live at `inspect --root branchpalindrome`) rests on this hypothesis: it uses the XXZ Δ ≠ 0 counterexample as its proof that on-the-line does NOT imply silent. Turning on Δ keeps the EP on Re λ = −4 (the dephasing half stays scalar at the midpoint) yet defects it. So the line is the palindrome's gift; the silence is integrability's, this hypothesis's. What is an interpretation here (the Δ ≠ 0 gate) is load-bearing there, inside a Tier-1 result. Plain-words sibling: `reflections/ON_WHO_WATCHES_WHOM.md`.

## Open / next

- **Generalize beyond N=4 path-3.** The gate is decisive at N=4; the "integrability protects the diabolic crossing" principle should be checked at path-5/6 and against the F90 (SE,DE) bridge: does every free-fermion (SE,DE) coalescence at the AT-midpoint come out diabolic, and does Δ always defect it?
- **Codimension argument: LANDED 2026-07-02.** Stated as a theorem in [the codimension-1 additivity proof](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §8 (Theorem A): two regimes, the automatic AT-locked crossings and the conditional twin-scalar residual regime, where additivity supplies the H-scalar half identically in q (this doc's degenerate-multiplet picture, now the "twin-scalar lemma") and the D-scalar half remains the genuine extra condition, proven here at the N=4 point. The honest general form was "codim-3 → codim-≤2, and codim-1 exactly where the D-half is supplied". RESOLVED 2026-07-02: the D-half check ran at N=5 (PROOF_CODIM1 §8, gate `TwinScalarDHalfTests`). The D-half is SUPPLIED (twin-scalar) at every genuinely-complex-q (Re q ≠ 0) N=5 residual diabolic, so the twin-scalar route extends from real q to complex q (codim-1); the pure-imaginary-q (Re q = 0, λ real) diabolics are semisimple by Hermiticity (the block is real-symmetric there), a separate mechanism. So the general form upgrades to codim-1 wherever the coalescence is twin-scalar. (Metric caveat: the test must run in the HS-orthonormal coherence basis; the ×2-cleared orbit basis is non-orthonormal at odd nBlock and gives a spurious non-scalar reading, which inverted a first result before it was caught.)
- **Markdown sweep: done.** The diabolic octic verdict and this WHY now live in `docs/ANALYTICAL_FORMULAS.md` (§ "F89 path-3 octic diabolic degeneracy", the "WHY diabolic" paragraph) and `compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md`, not only the typed C#. Gap closed.

## Anchors
- Typed *what*: `compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs` (Tier 1 derived, the Correction block); live `inspect --root f89octic` (`F89OcticCharacterWitness`); `compute/RCPsiSquared.Core/F89PathK/F89Path3OcticBlock.cs`.
- The AT-lock / overlap-no-overlap mirror: `F89PathKAtLockMechanismClaim`, `F89Path3SeDeFactorisationClaim`; `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` (§ Path-3 octic diabolic-degeneracy location).
- The *why* probes (gate-first, this hypothesis's evidence): `simulations/f89_zz_break_gate.py` (the decisive Δ-break), `simulations/f89_why_diabolic_probe.py` (the twin-scalar / free-fermion-additivity mechanism).
- The *line vs silence* downstream (uses the Δ-gate): `experiments/F89_BRANCH_LOCUS_PALINDROME.md` (Tier 1), `F89BranchLocusPalindromeClaim`, `inspect --root branchpalindrome`, `reflections/ON_WHO_WATCHES_WHOM.md`. Independent character corroboration: `inspect --root galoismonodromy` (a loop around q_EP returns the identity, the pair does not swap).
