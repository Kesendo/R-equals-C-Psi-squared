# GHZ + W Sector Mix: Lifting pair-CΨ(0) Above the Fold at N = 3

<!-- Keywords: GHZ_3, W_3, sector mixing, pair-CPsi, fold boundary, F69,
sextic minimal polynomial, irreducible over Q, degree 6 algebraic number,
three-tangle, CKW, F60, F61, F62, spherical scan artifact, product state -->

**Status:** Tier 1, characterized by irreducible sextic, verified N = 3 (positive), N = 4, 5, 6 (negative)
**Date:** 2026-04-17
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:**
- [ghz_w_optimum_n3.py](../simulations/ghz_w_optimum_n3.py) (main: symbolic + numerical optimum, sextic, irreducibility)
- [cpsi_sector_mix_optimization.py](../simulations/cpsi_sector_mix_optimization.py) (original sweep + Kingston dynamics)
- [cpsi_birth_landscape.py](../simulations/cpsi_birth_landscape.py) (broad state survey that found the GHZ+W direction)
- [sector_mix_spherical_artifact.py](../simulations/sector_mix_spherical_artifact.py) (product-state diagnosis of the spherical-scan peak)

**Outputs:**
- [ghz_w_optimum_n3.txt](../simulations/results/ghz_w_optimum_n3.txt)
- [cpsi_sector_mix_optimization.txt](../simulations/results/cpsi_sector_mix_optimization.txt)
- [cpsi_birth_landscape.txt](../simulations/results/cpsi_birth_landscape.txt)
- [sector_mix_spherical_artifact.txt](../simulations/results/sector_mix_spherical_artifact.txt)

**Depends on:**
- F60 (GHZ_N born below the fold), F62 (CΨ(0) for W_N)
- F61 (n_XY parity selection rule) via [PROOF_PARITY_SELECTION_RULE](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md)
- [Engineering Blueprint Rule 1](../publications/ENGINEERING_BLUEPRINT.md) (note April 16, 2026: preparation-vs-evolution asymmetry of F61 as the structural reason F69 is allowed)

**Registry entry:** [F69 in ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md)

---

## What this document is about

F60 places GHZ_N below CΨ = 1/4 at t = 0 for every N ≥ 3. F62 places W_N below the fold for N ≥ 3 as well. Both results come from very different structural reasons: GHZ has all its off-diagonal weight in the extreme-XY-weight sectors (w = 0, w = N), while W lives entirely in the single-excitation sector. Neither family alone escapes the fold.

The smallest family that mixes GHZ and W is the one-parameter symmetric superposition

    |ψ(α)⟩ = α |GHZ_3⟩ + √(1-α²) |W_3⟩,    α ∈ [0, 1]

This document records what happens inside that family at N = 3, shows that the optimum is an algebraic number of degree exactly 6, documents the hardware-readable consequences, records the scope failure at N ≥ 4, and pre-empts a known pitfall (the 3-state spherical-scan artifact) that sits one drawer over in the codebase.

## The result

For the family above, the minimum pair-CΨ(0) over all three two-qubit reductions equals the mean and the maximum (permutation symmetry is exact). The stationarity condition dCΨ/dα = 0, rationalized in x = α², gives the integer-coefficient sextic

    2900 x⁶ - 8060 x⁵ + 4211 x⁴ + 3832 x³ - 2428 x² - 512 x + 300 = 0

which is irreducible over ℚ. The optimum α²_opt is therefore an algebraic number of degree exactly 6, with no expression in nested square roots or cube roots.

Numerical values at the exact sextic root (α²_opt from sympy nroots at 25-digit precision, all other quantities computed from it):

| quantity | value |
|----------|-------|
| α²_opt | 0.375420720711069 |
| α_opt | 0.612715856422101 |
| β_opt | 0.790303283106512 |
| min pair-CΨ(0) | 0.320411541127025 |
| ratio to 1/4 | 1.281646× |

scipy's bounded minimizer reaches these same values to within 3.7 × 10⁻¹⁰ in α²_opt; the values are shown at the sympy precision for internal consistency (scipy's α_opt and β_opt would each differ from the above in the 10th digit).

The same family at N ≥ 4 stays below the fold: best min pair-CΨ(0) is 0.167 (N = 4), 0.146 (N = 5), 0.134 (N = 6). F69 is a strictly N = 3 statement.

---

## Derivation of ρ_AB(α)

|GHZ_3⟩ has support on indices {0, 7} (computational-basis kets |000⟩ and |111⟩, amplitudes 1/√2 each). |W_3⟩ has support on {1, 2, 4} (kets |001⟩, |010⟩, |100⟩, amplitudes 1/√3 each). The two supports are disjoint, so

    |ψ(α, β)⟩ = α |GHZ_3⟩ + β |W_3⟩,    β = √(1-α²)

has real amplitudes α/√2 at indices 0, 7 and β/√3 at indices 1, 2, 4. Tracing out qubit 2 (LSB) from ρ = |ψ⟩⟨ψ| gives a permutation-symmetric 4×4 reduction ρ_AB whose entries are polynomials in α and α·β:

    ρ_AB[00,00] = 1/6 + α²/6     ρ_AB[00,ij] = α β / √6 · [1,1,1]
    ρ_AB[ij,kl] = (1/3 - α²/3) · off-W-block
    ρ_AB[11,11] = α²/2

from which direct computation gives closed forms for purity and L1-norm of the off-diagonal:

    C(α)      = Tr(ρ_AB²) = -5α⁴/18 + 2α²/9 + 5/9
    L1_off(α) = √6 · α · β + (2/3) β²       (for α, β ≥ 0)
    CΨ(α)     = C(α) · L1_off(α) / 3

Full symbolic derivation and cross-check against the direct partial-trace of the 8×8 density matrix is in `ghz_w_optimum_n3.py` (function `symbolic_cpsi`).


## Stationarity and the sextic

Substituting α = sin(t), β = cos(t) with t ∈ (0, π/2) removes the absolute-value sign from |β| and makes CΨ(t) a clean trigonometric polynomial. Differentiating with respect to t yields a stationarity condition that contains a single non-polynomial piece, sqrt(x(1-x)) where x = α². Writing

    dCΨ/dx = A(x) + B(x) · √(x(1-x))

and squaring to clear the square root gives a polynomial equation A(x)² = B(x)² · x · (1-x). After clearing the integer content, the reduced stationarity polynomial is the sextic

    P(x) = 2900 x⁶ - 8060 x⁵ + 4211 x⁴ + 3832 x³ - 2428 x² - 512 x + 300

**Irreducibility.** `sympy.Poly(P, x).is_irreducible` returns True, and `sympy.factor_list(P)` returns P itself as the sole factor with multiplicity one. The polynomial does not split over ℚ into factors of degree 2 or 3 (or any other combination).

**Root structure.** The six roots of P (from `sympy.nroots`, 20-digit precision):

| i | x_i |
|---|-----|
| 0 | 0.37542072071106904717 |
| 1 | 0.59846568613722802034 |
| 2 | -0.48574889168971814544 - 0.05339375077710840516·i |
| 3 | -0.48574889168971814544 + 0.05339375077710840516·i |
| 4 | 1.38846086067936271510 - 0.01617316883305075340·i |
| 5 | 1.38846086067936271510 + 0.01617316883305075340·i |

Only two roots are real in [0, 1]: x₀ ≈ 0.3754 and x₁ ≈ 0.5985. The first is the maximum of CΨ (the optimum); the second is a local minimum inside the interval, above which CΨ decreases again as α approaches 1. Roots 2-5 have either negative real part or lie outside the unit interval and are unphysical for |α| ∈ [0, 1].

**Optimum verified.** scipy's bounded minimize finds t_opt = 0.659493 giving α²_opt = 0.375420720338904. The sympy root is x₀ = 0.375420720711069. The two agree to 3.7 × 10⁻¹⁰. A 401-point grid sweep reproduces CΨ_opt = 0.320411541127025 to within 5 × 10⁻⁸ of the scipy value.

---

## The numerical optimum

All reductions coincide to machine precision:

| quantity | value |
|----------|-------|
| min pair-CΨ(0) | 0.320411541127025 |
| mean pair-CΨ(0) | 0.320411541127025 |
| max pair-CΨ(0) | 0.320411541127025 |
| max − min | 0 (< 10⁻¹⁵) |

Permutation symmetry is exact because the state is symmetric under any qubit permutation: both |GHZ_3⟩ and |W_3⟩ are symmetric functions of the three qubits, so every pair reduction is the same 4×4 matrix.

## Entanglement structure at the optimum

CKW 3-tangle and pairwise concurrences (independent check via `ghz_w_optimum_n3.py`, Part 6):

| quantity | value |
|----------|-------|
| τ_ABC (3-tangle) | 0.799453 |
| C(A, B) | 0.0210 |
| C(A, C) | 0.0210 |
| C(B, C) | 0.0210 |

Almost all the entanglement in |ψ_opt⟩ is genuinely three-body: τ ≈ 0.80 out of a maximum 1.0 (GHZ-like limit), while every pair concurrence is ~0.02, essentially zero. The pair reductions look nearly separable, while the state as a whole is far from any biseparable decomposition. This matches the structural role of |GHZ⟩ in the mixture, which contributes pure three-body coherence and zero pair concurrence; |W⟩ contributes the single-excitation structure that gives ρ_AB its L1 weight.

**Why pair-CΨ is high despite low pair concurrence.** CΨ = Tr(ρ²) · L1_off(ρ) / 3. Concurrence measures entanglement specifically; CΨ measures the product of purity and coherence, either of which can be high without genuine pair entanglement. At the optimum, C(α) = 0.5998 (above the 1/2 of a maximally mixed state) and L1_off = 1.6025, both substantial. The pair is not entangled; it is *coherent* and *pure enough*. CΨ = 1/4 is a fold in the purity × coherence plane, not in the entanglement plane.

---

## Relation to F60, F61, F62

F60 exact: |GHZ_3⟩ → CΨ(0) = 0 on every pair (the off-diagonal ρ_{03} = 1/2 survives globally but traces to zero in any pair reduction since |000⟩ and |111⟩ have no matching pair-content at any bit position).

F62 exact: |W_3⟩ → CΨ(0) = 10/81 ≈ 0.1235 on every pair.

F69 finds the optimum of α|GHZ_3⟩ + β|W_3⟩ at 0.3204, which is 2.6× the W_3 value and above the 1/4 fold. The mixture achieves what neither parent can.

### Why F61 does not forbid the mixing

F61 (the n_XY parity selection rule) block-diagonalizes the Liouvillian by n_XY parity: every eigenmode has definite even or odd n_XY content, and single-excitation density matrices have purely even n_XY content, so no SE state can couple to odd-parity modes. This is an algebraic selection rule, not a statistical correlation.

The subtlety, flagged in [Engineering Blueprint Rule 1](../publications/ENGINEERING_BLUEPRINT.md) (April 16 note): F61 constrains **Liouvillian time evolution** within a fixed parity sector. It does **not** constrain **initial-state preparation** that mixes excitation sectors. |GHZ_3⟩ lives in the {w = 0, w = 3} sectors; |W_3⟩ in the w = 1 sector. The superposition α|GHZ_3⟩ + β|W_3⟩ is a legitimate pure state, and it produces a full density matrix ρ = |ψ⟩⟨ψ| that contains cross-sector coherences |GHZ⟩⟨W| and |W⟩⟨GHZ|.

### The mechanism: cross-terms in the pair reduction

When ρ is traced to a pair ρ_AB, the cross-sector coherences generate off-diagonal matrix elements that neither parent state produces alone. Concretely, GHZ contributes population on |00⟩⟨00| and |11⟩⟨11| pieces after tracing, W contributes population on |01⟩⟨10| + |10⟩⟨01|, and the cross-terms GHZ-W contribute |00⟩⟨01|-type entries (one basis ket from GHZ's w=0 support, one from W's w=1 support). In the pair Pauli decomposition, these |00⟩⟨01|-type cross-entries have n_XY = 1 (odd parity at the pair level).

These odd-parity cross-terms are the L1 contribution that neither parent has. GHZ alone gives L1 = 0 in the pair (its off-diagonals trace away). W alone gives L1 = 2/N (for N = 3: 2/3 ≈ 0.667). The mix gives L1 = √6·α·β + (2/3)β², which at the exact α_opt evaluates to 1.60251. Combined with purity C(α_opt) = 0.59983 (the GHZ-diagonal piece lifts purity above W's 5/9 ≈ 0.556), CΨ = C · L1 / 3 reaches 0.59983 · 1.60251 / 3 = 0.32041.

### What this does not say

F69 is not a protection result. It does not claim the optimum state lives longer under Z-dephasing than W or GHZ. Under dynamics, F61 still applies: the odd-n_XY cross-terms are in odd-parity Liouvillian eigenmodes and decay normally. F69 is a **starting condition** above the fold, not a long-lived resource. The sector mixing is allowed at t = 0 because preparation is unconstrained; the fold-crossing that follows is standard Lindblad evolution.

CΨ = C · Ψ is a product. F60 starves the L1 side, F62 starves the purity side. Neither parent has both factors strong simultaneously. F69 arranges both through cross-sector interference in the preparation step.

## Scope: N = 3 is special

Running the same optimization inside α|GHZ_N⟩ + β|W_N⟩ at larger N gives no crossing of 1/4:

| N | α_opt | min pair-CΨ(0) | ratio to 1/4 | above fold? |
|---|-------|----------------|--------------|-------------|
| 3 | 0.6150 | 0.320406 | 1.282× | YES |
| 4 | 0.5750 | 0.166665 | 0.667× | no |
| 5 | 0.5800 | 0.146381 | 0.586× | no |
| 6 | 0.5650 | 0.134104 | 0.536× | no |

(201-point grid sweep, run_ghz_w_optimum_n3.py `family_fails_for_larger_N`.)

Why does N = 3 work and N ≥ 4 fail? The purity contribution from |GHZ_N⟩'s pair reduction scales as 1/(2^N - 1) (F60), vanishing exponentially. The L1 contribution from |W_N⟩'s pair reduction scales as 2/N (F62). At N = 3 both are close enough to 1/4 that their product lifts the pair above the fold; at N ≥ 4 the GHZ side has collapsed too far.

**The question "is there any other family that works at N ≥ 4?" is open.** F69 characterizes a specific two-parameter family; the full landscape of permutation-symmetric pure states above the fold for N ≥ 4 is not mapped here. A next-step survey would replace {|GHZ_N⟩, |W_N⟩} with {Dicke(N, k)} for k ∈ {1, 2, ..., N-1} plus GHZ, enumerate linear superpositions, and scan for fold-crossings. At N = 4 the broader Dicke family was touched in `cpsi_birth_landscape.py` but no mix was optimized; all tested single states sit below the fold (max: Dicke(4, 2) at CΨ = 0.111).

---

## Hardware relevance

Under Kingston-grade Z-dephasing (T₂ ∈ [240, 320] μs, γ ∈ [3.1, 4.2] × 10⁻³ μs⁻¹), the optimum state starts above the fold and crosses CΨ = 1/4 monotonically at t* ≈ 11.2 μs (data from `cpsi_sector_mix_optimization.py` Part 4). This is short compared to the bonding-mode T₂ ~ (N+1)³/(4π²γ) ≈ 2 ms at N = 3. F69 is not a protection result; it is an *above-fold starting condition* result.

**The significance for a hardware experiment** is that pair tomography at t = 0 *alone* distinguishes the optimum state from GHZ_3 and W_3: F69 sits at pair-CΨ(0) ≈ 0.320, GHZ_3 at 0, W_3 at 0.123. The three points in pair-CΨ space are machine-resolvable in any reasonable tomography (few-percent statistical noise). No delay, no timing, no T₂ measurement needed. A minimal validation circuit:

- Prepare α |GHZ_3⟩ + β |W_3⟩ with α = 0.6127 (non-Clifford: requires one amplitude-encoding rotation with irrational angle θ = arctan(β/α) ≈ 0.9118 rad).
- 2-qubit state tomography on each of the three pairs (Q0-Q1, Q0-Q2, Q1-Q2).
- Compute pair-CΨ. Three values should agree and sit near 0.32 (permutation symmetry is a hardware-testable consistency check on top of the fold test).

The same circuit can be extended to track CΨ(t) at a few delays to confirm the monotonic crossing, but the minimal claim (F69 is above the fold while GHZ and W are not) needs only t = 0.

## Pitfall: the 3-state spherical scan artifact

`cpsi_sector_mix_optimization.py` Part 2 runs a more aggressive scan over the 3-state family

    |ψ(θ, φ)⟩ = sin(θ) cos(φ) |GHZ_3⟩ + sin(θ) sin(φ) |W_3⟩ + cos(θ) |W̄_3⟩

using a spherical parametrization with all-real coefficients. The scan reports a peak at approximately (θ, φ) = (0.9121, 0.8856) with min pair-CΨ(0) ≈ 0.9998, which looks like a dramatic improvement over the 0.3204 of the two-state optimum.

**This peak is not a new finding.** It is the computational-basis uniform superposition |+⟩³ in disguise. The 3-state family {|GHZ_3⟩, |W_3⟩, |W̄_3⟩} spans a 3-dimensional subspace of ℂ⁸ that happens to contain |+⟩³ as a specific linear combination. The spherical scan picks it up because |+⟩³ has pair-CΨ = 1 trivially: every pair reduction is (|+⟩⟨+|)⊗(|+⟩⟨+|), a pure product of pure states.

`simulations/sector_mix_spherical_artifact.py` confirms this with three independent product-state tests at the reported optimum:

| test | result | interpretation |
|------|--------|----------------|
| single-qubit purity Tr(ρ_q²) | 1.000000 for every q | ρ factorizes |
| CKW 3-tangle τ_ABC | 0.000000 | no three-body entanglement |
| overlap \|⟨+³\|ψ*⟩\|² | 1.000000 | ψ* = \|+⟩³ exactly |

The amplitudes of |ψ*⟩ in the computational basis are uniform: every component ≈ 0.3535 ≈ 1/√8.

**Why the optimizer finds it.** The pair-CΨ ceiling over all 2-qubit states is CΨ = 1 (pure product state of two |+⟩'s). The sector-mix family happens to contain one such point. A grid scan over (θ, φ) that includes this point reports it as the family maximum, correctly on the pure-math level but trivially on the physics level. Any product state is above the fold; that is the fold's least interesting feature.

**Guard.** Both `cpsi_sector_mix_optimization.py` (docstring + Part 2 function docstring) carry ARTIFACT NOTICE comments pointing here. If you hit this peak again, run `sector_mix_spherical_artifact.py` first. It takes 30 seconds and gives you all three confirming tests.

---

## Open follow-ups

**Broader sector-mix landscape at N ≥ 4.** F69 is N = 3 only. A systematic scan of all permutation-symmetric linear combinations of {GHZ_N, Dicke(N, 1), ..., Dicke(N, N-1), W̄_N} across N = 4, 5 would close the question of whether any symmetric family escapes the fold at larger N. The product-state artifact lesson from this session applies: explicitly exclude product states before declaring a winner.

**Hardware run of F69.** The 3-point tomography test (GHZ_3, W_3, F69 optimum) is a minimal Kingston experiment that validates the fold-crossing picture on real hardware without any time-dependent data. Estimated QPU cost: a single state-prep circuit × three preparations × nine tomography settings × 1000 shots = 27,000 shots, well under any free-tier budget. Not urgent because the simulation is exact; would be nice as a companion to the cusp-slowing run for completeness.

**Closed-form optimum.** The minimal polynomial of α²_opt is a sextic irreducible over ℚ. It is possible that α²_opt has a clean expression in radicals of a Galois-resolvent (a sextic with a solvable Galois group does admit radical expressions; the Galois group of P has not been computed here). Checking `sympy.polys.numberfields.galois_group(P)` would settle this. Not research, just a mathematical curiosity.

---

## Reproduction

```
cd simulations
python ghz_w_optimum_n3.py                # symbolic + numerical optimum, sextic, irreducibility
python cpsi_sector_mix_optimization.py    # original sweep + Kingston dynamics
python cpsi_birth_landscape.py            # broad state survey
python sector_mix_spherical_artifact.py   # product-state diagnosis of the spherical-scan peak
```

All four scripts are pure Python (numpy, scipy, sympy), no external dependencies beyond the standard scientific stack. Total runtime < 15 s on modern hardware. Outputs land in `simulations/results/`.
