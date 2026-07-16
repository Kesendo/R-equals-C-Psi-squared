# The mirror's order-sorting law: one principle, two theorems

**Date:** 2026-07-16
**Status:** assembly proof. One new identity family; every cell witness is an already-proven Tier-1 result or the committed gate. Gate: [`simulations/mirror_order_sorting.py`](../../simulations/mirror_order_sorting.py) (N = 4, 5 default, N = 6 via `--deep`; all PASS, exit 0). The antiunitary column's own gate is [`simulations/zeta2_anti_protection.py`](../../simulations/zeta2_anti_protection.py).
**Registry:** minted **F131** (2026-07-16, Tom's call after the review rounds); entry at the end of [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md).

Before the algebra, the shape: a mirror does not throw information away, it
sorts it. What the noise does to a mirror-respecting system arrives in two
channels, and the mirror decides by pure parity bookkeeping which orders of
the disturbance may enter which channel. Protection (odd orders cancel) and
anti-protection (even orders double) are the two visible corners of one
four-cell table; [PROOF_ZETA2_ANTI_PROTECTION](PROOF_ZETA2_ANTI_PROTECTION.md)
said it for its one pair: "Protection and anti-protection are not two
symmetries. They are one symmetry read at odd and at even order." This
document is that sentence made general.

## 1. The principle

Let M be a mirror (a unitary or antiunitary involution) and G(x) a generator
family (Hamiltonian, Liouvillian, or Floquet step) over parameters x. Suppose
mirror conjugation reflects a scan: for a base point x₀ and a direction δ,

M · G(x₀ + s·δ) · M⁻¹ = G(x₀ + σ_eff·s·δ),

where **σ_eff is the sign mirror-conjugation puts on the scan parameter**:
σ_eff = σ_op · χ_M, with σ_op the operator parity of the direction δ under M
and χ_M = +1 for a linear (unitary) mirror, χ_M = −1 for an antilinear
(antiunitary) one (the antiunitary flips the i in the exponent; that is how a
mirror-FIXED perturbation can still scan with σ_eff = −1). For a readout of
definite mirror parity q, the response orders sort by the product q·σ_eff:

| | σ_eff = +1 | σ_eff = −1 |
|---|---|---|
| **q = +1** | generic (no constraint) | EVEN response: only even orders enter |
| **q = −1** | IDENTICALLY ZERO (all t, all times) | ODD response: only odd orders enter |

Two honesty notes on the table itself. The σ_eff = +1 column is a static
selection rule (the generator does not move along the scan as seen by the
mirror), not a response-order statement. And the law asserts PARITIES only,
never magnitudes: which orders are absent is a theorem; how large the
surviving orders are is instance physics (F91's κ is Tier-2 empirical, the
ζ² factor 2 is Floquet-pair-specific).

## 2. Theorem A: the unitary column (unconditional)

**Setting.** The open N-site XX chain with local Z dephasing; R = the F71 site
reversal; parameters J (bonds), γ (dephasing profile), h (longitudinal field)
split into R-even base + t times an R-odd direction (σ_op = −1, χ_R = +1, so
σ_eff = −1). The conjugation identity is entry-wise and exact,

(R⊗R) · L(base + t·dir) · (R⊗R) = L(base − t·dir),

(gate T0: residual exactly 0.0 for J-scans and γ-scans; h-scans by the same
one-line algebra; the general vectorization condition is U⊗Ū, and kron(R,R)
is correct here because R is a real permutation). Fence: the identity reflects
ALL scanned parameters at once; a mixed J/γ scan direction is R-odd only if
every component negates.

**Theorem A.** For an OPERATOR-R-even preparation (RρR = ρ; expectation-even
is provably insufficient, see §4) and a readout O with ROR = qO:

⟨O⟩(t; time) = q · ⟨O⟩(−t; time), for every evolution time.

*Proof.* Write 𝓡(X) = RXR, so 𝓡² = id and e^{L(−t)·T} = 𝓡 e^{L(t)·T} 𝓡.
Then ⟨O⟩(−t) = Tr[O · 𝓡(e^{L(t)T}(𝓡 ρ))] = Tr[(ROR) · e^{L(t)T} ρ] = q·⟨O⟩(t). ∎

**The cells, with their witnesses:**

- **(−,−) odd response:** the bond-mirror deviation D(b) = c₁(b) − c₁(N−2−b)
  is exactly odd in the anti-palindromic direction. Owned:
  [PROOF_F100](PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md) (J axis) and F101
  (γ axis), machine-zero even-power fits. Gate T2 re-pins the cell on plain
  ⟨Z⟩ readouts (odd to 4·10⁻¹⁶, and non-vacuous, |⟨O_odd⟩| up to 0.72).
- **(+,−) even response, spectral face:** the FULL Liouville spectrum is an
  even function of t (a one-line corollary: L(t) and L(−t) are conjugate, so
  their spectra coincide as multisets; gate T1, optimal-assignment matching
  4·10⁻¹³ worst; a naive complex sort mispairs the near-degenerate −2γ groups
  and must not be used). The owned F91/F92/F93 statements (diagonal-block
  spectra depend only on pair-sums, [PROOF_F92](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)
  "δλ = 0 to first order") are the STRONGER all-orders invariance of a
  sub-object; the full spectrum is merely even. Spectral evenness and
  trajectory parity are two faces of the one identity, not one object.
- **(+,−) even response, trajectory face:** gate T2, ⟨O_even⟩(t) − ⟨O_even⟩(−t)
  at 5·10⁻¹⁶ across scan scales and evolution times.
- **(−,+) the zero cell:** R-even generator, operator-R-even preparation,
  R-odd readout: ⟨O⟩ ≡ 0 for all times (gate T3, 3·10⁻¹⁶ at N = 4, 5). The
  witness for this cell is the gate itself; see §5 for what this cell is NOT.
- **(+,+):** no claim. Generic response. And the reader's converse fence: no
  mirror, no theorem. A pair that is not a mirror pair (the n = 12 collision
  pair, the A1/A2 arms) has no definite-parity readout in this sense and sits
  outside the table entirely; expect generic (first-order) response there.

## 3. Theorem B: the antiunitary column (the ζ² law)

The same principle entered through an antiunitary mirror is
[PROOF_ZETA2_ANTI_PROTECTION](PROOF_ZETA2_ANTI_PROTECTION.md): Θ = T·K (class
BDI, owner `ChiralKClaim`) commutes with the brickwork Floquet step and FIXES
the occupation-diagonal ZZ perturbation (σ_op = +1), but antiunitarity flips
the exponent's i (χ_Θ = −1), so σ_eff = −1 and the scan reflects:
Θ·W(ζ)·Θ⁻¹ = W(−ζ). For the Θ-mirror branch pair, the pair DIFFERENCE is the
q = +1 readout of the table (even response: the ζ² law, exact factor 2) and
the pair SUM is the q = −1 readout (odd response: common drift, invisible to
the fringe).

**The extra hypotheses are real and belong to this column only:** the
branch-resolved identity θ_τ(ζ) = −θ_ν(−ζ) needs the branch to be a simple,
isolated eigenphase tracked continuously (no crossing relabels it); the
MULTISET version of the eigenphase identity is unconditional, exactly as in
Theorem A. One principle, two theorems: Theorem A never tracks anything;
Theorem B pays tracking hypotheses to speak about a single pair.

## 4. The boundary: the hypotheses are the physics

- **Operator-even preparation.** The proof consumes RρR = ρ as an operator
  identity. A state that is merely expectation-even (⟨O_odd⟩(0) = 0 with a
  nonzero R-odd operator component) leaks: the R-odd component evolves inside
  the R-odd subspace and reads out against an R-odd observable.
- **The leak is exactly affine.** For ρ = ρ_even + ε·ρ_odd the even-cell
  readout splits as A(t) + ε·B(t) with A even and B odd in t, EXACTLY (the
  master equation is linear in ρ; there is no higher-order ε tail). The
  forbidden channel opens at leading order O(ε·t). Gate T4: ratio 2.000000 on
  halving ε, slope drift 0.1 % on halving t, coefficient of the same order
  across N with a mild drift (2.5/2.4/2.3·10⁻² at N = 4/5/6 in the test
  frame: local end-site physics).
- **The same shape on the antiunitary side (a rhyme, not an identification):**
  the flown b_qs budget (quasi-static disorder pulling the A0 slope through
  prep-basis mismatch, [IBM_F129_RAMSEY_FRINGE](../../experiments/IBM_F129_RAMSEY_FRINGE.md)
  §5) is the same STRUCTURAL failure, a preparation without definite mirror
  parity opening the forbidden linear order. It lives on the Θ side with its
  own mechanism; we do not equate the two instances, we note the recurring
  shape: **in both linear leaks we have examined (the ε-leak here, the b_qs
  budget on the Θ side) the failure was a broken parity hypothesis, not a
  broken mirror.**

## 5. What the zero cell is NOT: Π-protected observables

The repo's Π-protected observables (`observables.py`,
`PiProtectedObservables.cs`, hardware entry `pi_protected_xiz_yzzy`) are a
DIFFERENT zero mechanism and are deliberately NOT cited as the (−,+) cell's
witness: their zeros come from degenerate-cluster cancellations in the
eigenmode expansion (the code's own docstring calls this strictly weaker than
per-mode vanishing), they are computed per initial state, and the flagship
⟨X₀IZ₂⟩ does not even have definite R-parity (site reversal maps it to
Z₀IX₂). A cluster cancellation and a mirror selection rule can both print
zero; they are not the same object. The cell's witness is gate T3, a pure
selection rule.

## 6. Hardware, honestly sized

Confirmation 24 (ibm_kingston, the standing fringe) is a hardware sighting of
the FIRST-order protection: the odd orders cancel in the pair difference and
the fringe stands. The even-order anti-protection term stayed inside its
pre-registered budget (b_zz2 ≈ 38 % of the observed A0 excess, the remainder
under the unknown-sign b_qs draw) and is NOT independently measured. No cell
of this table beyond the odd-order cancellation has a hardware witness yet;
the (+,−) coefficient measurement (a ζ²-meter on a mirror pair) is a
candidate future instrument, not a result.

## 7. Where the pieces live

- This proof + gate [`mirror_order_sorting.py`](../../simulations/mirror_order_sorting.py)
  (T0 identity + γ-corollary, T1 spectral evenness, T2 both response cells,
  T3 the zero cell, T4 the leak; N = 4, 5, `--deep` 6).
- Theorem B: [PROOF_ZETA2_ANTI_PROTECTION](PROOF_ZETA2_ANTI_PROTECTION.md) +
  its gate; owner of Θ: `ChiralKClaim` + [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md).
- The owned cell witnesses: [PROOF_F92](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)
  (+ F91/F93), [PROOF_F100](PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md) (+ F101).
- The engines: F118 ([PROOF_PI_FACTORS_AS_R_TIMES_D](PROOF_PI_FACTORS_AS_R_TIMES_D.md)),
  F119 ([PROOF_ANTILINEAR_TRIANGLE](PROOF_ANTILINEAR_TRIANGLE.md)).
- MirrorWorld adoption (2026-07-16): `compute/MirrorWorld/OrderSorting.cs`, run
  mode `sorting N` -- Theorem A's trajectory face by twin RK4 on the γ axis, the
  pencil face as `ParameterKlein.MirrorConjugationResidual`; guarded from below
  by `OrderSortingTests` (6). Theorem B and the spectral evenness stay in the
  main repo (the drawn boundary, same style as `mirror`/`seed`).
- Typed surface (2026-07-16): `MirrorOrderSortingClaim`
  (`compute/RCPsiSquared.Core/Symmetry/MirrorOrderSortingClaim.cs`), parents
  `ChiralKClaim` + `AntilinearTriangleClaim` + the F91 family (F71/F92/F93);
  self-check battery at N = 3 (moment-level, ten cases), guarded by
  `MirrorOrderSortingClaimTests` (Core.Tests) and the registration + wiring
  audit (Runtime.Tests).
