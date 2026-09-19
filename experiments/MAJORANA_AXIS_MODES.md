<!-- QUARTER-CURRENT -->
# N=4 axis-mode coefficient census under R parity

Current reading: the script reports a finite operator-basis coefficient census,
including `94+2` leakage rather than the predicted 96-row layer.  R-parity of an
operator-space subspace does not make each Liouvillian eigenvector individually
self-conjugate or Pi-fixed. A separate full-spectrum match has 14 F1-fixed eigenvalues
and 121 two-member F1 orbits; that orbit census is not the 18-entry silent bin.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the Majorana and self-conjugacy narrative below records
the prediction that the finite census tested and narrowed.

# Majorana lens on N=4 axis modes: a finite R-parity coefficient census

**Status:** Finite N=4 coefficient census; individual-mode self-conjugacy interpretation withdrawn
**Date:** 2026-05-15
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [the F80 Bloch sign-walk proof](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md), "Majorana bridge" passage and Spec(M) = ±2i·Spec(H) structural identity
- [the F86b obstruction proof](../docs/proofs/PROOF_F86B_OBSTRUCTION.md), "What survives from the N=4 golden-ratio structure" paragraph
- Frozen N=4 script: [`simulations/axis_modes_n4.py`](../simulations/axis_modes_n4.py)
- Parameterized N script: [`simulations/axis_modes.py`](../simulations/axis_modes.py) (run as `python simulations/axis_modes.py 6` for the N=6 witness probe)

---

## Abstract

Under uniform Z-dephasing of the XY-summed chain H = (J/2)·Σ(X_l X_{l+1} + Y_l Y_{l+1}), the axis subspace at N=4 (Re(λ) = −Σγ, the n_XY = 2 Pauli layer) has 94 modes (96 predicted by layer count, 2 numerically leaked off-axis). Two empirical findings sharpen the existing N=4 golden-ratio observation in F86B:

1. **Site-reflection R sorts the axis modes:** 58 R-even and 36 R-odd. All 18 Im(λ) = 0 (silent) modes are R-even with eigenvalue +1 exactly.
2. **Im(λ) decomposes into integer combinations of the single-particle Majorana dispersion {±φ, ±1/φ}** where φ = (1+√5)/2 is the golden ratio. Three Im clusters are R-parity-protected (single signature only): ±√5 = ε₁+ε₂ is R-even-only (18 modes), ±2√5 is R-even-only (2 modes), ±1 = ε₁−ε₂ is R-odd-only (20 modes).

The earlier F80 "Majorana bridge" passage suggested reading the framework's
Π conjugation through Majorana particle-hole language. The 18 silent entries
reported here are centered eigenvalue occurrences in one numerical census. The
calculation does not establish that their individual eigenvectors are
self-conjugate or fixed by Π.

---

## Background

### F80's Majorana bridge

The Π²-odd chain bilinear under Jordan-Wigner becomes a pure Majorana bilinear in the γ' = i(c†−c) modes with single-particle spectrum ε(k) = 2J·cos(πk/(N+1)), paired as ±ε(k). The F80 "Majorana bridge" passage notes: *"The framework's Π conjugation is, in disguise, the operator-space realization of Majorana's particle-hole self-conjugacy. He had it right; we just have a richer vocabulary to express it now."*

This experiment applies the same lens to the Π²-EVEN XY-summed Hamiltonian H = (J/2)·Σ(XX+YY). Its JW transformation yields the standard Dirac-fermion tight-binding model with the SAME single-particle dispersion ε(k). At N=4: ε(k) = {2cos(π/5), 2cos(2π/5), 2cos(3π/5), 2cos(4π/5)} = {φ, 1/φ, −1/φ, −φ}. The many-body spectrum is richer than F80's Π²-odd case (9 distinct eigenvalues with multiplicities {4, 2, 2, 2, 2, 1, 1, 1, 1} on the 16-dim Hilbert space), reflecting all 2⁴ = 16 fermion-occupation patterns rather than the constrained Π²-odd subspace.

### F86B's golden-ratio observation

The PROOF_F86B_OBSTRUCTION "What survives from the N=4 golden-ratio structure" paragraph notes that the N=4 OBC dispersion is {±φ, ±1/φ} and that the Liouvillian Im(λ) values in the n_XY=2 sector decompose into integer combinations of {φ, 1/φ, 1, √5}, but this N=4 special case does not propagate to a closed form for Q_peak at general (c, N). This experiment adds the R-parity sorting and the explicit Majorana cluster identification to that observation.

---

## Empirical findings at N=4 (J = 1.0, γ = 0.05; reproduce with `python simulations/axis_modes_n4.py`)

### Site-reflection R on the axis subspace

R reverses qubit ordering on basis states; the operator-space action is R_op = R_h ⊗ R_h.

| | count |
|---|---|
| R-even axis modes | 58 |
| R-odd axis modes | 36 |
| total | 94 (of 96 layer-predicted; 2 leaked off-axis) |

All 18 silent (Im = 0) modes are R-even with eigenvalue +1 exactly (zero R-odd, zero ambiguous).

### Im(λ) by R-parity

| \|Im(λ)\| | exact identification | R-even | R-odd |
|---:|---|---:|---:|
| 0 | 0 (silent) | 18 | 0 |
| 0.386 | 1 − 1/φ = (3−√5)/2 | 4 | 4 |
| 1.000 | ε₁ − ε₂ = φ − 1/φ | 0 | 20 |
| ≈2.003 | 2 | 2 | 0 |
| ≈2.236 | √5 = ε₁ + ε₂ | 18 | 0 |
| ≈2.617 | φ² = φ + 1 | 4 | 4 |
| ≈2.849 | √5 + 1/φ | 4 | 4 |
| ≈3.852 | √5 + φ | 4 | 4 |
| ≈4.468 | 2√5 | 2 | 0 |

(One additional cluster at \|Im\| ≈ 2.224 with 2 R-even modes is a hybridized neighbour of the ±√5 line; small numerical drift from the eigendecomposition.)

Three Im clusters carry a single R-parity signature exclusively:
- **±√5 R-even-only (18 modes).** Sum-type pair bilinear ε₁ + ε₂.
- **±2√5 R-even-only (2 modes).** Doubled-pair type.
- **±1 R-odd-only (20 modes).** Difference-type bilinear ε₁ − ε₂.

The R-even-only clusters (18 + 2 = 20) plus the silent kernel (18) account for 38 of the 58 R-even modes; the R-odd-only cluster (20) accounts for 20 of the 36 R-odd modes.

---

## Candidate interpretation: R and momentum reversal

In the JW-Bogoliubov basis, single-particle modes are labelled by k ∈ {1, ..., N}. The single-particle dispersion ε(k) = 2cos(πk/(N+1)) satisfies ε(N+1−k) = −ε(k), i.e. the momentum mirror k ↔ N+1−k is the particle-hole map on Bogoliubov modes.

Site-reflection R on the spin chain corresponds to this momentum-reversal on the JW fermion modes. Under R:
- Sum-type bilinears (operators on Bogoliubov modes (p, q) with both same particle-hole sign) are R-EVEN.
- Difference-type bilinears (one particle, one hole) are R-ODD.

The 18 silent entries have equal-energy bra/ket labels in the chosen
reconstruction. Equal energy and an R-even coefficient do not by themselves
show that each operator |a⟩⟨b| is invariant under momentum reversal; a
degenerate eigenspace permits basis rotations. Thus "self-conjugate" is a
candidate lens, not a property established for the individual vectors.

F80 supplies an operator-level conjugation statement. This experiment has not
proved that it projects onto the 18 selected vectors or that
σ_(a,b)∝σ_(a,b̄) for each of them. A projector-level comparison across the
degenerate block remains the appropriate test.

---

## F1 rotation in the complex λ-plane: fixed eigenvalues versus fixed vectors

The F1 palindrome Π · L · Π⁻¹ = −L − 2σ·I is, geometrically, a **180° rotation in the complex λ-plane around the point (−σ, 0)** with σ = Nγ₀. For λ = a + bi, Π takes λ to −λ − 2σ = (−a − 2σ) + (−b)·i: BOTH real and imaginary parts are reflected through the rotation centre. This is point reflection, not line reflection.

**Fixed points of Π-rotation.** λ = −λ − 2σ has the unique solution λ = −σ. The only fixed point of the F1 rotation in the complex λ-plane is the single complex number (−σ, 0).

A centered eigenvalue satisfies Re(λ)=−σ and Im(λ)=0. That scalar
condition says only that the eigenvalue is fixed by the F1 rate map. It does not
say an individual eigenvector is fixed by the implementing operator.

**Current distinction:** the silent entries are eigenvalues at the scalar center
within this selected axis census. Individual R parity, Π invariance, and
Majorana self-conjugacy are separate vector claims and are not proved here.

**Silent count = multiplicity of the eigenvalue (−σ, 0) in L:**

- N=4: this axis selection reports 18 silent entries. The independent full
  F1 matching reports 14 fixed eigenvalues and 121 two-member orbits, so the
  18-entry basis census must not be identified with the fixed-orbit count.
- N=6: 0 silent modes. The fixed point exists geometrically at (−σ, 0), but no L-eigenvalue happens to coincide with it. The richer many-body spectrum at N=6 (three positive single-particle energies {±1.802, ±1.247, ±0.445}) spreads eigenmodes away from exact Im=0 inside the axis layer.

The scalar center is present at every N. Whether the spectrum contains that
eigenvalue, with what multiplicity, and how Π acts inside its eigenspace are
separate questions.

**Connection to the γ-ramping (axis shift) theme.**

At γ = 0: the F1 rotation centre is at (0, 0). L = L_H is anti-Hermitian, its spectrum lies on the Im-axis through the origin. The rotation centre sits at the *edge* of where L has eigenvalues (touching the Im-axis at one point).

At γ > 0 the F1 rate center moves linearly along the negative real axis to
(−σ,0)=(−Nγ,0). The n_XY=N/2 diagonal layer is centered there in the
uniform-rate book; a finite numerical eigenvalue selection near that line is
not automatically an invariant subspace.

The scalar center therefore translates along the real axis as γ grows. Calling
selected eigenvectors an interior pivot would require the missing eigenspace
action test.

The earlier notebook tried to unify four threads in one geometric picture:
- F1 palindrome (Π as 180° point reflection in complex λ-plane)
- The γ-driven shift from 0 to −σ as motion of the rotation centre
- centered eigenvalues in the finite axis census
- a proposed, still unproved Majorana/self-conjugacy interpretation

---

## Connection to existing results

- **F80 Spec(M) = ±2i·Spec(H):** At N=4 Π²-odd, M's nontrivial Im eigenvalues are exactly {±2√5, ±2} (F80 numerical verification). Our axis-modes Im distribution contains both values: ±√5 (R-even-only, 18 modes) maps to the M ±2√5 via the 2i factor, and ±2 (R-even-only, 2 modes) is direct. The richer Im set of the Π²-even XY-summed Hamiltonian (±1, ±φ², ±√5+φ, etc.) reflects its richer many-body spectrum, which is no longer Π-projected.
- **F86B golden-ratio observation:** Integer combinations of {φ, 1/φ, 1, √5} are concretely identified as Majorana-bilinear sums on the single-particle dispersion {±φ, ±1/φ}.
- **Slow-mode end of the spectrum (Re(λ) ≈ 0):** the same R-parity decomposition sorts the slow-mode landscape too. The stationary subspace is exclusively R-even at every tested N; the first slow band (Re=−2γ₀) is R-balanced; F86's L_eff lives entirely in R-even, with a parallel R-odd channel invisible to standard F86 probes. See [the slow-mode R-parity decomposition](SLOW_MODE_R_PARITY.md).

---

## Open questions (status update 2026-05-15)

1. **Even N=6 axis modes: partial closure.** Verified via `python simulations/axis_modes.py 6`. **Predicted** Im(λ) clusters at integer combinations of single-particle Bloch dispersion {±1.8019, ±1.2470, ±0.4450} (= 2cos(πk/7), k=1,2,3): **confirmed**, with prominent clusters at ε_k themselves (Im = ±0.445, ±1.247, ±1.802, each count 32) plus integer combinations (Im = ±1.0 (24) ≈ ε_1−ε_2+ε_3; Im = ±2.604 (16) = ε_1+ε_2−ε_3; Im = ±3.494 (16) = ε_1+ε_2+ε_3; etc.). **Surprise** at N=6:
   - Silent-mode count is **zero** (vs 18 at N=4). No Im=0 mode survives in the n_XY=3 axis layer.
   - R-decomposition: 360 R-even / 368 R-odd of 728 axis modes (R-odd-majority, opposite to N=4's R-even majority).
   - 552 of 1280 layer-predicted modes leaked off-axis at the 10⁻⁹ threshold (43%), vs 2/96 at N=4 (2%); substantial N-scaling leakage.

   The qualitative pattern (Im-cluster decomposition into integer combinations of single-particle dispersion) holds at N=6. The quantitative details (silent-mode count, R-even/R-odd balance, axis-layer protection) do NOT generalize directly from N=4. N=4's specific silent-mode richness appears to be a small-N coincidence tied to the golden-ratio degeneracy structure.

2. **Odd N (N=3, 5): ill-posed as written.** For odd N, n_XY = N/2 is not integer, so there is no exact axis layer at Re(λ) = −Nγ₀ on the absorption grid. The Majorana zero mode at k = (N+1)/2 (which motivated the original question) instead manifests in the V_inter SVD R-parity split at odd N (σ_0⁺ ≠ σ_0⁻ at N=3, 5; σ_0⁺ = σ_0⁻ at N=4, 6), documented in [Slow-Mode R-Parity](SLOW_MODE_R_PARITY.md) "V_inter SVD R-parity decomposition" section. The odd-N structural effect is real and characterised, just not in the axis-layer language.

3. **Formal F80 Step 5 (Π action on Bogoliubov modes): open.** The finite
   R-parity census is compatible with a momentum-reversal lens but does not
   demonstrate operator equivalence. A projector-level proof or counterexample
   across degenerate blocks remains the analytical task.
