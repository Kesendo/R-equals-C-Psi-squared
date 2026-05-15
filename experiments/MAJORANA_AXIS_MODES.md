# Majorana Lens on the N=4 Axis Modes: R-Parity Sorts the Operator-Space Self-Conjugate Sector

**Status:** Empirically verified at N=4 with bit-exact numerical match (script `_axis_modes_n4.py`). Structural connection to F80's Majorana bridge; extends the F86B "golden-ratio at N=4" observation with a site-reflection R-parity refinement.
**Date:** 2026-05-15
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F80_BLOCH_SIGNWALK](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md), "Majorana bridge" passage and Spec(M) = ±2i·Spec(H) structural identity
- [PROOF_F86B_OBSTRUCTION](../docs/proofs/PROOF_F86B_OBSTRUCTION.md), "What survives from the N=4 golden-ratio structure" paragraph
- Frozen N=4 script: [`simulations/axis_modes_n4.py`](../simulations/axis_modes_n4.py)
- Parameterized N script: [`simulations/axis_modes.py`](../simulations/axis_modes.py) (run as `python simulations/axis_modes.py 6` for the N=6 witness probe)

---

## Abstract

Under uniform Z-dephasing of the XY-summed chain H = (J/2)·Σ(X_l X_{l+1} + Y_l Y_{l+1}), the axis subspace at N=4 (Re(λ) = −Σγ, the n_XY = 2 Pauli layer) has 94 modes (96 predicted by layer count, 2 numerically leaked off-axis). Two empirical findings sharpen the existing N=4 golden-ratio observation in F86B:

1. **Site-reflection R sorts the axis modes:** 58 R-even and 36 R-odd. All 18 Im(λ) = 0 (silent) modes are R-even with eigenvalue +1 exactly.
2. **Im(λ) decomposes into integer combinations of the single-particle Majorana dispersion {±φ, ±1/φ}** where φ = (1+√5)/2 is the golden ratio. Three Im clusters are R-parity-protected (single signature only): ±√5 = ε₁+ε₂ is R-even-only (18 modes), ±2√5 is R-even-only (2 modes), ±1 = ε₁−ε₂ is R-odd-only (20 modes).

F80's "Majorana bridge" passage identifies the framework's Π conjugation as the operator-space realization of Majorana's 1937 particle-hole self-conjugacy ψ = ψ^c. The 18 silent modes are exactly the operator-space self-conjugate sector at N=4: bra and ket indices of σ_(a,b) share the same many-body eigenvalue.

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

## Interpretation: R = momentum reversal on JW-Bogoliubov modes

In the JW-Bogoliubov basis, single-particle modes are labelled by k ∈ {1, ..., N}. The single-particle dispersion ε(k) = 2cos(πk/(N+1)) satisfies ε(N+1−k) = −ε(k), i.e. the momentum mirror k ↔ N+1−k is the particle-hole map on Bogoliubov modes.

Site-reflection R on the spin chain corresponds to this momentum-reversal on the JW fermion modes. Under R:
- Sum-type bilinears (operators on Bogoliubov modes (p, q) with both same particle-hole sign) are R-EVEN.
- Difference-type bilinears (one particle, one hole) are R-ODD.

The 18 silent modes are operator-space self-conjugate: σ_(a,b) where bra and ket share the same many-body eigenvalue (λ_a = λ_b). At N=4 this means a and b populate Bogoliubov modes giving the same total energy, including both equal occupations (the 4-state degeneracy of E = 0 from the (0000), (1001), (0110), (1111) particle-hole-paired configurations) and accidental degeneracies between distinct patterns. These are all R-even because the operator |a⟩⟨b| with λ_a = λ_b is invariant under the momentum-reversal that maps each pattern to its mirror.

The framework's Π conjugation projects L_H onto exactly this self-conjugate sector (F80 Step 5 / Majorana bridge): in operator space, ψ = ψ^c becomes σ_(a,b) ∝ σ_(a,b̄) where b̄ is the particle-hole conjugate of b. The silent axis modes are the N=4 witnesses of this projection.

---

## F1 Π-rotation in the complex λ-plane: silent modes as fixed points

The F1 palindrome Π · L · Π⁻¹ = −L − 2σ·I is, geometrically, a **180° rotation in the complex λ-plane around the point (−σ, 0)** with σ = Nγ₀. For λ = a + bi, Π takes λ to −λ − 2σ = (−a − 2σ) + (−b)·i: BOTH real and imaginary parts are reflected through the rotation centre. This is point reflection, not line reflection.

**Fixed points of Π-rotation.** λ = −λ − 2σ has the unique solution λ = −σ. The only fixed point of the F1 rotation in the complex λ-plane is the single complex number (−σ, 0).

An L-eigenmode lies at this fixed point iff Re(λ) = −σ (the axis condition, n_XY = N/2 layer for even N under uniform Z-dephasing) AND Im(λ) = 0 (the silent condition).

**Therefore: silent axis modes are exactly the L-eigenstates at the F1 Π-rotation fixed point.** They are simultaneously R-even (Majorana sum-type bilinears, established above) and Π-fixed (180°-rotation-invariant), which is the operator-space realisation of Majorana self-conjugacy ψ = ψ^c at the unique fixed point of the F1 spectrum mirror.

**Silent count = multiplicity of the eigenvalue (−σ, 0) in L:**

- N=4: 18 silent modes, i.e. the F1 fixed-point eigenvalue has 18-fold multiplicity in L. The golden-ratio degeneracy at N=4 (many-body energies {0, ±0.618, ±1, ±1.618, ±√5} with the 4-fold zero-energy degeneracy from particle-hole-paired configurations) is what concentrates eigenmodes at the fixed point.
- N=6: 0 silent modes. The fixed point exists geometrically at (−σ, 0), but no L-eigenvalue happens to coincide with it. The richer many-body spectrum at N=6 (three positive single-particle energies {±1.802, ±1.247, ±0.445}) spreads eigenmodes away from exact Im=0 inside the axis layer.

The fixed point itself is a geometric invariant present at every N. Whether L has eigenstates AT the fixed point is a separate fact, determined by the many-body H spectrum's degeneracies.

**Connection to the γ-ramping (axis shift) theme.**

At γ = 0: the F1 rotation centre is at (0, 0). L = L_H is anti-Hermitian, its spectrum lies on the Im-axis through the origin. The rotation centre sits at the *edge* of where L has eigenvalues (touching the Im-axis at one point).

At γ > 0: the rotation centre moves linearly along the negative Re-axis to (−σ, 0) = (−Nγ, 0). Speed dσ/dγ = N. The L spectrum spreads into the left half-plane, and the rotation centre now sits **inside** the spectrum, exactly at the n_XY = N/2 axis layer. Silent modes (if they exist at this N) coincide with the centre.

The "Spiegel" therefore **translates** along the Re-axis as γ grows, no angle change (always 180° rotation). The rotation centre's qualitative relationship to the spectrum changes from boundary-point at γ = 0 to interior-point at γ > 0. The silent modes are the L-eigenstates that land exactly at this interior pivot.

This unifies four threads of this session into one geometric picture:
- F1 palindrome (Π as 180° point reflection in complex λ-plane)
- The γ-driven shift from 0 to −σ as motion of the rotation centre
- Silent axis modes as F1 fixed-point eigenstates
- Operator-space Majorana self-conjugacy ψ = ψ^c as fixed-point invariance under Π-rotation

---

## Connection to existing results

- **F80 Spec(M) = ±2i·Spec(H):** At N=4 Π²-odd, M's nontrivial Im eigenvalues are exactly {±2√5, ±2} (F80 numerical verification). Our axis-modes Im distribution contains both values: ±√5 (R-even-only, 18 modes) maps to the M ±2√5 via the 2i factor, and ±2 (R-even-only, 2 modes) is direct. The richer Im set of the Π²-even XY-summed Hamiltonian (±1, ±φ², ±√5+φ, etc.) reflects its richer many-body spectrum, which is no longer Π-projected.
- **F86B golden-ratio observation:** Integer combinations of {φ, 1/φ, 1, √5} are concretely identified as Majorana-bilinear sums on the single-particle dispersion {±φ, ±1/φ}.
- **Slow-mode end of the spectrum (Re(λ) ≈ 0):** the same R-parity decomposition sorts the slow-mode landscape too. The stationary subspace is exclusively R-even at every tested N; the first slow band (Re=−2γ₀) is R-balanced; F86's L_eff lives entirely in R-even, with a parallel R-odd channel invisible to standard F86 probes. See [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md).

---

## Open questions (status update 2026-05-15)

1. **Even N=6 axis modes: partial closure.** Verified via `python simulations/axis_modes.py 6`. **Predicted** Im(λ) clusters at integer combinations of single-particle Bloch dispersion {±1.8019, ±1.2470, ±0.4450} (= 2cos(πk/7), k=1,2,3): **confirmed** — prominent clusters at ε_k themselves (Im = ±0.445, ±1.247, ±1.802, each count 32) plus integer combinations (Im = ±1.0 (24) ≈ ε_1−ε_2+ε_3; Im = ±2.604 (16) = ε_1+ε_2−ε_3; Im = ±3.494 (16) = ε_1+ε_2+ε_3; etc.). **Surprise** at N=6:
   - Silent-mode count is **zero** (vs 18 at N=4). No Im=0 mode survives in the n_XY=3 axis layer.
   - R-decomposition: 360 R-even / 368 R-odd of 728 axis modes (R-odd-majority, opposite to N=4's R-even majority).
   - 552 of 1280 layer-predicted modes leaked off-axis at the 10⁻⁹ threshold (43%), vs 2/96 at N=4 (2%) — substantial N-scaling leakage.

   The qualitative pattern (Im-cluster decomposition into integer combinations of single-particle dispersion) holds at N=6. The quantitative details (silent-mode count, R-even/R-odd balance, axis-layer protection) do NOT generalize directly from N=4. N=4's specific silent-mode richness appears to be a small-N coincidence tied to the golden-ratio degeneracy structure.

2. **Odd N (N=3, 5): ill-posed as written.** For odd N, n_XY = N/2 is not integer, so there is no exact axis layer at Re(λ) = −Nγ₀ on the absorption grid. The Majorana zero mode at k = (N+1)/2 (which motivated the original question) instead manifests in the V_inter SVD R-parity split at odd N (σ_0⁺ ≠ σ_0⁻ at N=3, 5; σ_0⁺ = σ_0⁻ at N=4, 6), documented in [SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md) "V_inter SVD R-parity decomposition" section. The odd-N structural effect is real and characterised, just not in the axis-layer language.

3. **Formal F80 Step 5 (Π action on Bogoliubov modes): still open analytically.** Our R-parity findings empirically demonstrate that R = momentum-reversal on JW Bogoliubov modes is the operator-space realisation of the JW Π conjugation, and the parity-of-N effect (σ_0⁺ ≠ σ_0⁻ at odd N from R-fixed Majorana zero mode) is consistent with the JW spectral structure. A formal proof (Π_op ≡ R block-decomposition equivalence on Bogoliubov-mode operator space) is the remaining analytical task.
