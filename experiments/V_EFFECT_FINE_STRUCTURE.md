# V-Effect Fine Structure: Hard vs Soft Breaks via Framework

**Status:** Computational (Tier 1-2). Framework-based re-examination of the V-Effect's 14-of-36 result. Reproduces March 2026 finding via stricter operator-equation test, plus reveals a 19-case "soft break" intermediate category that the spectrum-pairing test missed.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `simulations/_veffect_36_combos_via_framework.py` (uses `framework.py`)
**See also:** [V_EFFECT_PALINDROME](V_EFFECT_PALINDROME.md), [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md), [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md)

---

## What this clarifies

The V-Effect (March 2026) found that of 36 two-term Pauli-pair Hamiltonians H = J(term1 + term2) at N=3 with two bonds, exactly **14 break** the palindromic structure and **22 do not**. The test used was eigenvalue pairing λ ↔ −λ − 2Σγ.

The framework's `palindrome_residual` function tests the stricter operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0. Operator equation implies eigenvalue pairing; the converse does not hold. Re-running all 36 combos with both criteria gives:

| Test | Broken | Unbroken |
|------|--------|----------|
| Operator equation (strict) | **33** | **3** |
| Eigenvalue pairing (V-Effect) | **14** | **22** |

The 22 V-Effect-unbroken cases decompose into:
- **3 truly unbroken**: {XX+YY, XX+ZZ, YY+ZZ}, both terms in the both-parity-even set.
- **19 soft-broken**: operator equation residual ~22-45 (substantial), but eigenvalue pairing accurate to ~10⁻¹⁴ (machine precision).

The 14 hard-broken cases break both criteria simultaneously. These are the V-Effect's original 14.

## Numerical results

Pipeline: `_veffect_36_combos_via_framework.py`. N=3, γ_l = 0.1 per site, two bonds (0,1), (1,2), threshold 10⁻¹⁰ for OP, 10⁻⁶ for SPEC.

**Truly unbroken (3):**

| Combo | OP residual | SPEC error | parities |
|-------|-------------|------------|----------|
| XX+YY | 7.8e-16 | 1.2e-14 | (ab)+(ab) |
| XX+ZZ | 7.8e-16 | 3.7e-14 | (ab)+(ab) |
| YY+ZZ | 7.8e-16 | 3.7e-14 | (ab)+(ab) |

All three have both terms in {XX, YY, ZZ}: the framework's both-parity-even selection. These are the Heisenberg/XXZ-form Hamiltonians (modulo identity).

**Soft-break sample (19 total):**

| Combo | OP residual | SPEC error | parities |
|-------|-------------|------------|----------|
| XX+XZ | 22.6 | 1.7e-14 | (ab)+(──) |
| XX+YZ | 32.0 | 1.3e-14 | (ab)+(─b) |
| XY+YX | 32.0 | 2.7e-14 | (a─)+(a─) |
| YZ+ZY | 45.3 | 9.3e-15 | (─b)+(─b) |
| ... | ... | ... | ... |

Operator equation broken at order 10¹, eigenvalue pairing intact at machine precision. The breaks in the operator sit in matrix elements that connect Π-paired sectors but happen not to shift the eigenvalue spectrum.

**Hard breaks (14):**

| Combo | OP residual | SPEC error | parities |
|-------|-------------|------------|----------|
| XX+XY | 22.6 | 0.198 | (ab)+(a─) |
| XX+YX | 22.6 | 0.200 | (ab)+(a─) |
| XY+YY | 22.6 | 0.200 | (a─)+(ab) |
| YX+YY | 22.6 | 0.200 | (a─)+(ab) |
| XY+XZ | 32.0 | 0.142 | (a─)+(──) |
| XY+ZX | 32.0 | 0.184 | (a─)+(──) |
| XZ+YX | 32.0 | 0.184 | (──)+(a─) |
| YX+ZX | 32.0 | 0.142 | (a─)+(──) |
| XY+YZ | 39.2 | 0.184 | (a─)+(─b) |
| XY+ZY | 39.2 | 0.142 | (a─)+(─b) |
| YX+YZ | 39.2 | 0.142 | (a─)+(─b) |
| YX+ZY | 39.2 | 0.184 | (a─)+(─b) |
| YZ+ZX | 39.2 | 0.133 | (─b)+(──) |
| XZ+ZY | 39.2 | 0.133 | (──)+(─b) |

These are the V-Effect's 14 broken cases.

## Structural pattern

Reading parity labels: 'a' = bit_a parity even, 'b' = bit_b parity even, '─' = parity violated. So `(a─)` means term has bit_a even but bit_b odd (a Pauli pair with one Y, like XY or YX).

**Hard-break signature:** every hard-break case contains at least one term of class `(a─)` (only bit_b violated, like XY or YX) OR the combination `(─b)+(──)` / `(──)+(─b)`. The bit_b violation propagates to spectral asymmetry without compensating bit_a structure to absorb it.

**Soft-break signature:** every soft-break case has either:
- One term in {XX, YY, ZZ} = `(ab)` (fully good) plus another that violates both parities `(──)` or only bit_a `(─b)`.
- OR matched parity violations between the two terms (`(a─)+(a─)`, `(─b)+(─b)`, etc.).
- OR mixed cases where the two violations cancel in the spectrum despite the operator equation showing residue.

The cleanest soft-break families are:
- `(ab)+(─b)` and `(ab)+(──)`: one good term provides palindromic structure, the other's violations stay sub-spectral.
- `(a─)+(a─)`: both terms violate only bit_b, identically. The matched violation gives matched eigenvalue shifts that pair.

## What the framework adds beyond the V-Effect

The V-Effect tested EIGENVALUE pairing as the criterion for "palindrome holds". This is the SPECTRAL consequence. The framework's OPERATOR-equation test is finer: it catches breaks in the matrix structure that don't propagate to eigenvalues.

The 19 soft-break cases are interesting: they have **non-zero** operator residual (10¹-magnitude) but **zero** eigenvalue pairing error (machine precision). The breaks are real, but they live entirely in the OFF-DIAGONAL matrix elements between Π-paired sectors, not in the diagonal eigenvalue structure.

For physical observation:
- Eigenvalue pairing is what shows up in spectroscopy (transition frequencies pair).
- Operator equation is what governs the full propagator e^(L·t) and operator-level relations.

The soft-break category corresponds to systems where spectroscopy looks symmetric but the dynamics has off-diagonal asymmetric mixing. These would only show up in non-spectral measurements (e.g., specific operator expectations, time-resolved cross-correlations).

## Cross-reference to PROOF_ZERO_IMMUNITY

[PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md) showed that the (w=0, w=N) extreme sectors of any 2-body Hamiltonian's Liouvillian preserve palindrome exactly. The break, when present, lives in boundary sectors 0 < w < N.

For the 33 OP-broken cases here, the break lives in boundary sectors at N=3 (w=1, w=2). Among these:
- 14 hard-breaks: the boundary break propagates to off-diagonal eigenvalue structure ⇒ eigenvalue pairing fails.
- 19 soft-breaks: the boundary break is structured such that off-diagonal eigenvalue parity-symmetry still holds ⇒ eigenvalue pairs intact.

The two categories share the boundary-residing break but differ in whether the boundary-mode coupling structure is "spectrally symmetric" or not.

## Hardware verification (2026-04-26)

The 3/19/14 distinction was tested live on `ibm_marrakesh` (Heron r2) via three representative Hamiltonians (XX+YY truly, XY+YX soft, XX+XY hard) on a 3-qubit chain at [48, 49, 50] from |+−+⟩ initial state, t=0.8, n_trotter=3, 4096 shots/basis. Pipeline: `D:\...\ibm_quantum_tomography\run_soft_break.py` (first script in that directory to import `framework.py`).

| Observable | Framework idealised | Aer w/ Marrakesh noise | **Hardware Marrakesh** |
|------------|---------------------|------------------------|--------------------------|
| ⟨X₀Z₂⟩ truly_unbroken (XX+YY) | 0.000 | -0.020 | **+0.011** |
| ⟨X₀Z₂⟩ soft_broken (XY+YX) | -0.623 | -0.660 | **-0.711** |
| ⟨X₀Z₂⟩ hard_broken (XX+XY) | +0.195 | +0.230 | **+0.205** |
| Δ(soft − truly) | -0.62 | -0.64 | **-0.72** |

All three categories separately resolved at SNR ~13-47σ. Job `d7mjnjjaq2pc73a1pk4g`. Raw JSON + breakdown in [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/README.md).

The hardware result is slightly STRONGER than the idealized framework predicts. T1 thermal relaxation and ZZ crosstalk compound the soft-break signal in operators not purely diagonal in L's eigenbasis. Same mechanism as the γ-profile amplification observed in CMRR_BREAK and GAMMA_AS_SIGNAL: hardware non-uniformity amplifies symmetry-breaking signatures.

## What this does and does not establish

**Establishes:**

- Framework primitives (palindrome_residual, lindbladian_z_dephasing, _build_bilinear) reproduce the V-Effect's 14/22 split when filtered through the SPEC criterion.
- Framework reveals a 3 / 19 / 14 fine structure that the V-Effect's two-category 22 / 14 partition obscured.
- The 3 truly-unbroken cases are predictable structurally: both terms in the both-parity-even set {XX, YY, ZZ}.
- The 3/19/14 distinction is operationally observable on real IBM Heron r2 hardware via specific 2-qubit Pauli expectations (above).

**Does not establish:**

- An analytical predictor for which combos are soft vs hard. The structural patterns (matched parities, one-good-term, etc.) are descriptive but not yet derived.
- Generalization to N > 3. Framework supports it; not tested here.
- The full 19 soft-broken cases on hardware. We tested 1 representative each from {3 truly, 19 soft, 14 hard}.

## Open questions

- Soft-break analytical structure: derive when off-diagonal Π-paired matrix elements give zero shift to eigenvalue pairing despite non-zero operator residual.
- Multi-bond V-Effect at N≥4: same combos extended to longer chains; does the soft/hard split depend on chain length?
- Hardware signature: is there an IBM tomography observable that distinguishes the 19 soft breaks from the 3 truly-unbroken? If yes, this would be the first hardware-level test of the framework's strict-vs-spectral distinction.

## References

- [V_EFFECT_PALINDROME](V_EFFECT_PALINDROME.md): the original 14-of-36 result, March 2026.
- [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md): structural location of the break (w=1, 2 sectors at N=3).
- [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md): analytical proof of (w=0, w=N) extreme-sector immunity.
- `simulations/framework.py`: framework primitives.
- `simulations/_veffect_36_combos_via_framework.py`: this calculation's pipeline.
