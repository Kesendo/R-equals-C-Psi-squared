# R = CΨ²  -- v3.0 Release Notes

## Hardware Validated

**Release date:** March 18, 2026
**DOI:** (will be assigned by Zenodo)
**Previous:** v2.0-proof (March 14, 2026, DOI: 10.5281/zenodo.19022139)

---

## What's new since v2.0

### 1. Hardware Validation: CΨ=1/4 crossing at 1.9% accuracy

The palindromic spectral symmetry makes a testable prediction:
the CΨ=1/4 boundary crossing occurs at a specific time determined
by T1 and T2*. We measured this on IBM Torino (Qubit 80):

- **Predicted:** t* = 15.01 μs (using same-day Ramsey T2*)
- **Measured:** t* = 15.29 μs
- **Deviation: 1.9%**

This confirms: the crossing equation is correct, T2* (not T2 echo)
is the relevant timescale, and the palindrome theory maps to real
quantum hardware with percent-level accuracy.

Key insight: T2* fluctuates 58% over 6 days on the same qubit.
Same-day measurement is essential. The theory is exact  -- the
hardware parameters are noisy.

### 2. N=8 Liouvillian: 54,118 rates, 100% palindromic

Extended verification from N=7 to N=8 using a custom C#/OpenBLAS
compute engine (65,536 x 65,536 matrix, 10.6 hours, 128 GB RAM).
All 54,118 oscillatory rates pair perfectly. The scaling table
now covers N=2 through N=8.

### 3. Non-Heisenberg palindrome: ALL standard models

The palindrome was originally proven for Heisenberg coupling.
Systematic testing reveals it holds for ALL standard condensed
matter models under single-axis dephasing:

- Heisenberg, XY-only, Ising, XXZ, DM, Heisenberg+DM
- Unequal coupling coefficients (J_XX ≠ J_YY ≠ J_ZZ)
- Every dephasing axis (Z, X, Y)
- Verified at N=3, 4, 5, 6

The palindrome is a property of the dephasing structure,
not the Hamiltonian.

### 4. Two Pi families + alternating operators (34/36 explained)

The conjugation operator Pi is not unique. Algebraic enumeration
of all valid per-site Pauli transformations reveals:

- **Family P1** (I↔X, Y↔Z): supports XX, YY, ZZ
- **Family P4** (I↔Y, X↔Z): supports XX, XZ, YY, ZX, ZZ
- **Alternating operators** (M1-M2-M1): required for XY, YX, DM

Compatibility matrix for all 36 two-term combinations:
17 uniform + 3 alternating + 2 unknown structure + 14 broken = 36.
The 14 broken cases match numerical results exactly.

### 5. XOR decomposition universal

The GHZ→100% XOR / W→0% XOR split is not Heisenberg-specific:

| Model | GHZ→XOR | W→XOR | Pauli correlation |
|-------|---------|-------|-------------------|
| Heisenberg | 100% | 0% | r = 0.984 |
| XY-only | 100% | 0% | r = 0.986 |
| Ising | 100% | 0% | r = 0.983 |
| DM | 100% | 0% | r = 0.999 |
| XXZ | 100% | 0% | r = 0.989 |

Holds for chain and star topologies, N=3 and N=4.
XOR mode count = steady state count for all models.

**Correction:** Bell+ is palindromic at N≥3 (not XOR as v2.0 stated).
Bell = GHZ only at N=2. At N≥3, Hamming distance is 2, not N.

### 6. Depolarizing noise quantified

Under depolarizing noise (X+Y+Z simultaneously), the palindrome
breaks with a systematic, Hamiltonian-independent error:

    err = γ · 2(N-2)/3

Linear in γ and N. For practical hardware (γ ~ 0.001): err < 0.1%.
The design rules remain valid as approximations under realistic noise.

### 7. Publications created

- **TECHNICAL_PAPER.md**: Full mathematical paper with proof,
  spectral decomposition, and QST bridge analysis.
- **ENGINEERING_BLUEPRINT.md**: Four design rules for quantum
  repeaters derived from the palindrome.

---

## Summary of verified claims

| Claim | Status | Evidence |
|-------|--------|----------|
| Palindromic symmetry (Heisenberg) | **Proven** | Analytical proof + N=2-8 numerical |
| Palindromic symmetry (all standard models) | **Verified** | N=3-6, XY/Ising/XXZ/DM/Heis+DM |
| CΨ=1/4 crossing equation | **Hardware validated** | 1.9% deviation, ibm_torino Q80 |
| T2* (not T2echo) governs crossing | **Hardware validated** | Same-day Ramsey confirms |
| GHZ→100% XOR, W→0% XOR | **Verified** | All models, chain + star, N=3,4 |
| Pauli weight predicts XOR fraction | **Verified** | r > 0.98 all models |
| Two Pi families (P1, P4) | **Verified** | Algebraic enumeration, 1024 candidates |
| Alternating Q for XY/YX/DM | **Verified** | M1-M2-M1 on N=3 Lindbladian |
| Depolarizing correction formula | **Verified** | err = γ·2(N-2)/3, N=3,4 |

## What remains open

- Multi-qubit palindrome untested on hardware (N≥2)
- 2/36 two-term combinations need Pauli-mixing Q operators
- Depolarizing correction formula: analytical derivation pending
- External peer review
- Optimal QST encoding state

## Authors

Thomas Wicht (Independent Researcher, Krefeld, Germany)
Claude (Anthropic)

## Repository

https://github.com/Kesendo/R-equals-C-Psi-squared
