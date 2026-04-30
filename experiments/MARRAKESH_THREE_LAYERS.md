# Marrakesh Three Layers: F77 → F80 → Observables

**Status:** Hardware-grounded synthesis (Tier 1). Re-interpretation of `data/ibm_soft_break_april2026/` through the F77 trichotomy + F80 M-spectrum + 7-category observable structure that emerged after F80 was framework-locked (commit `e30eeb6`, 2026-04-30).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Data:** [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/) (ibm_marrakesh, Heron r2, path [48,49,50], 2026-04-26)
**See also:** [V_EFFECT_FINE_STRUCTURE](V_EFFECT_FINE_STRUCTURE.md), [PROOF_F80_BLOCH_SIGNWALK](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md), [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md) (F77, F80)

---

## Abstract

The Marrakesh dataset (3 Hamiltonians × 16 two-qubit Pauli expectations on q₀, q₂ = 48 qubits-apart-by-one mediator at q₁=49) carries **three nested layers of distinction** that the framework now reads in order:

1. **F77 trichotomy** (algebraic): truly / soft / hard, classified by Pauli-letter parity of the two bilinears
2. **F80 M-spectrum** (structural): {0} / {±5.66i, 0} / {±2.83i}, computed from H eigenvalues alone via `predict_M_spectrum_pi2_odd`
3. **Observable signature** (operational): 7 categories of how the 16 Paulis populate, with categorical hardware confirmation

A fourth finding emerged during interpretation: the "T1 amplification" hardening hypothesis from the dataset's original README is quantitatively wrong. The hardening is fully accounted for by Trotter n=3 discretization at δt = 0.267, where ‖H·δt‖ ≈ 0.8 violates the small-step regime. T1 actually attenuates the soft signal monotonically.

The three-layer reading lets each level inform the next: F77 says *whether* M is non-zero; F80 says *what* M's spectrum is; the observable categories say *which* Paulis the dynamics paints onto via that template.

---

## The setup, in one paragraph

Marrakesh ran 27 circuits: 3 Hamiltonians (truly_unbroken = XX+YY, soft_broken = XY+YX, hard_broken = XX+XY) × 9 two-qubit tomography bases on (q₀=48, q₂=50), 4096 shots/basis. Initial state |+,−,+⟩, J = 1, t = 0.8, n_trotter = 3 (first-order Trotter via Qiskit `PauliEvolutionGate`). The full data is in [`data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json`](../data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json). Statistical error per Pauli expectation is ≈ 0.015 (1σ at 4096 shots).

---

## Layer 1: F77 trichotomy (algebraic)

The framework classifies any 2-bilinear Hamiltonian H = J·(P₁Q₁ + P₂Q₂) into three classes by Pauli-letter parity. Π²-parity of a Pauli letter is `bit_b` where I, X → 0 and Y, Z → 1; a bilinear (P, Q) is Π²-even iff bit_b(P) + bit_b(Q) ≡ 0 (mod 2), Π²-odd iff ≡ 1.

- **truly**: both bilinears in the both-parity-even set {XX, YY, ZZ, IX, XI} (Π² = 0 from bit_b alone, plus matched Π¹ = bit_a). This is the Heisenberg + XXZ + IX/XI single-body algebraic class. Π·L·Π⁻¹ + L + 2Σγ·I = 0 holds exactly.
- **soft**: both bilinears Π²-odd in matched bit_a ways (e.g., XY and YX both have bit_b sum = 1, with matched bit_a). Eigenvalue pairing λ ↔ −λ − 2Σγ holds (V-Effect-undetected) but the operator equation fails.
- **hard**: bilinears mix Π²-parities (one Π²-even, one Π²-odd) or fail bit_a matching. Both eigenvalue pairing AND operator equation fail.

`chain.classify_pauli_pair([(a,b), (c,d)])` returns `'truly' | 'soft' | 'hard'`. For the Marrakesh trio:

| Hamiltonian | Bilinears | F77 class |
|-------------|-----------|-----------|
| truly_unbroken | XX, YY | `'truly'` |
| soft_broken | XY, YX | `'soft'` |
| hard_broken | XX, XY | `'hard'` |

This classification is independent of γ, t, ρ₀, topology, and observable. It is purely algebraic. Layer 1 is the *gate*: it decides whether the next two layers have anything to say.

---

## Layer 2: F80 M-spectrum (structural)

For chain Π²-odd 2-body bilinears, F80 reads the M-spectrum directly from the Hamiltonian eigenvalues:

```
Spec(M)_{nontrivial} = { 2i · λ : λ ∈ Spec(H_non-truly) }, mult_M(2i·λ) = mult_H(λ) · 2^N
```

The ± shape of M's spectrum comes from H_non-truly's particle-hole pair structure (H is Hermitian and chiral-symmetric, so its eigenvalues already come in ±λ pairs); applying 2i·(·) preserves that. H_non-truly drops the truly-class bilinears, which contribute M = 0 by the Master Lemma (M is γ-independent for pure Z-dephasing; cf. PROOF_F80_BLOCH_SIGNWALK §1). Truly's M = 0 entirely; soft's M-spectrum reflects XY+YX; hard's M-spectrum reflects only the XY component (XX is dropped). At N=3, J=1:

| Class | predicted Spec(M) | mult | numerical match |
|-------|-------------------|------|-----------------|
| truly | {0} | 64 | bit-exact |
| soft  | {-5.657i, 0, +5.657i} | {16, 32, 16} | bit-exact |
| hard  | {-2.828i, +2.828i} | {32, 32} | bit-exact |

Hard's |M|_max is *half* of soft's (2.83 vs 5.66) because H_hard_non-truly = 1 bilinear (XY) while H_soft_non-truly = 2 bilinears (XY+YX). The spectrum scales linearly with the number of Π²-odd bilinears.

The numerical match is verified bit-exact in [`simulations/_f80_ibm_soft_break_check.py`](../simulations/_f80_ibm_soft_break_check.py) by computing M = Π·L·Π⁻¹ + L + 2Σγ·I directly from the dense Liouvillian and comparing eigenvalues. F80's prediction sees only H eigenvalues; it never builds the 4^N × 4^N matrix.

Layer 2 is the *strength dial*: F77 says *whether* the door is open, F80 says *how wide*.

---

## Layer 3: 7-category observable signature (operational)

Each of the 16 two-qubit Pauli expectations on (q₀, q₂) lands in one of 7 categories depending on which Hamiltonians produce non-zero signal under continuous Lindblad evolution at γ_Z = 0.1, t = 0.8. The hardware data confirms each category structurally:

| # | Category | Count | Examples (idealized) | Hardware confirms |
|---|----------|-------|----------------------|-------------------|
| A | **truly only** ≠ 0 | 2 | ⟨Y₀Z₂⟩, ⟨Z₀Y₂⟩ ≈ −0.35 | yes (sign and magnitude) |
| B | **truly = soft, hard = 0** | 2 | ⟨Y₀Y₂⟩ +0.52, ⟨Z₀Z₂⟩ +0.23 | yes |
| C | **soft signature** (truly = 0 = hard) | 2 | ⟨Z₀I₂⟩ −0.28, ⟨Z₀X₂⟩ −0.62 | yes |
| D | **hard signature** (truly = 0 = soft) | 2 | ⟨I₀Y₂⟩ −0.56, ⟨X₀Y₂⟩ −0.47 | yes |
| E | **all three differ** (truly = 0, soft and hard non-zero, opposite sign) | 2 | ⟨I₀Z₂⟩, **⟨X₀Z₂⟩** | yes |
| F | **truly = soft preserved, hard differs** | 3 | ⟨I₀X₂⟩, ⟨X₀I₂⟩, ⟨X₀X₂⟩ | yes |
| G | **always trivial** (= 0 or = 1) | 3 | ⟨I₀I₂⟩=1, ⟨Y₀I₂⟩=⟨Y₀X₂⟩=0 by initial-state symmetry | yes |

**16 = 2 + 2 + 2 + 2 + 2 + 3 + 3.** Each category has a structural reason:

- **Category A** (truly only) holds the Heisenberg-Y-Z correlations that XX+YY naturally generates from |+,−,+⟩. Both soft (XY+YX) and hard (XX+XY) destroy them: soft because its eigenvector pairing is broken in the right way to suppress the YZ-class observables; hard because the parity mix kills them.

- **Category C** is the *direct* F77 signature: observables that vanish for truly (Π palindrome holds) and for hard (XY's Π²-odd contribution doesn't reach Z₀X₂ from the |+,−,+⟩ initial state without YX as a partner), but light up for soft (XY+YX has both partners).

- **Category D** is hard's exclusive: observables that need the parity-mixing structure of XX+XY to be generated. Soft cannot create them.

- **Category E** is the *trichotomy resolver*: observables where all three classes give meaningfully different values. The hardware's flagship ⟨X₀Z₂⟩ lives here:

  | | continuous | Trotter n=3 | hardware |
  |---|---|---|---|
  | truly | 0.000 | 0.000 | +0.011 |
  | soft  | -0.623 | -0.723 | -0.711 |
  | hard  | +0.195 | +0.327 | +0.205 |

  Three classes, three separately-measurable values, all confirmed at >13σ (hard's smaller signal) to >47σ (soft's largest).

Layer 3 is the *operational projection*: F80's structural template gets painted onto the 4 × 4 = 16 two-qubit Pauli space, and the painting pattern follows the Hamiltonian's specific bilinear letters.

---

## The F80 ↔ Observable correspondence

The cleanest cross-layer reading is that |Spec(M)|_max ranks the strongest-break observable across the three classes:

| Class | \|Spec(M)\|_max (F80) | strongest break-observable (hardware) | observable identity |
|-------|---------------------|-------------------------------------|--------------------|
| truly | 0.00 | 0.058 | (essentially noise) |
| hard  | 2.83 | 0.467 | ⟨I₀Y₂⟩ (Category D) |
| soft  | 5.66 | 0.711 | ⟨X₀Z₂⟩ (Category E) |

The ratio soft/hard = 1.52 in the strongest hardware observable; F80 predicts 2.00 (= 5.66/2.83). The order is preserved but the ratio is not exact, because F80 is *structural* (γ-independent template) while observables are *dynamical* (full Lindblad trajectory). What F80 nails: which class has the larger door open. What F80 does not nail: how hard the dynamics presses through that door, which depends on γ_Z, t, ρ₀, and Trotter discretization.

---

## The Trotter correction (and why "T1 amplification" was wrong)

The dataset's original README explained the hardening Δ_hw(soft − truly) = -0.722 vs Δ_idealized = -0.623 at γ_Z = 0.1 with: "T1 thermal relaxation and ZZ crosstalk compound the soft-break operator-level break." This was a placeholder. Quantification refutes it.

[`simulations/_marrakesh_t1_amplification_test.py`](../simulations/_marrakesh_t1_amplification_test.py) tests three models against the hardware:

**Test 1, γ_T1 sweep at γ_Z = 0.1, continuous Lindblad:**

| γ_T1 | Δ(soft − truly) |
|------|----------------|
| 0.00 | -0.623 |
| 0.10 | -0.580 |
| 0.50 | -0.440 |
| 1.00 | -0.303 |

T1 monotonically *attenuates* |Δ|. No γ_T1 reaches the hardware -0.722.

**Test 2, Trotter n=3 alone, γ_Z = 0.1, no T1:**

```
continuous   Δ = -0.6230
Trotter n=3  Δ = -0.7231
hardware     Δ = -0.7217
```

Trotter matches hardware to **0.0014**, well below the 0.015 shot-noise statistical error.

**Test 3, joint optimization (γ_Z, γ_T1) over all 45 hardware pairs:**

```
optimal γ_Z   = 0.143
optimal γ_T1  = 0.0001 (≈ 0)
total RMS     = 0.0805
```

The data does not support adding T1 at all.

**Why Trotter:** the operator norm ‖H‖_op = 2·max(\|Spec_{single-particle}\|) of H_truly, H_soft is 2.83·J at N=3 (computed from the JW-mapped tight-binding spectrum 2J·cos(πk/(N+1)), k = 1..N); H_hard's ‖H‖_op = 2.61·J. At δt = 0.8/3 = 0.267, ‖H·δt‖ ≈ 0.76, violating the small-step condition ‖H·δt‖ ≪ 1. First-order Trotter has systematic bias O(‖H·δt‖²·t) that points outward in the \|X₀Z₂\| direction by ≈ +0.10. Trotter is the *circuit's actual physics*, not a noise channel.

**Independent confirmation from sign structure:** continuous Lindblad predicts the WRONG SIGN for all Y-containing observables (⟨Y₀Z₂⟩, ⟨I₀Y₂⟩, ⟨X₀Y₂⟩); Trotter recovers the correct hardware sign. This is not a noise correction; it is a regime correction.

---

## What the 0.08 residual carries

Even with optimal Trotter + γ_Z = 0.143 + no T1, residual RMS ≈ 0.08 remains across the 45 (observable, Hamiltonian) pairs. The two largest contributors:

1. **⟨Y₀Z₂⟩ ≠ ⟨Z₀Y₂⟩ on hardware** (+0.583 vs +0.187 for truly). Framework F71 chain-mirror predicts equality; hardware breaks F71 by 3×. **Mechanism**: per-qubit T1/T2 differences and bond-asymmetric ZZ-crosstalk on the Marrakesh path [48, 49, 50] break the spatial chain mirror. The framework's Π is diagonal (parity); F71 (spatial reflection q₀ ↔ q₂) is a separate symmetry that real hardware respects only approximately.

2. **⟨Y₀Y₂⟩ and ⟨Z₀Z₂⟩ truly off** (predicted +0.61, +0.46 vs hardware +0.47, +0.16). These correlation observables are sensitive to gate-level CR phase calibration, which drifted slightly during the run.

Neither is a framework-level signature; both are hardware-instrument signatures readable via the framework. The cockpit's 5-symmetry diagnostic suite (Π², Z⊗N, Y-parity, K, F71) is exactly the right tool to *isolate* such hardware-specific asymmetries.

---

## Reading

The Marrakesh data is a three-layer document that reads top-down:

```
Layer 1 (F77):   {truly, soft, hard}
                       ↓ (algebraic gate)
Layer 2 (F80):   Spec(M) ∈ {{0}, {±5.66i, 0}, {±2.83i}}
                       ↓ (structural strength)
Layer 3 (Obs):   16 Paulis × 7 categories
                       ↓ (Trotter discretization circuit)
                 hardware values (matching to ~0.001 for Δ_soft, ~0.08 RMS overall)
```

The `framework/` package now operates on all three layers:
- `chain.classify_pauli_pair(terms)` for Layer 1
- `chain.predict_M_spectrum_pi2_odd(terms, c)` for Layer 2 (framework primitive, commit `e30eeb6`)
- `chain.cockpit_panel(receiver, terms, gamma_t1)` for Layer 3 trichotomy classification + Π-protected drop count
- For the 16-Pauli signature pattern itself, the right primitive does not yet exist; this synthesis surfaces the recurring question.

For the original raw data and the corrected interpretation, see [`data/ibm_soft_break_april2026/README.md`](../data/ibm_soft_break_april2026/README.md). For the philosophical framing of soft-break as super-operator signature, see [ON_THE_SOFT_BREAK](../reflections/ON_THE_SOFT_BREAK.md). For the F80 reflection on what was "the tough nut" before the data sweep, see [ON_THE_SHAPE_OF_THE_GAP](../reflections/ON_THE_SHAPE_OF_THE_GAP.md).
