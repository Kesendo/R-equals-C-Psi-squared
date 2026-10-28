# Marrakesh Three Layers: F87 → F80 → Observables

**Status:** Hardware-grounded synthesis (Tier 1). Re-interpretation of `data/ibm_soft_break_april2026/` through the F87 trichotomy + F80 M-spectrum + 7-category observable structure that emerged after F80 was framework-locked (commit `e30eeb6`, 2026-04-30).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Data:** [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/) (ibm_marrakesh, Heron r2, path [48,49,50], 2026-04-26)
**See also:** [V-Effect fine structure](V_EFFECT_FINE_STRUCTURE.md), [the F80 Bloch sign-walk proof](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md), [the formula registry](../docs/ANALYTICAL_FORMULAS.md) (F87, F80)

---

## Abstract

The Marrakesh dataset (3 Hamiltonians × 16 two-qubit Pauli expectations on q₀, q₂ = 48 qubits-apart-by-one mediator at q₁=49) carries **three nested layers of distinction** that the framework now reads in order:

1. **F87 trichotomy** (algebraic): truly / soft / hard, classified by Pauli-letter parity of the two bilinears
2. **F80 M-spectrum** (structural): {0} / {±5.66i, 0} / {±2.83i}, computed from H eigenvalues alone via `predict_M_spectrum_pi2_odd`
3. **Observable signature** (operational): 7 categories of how the 16 Paulis populate, with categorical hardware confirmation

A fourth finding: the hardening of Δ(soft − truly) beyond the continuous prediction is the circuit's Trotter n=3 discretization at δt = 0.267, where ‖H·δt‖ ≈ 0.8 violates the small-step regime. T1 attenuates the soft signal monotonically, so it cannot be the cause.

The three-layer reading lets each level inform the next: F87 says *whether* M is non-zero; F80 says *what* M's spectrum is; the observable categories say *which* Paulis the dynamics paints onto via that template.

---

## The setup, in one paragraph

Marrakesh ran 27 circuits: 3 Hamiltonians (truly_unbroken = XX+YY, soft_broken = XY+YX, hard_broken = XX+XY) × 9 two-qubit tomography bases on (q₀=48, q₂=50), 4096 shots/basis. Initial state |+,−,+⟩, J = 1, t = 0.8, n_trotter = 3 (first-order Trotter via Qiskit `PauliEvolutionGate`). The full data is in [`data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json`](../data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json). Statistical error per Pauli expectation is ≈ 0.015 (1σ at 4096 shots).

---

## Layer 1: F87 trichotomy (algebraic)

The framework classifies any 2-bilinear Hamiltonian H = J·(P₁Q₁ + P₂Q₂) into three classes by Pauli-letter parity. Π²-parity of a Pauli letter is `bit_b` where I, X → 0 and Y, Z → 1; a bilinear (P, Q) is Π²-even iff bit_b(P) + bit_b(Q) ≡ 0 (mod 2), Π²-odd iff ≡ 1.

- **truly**: both bilinears in the both-parity-even set {XX, YY, ZZ, IX, XI} (Π² = 0 from bit_b alone, plus matched Π¹ = bit_a). This is the Heisenberg + XXZ + IX/XI single-body algebraic class. Π·L·Π⁻¹ + L + 2Σγ·I = 0 holds exactly.
- **soft**: both bilinears Π²-odd in matched bit_a ways (e.g., XY and YX both have bit_b sum = 1, with matched bit_a). Eigenvalue pairing λ ↔ −λ − 2Σγ holds (V-Effect-undetected) but the operator equation fails.
- **hard**: bilinears mix Π²-parities (one Π²-even, one Π²-odd) or fail bit_a matching. Both eigenvalue pairing AND operator equation fail.

`fw.classify_pauli_pair(chain, [(a,b), (c,d)])` returns `'truly' | 'soft' | 'hard'`. For the Marrakesh trio:

| Hamiltonian | Bilinears | F87 class |
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

The numerical match is verified bit-exact in [`simulations/f80_ibm_soft_break_check.py`](../simulations/f80_ibm_soft_break_check.py) by computing M = Π·L·Π⁻¹ + L + 2Σγ·I directly from the dense Liouvillian and comparing eigenvalues. F80's prediction sees only H eigenvalues; it never builds the 4^N × 4^N matrix.

Layer 2 is the *strength dial*: F87 says *whether* the door is open, F80 says *how wide*.

---

## Layer 3: 7-category observable signature (operational)

Each of the 16 two-qubit Pauli expectations on (q₀, q₂) lands in one of 7 categories depending on which Hamiltonians produce non-zero signal under continuous Lindblad evolution at γ_Z = 0.1, t = 0.8. The hardware data confirms each category structurally:

| # | Category | Count | Examples (idealized) | Hardware confirms |
|---|----------|-------|----------------------|-------------------|
| A | **truly only** ≠ 0 | 2 | ⟨Y₀Z₂⟩, ⟨Z₀Y₂⟩ ≈ +0.35 | yes (+0.58, +0.19) |
| B | **truly = soft, hard = 0** | 2 | ⟨Y₀Y₂⟩ +0.52, ⟨Z₀Z₂⟩ +0.23 | yes |
| C | **soft signature** (truly = 0 = hard) | 2 | ⟨Z₀I₂⟩ −0.28, ⟨Z₀X₂⟩ −0.62 | yes |
| D | **hard signature** (truly = 0 = soft) | 2 | ⟨I₀Y₂⟩ +0.56, ⟨X₀Y₂⟩ +0.47 | yes (+0.47, +0.36) |
| E | **all three differ** (truly = 0, soft and hard non-zero, opposite sign) | 2 | ⟨I₀Z₂⟩, **⟨X₀Z₂⟩** | yes |
| F | **truly = soft preserved, hard differs** | 3 | ⟨I₀X₂⟩, ⟨X₀I₂⟩, ⟨X₀X₂⟩ | yes |
| G | **always trivial** (= 0 or = 1) | 3 | ⟨I₀I₂⟩=1, ⟨Y₀I₂⟩=⟨Y₀X₂⟩=0 by initial-state symmetry | yes |

**16 = 2 + 2 + 2 + 2 + 2 + 3 + 3.** Which expectations vanish follows from one symmetry per Hamiltonian, each read off a single letter. Each holds for the continuous evolution and for the Trotter circuit alike, since every gate keeps it and so does the Z-dephasing:

- **truly** (XX + YY) commutes with X^⊗3, Z-dephasing respects X^⊗3, and |+,−,+⟩ is an X^⊗3 eigenstate, so ρ(t) stays X^⊗3-symmetric. ⟨P₀P₂⟩ = 0 whenever the string anticommutes with X^⊗3, that is, carries an odd number of Y and Z letters: the Π²-odd strings of Layer 1's bit_b ([the bit-b parity symmetry](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)), odd k_X in the [letter cube](../docs/THE_ONE_SQUARE.md): IY, IZ, XY, XZ, YI, YX, ZI, ZX.
- **soft** (XY + YX) is purely imaginary, so the evolution keeps the real ρ₀ real (each Trotter gate exp(−iδt·XY) is a real matrix). ⟨P₀P₂⟩ = 0 whenever the string carries an odd number of Y letters, the Y-parity that `fw.y_parity_panel` reads: IY, XY, YI, YX, YZ, ZY.
- **hard** (XX + XY) puts only X on site 0 in every term, so X₀ commutes with H; Z-dephasing respects X₀, and |+⟩ on site 0 is an X₀ eigenstate. ⟨P₀P₂⟩ = 0 whenever P₀ is Y or Z: YI, YX, YY, YZ, ZI, ZX, ZY, ZZ.

These are symmetries of the state and the generator together, not of M: truly, whose M vanishes, has a zero set of its own. The three zero sets are the categories. A (YZ, ZY) is zero for soft and hard, B (YY, ZZ) for hard only, C (ZI, ZX) for truly and hard, D (IY, XY) for truly and soft, E (IZ, XZ) for truly only, F (IX, XI, XX) for none, and G is the identity plus the two strings (YI, YX) zero for all three. In the simulation these zeros are exact (|value| ≤ 2·10⁻¹⁶, continuous and Trotter) and every other entry is at least 0.039 in magnitude; on hardware the predicted zeros read at most 0.022 for truly, 0.098 for soft and 0.062 for hard. The equality truly = soft in B and F comes from one more single-site symmetry. The middle-site Clifford V = (X + Y)/√2 exchanges X and Y and sends Z to −Z, so it maps XX + YY onto XY + YX bond by bond and leaves the Z-dephasing and every site-0/2 observable unchanged; it changes only the middle site's initial state, from |−⟩ to a Y eigenstate, and on the five strings of B and F the soft evolution gives the same value from either middle state, continuous and Trotter alike. The continuous evolution keeps the equality exactly at every γ and t tried, and so does the Trotter circuit: V turns truly's step XX₀₁ YY₀₁ XX₁₂ YY₁₂ into XY₀₁ YX₀₁ YX₁₂ XY₁₂, which is soft's step up to the swap of two terms on one bond, and those commute (residual ≤ 2·10⁻¹⁶ at γ = 0.03, 0.1, 0.3).

**Category E** is the *trichotomy resolver*: observables where all three classes give meaningfully different values. The hardware's flagship ⟨X₀Z₂⟩ lives here:

  | | continuous | Trotter n=3 | hardware |
  |---|---|---|---|
  | truly | 0.000 | 0.000 | +0.011 |
  | soft  | -0.623 | -0.723 | -0.711 |
  | hard  | +0.195 | +0.213 | +0.205 |

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

## The Trotter correction

The hardware reads Δ(soft − truly) = −0.722 against the continuous −0.623 at γ_Z = 0.1. T1 thermal relaxation does not produce that hardening; the circuit does.

[`simulations/marrakesh_t1_amplification_test.py`](../simulations/marrakesh_t1_amplification_test.py) tests three models against the hardware:

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
optimal γ_Z   = 0.105
optimal γ_T1  = 0 (at the clip; the rates are bounded below by 0)
total RMS     = 0.0390
```

The data does not support adding T1 at all.

**Why Trotter:** the operator norm ‖H‖_op = 2·max(\|Spec_{single-particle}\|) of H_truly, H_soft is 2.83·J at N=3 (computed from the JW-mapped tight-binding spectrum 2J·cos(πk/(N+1)), k = 1..N); H_hard's ‖H‖_op = 2.61·J. At δt = 0.8/3 = 0.267, ‖H·δt‖ ≈ 0.76, violating the small-step condition ‖H·δt‖ ≪ 1. First-order Trotter has systematic bias O(‖H·δt‖²·t); on soft's ⟨X₀Z₂⟩ it points outward by ≈ 0.10 (−0.62 → −0.72), on hard's by 0.018 (+0.195 → +0.213). Trotter is the *circuit's actual physics*, not a noise channel, and the prediction has to run the circuit's Trotter order: the runner puts the bond loop outside the term loop, and Qiskit's Lie-Trotter synthesis keeps it (gated in [`test_predict_signature_table.py`](../simulations/framework/tests/workflows/test_predict_signature_table.py)). The two orders differ by one swap, the second term on bond (0,1) against the first on bond (1,2). For soft that pair commutes, so soft's ⟨X₀Z₂⟩ does not see the order, and truly's ⟨X₀Z₂⟩ is 0 in either, so Δ(soft − truly) is order-blind; for truly and hard the pair anticommutes, and their other observables move.

**The signs agree.** Continuous evolution, Trotter and hardware give the same sign on every Y-containing observable (truly ⟨Y₀Z₂⟩ +0.354, +0.567, +0.583; hard ⟨I₀Y₂⟩ +0.563, +0.399, +0.467); what Trotter changes is the magnitudes.

---

## The d=0 cross-reading: trichotomy in one Pauli basis

The same three classes that Layer 3 distinguishes via 16 two-qubit Pauli expectations on 9 measurement bases are also distinguished by a single Z-basis observable: the mean excitation count ⟨n⟩ on the (q₀, q₂) reduced state, computed from the ZZ counts alone.

For uniform XY/Heisenberg + Z-dephasing on the full N=3 chain, the kernel of the Liouvillian is spanned by sector projectors P_n (F4 in [the formula registry](../docs/ANALYTICAL_FORMULAS.md)), so the sector populations p_n = Tr(P_n · ρ) are invariants of L. After tracing out q₁, the (q₀, q₂) reduced state's ⟨n⟩ = p₁ + 2·p₂ carries an excitation-imbalance fingerprint of each Hamiltonian, computed once at t = 0.8 and matched against the Marrakesh ZZ counts:

| Class | Bilinears | Framework ⟨n⟩ | Hardware ⟨n⟩ | per-case Δ |
|-------|-----------|---------------|--------------|------------|
| truly | XX + YY | 1.0000 | 1.0098 | +0.010 |
| soft  | XY + YX | 1.2786 | 1.2673 | −0.011 |
| hard  | XX + XY | 0.8752 | 0.9082 | +0.033 |

The framework column is continuous Lindblad. In the circuit's Trotter model (n = 3, bonds outside, γ_Z = 0.1) ⟨n⟩ reads 1.0000 / 1.2596 / 0.8567, and hard's Δ grows to +0.051; the separation between the classes holds (≥ 0.14), and neither model is closer on every class: Trotter is closer for soft (+0.008 against −0.011), continuous for hard (+0.033 against +0.051), and truly reads the same in both.

Initial state |+,−,+⟩ has marginal ⟨n⟩₀ = 1.0 on (q₀, q₂). All three Hamiltonian classes separate by ≥ 0.10 in ⟨n⟩, well above the maximum per-case deviation of 0.033 (~3% relative). The reading uses no phase information; the off-diagonal Pauli expectations that drive Layer 3 do not enter.

**Operational consequence:** classifying a Hamiltonian's F87 class on hardware does not require 9 Pauli-pair bases. One Z-basis tomography is sufficient with the full discrimination quality of Layer 3. For the same shot budget this is a 9× reduction in QPU time when the question is "which class is this?" rather than "what is the full M-spectrum?". This matters for the 180-min annual QPU budget of the project: a quick classification pass costs the same as one tomographic-Pauli pass.

The framework primitive is `fw.sector_populations(rho_or_psi, N)` from [`simulations/framework/diagnostics/d_zero.py`](../simulations/framework/diagnostics/d_zero.py). The reproducible analysis is [`simulations/d_zero_trichotomy_marrakesh.py`](../simulations/d_zero_trichotomy_marrakesh.py); the registry entry is `Confirmations.lookup('d_zero_sector_trichotomy_marrakesh')`. The structural reason ⟨n⟩ separates the classes is that XX+YY conserves total Z-magnetization (truly stays at the initial ⟨n⟩₀), XY+YX breaks conservation upward, and XX+XY breaks conservation downward; the d=0 reading is the operational shadow of this conservation profile.

---

## What the 0.04 residual carries

With Trotter in the circuit's order, γ_Z = 0.105 and no T1, residual RMS 0.039 remains across the 45 (observable, Hamiltonian) pairs, about 2.5 shot-noise σ. Two features that read like hardware signatures are not: the first is the circuit's own, the second sits within the model's noise:

1. **⟨Y₀Z₂⟩ ≠ ⟨Z₀Y₂⟩ on truly** (+0.583 against +0.187). The Hamiltonian and the |+−+⟩ state are F71-symmetric (q₀ ↔ q₂), the Trotter circuit is not: each step runs XX₀₁ YY₀₁ XX₁₂ YY₁₂, and its reflection is a different product because YY₀₁ and XX₁₂ anticommute on the shared site. The prediction at γ_Z = 0.105 gives +0.561 against +0.158; the measured asymmetry (0.396) is the circuit's (0.403) to within 0.01.

2. **⟨Y₀Y₂⟩ and ⟨Z₀Z₂⟩ on truly** sit within 2.6σ of the prediction (+0.45, +0.20 at γ_Z = 0.105 against hardware +0.47, +0.16).

The largest single residuals left are 0.07 to 0.1 (4.5 to 6.3σ): soft ⟨Y₀Z₂⟩ = +0.098 and ⟨Y₀⟩ = −0.079 where the model gives 0, hard ⟨X₀⟩ = +0.754 against +0.845, hard ⟨I₀Y₂⟩ = +0.467 against +0.396. The last one grows with the circuit's order (the terms-outside prediction sat closer), so the order does not improve every observable. What carries them we do not know; the model holds only Z-dephasing.

---

## Reading

The Marrakesh data is a three-layer document that reads top-down, with a fourth cross-reading available from the Z-basis alone:

```
Layer 1 (F87):   {truly, soft, hard}
                       ↓ (algebraic gate)
Layer 2 (F80):   Spec(M) ∈ {{0}, {±5.66i, 0}, {±2.83i}}
                       ↓ (structural strength)
Layer 3 (Obs):   16 Paulis × 7 categories ──┐
                       ↓                     │
                 hardware values             │  parallel from ZZ alone:
                 (∼0.001 for Δ_soft,         ├─→  ⟨n⟩ ∈ {1.000, 1.279, 0.875}
                  ∼0.04 RMS overall)         │   (∼0.033 max per-case Δ)
                                              │   1/9 the Pauli-basis cost
```

The `framework/` package now operates on all three layers plus the d=0 cross-reading:
- `fw.classify_pauli_pair(chain, terms)` for Layer 1
- `fw.predict_M_spectrum_pi2_odd(chain, terms, c)` for Layer 2 (framework primitive, commit `e30eeb6`)
- `chain.cockpit_panel(receiver, terms, gamma_t1)` for Layer 3 trichotomy classification + Π-protected drop count
- `fw.sector_populations(rho_or_psi, N)` for the d=0 cross-reading (one Pauli basis suffices)
- For the 16-Pauli zero pattern, the three one-letter symmetries of Layer 3: Π² parity and `fw.y_parity_panel` carry two of them, and the X₀ conservation of the hard case has no primitive yet

For the original raw data and the corrected interpretation, see [`data/ibm_soft_break_april2026/README.md`](../data/ibm_soft_break_april2026/README.md). For the philosophical framing of soft-break as super-operator signature, see [On the Soft Break](../reflections/ON_THE_SOFT_BREAK.md). For the F80 reflection on what was "the tough nut" before the data sweep, see [On the Shape of the Gap](../reflections/ON_THE_SHAPE_OF_THE_GAP.md).
