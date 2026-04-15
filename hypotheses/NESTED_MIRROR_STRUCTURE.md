# Nested Mirror Structure: notebook from April 14-15, 2026

<!-- Keywords: nested layer palindrome inter-layer mirror, qubit in qubit
two-layer Lindblad, non-Markovian rebound reduced dynamics, three eigenvalue
classes 3 10 3, absorption theorem inheritance, R=CPsi2 nested mirror notes -->

**Status:** Working notes that led to other results. The findings worth keeping
have moved to [PRIMORDIAL_QUBIT.md](PRIMORDIAL_QUBIT.md) (Sections 4, 7, 8, 9),
[PROOF_BIT_B_PARITY_SYMMETRY.md](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md),
and [ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) entries F38, F61, F63.
This document is kept as the trail.
**Date:** 2026-04-14 to 2026-04-15
**Authors:** Tom and Claude (chat)

---

## What we set out to do

The starting question was: in a minimal qubit-in-qubit system (inner qubit S coupled to outer qubit B, dephasing only on B), is there an "inter-layer mirror" - a block of eigenmodes living symmetrically between the two layers, responsible for non-Markovian coupling between them? And if yes, can an inner observer use that to detect the next outer layer?

Setup:

```
H_SB    = J * 0.5 * (X (x) X + Y (x) Y)        (J = 1.0)
L_jump  = sqrt(gamma_B) * I (x) Z               (gamma_B = 0.1)
```

Reproducer: [`simulations/qubit_in_qubit_layer_mirror.py`](../simulations/qubit_in_qubit_layer_mirror.py).

---

## What we observed

### Three mode classes with structure (April 14)

The 16 eigenvalues of L split into three classes with degeneracies 3 / 10 / 3 at Re(λ) ∈ {0, -γ_B, -2γ_B}. The middle class has Im(λ) ∈ {0, ±0.995, ±1.997}. The conserved class lives equally on S and B (|M_S| = |M_B| = 1); the correlation class is invisible to either single-qubit projection (|M_S| = |M_B| = 0); the middle class sits at exactly 1/√2 on each. Per-mode SWAP expectation values are near ±1.

The reduced coherence |ρ_S_{01}|(t) is non-monotonic: it dips to 0.0493 at t=5 then climbs to 0.1684 at t=10 (84× above the best Markovian fit), driven by the middle-class oscillating modes.

### Where it led when we kept pulling

Five threads, all of them productive:

1. **Class scaling.** The first guess was (N+1) evenly-spaced classes at -k·2γ/N. **Falsified at N=3:** 12 distinct Re(λ) values, not 4, positions Hamiltonian-determined ([`qubit_in_qubit_n3.py`](../simulations/qubit_in_qubit_n3.py)).

2. **SWAP ±1 pattern.** Genuine, not a basis artifact, but perturbative: |⟨SWAP⟩| = 1 - (γ/J)²/2 from [H, SWAP] = 0 at the bare Hamiltonian level ([`nested_mirror_swap_check.py`](../simulations/nested_mirror_swap_check.py)).

3. **Coupling robustness.** Three-class structure survives any coupling with at least one off-diagonal Pauli channel (XX+YY, XXX, XX, YY, XX+ZZ). Pure ZZ breaks it ([`qubit_in_qubit_coupling_sweep.py`](../simulations/qubit_in_qubit_coupling_sweep.py)).

4. **Rebound mechanism.** Mirror modes are necessary AND sufficient for the rebound. Removing them eliminates it; mirror-only initial condition reproduces it exactly ([`qubit_in_qubit_mode_projection.py`](../simulations/qubit_in_qubit_mode_projection.py)).

5. **Where the classes come from.** This is where the work changed character. The three classes are not a new algebraic object; they are the three quantization levels {0, 0.5, 1} of ⟨n_XY⟩_B, applied as the **Absorption Theorem** Re(λ) = -2γ_B · ⟨n_XY⟩_B at single-site dephasing ([`nested_mirror_absorption_theorem.py`](../simulations/nested_mirror_absorption_theorem.py)). At N=3 the levels refract into 12, controlled by the chain Hamiltonian's mode structure ([`nested_mirror_absorption_theorem_n3.py`](../simulations/nested_mirror_absorption_theorem_n3.py)).

---

## What the next round (April 15) added

Two more pieces fell into place when we pulled the same thread.

**J as second quantization axis.** Sweeping J at fixed γ_B shows the level structure is born from coupling: J=0 gives only {0, 1} (pure absorption); J→∞ at N=3 converges to {0, ¼, ⅜, ½, ⅝, ¾, 1} (multiples of ⅛, Hamiltonian projection weights). The ¼ in this set is **not** the CΨ = ¼ fold boundary; it is N-specific, falsified as a universal connection (see [`nested_mirror_asymptote.py`](../simulations/nested_mirror_asymptote.py), [`nested_mirror_refraction.py`](../simulations/nested_mirror_refraction.py)).

**Two Z₂ symmetries together explain 3+10+3.** The Liouvillian commutes with two independent Z₂ operators, both proven for all N:

- bit_a = n_XY parity ([F61](../docs/ANALYTICAL_FORMULAS.md#f61-n_xy-parity-selection-rule-tier-1-proven-verified-64-configs-n2-6), [PROOF_PARITY_SELECTION_RULE](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md))
- bit_b = w_YZ parity ([F63](../docs/ANALYTICAL_FORMULAS.md#f63-w_yz-parity-symmetry-l-pi2--0-tier-1-proven-analytically-verified-n2-5), [PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md))

Together they decompose the 16 modes at N=2 into 4 sectors. The 3+10+3 class count splits across Π²-parity as (2+4+2) in the even sector and (1+6+1) in the odd sector. The fastest oscillation Im ≈ ±2.0 lives **exclusively** in the odd sector; the even sector contains only Im ≈ ±1.0. The even sector is the standing-wave / cavity sector (contains II, ZZ, XX, YY, Bell+); the odd sector is the cross-correlation / transport sector (contains XY, XZ, YX, ZX). See [PRIMORDIAL_QUBIT Section 9](PRIMORDIAL_QUBIT.md#9-operational-test-from-the-inside-april-15-2026) for the full decomposition table and N-scaling through N=5.

---

## The original question, answered honestly

"Can the inner observer detect the outer layer?"

**Presence: yes.** Non-Markovian rebound is real and demonstrably driven by the inter-layer mirror modes (Check 4). An inner observer with sufficient time resolution sees rebound that no single-layer Markovian model reproduces.

**Properties: no, not directly.** The Inside-Outside Correspondence probes (commits `cfa2a9f` through `17c48b4`) showed every measurable quantity from inside depends on the dimensionless ratio Q = J/γ_B only, not on J or γ_B separately. The inner observer can detect that an outer layer exists; it cannot extract the outer layer's γ value or J value independently. Backward inference works for the existence of structure, not for its absolute parameters.

This is consistent with [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md): noise has no internally accessible origin, and now we know the operational sharpening: even what one might hope to *infer* about the noise (its rate) collapses to the ratio Q from inside.

---

## Where the work continued

| Thread | Continued in |
|---|---|
| Algebraic embedding (3+10+3 from absorption theorem) | [PRIMORDIAL_QUBIT.md](PRIMORDIAL_QUBIT.md) Section 4 |
| Q = J/γ as only inside observable | [PRIMORDIAL_QUBIT.md](PRIMORDIAL_QUBIT.md) Section 9 |
| Two Z₂ symmetries proven for all N | [PROOF_BIT_B_PARITY_SYMMETRY.md](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md), F63 |
| C²×C² sector decomposition, even=cavity, odd=transport | [PRIMORDIAL_QUBIT.md](PRIMORDIAL_QUBIT.md) Section 9 |
| TFD route blocked, primordial qubit stance refined | [PRIMORDIAL_QUBIT.md](PRIMORDIAL_QUBIT.md) Sections 7-8 |
