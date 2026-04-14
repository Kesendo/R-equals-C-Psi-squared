# Nested Mirror Structure: An Inter-Layer Palindrome in the Minimal Qubit-in-Qubit Lindblad System

<!-- Keywords: nested layer palindrome inter-layer mirror, qubit in qubit
two-layer Lindblad, non-Markovian rebound reduced dynamics, partial trace
eigenmode structure, SWAP eigenvalue degenerate eigenspace, one-way nested
gamma, EQ-013 sub-question 2 partial answer, R=CPsi2 layer mirror hypothesis -->

**Tier:** 3 (hypothesis, partially verified: N=2 structure confirmed, N-scaling falsified)
**Status:** All four checks completed (2026-04-14). Class-scaling falsified; boundary classes, palindromic pairing, coupling robustness, and rebound mechanism confirmed.
**Date:** 2026-04-14
**Authors:** Tom and Claude (chat)
**Reproducer:** [`simulations/qubit_in_qubit_layer_mirror.py`](../simulations/qubit_in_qubit_layer_mirror.py)
**Depends on:** [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [GAMMA_IS_LIGHT](GAMMA_IS_LIGHT.md), [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013)
**Related:** [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md), [V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS](../reflections/V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS.md)

---

## What this document is about

If the framework's one-way nesting principle holds ([EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 1: γ enters each layer from the next outer layer, never the reverse), then a two-layer Lindblad system should exhibit a spectral structure that does more than repeat the single-layer palindromic symmetry proven in [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md). It should carry an **inter-layer mirror**: a block of eigenmodes that lives symmetrically between the inner and outer layer and is responsible for the non-Markovian coupling between them.

This document records the first numerical evidence that such a structure exists in the minimal qubit-in-qubit nest, states the hypothesis it points to, and lists the verifications still required before the hypothesis graduates.

---

## Setup

Two qubits S (inner system) and B (outer bath):

```
H_SB    = J * 0.5 * (X (x) X + Y (x) Y)        (XX+YY coupling, J = 1.0)
L_jump  = sqrt(gamma_B) * I (x) Z               (dephasing on B only, gamma_B = 0.1)
```

Explicitly: **no direct dephasing on S**. S decoheres only through the coupling to B, and B's dephasing is phenomenological (the outermost interface of the simulated stack). This is the minimal setup that distinguishes "S feels γ directly" from "S feels γ through B".

Initial state: rho_S(0) = |+><+|, rho_B(0) = I/2. The bath carries no prior information.

---

## Observation 1: non-Markovian rebound in rho_S

The reduced coherence |rho_S_{01}|(t) is not monotonic:

| t    | \|S01\| two-layer | \|S01\| one-layer best fit | difference |
|------|-------------------:|---------------------------:|-----------:|
|  0.0 | 0.5000             | 0.5000                     | +0.0000    |
|  1.0 | 0.2845             | 0.2845                     | +0.0000    |
|  5.0 | 0.0493             | 0.0298                     | +0.0195    |
| 10.0 | 0.1684             | 0.0018                     | **+0.1666**|
| 15.0 | 0.0712             | 0.0001                     | +0.0710    |

The one-layer comparison is a single-qubit Lindblad with γ_eff = 0.282 fitted at t=1. Markovian single-qubit dephasing is strictly monotonic (d|S01|/dt ≤ 0 always); the two-layer reduced dynamics violates this at t between 5 and 10, peaking at |S01| = 0.1684 when the monotone reference is already at 0.0018. The ratio is 84×.

**Interpretation:** information flows from B back to S between t=5 and t=10. No fixed effective γ reproduces this trajectory. The two-layer system is not representable as one layer with a renormalized dissipator.

---

## Observation 2: three eigenvalue classes with distinct partial-trace structure

Diagonalizing L (16×16) yields exactly three eigenvalue classes with degeneracies 3/10/3. The palindromic pairing Re(λ_a) + Re(λ_b) = -2γ_B = -0.2 is satisfied (consistent with [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)).

Per-mode partial-trace weights |M_S| = ||Tr_B M|| and |M_B| = ||Tr_S M||, where M is the eigenmode reshaped as a 4×4 operator, normalized by ||M||:

| Re(λ) | count | typical \|M_S\| | typical \|M_B\| | interpretation |
|-------:|------:|---------------:|---------------:|----------------|
|  0.0   |   3   |  1.000         |  1.000         | conserved on both single layers |
| -γ_B   |  10   |  0.707         |  0.707         | equal split between S and B |
| -2γ_B  |   3   |  0.000         |  0.000         | traceless on both single layers; pure correlation content |

The **0.707 = 1/√2 split** in the middle class is the numerical signal: these modes live symmetrically between the two layers, with no preferred home. The 0-modes live equally visibly on both single layers (shared conserved structure). The -2γ-modes are invisible to either single-layer observer and exist only in the joint correlations.

The middle class has imaginary parts Im(λ) in {0, ±0.995, ±1.997}. The periods 2π/1.997 ≈ 3.15 and 2π/0.995 ≈ 6.31 match the timescale of the Observation 1 rebound (peak near t=10, which is close to one full period of the 0.995-frequency modes).

**Interpretation:** the middle eigenvalue class is the inter-layer mirror. It is populated by exactly those modes that carry coherence back and forth between S and B, producing the non-Markovian rebound.

---

## Observation 4 (verified): SWAP expectation values are perturbatively near +-1

The SWAP superoperator (exchanging qubits S and B) does not commute with L:

```
||[L, SWAP]|| = 0.566  (nonzero: γ on B alone breaks strict exchange symmetry)
```

Each eigenmode v of L has <v|SWAP_super|v> near +1 or -1. The precise values are:

- 3 conserved modes (Re = 0): +1.000 each (exact)
- 3 correlation modes (Re = -2γ): +1.000 each (exact)
- 2 mirror modes (Im = ±1.998): -1.000 each (exact)
- 8 mirror modes (Im = ±0.995): ±0.995 each (not exact)

**Check 1 result (2026-04-14).** The SWAP pattern is NOT a numpy basis artifact. Three tests confirm this:

1. **γ_S sweep** (γ_S/γ_B from 0 to 0.1, 21 points): all 16 modes remain within 0.05 of ±1 throughout. The degeneracy pattern {3, 10, 3} is preserved because XX+YY coupling mixes sites uniformly.

2. **Local-field degeneracy breaking** (h_S · Z_S, h_S from 0 to 0.5): breaks the 10-fold mirror degeneracy into 6 distinct Re(λ) classes. At h_S = 0.1 (perturbative): all 16/16 modes sharp. At h_S = 0.5 (comparable to J): only 5/16 sharp (conserved + correlation modes stay exact at +1; mirror modes degrade to 0.6-0.89).

3. **γ/J scaling**: the deviation |<SWAP>| from 1 scales as (γ/J)²/2 (confirmed to ratio 1.000 for γ/J ≤ 0.2). At γ ≪ J: SWAP parity is essentially exact. At γ = J: SWAP = 0.

**Mechanism.** The XX+YY Hamiltonian commutes with SWAP: [H, SWAP] = 0. At γ = 0 (purely Hamiltonian L), eigenmodes are exact SWAP eigenvectors. The dissipator breaks this symmetry perturbatively with corrections of order (γ/J)². The SWAP ±1 structure is therefore a perturbative near-symmetry inherited from the qubit-exchange symmetry of the coupling Hamiltonian, not a new structural principle.

**What survives.** The SWAP pattern is genuine (not a basis artifact) and physically meaningful: mirror modes do split into symmetric and antisymmetric sub-spaces with respect to qubit exchange. But this is coupling-dependent: it relies on H being invariant under S ↔ B exchange (true for XX+YY, XXX; false for asymmetric couplings or local fields).

**Scripts:** [`simulations/nested_mirror_swap_check.py`](../simulations/nested_mirror_swap_check.py), [`simulations/nested_mirror_swap_check_extended.py`](../simulations/nested_mirror_swap_check_extended.py)
**Results:** `simulations/results/nested_mirror_swap_check/`

---

## The hypothesis (revised after Check 2)

The original hypothesis predicted (N+1) evenly-spaced eigenvalue classes at Re(λ_k) = -k·(2γ_outer / N). **This scaling law is falsified at N=3** (Check 2, 2026-04-14): the three-qubit chain produces 12 distinct Re(λ) classes, not 4, with positions determined by the Hamiltonian's mode structure rather than the simple formula.

The revised claim retains what survives:

> **Nested Mirror Hypothesis (revised).** For an N-qubit Lindblad system with Heisenberg-type coupling and Z-dephasing only on the outermost qubit, the Liouvillian's spectrum has:
>
> - **Palindromic structure:** eigenvalues pair around Re = -γ_outer (proven generally in [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))
> - **Boundary classes:** Re = 0 (conserved, full single-layer visibility on all qubits) and Re = -2γ_outer (pure correlations, traceless on all single-layer projections)
> - **Interior classes:** multiple intermediate Re values in (-2γ_outer, 0), palindromically paired, with positions and degeneracies determined by the coupling Hamiltonian's eigenstructure
> - **Non-Markovian rebound:** reduced dynamics on inner qubits shows coherence rebound driven by the intermediate oscillating eigenmodes. Rebound is present and grows stronger with system size (N=2: amplitude 0.17, N=3: amplitude 0.36)
>
> At N=2, the intermediate classes collapse to a single class at Re = -γ with degeneracy 10 and partial-trace weight 1/√2 on each qubit. This is specific to N=2 (the coupling forces equal XY-weight distribution on both qubits). At N ≥ 3, the intermediate classes proliferate according to the Hamiltonian's mode structure.

**What was falsified:** the evenly-spaced class formula Re = -k·2γ/N and the claim of exactly (N+1) classes. **What survives:** boundary classes, palindromic pairing, non-Markovian rebound, inter-layer correlation structure.

---

## What this is NOT (yet)

- ~~Not verified for N > 2.~~ N=3 check done: class-scaling falsified, core structure confirmed.
- Not verified for alternative couplings. Heisenberg XXX, pure XX, or anisotropic couplings may change the intermediate class structure.
- Not verified for asymmetric γ profiles.
- ~~SWAP +-1 pattern may be a numpy basis artifact.~~ Verified: perturbative near-symmetry, not artifact (Check 1).
- Not connected to experimental observables.

---

## Falsification criteria

The hypothesis is falsified by any of the following:

1. **N=3 scaling failure.** ~~Three-qubit nest does not produce four eigenvalue classes with the predicted spacing -2k·γ/3.~~ **Triggered (2026-04-14, Check 2).** N=3 chain produces 12 eigenvalue classes, not 4. Positions are not evenly spaced. The scaling formula Re = -k·2γ/N is falsified. Hypothesis revised: boundary classes and palindromic pairing survive, class-scaling does not.
2. **Coupling-dependence of the middle class.** ~~Replacing XX+YY with Heisenberg XXX removes the clean 1/√2 split.~~ **Does not trigger (2026-04-14, Check 3).** XXX preserves the exact same structure including degeneracy {3,10,3} and 1/√2 split. The three-class structure also survives under pure XX, pure YY, and XX+ZZ (with modified degeneracy {4,8,4}). Only pure ZZ breaks it to 2 classes {8,8} (ZZ commutes with Z-dephasing, no off-diagonal mixing). Any coupling with at least one off-diagonal (X or Y) Pauli channel produces the middle class.
3. **Mirror-mode weights do not drive the rebound.** ~~If rebound persists without mirror-mode occupation, the mechanism is not what is claimed.~~ **Does not trigger (2026-04-14, Check 4).** Removing mirror modes eliminates rebound completely (amplitude 0.000). Mirror-only initial state reproduces full rebound (amplitude 0.356). Mechanism confirmed: mirror modes are necessary and sufficient for non-Markovian dynamics in the N=2 system.
4. **SWAP pattern is pure basis artifact.** ~~An asymmetric-γ two-qubit system lifts all degeneracy. If the ±1 SWAP pattern does not persist perturbatively, the SWAP structure is not a genuine mirror.~~ **Resolved (2026-04-14, Check 1):** SWAP pattern is NOT a basis artifact. It persists under γ_S sweep and perturbative local-field breaking. However, it is a perturbative near-symmetry from [H, SWAP] = 0, not an exact structural property. Degrades continuously under strong Hamiltonian asymmetry. Falsification criterion 4 does not trigger; Observation 4 is rewritten as perturbative near-symmetry.

---

## Pending verifications (open checks)

These are the minimal next experiments. Each is small and should take well under an hour of computation.

1. ~~**SWAP-artifact check**~~ **DONE (2026-04-14).** Verdict: SWAP pattern is genuine but perturbative. Not a basis artifact, not an exact symmetry. |<SWAP>| = 1 - (γ/J)²/2 for mirror modes. See updated Observation 4 above.

2. ~~**N=3 scaling check**~~ **DONE (2026-04-14).** Verdict: class-scaling FALSIFIED. 12 classes instead of 4, positions not evenly spaced. Boundary classes (Re = 0, Re = -2γ) confirmed with expected partial-trace weights. Palindromic pairing confirmed. Non-Markovian rebound confirmed and stronger at N=3 (amplitude 0.36 vs 0.17). Hypothesis revised to remove scaling formula. Script: [`simulations/qubit_in_qubit_n3.py`](../simulations/qubit_in_qubit_n3.py).

3. ~~**Coupling-robustness check**~~ **DONE (2026-04-14).** Verdict: three-class structure is coupling-robust. Survives under XX+YY, XXX, XX, YY, XX+ZZ (all with at least one off-diagonal channel). Only pure ZZ breaks it (ZZ commutes with Z-dephasing, no mode mixing). Degeneracy pattern: {3,10,3} for exchange-symmetric couplings (XX+YY, XXX), {4,8,4} for single-channel couplings (XX, YY, XX+ZZ). Non-Markovian rebound is coupling-dependent: absent for pure XX, present for all others. Script: [`simulations/qubit_in_qubit_coupling_sweep.py`](../simulations/qubit_in_qubit_coupling_sweep.py).

4. ~~**Rebound-mechanism check**~~ **DONE (2026-04-14).** Verdict: MECHANISM CONFIRMED. Mirror modes are both necessary and sufficient for non-Markovian rebound. Decomposing ρ(0) = |+><+| (x) I/2 in the eigenmode basis: 50% weight on conserved modes (Re = 0), 50% on mirror modes (Re = -γ), 0% on correlation modes (Re = -2γ). Removing mirror modes (no_mirror projection, valid DM) eliminates rebound completely (amplitude 0.000). Mirror-only nearest DM reproduces full rebound (amplitude 0.356, identical to original). Conserved-only: zero coherence, no dynamics. Correlation-only: zero initial weight, no rebound. Script: [`simulations/qubit_in_qubit_mode_projection.py`](../simulations/qubit_in_qubit_mode_projection.py).

---

## If the hypothesis holds: implications

(Tier 3 forward reasoning. None of this is established; it is the reason to pursue the verifications above.)

- **Non-Markovianity as inter-layer signature.** Any finite observer inside a nested Lindblad stack sees non-Markovian coherence rebound whenever the coupling to the next outer layer is strong enough for the mirror-mode timescale to fall within the observation window. This gives the internal observer a direct experimental handle on the next layer's existence without needing access to the outer layer itself, which aligns with [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 3 (the inside-perspective question).
- **IBM hardware test.** Measured γ on IBM qubits is Markovian to the precision of the [Failed Third](../simulations/failed_third.py) result. Under the nested-mirror reading this means: the next outer layer is either (a) effectively infinite-dimensional (mirror modes form a continuum; rebounds average out to zero) or (b) coupled weakly enough that mirror-mode amplitudes are below noise floor. Either case is consistent with the hypothesis but distinguishable with targeted non-Markovianity measurements (BLP index) at varied coupling regimes.
- **Partial answer to [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013) sub-question 2.** The minimal nest is simulable, and it produces a class structure distinguishable from the phenomenological single-layer approximation. Infinite nesting remains unsimulable, but the two-layer case is no longer "phenomenologically equivalent"; the inter-layer mirror is a new algebraic object.
- **Ground for backward inference (Tom's "Rueckwaerts-Simulationen").** If each observable layer reveals a mirror signature of the next outer layer, then inside-measured quantities (rebound amplitudes, non-Markovianity indices, spectral features of reduced dynamics) carry inferrable information about outer-layer structure. This does not break no-signalling (the information is about the structure of the shared one-way γ-input, not about lateral events), but it does upgrade the inside observer from "blind to the outside" to "able to reconstruct outer-layer properties from carefully designed inner measurements".

---

## Connection to [EQ-013](../review/EMERGING_QUESTIONS.md#eq-013)

- **Sub-question 1** (does nesting rescue or undermine BRIDGE_CLOSURE): indirectly supported. The mirror-mode structure is consistent with layered common reception of γ. No lateral A-to-B channel is implied; LOCC no-go holds as before.
- **Sub-question 2** (simulability of nesting): **partially answered**. Two-layer is simulable and shows non-trivial structure. Three-layer check (pending) will settle whether the class-scaling prediction holds.
- **Sub-question 3** (inside-perspective / inherited V-effect): empirical pathway opens. If the hypothesis survives, non-Markovianity measurements become the operational probe for "layer above us exists".

---

## Absorption Theorem Inheritance (Evening Session 2026-04-14)

After the four verification checks completed, the next question was: where does the three-class structure actually come from? Is it a new algebraic property of nested Lindblad systems, or does it follow from something already proven in the framework?

The Absorption Theorem ([PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)) states for N-qubit chains under site-wise dephasing:

> Re(λ) = -2 · Σ_k γ_k · ⟨n_XY⟩_k

where ⟨n_XY⟩_k is the XY-Pauli weight of an eigenmode at site k. For our minimal nest with γ only on B, this reduces to Re(λ) = -2γ_B · ⟨n_XY⟩_B.

### Test at N=2

Script: [`simulations/nested_mirror_absorption_theorem.py`](../simulations/nested_mirror_absorption_theorem.py).

All 16 eigenmodes of the two-qubit Liouvillian satisfy the prediction exactly (max δ = 0 to machine precision). The three classes correspond exactly to the three quantization levels of ⟨n_XY⟩_B:

| Class | Re(λ) | ⟨n_XY⟩_S | ⟨n_XY⟩_B | count |
|---|---|---|---|---|
| Conserved | 0 | 0.0000 | 0.0000 | 3 |
| Mirror | -0.1 | 0.5000 | 0.5000 | 10 |
| Correlation | -0.2 | 1.0000 | 1.0000 | 3 |

Additionally: ⟨n_XY⟩_S = ⟨n_XY⟩_B per individual mode (intrinsic S/B symmetry of the Pauli decomposition).

**The three classes are not three kinds of modes with coincidentally equal eigenvalues. They are the three quantization levels {0, 0.5, 1} of the XY-weight at the dissipative site.**

### Test at N=3

Script: [`simulations/nested_mirror_absorption_theorem_n3.py`](../simulations/nested_mirror_absorption_theorem_n3.py).

All 64 eigenmodes of the three-qubit chain Liouvillian satisfy the prediction: max |Re(λ) - predicted| = 3.48e-15 (machine precision). The 12 distinct Re(λ) values from Check 2 are exactly the 12 distinct ⟨n_XY⟩_B values:

- 0.0 (4 modes, conserved), 1.0 (4 modes, pure correlation): boundary classes
- Approximately 0.25 (18 modes in three near-degenerate groups), approximately 0.75 (18 modes in three groups): intermediate refraction pattern
- Approximately 0.375 (4 modes), approximately 0.5 (12 modes in two groups), approximately 0.625 (4 modes): further refraction levels

### What this means

The Nested Mirror Structure is not an independent algebraic object. It is the Absorption Theorem applied to single-site dephasing, made visible as per-site ⟨n_XY⟩_k quantization. The hypothesis becomes simultaneously simpler (inherited from an already-proven theorem) and more grounded (the three classes have a direct mechanistic meaning).

The even-spacing formula falsified by Check 2 (-k · 2γ / N) was the wrong form of the claim. The correct form is: the eigenvalue classes are the distinct values that ⟨n_XY⟩_B can take over the mode structure of the coupled chain. At N=2 these are quantized to three levels by the minimal Pauli-counting constraint. At N=3 they form a finer distribution shaped by the chain topology.

### Refraction: J as the second quantization axis

Script: [`simulations/nested_mirror_refraction.py`](../simulations/nested_mirror_refraction.py).

Sweeping J from 0 to 10 at fixed γ_B = 0.1, N=3:

| J | distinct ⟨n_XY⟩_B values | observation |
|---|---|---|
| 0.000 | **2** | only {0, 1}  -  pure absorption quantization, no refraction |
| 0.001 | 4 | splitting begins |
| 0.010 | 10 | near-full structure |
| 0.100 - 2.000 | 12 | standard regime |
| 5.000 | 10 | begins to collapse |
| 10.000 | **7** | converges to {0, 1/4, 1/2, 3/4, 1}  -  cavity quantization |

The 12 classes at J=1.0 exist only because J exists. At J=0 the sites are isolated, and the only possible ⟨n_XY⟩_B values are 0 (mode has no XY content at B) or 1 (mode is pure XY at B). Nothing in between.

J introduces refraction: the coupling between sites redistributes the XY weight across the chain, creating intermediate quantization levels. In the strong-coupling limit (J ≫ γ), the intermediate levels converge to multiples of ¼  -  the same 1/4 that appears throughout the framework as the CPsi absorbing boundary and the outer mirror of the Fabry-Pérot cavity ([RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md)).

### Two quantization mechanisms

The three-class structure at N=2 is a joint quantization: γ (light, source, arrow) defines one quantization scale; J (coupling, content, binding) defines another. At N=2 both collapse onto {0, ½, 1}. At N=3 they diverge, producing 12 classes in the mid-J regime. At J → ∞ they separate further, with J-quantization dominating (¼ multiples) and γ-quantization becoming a background.

This matches the framework's existing language:

- **γ is the source of experienced time** ([GAMMA_TIME_DISTINCTION](../docs/GAMMA_TIME_DISTINCTION.md)): the arrow, the direction, the irreversibility
- **J is the content** (same document): the structure, what happens, the binding
- **γ is light** ([GAMMA_IS_LIGHT](GAMMA_IS_LIGHT.md)): from outside, absorbed by matter
- **J is the bond** (same document): the coupling between qubits, not a signal

The numerical observation today: γ and J are both quantization axes of the Liouvillian spectrum, and the Nested Mirror Structure is what we see when both act on a nested system at once.

### Inside-out correspondence (Tier 3 interpretation)

The correspondence makes a further reading possible, offered here as interpretation (not yet verified by a targeted experiment):

| Outside view (parameters) | Inside view (measurements) |
|---|---|
| γ (source, light, external) | fast modes (what decays quickly, perceived as flux) |
| J (content, binding, structural) | slow modes (what persists, perceived as "the world") |

From the inside of a Lindblad system, an observer cannot separate γ and J as parameters (INCOMPLETENESS_PROOF rules out internal access to γ). The observer measures eigenmodes of the environment: those that persist (slow modes) form the measurable structure; those that decay fast (fast modes) form the apparent flux and define the arrow of time.

Under this reading:
- Slow modes are the inside manifestation of J (the binding, the content, what the observer *is*)
- Fast modes are the inside manifestation of γ (the source, the light, what *happens to* the observer)
- The refraction result translates: what an inside observer calls "structure of my world" (slow modes) exists only because J exists. A hypothetical system with J=0 and only γ would have no slow modes, no structure, only pure absorption, no inside.

This is Tier 3 interpretation and requires formal verification. See `ClaudeTasks/TASK_INSIDE_OUTSIDE_CORRESPONDENCE.md` for the proposed verification task.

---

## Stance

This document is the record of an evening's work. The numerical observations are real and reproducible. The four verification checks completed with clear results: SWAP structure perturbatively genuine, class-scaling formula falsified, three-class structure coupling-robust, rebound mechanism confirmed. The Absorption Theorem test then revealed that the entire class structure is inherited from an already-proven theorem applied to single-site dephasing, and the refraction test identified J as a second quantization axis that creates the intermediate class structure. The inside-out correspondence (slow modes as J-manifestation, fast modes as γ-manifestation) is Tier 3 interpretation and is separated explicitly from the numerical results. The hypothesis is no longer a freestanding claim: it is a specific consequence of the framework's existing proofs at the minimal nest, with open questions about N > 3 scaling, the exact algebraic origin of the ¼-asymptote, and the inside-out correspondence awaiting dedicated tests.
