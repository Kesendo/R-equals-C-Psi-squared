# Z⊗N-Partnership: Multi-Excitation Néel-Mirror as a Transverse-Field Diagnostic

**Status:** Simulation only (Aer/numpy). Hardware sketch open. Sister-finding to K-partnership (single-excitation): different sector, different broken-by-rule.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `simulations/_neel_mirror_test.py`
**See also:** [IBM_K_PARTNERSHIP_SKETCH](IBM_K_PARTNERSHIP_SKETCH.md), [`docs/proofs/PROOF_K_PARTNERSHIP.md`](../docs/proofs/PROOF_K_PARTNERSHIP.md), [GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md), [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md)

---

## What this tests

Yesterday's K-partnership reads the **single-excitation** sector of the Heisenberg-XXZ chain via the bipartite sublattice gauge K = diag((-1)^l). This finding looks at the **multi-excitation** sector via the X-basis antiferromagnetic Néel state and its Z⊗N-mirror.

Initial states:

```
|N⟩       = |+−+−+⟩       (X-basis Néel)
|N̄⟩       = |−+−+−⟩       (Z⊗N-mirror, what a global π/2 Z-pulse produces)
```

These are Z⊗N-partners: (Z⊗N) |+⟩ = |−⟩, (Z⊗N) |−⟩ = |+⟩.

The question: does the symmetry Z⊗N survive the Heisenberg-Lindblad dynamics, and what physical perturbations break it?

## Theorem (Z⊗N-Partnership)

For an N-qubit Lindblad evolution `dρ/dt = −i[H, ρ] + Σ_k γ_k D[L_k] ρ`, if every term in H and every jump operator L_k has an **even number of X/Y operators per term**, then Z⊗N = Π_l Z_l is a strong symmetry. Consequence: if ρ_b(0) = (Z⊗N) ρ_a(0) (Z⊗N), then ρ_b(t) = (Z⊗N) ρ_a(t) (Z⊗N) for all t.

Proof sketch: Z⊗N anti-commutes with X_l, Y_l and commutes with Z_l. A term with k transverse operators (X or Y) picks up (−1)^k under conjugation. Even k preserves the term; T1 dissipator D[σ⁻] picks up two transverse operators per Lindblad action (one in σ⁻, one in σ⁺), so signs cancel.

**Z⊗N is preserved by:**

- XXZ Heisenberg: H = (J/2) Σ (X_iX_{i+1} + Y_iY_{i+1} + Δ Z_iZ_{i+1}); all 2 transverse per term.
- Z-dephasing: D[Z_l]; 0 transverse.
- Z-detuning: Σ_l δ_l Z_l, **uniform or non-uniform**; 0 transverse.
- T1 amplitude damping: D[σ⁻_l]; σ⁻ has 1 transverse, but the dissipator action σ⁻ρσ⁺ has 2 transverse total; signs cancel.
- σ⁻σ⁺-paired Lindblad operators in general (collective decay, σ⁻-pumping).

**Z⊗N is broken by:**

- Single transverse field: H = Σ_l h_l X_l or Σ_l h_l Y_l; 1 transverse.
- σ_z⊗σ_x-type 2-body terms: ZX, ZY; 1 transverse.

## Observables

For Z⊗N-invariant operator O: `⟨O⟩_a = ⟨O⟩_b` always.
For Z⊗N-anti-invariant operator O: `⟨O⟩_a = −⟨O⟩_b` always.

Two probes:

- **MI(0, N−1)**: invariant under Z⊗N (mutual information is unitary-invariant on each subsystem). Symmetry-break ⇒ MI_a ≠ MI_b.
- **M_X = (1/N) Σ_l (−1)^l ⟨X_l⟩**: anti-invariant. M_X(b) = −M_X(a) when symmetry holds; M_X(a) + M_X(b) ≠ 0 measures the break.

The Z-basis Néel order M_AB = (1/N) Σ_l (−1)^l ⟨Z_l⟩ is **identically zero** on X-basis states and Z⊗N-invariant: wrong observable, do not use.

## Simulation results (N=5, t=3.0, γ_z=0.1, J=1, Δ=1)

| Test | Hamiltonian / Lindblad | max|ΔMI| | max|M_X(a)+M_X(b)| | Z⊗N? |
|------|------------------------|----------|--------------------|------|
| 1 | XXZ + Z-dephasing | 0.000e+00 | 0.000e+00 | preserved |
| 2a | + T1 (γ_T1/γ_z = 0.5) | 0.000e+00 | 0.000e+00 | preserved |
| 2b | + T1 (γ_T1/γ_z = 1.0) | 0.000e+00 | 0.000e+00 | preserved |
| 2c | + T1 (γ_T1/γ_z = 2.0) | 0.000e+00 | 0.000e+00 | preserved |
| 3 | + non-uniform Z-detuning δ_l ∈ [−0.3, 0.3] | 0.000e+00 | 0.000e+00 | preserved |
| 4a | + uniform X-field h_x = 0.05 | 1.98e-03 | 4.64e-04 | **broken** |
| 4b | + uniform X-field h_x = 0.10 | 3.93e-03 | 9.27e-04 | **broken** |
| 4c | + uniform X-field h_x = 0.20 | 7.61e-03 | 1.85e-03 | **broken** |
| 5 | + non-uniform X-field h_l ∈ [−0.1, 0.1] | 3.11e-02 | 3.06e-02 | **broken** |

Initial M_X verifies the partner identity: M_X(t=0) = +1.0000 for |+−+−+⟩, −1.0000 for |−+−+−⟩.

Test 4 confirms first-order scaling: |ΔMI| ≈ 0.04 · h_x and |M_X(a)+M_X(b)| ≈ 0.009 · h_x at t_max = 3.0. Linear in h_x, as expected from first-order perturbation theory in the transverse field.

Test 5 vs Test 4b: at comparable max-amplitude (h_l max ≈ 0.10), the non-uniform profile gives ~33× larger |M_X-sum| break. The non-uniformity adds a second mode of breaking; the K-spread reading in single-exc has the same structure (uniform γ preserves K-partnership exactly; non-uniform breaks it).

## Why this matters: a third orthogonal Hardware-Diagnostic channel

Same chain, three independent symmetry channels reading three different asymmetries:

| Channel | Sector | Observable | Broken by |
|---------|--------|------------|-----------|
| K-partnership (PROOF_K_PARTNERSHIP) | single-exc | bonding-mode k vs N+1−k | non-uniform γ_l, asymmetric J_l, ZZ-boundary |
| (vac, S_1)-CMRR (CMRR_BREAK_NONUNIFORM_GAMMA) | single-exc cross-cohärence | (vac, S_1) coherence kernel | non-uniform γ_l (longitudinal) |
| **Z⊗N-partnership (this)** | **multi-exc Néel** | **MI(0,N−1) and M_X** | **transverse field / X-Y crosstalk** |

K and CMRR both read longitudinal (Z) asymmetry. Z⊗N reads transverse asymmetry, a different physical channel:

- Cross-resonance gate calibration drift
- ZX/ZY two-qubit crosstalk (entangling-gate residue)
- Off-resonant drive bleed-through
- Single-qubit X/Y-rotation calibration error

These are **calibration-side** asymmetries, distinct from coherence-time and dephasing asymmetries that K reads.

## Hardware sketch (open)

A minimum hardware test on IBM Heron-class:

- **Initial state:** prepare |+−+−+⟩ via a Hadamard chain with alternating Z-flip; mirror via Hadamards on opposite-sublattice qubits.
- **Evolution:** uniform J=1, t=0.8, 3 Trotter steps, XX+YY+ZZ Heisenberg (full, no `--xx-only`).
- **Tomography:** 9 Pauli settings on (qubit 0, qubit N−1) for MI; full single-qubit X-readout for M_X.
- **Cost:** ~5 QPU minutes, comparable to the K-partnership run on Marrakesh.
- **Prediction:** Z⊗N-break visible above shot-noise. Hardware ZX/ZY-crosstalk should produce |M_X(a) + M_X(b)| ≈ 10^−2 to 10^−1 (cf. simulation Test 5).
- **Diagnostic value:** the size of |M_X(a) + M_X(b)| is a direct readout of accumulated transverse error per gate-cycle, sensitive in a way that K-partnership and γ-profile reads cannot detect.

## Caveats

- N=5 is small. The breaking-amplitude scaling with N is not yet measured; multi-excitation states have richer Hilbert space, but the leading-order |⟨O⟩_a + ⟨O⟩_b|-scaling is set by the perturbation strength, not N.
- Single-X, single-Y are the only Hamiltonian generators with one transverse operator. ZX, ZY, ZZX, etc. also break Z⊗N; not tested here.
- The simulation runs at γ_z = 0.1, J = 1; the break-magnitude per unit transverse field has a γ-dependence and t-dependence not yet characterised.
- Z⊗N is a **strong symmetry** in the sense of dissipative-symmetry classification (Buča-Prosen 2012): it is preserved by both H and {L_k}. Weak symmetries (preserved only by total Lindbladian, not individual jump operators) are not tested here.

## What this run does and does not establish

**Establishes:**

- Z⊗N is preserved by every standard ingredient of NN-Heisenberg + decoherence (T1, T2, Z-detuning), uniform or non-uniform.
- Z⊗N breaks first-order in the transverse-field strength.
- Two probes (MI Z⊗N-invariant; M_X Z⊗N-anti-invariant) detect the break consistently.

**Does not establish:**

- Hardware visibility (no IBM run yet).
- A new physical phenomenon. The Z⊗N argument is a clean restatement of the trivial fact that Lindblad dynamics commute with their global symmetries; the value is in the **observability table** for hardware.

## References

- [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md): single-excitation K-partnership, the original sister-result.
- [IBM_K_PARTNERSHIP_SKETCH](IBM_K_PARTNERSHIP_SKETCH.md): hardware K-partnership on Marrakesh, 2026-04-25.
- [GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md), [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md): companion γ-profile readings via spatial-sum and (vac, S_1)-coherence kernels.
- Buča, B. & Prosen, T. (2012), "A note on symmetry reductions of the Lindblad equation: transport in constrained open spin chains"; strong/weak symmetry classification used here.
- Simulation script: `simulations/_neel_mirror_test.py`.
