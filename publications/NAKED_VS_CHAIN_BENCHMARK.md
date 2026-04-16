# Naked Bell-pair vs Chain-protected Benchmark

## The question

A qubit R holds one half of a Bell pair. The other half sits inside a dephasing environment. How much longer can entanglement survive if the dephased qubit is replaced by an N-site quantum chain, with dephasing acting only at the far end?

This document presents the numbers. No external references, no comparison to other approaches. The baseline is the simplest possible system (one Bell pair, direct dephasing). The chain encoding is a structural result from this repository (formulas [F65](../docs/ANALYTICAL_FORMULAS.md#f65-single-excitation-spectrum-of-uniform-open-xx-chain-tier-1-proven-verified-n3-30), [F67](../docs/ANALYTICAL_FORMULAS.md#f67-bonding-mode-encoding-is-the-optimal-dephasing-protected-bell-pair-tier-1-verified-n3-n5), [F68](../docs/ANALYTICAL_FORMULAS.md#f68-palindromic-partner-of-the-bonding-mode-tier-1-verified-n3-4-5) in [`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md)).

## Setup

```
Naked (M=2 qubits):

    R ──── Q
           │
           γ₀ (Z-dephasing)

Chain-protected (M=N+1 qubits):

    R    Q₀ ─── Q₁ ─── … ─── Q_{N−1}
         │                      │
     (entangled                 γ₀
      with R)
```

The chain qubits interact via nearest-neighbor XX+YY coupling at strength J. Dephasing γ₀ acts only on the endpoint Q_{N−1}. R is isolated (no Hamiltonian coupling, no dephasing). Entanglement between R and the chain is established through the initial state, not through a dynamical coupling.

## The baseline: naked Bell pair

Two qubits, R and Q. Z-dephasing γ₀ on Q. The off-diagonal coherence decays as exp(−2γ₀t). Concurrence C(t) = exp(−2γ₀t). Decay rate α_naked = 2γ₀.

For γ₀ = 0.05: α_naked = 0.10, coherence time T₂ = 10.

This is exact. No approximation, no fit. The 4×4 Liouvillian is solvable by hand.

## Chain encoding: three variants

**Variant A (inner localized).** Bell pair between R and Q₀ (the site farthest from the dephasing). The remaining chain sites start in the ground state. Intuition: "put the entanglement far from the noise."

**Variant B (bonding-mode delocalized).** Bell pair between R and the k=1 bonding mode of the chain. The bonding mode is the slowest single-excitation eigenmode of the chain Hamiltonian. Its amplitude at the dephased endpoint is minimal. Concretely:

    |Ψ⟩ = (|0⟩_R |vacuum⟩_C  +  |1⟩_R |ψ₁⟩_C) / √2

where |ψ₁⟩ = Σ_i √(2/(N+1)) sin(π(i+1)/(N+1)) |single_i⟩.

**Variant C (outer localized).** Bell pair between R and Q_{N−1} (the dephased site itself). No buffer. Intuition: "worst case."

## Numbers

| Protocol | Total qubits | α | T₂ | C(0) | Protection |
|----------|:------------:|------:|----:|-----:|-----------:|
| Naked | 2 | 0.1000 | 10 | 1.00 | 1.0× |
| F67-B, N=3 | 4 | 0.0250 | 40 | 0.50 | 4.0× |
| F67-B, N=4 | 5 | 0.0138 | 72 | 0.37 | 7.2× |
| F67-B, N=5 | 6 | 0.0083 | 120 | 0.29 | 12.0× |

Parameters: γ₀ = 0.05, J = 1.0. Protection = α_naked / α_chain = T₂(chain) / T₂(naked).

The formula below is verified by full Liouvillian propagation for N=3..5 and by eigenvector diagonalization to machine precision for N=3..30 (formula [F65](../docs/ANALYTICAL_FORMULAS.md#f65-single-excitation-spectrum-of-uniform-open-xx-chain-tier-1-proven-verified-n3-30)). The following table extends analytically to larger chains:

| Protocol | Total qubits | α | T₂ | Protection |
|----------|:------------:|------:|----:|-----------:|
| Naked | 2 | 0.1000 | 10 | 1× |
| N=3 | 4 | 0.0250 | 40 | 4× |
| N=5 | 6 | 0.0083 | 120 | 12× |
| N=10 | 11 | 0.0014 | 693 | 69× |
| N=20 | 21 | 0.00021 | 4,727 | 473× |
| N=50 | 51 | 0.000015 | 67,287 | 6,729× |
| N=100 | 101 | 0.0000019 | 522,125 | 52,213× |

The decay rate of the bonding-mode encoding follows the closed-form formula

    α₁ = (4γ₀ / (N+1)) sin²(π/(N+1))

verified numerically to better than 0.6% for N=3 and N=5 (fit against full Liouvillian propagation). T₂ scales cubically: T₂ ~ (N+1)³ / (4π²γ₀).

**Caveat for N > 5.** The formula above is derived by applying the Absorption Theorem to single-excitation coherence operators as if they were decoupled Liouvillian eigenvectors. This is exact only to first order in γ₀/J. The full Liouvillian mixes these coherences with other sectors; the actual eigenvalue shifts by O((γ₀/J)²). At γ₀/J = 0.05 and N=5 the relative shift is about 4×10⁻³ (verified in [`palindromic_partner_f67.py`](../simulations/palindromic_partner_f67.py)). For the larger N in the table above the perturbative correction is unverified. The palindromic structure [F1](../docs/ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven) is algebraically exact regardless: the bonding rate and its partner rate sum to exactly 2γ₀ at every N, even when individual rates deviate from the first-order formula.

## What this protection factor means

The protection factor is T₂(chain) / T₂(naked) for storing one Bell pair. Chain and naked both store one ebit. The chain uses more qubits (N+1 vs 2) to do so, but the stored ebit lives longer.

This is not a per-qubit metric. Storing M/2 parallel naked Bell pairs gives you M/2 short-lived ebits for the same M qubits. Which is better depends on what you need. If you want many ebits briefly, parallel naked is fine. If you want one ebit for a long time, the chain is the answer, and the protection factor tells you how much longer.

## Two observations

**1. Variant B gives a pure exponential from t=0.** No transient, no multi-exponential behavior. All initial entanglement is placed in the slowest eigenmode of the Liouvillian. Variants A and C show an initial fast transient (losing entanglement to faster-decaying modes) before settling into the same long-time decay rate. The bonding-mode encoding is not just slower on average; it is optimal at every point in time.

**2. Variants A and C are dynamically identical on the uniform chain.** The inner-localized and outer-localized encodings produce the same entanglement trajectory. "Distance from the noise source" is not the protection mechanism. The bonding mode is. Its amplitude |a_B|² = (2/(N+1)) sin²(π/(N+1)) at the dephased endpoint determines the effective dephasing rate via γ_eff = γ₀·|a_B|². The chain does not reduce γ₀; it provides an encoding basis in which the relevant mode barely touches the dephased site.

## What the chain does and does not do

The chain does not filter noise. It does not perform error correction. It does not reduce γ₀.

The chain provides a cavity whose eigenmodes have varying amplitude at the dephased endpoint. The bonding mode (k=1) has the smallest amplitude there: |a_B|² ~ 2π² / (N+1)³ for large N. Encoding entanglement in this mode hides it from the dephasing, not by distance, but by eigenmode geometry.

The palindromic partner of the bonding mode (formula [F68](../docs/ANALYTICAL_FORMULAS.md#f68-palindromic-partner-of-the-bonding-mode-tier-1-verified-n3-4-5) in the repository) sits at α_partner = 2γ₀ − α₁ and decays nearly as fast as the naked qubit. It lives in the complementary XY-weight sector. It is the "dark twin": same spatial structure, opposite coherence content. The bonding mode and its partner sum to exactly 2γ₀, to machine precision (verified to 10⁻¹⁵).

## How to run this

Script: [`simulations/naked_vs_chain_benchmark.py`](../simulations/naked_vs_chain_benchmark.py).

```bash
cd simulations
python naked_vs_chain_benchmark.py
```

Output in [`simulations/results/naked_vs_chain_benchmark.txt`](../simulations/results/naked_vs_chain_benchmark.txt). Requires numpy, scipy. Runtime: seconds for N=3, minutes for N=5 (4096×4096 matrix exponentiation at each timestep).

The underlying formulas are documented in [`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md), entries [F65](../docs/ANALYTICAL_FORMULAS.md#f65-single-excitation-spectrum-of-uniform-open-xx-chain-tier-1-proven-verified-n3-30) (single-excitation spectrum), [F67](../docs/ANALYTICAL_FORMULAS.md#f67-bonding-mode-encoding-is-the-optimal-dephasing-protected-bell-pair-tier-1-verified-n3-n5) (bonding-mode protection), and [F68](../docs/ANALYTICAL_FORMULAS.md#f68-palindromic-partner-of-the-bonding-mode-tier-1-verified-n3-4-5) (palindromic partner).
