# Naked Bell-pair vs Chain-protected Benchmark

## The question

A qubit R holds one half of a Bell pair. The other half sits inside a dephasing environment. How much longer can entanglement survive if the dephased qubit is replaced by an N-site quantum chain, with dephasing acting only at the far end?

This document presents the numbers. No external references, no comparison to other approaches. The baseline is the simplest possible system (one Bell pair, direct dephasing). The chain encoding is a structural result from this repository (formulas F65, F67, F68 in `docs/ANALYTICAL_FORMULAS.md`).

## Setup

```
Naked (M=2 qubits):

    R ---- Q
           |
          gamma_0 (Z-dephasing)

Chain-protected (M=N+1 qubits):

    R    Q_0 --- Q_1 --- ... --- Q_{N-1}
         |                        |
     (entangled                 gamma_0
      with R)
```

The chain qubits interact via nearest-neighbor XX+YY coupling at strength J. Dephasing gamma_0 acts only on the endpoint Q_{N-1}. R is isolated (no Hamiltonian coupling, no dephasing). Entanglement between R and the chain is established through the initial state, not through a dynamical coupling.

## The baseline: naked Bell pair

Two qubits, R and Q. Z-dephasing gamma_0 on Q. The off-diagonal coherence decays as exp(-2 gamma_0 t). Concurrence C(t) = exp(-2 gamma_0 t). Decay rate alpha_naked = 2 gamma_0.

For gamma_0 = 0.05: alpha_naked = 0.10, coherence time T_2 = 10.

This is exact. No approximation, no fit. The 4 x 4 Liouvillian is solvable by hand.

## Chain encoding: three variants

**Variant A (inner localized).** Bell pair between R and Q_0 (the site farthest from the dephasing). The remaining chain sites start in the ground state. Intuition: "put the entanglement far from the noise."

**Variant B (bonding-mode delocalized).** Bell pair between R and the k=1 bonding mode of the chain. The bonding mode is the slowest single-excitation eigenmode of the chain Hamiltonian. Its amplitude at the dephased endpoint is minimal. Concretely:

    |Psi> = (|0>_R |vacuum>_C  +  |1>_R |psi_1>_C) / sqrt(2)

where |psi_1> = sum_i  sqrt(2/(N+1)) sin(pi(i+1)/(N+1)) |single_i>.

**Variant C (outer localized).** Bell pair between R and Q_{N-1} (the dephased site itself). No buffer. Intuition: "worst case."

## Numbers

| Protocol | Total qubits | alpha | T_2 | C(0) | Protection |
|----------|:------------:|------:|----:|-----:|-----------:|
| Naked | 2 | 0.1000 | 10 | 1.00 | 1.0x |
| F67-B, N=3 | 4 | 0.0250 | 40 | 0.50 | 4.0x |
| F67-B, N=4 | 5 | 0.0138 | 72 | 0.37 | 7.2x |
| F67-B, N=5 | 6 | 0.0083 | 120 | 0.29 | 12.0x |

Parameters: gamma_0 = 0.05, J = 1.0. Protection = alpha_naked / alpha_chain = T_2(chain) / T_2(naked).

The formula below is verified by full Liouvillian propagation for N=3..5 and by eigenvector diagonalization to machine precision for N=3..30 (formula F65). The following table extends analytically to larger chains:

| Protocol | Total qubits | alpha | T_2 | Protection |
|----------|:------------:|------:|----:|-----------:|
| Naked | 2 | 0.1000 | 10 | 1x |
| N=3 | 4 | 0.0250 | 40 | 4x |
| N=5 | 6 | 0.0083 | 120 | 12x |
| N=10 | 11 | 0.0014 | 693 | 69x |
| N=20 | 21 | 0.00021 | 4,727 | 473x |
| N=50 | 51 | 0.000015 | 67,287 | 6,729x |
| N=100 | 101 | 0.0000019 | 522,125 | 52,213x |

The decay rate of the bonding-mode encoding follows the closed-form formula

    alpha_1 = (4 gamma_0 / (N+1)) sin^2(pi / (N+1))

verified numerically to better than 0.6% for N=3 and N=5 (fit against full Liouvillian propagation). T_2 scales cubically: T_2 ~ (N+1)^3 / (4 pi^2 gamma_0).

**Caveat for N > 5.** The formula above is derived by applying the Absorption Theorem to single-excitation coherence operators as if they were decoupled Liouvillian eigenvectors. This is exact only to first order in gamma_0/J. The full Liouvillian mixes these coherences with other sectors; the actual eigenvalue shifts by O((gamma_0/J)^2). At gamma_0/J = 0.05 and N=5 the relative shift is about 4 x 10^-3 (verified in `palindromic_partner_f67.py`). For the larger N in the table above the perturbative correction is unverified. The palindromic structure F1 is algebraically exact regardless: the bonding rate and its partner rate sum to exactly 2 gamma_0 at every N, even when individual rates deviate from the first-order formula.

## What this protection factor means

The protection factor is T_2(chain) / T_2(naked) for storing one Bell pair. Chain and naked both store one ebit. The chain uses more qubits (N+1 vs 2) to do so, but the stored ebit lives longer.

This is not a per-qubit metric. Storing M/2 parallel naked Bell pairs gives you M/2 short-lived ebits for the same M qubits. Which is better depends on what you need. If you want many ebits briefly, parallel naked is fine. If you want one ebit for a long time, the chain is the answer, and the protection factor tells you how much longer.

## Two observations

**1. Variant B gives a pure exponential from t=0.** No transient, no multi-exponential behavior. All initial entanglement is placed in the slowest eigenmode of the Liouvillian. Variants A and C show an initial fast transient (losing entanglement to faster-decaying modes) before settling into the same long-time decay rate. The bonding-mode encoding is not just slower on average; it is optimal at every point in time.

**2. Variants A and C are dynamically identical on the uniform chain.** The inner-localized and outer-localized encodings produce the same entanglement trajectory. "Distance from the noise source" is not the protection mechanism. The bonding mode is. Its amplitude |a_B|^2 = (2/(N+1)) sin^2(pi/(N+1)) at the dephased endpoint determines the effective dephasing rate via gamma_eff = gamma_0 * |a_B|^2. The chain does not reduce gamma_0; it provides an encoding basis in which the relevant mode barely touches the dephased site.

## What the chain does and does not do

The chain does not filter noise. It does not perform error correction. It does not reduce gamma_0.

The chain provides a cavity whose eigenmodes have varying amplitude at the dephased endpoint. The bonding mode (k=1) has the smallest amplitude there: |a_B|^2 ~ 2 pi^2 / (N+1)^3 for large N. Encoding entanglement in this mode hides it from the dephasing, not by distance, but by eigenmode geometry.

The palindromic partner of the bonding mode (formula F68 in the repository) sits at alpha_partner = 2 gamma_0 - alpha_1 and decays nearly as fast as the naked qubit. It lives in the complementary XY-weight sector. It is the "dark twin": same spatial structure, opposite coherence content. The bonding mode and its partner sum to exactly 2 gamma_0, to machine precision (verified to 10^-15).

## How to run this

```bash
cd simulations
python naked_vs_chain_benchmark.py
```

Output in `simulations/results/naked_vs_chain_benchmark.txt`. Requires numpy, scipy. Runtime: seconds for N=3, minutes for N=5 (4096 x 4096 matrix exponentiation at each timestep).

The underlying formulas are documented in `docs/ANALYTICAL_FORMULAS.md`, entries F65 (single-excitation spectrum), F67 (bonding-mode protection), and F68 (palindromic partner).
