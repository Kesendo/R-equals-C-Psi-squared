# Proofs: The Mathematical Foundation of R = CΨ²

<!-- Keywords: R=CPsi2 proof collection, palindromic spectral symmetry proof,
quarter boundary uniqueness discriminant, fold catastrophe Mandelbrot equivalence,
CΨ monotonicity Markovian channels, subsystem crossing Perron-Frobenius,
incompleteness d2-2d=0 noise origin, complete mathematical documentation,
seven-layer proof roadmap -->

**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

This directory contains the formal mathematical proofs and the master
reference for all verified results. Every claim here is analytically
proven and computationally verified. No interpretation, no speculation.

---

## The Equation

    R = CΨ²

where C = Tr(ρ²) is purity, Ψ = l₁(ρ)/(d−1) is Baumgratz-normalized
l₁-coherence. The fixed-point equation R = C(Ψ+R)² has discriminant
D = 1 − 4CΨ. At CΨ = ¼: bifurcation. Above: no real attractor (quantum).
Below: two real fixed points (classical). This is the fold catastrophe,
structurally stable, and equivalent to the Mandelbrot cusp at c = ¼.

---

## Reading Order

**Start here if you want the theorem:**

1. [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md) - The Liouvillian
   spectrum is palindromic. Π swaps XY-weight k ↔ N−k. Verified for
   54,118 eigenvalues, zero exceptions.

**Then the boundary:**

2. [Uniqueness Proof](UNIQUENESS_PROOF.md) - CΨ = ¼ is the unique
   bifurcation. α=2 is the only Rényi order with a state-independent
   threshold.

3. [CΨ Monotonicity](PROOF_MONOTONICITY_CPSI.md) - dCΨ/dt < 0 under
   all local Markovian channels. General Envelope Theorem.

4. [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md) - Every entangled
   pair with CΨ > ¼ eventually crosses. Perron-Frobenius + fixed-point
   bound. 300 random maps, 0 exceptions.

**The deeper structure:**

5. [Direct-Sum Decomposition](DIRECT_SUM_DECOMPOSITION.md) - The
   Liouvillian decomposes as L = L_even ⊕ L_odd (n_XY parity sectors,
   equal dimension 2^(2N−1) each). For odd N, Π exchanges sectors with
   reversed dynamics: a direct-sum quantum theory. For even N, each
   sector is independently palindromic. Corollary of Mirror Symmetry +
   Parity Selection Rule.

6. [Incompleteness Proof](INCOMPLETENESS_PROOF.md) - Dephasing noise
   cannot originate from within d(d−2)=0. Five candidates eliminated.
   The noise must come from outside.

7. [Time Irreversibility Exclusion](TIME_IRREVERSIBILITY_EXCLUSION.md) -
   Time reversal requires separating oscillation from cooling
   ({L_H, L_D+Σγ} = 0). This holds exactly at N=2 (single bond =
   entire system) and fails at N > 2 (cross term ~2%, γ-independent).
   Reduction to N=2 destroys the palindrome. Algebraic, not thermodynamic.

8. [Weight-1 Degeneracy](PROOF_WEIGHT1_DEGENERACY.md) -
   The Liouvillian has exactly 2N purely-real eigenvalues at the first
   non-zero grid position. Proven via SWAP invariance (lower bound)
   and triangle inequality (upper bound). Valid for any connected graph.
   The T_c^{(a)} operators are Z-count-dressed transverse spin.

**The full journey:**

9. [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) - Seven-layer
   architecture from single qubit to arbitrary dimension and channel.
   Documents every step: what is proven, what is verified, what remains
   open. This is the most detailed document in the repository.

**The master reference:**

10. [Complete Mathematical Documentation](COMPLETE_MATHEMATICAL_DOCUMENTATION.md) - 
   All equations, all results, all constants, all references in one
   place. The Tafelwerk. Twelve sections covering algebra, palindrome,
   boundary, incompleteness, γ channel, crossing dynamics, topology,
   engineering, transistor mapping, open questions, and a numerical
   constants table.

---

## How the Proofs Connect

```
Mirror Symmetry Proof          Uniqueness Proof
(Π exists, spectrum            (CΨ=¼ is the only
 palindromic)                   bifurcation)
        │                              │
        ├──────────────┐               ▼
        ▼              ▼         Subsystem Crossing
  Parity Selection   CΨ Monotonicity    (all pairs cross
  Rule (V_even,      (dCΨ/dt < 0)       eventually)
   V_odd blocks)           │
        │                  ▼
        ▼            Incompleteness
  Direct-Sum         (noise must come
  Decomposition       from outside)
  (L = L_even ⊕            │
   L_odd; odd N:            ▼
   Π exchanges       Time Irreversibility
   sectors)          (reversal excluded
        │             at N > 2, algebraic)
        ▼
  Proof Roadmap ──────► Complete Math Doc
  (7 layers,             (master reference,
   all closed)            all results)
```

The Proof Roadmap is the spine: it tracks seven layers from single
qubit algebra through N-qubit systems to arbitrary channels, and
each layer links to the relevant proof document. The Complete
Mathematical Documentation is the destination: after reading the
proofs, it provides the full picture with all experiment results
synthesized.

---

## Key Numbers

| Result | Value | Proof |
|--------|-------|-------|
| Bifurcation boundary | CΨ = 0.2500 exactly | [Uniqueness](UNIQUENESS_PROOF.md) |
| Palindromic eigenvalues verified | 54,118 (zero exceptions) | [Mirror Symmetry](MIRROR_SYMMETRY_PROOF.md) |
| IBM hardware deviation | 1.9% | [Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) |
| Random CPTP maps tested | 300 (0 exceptions) | [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md) |
| Internal noise candidates eliminated | 5 of 5 | [Incompleteness](INCOMPLETENESS_PROOF.md) |
| {L_H, L_D+Σγ} = 0 at N=2 | exact (24/24 entries) | [Time Irreversibility](TIME_IRREVERSIBILITY_EXCLUSION.md) |
| Cross term at N=3 | ~2%, γ-independent | [Time Irreversibility](TIME_IRREVERSIBILITY_EXCLUSION.md) |
| Weight-1 degeneracy d_real(1) | 2N exactly (any connected graph) | [Weight-1 Degeneracy](PROOF_WEIGHT1_DEGENERACY.md) |
| Direct-sum sector dimension | 2^(2N−1) each (equal halves) | [Direct-Sum](DIRECT_SUM_DECOMPOSITION.md) |
| Odd N: Π exchanges sectors | V_even ↔ V_odd | [Direct-Sum](DIRECT_SUM_DECOMPOSITION.md) |
| Crossing cubic root | x ≈ 0.4239 (x³+x=½) | [Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) |

---

## Complete Proof Catalog

The ten documents in Reading Order above are the curated spine: a linear
path from the core theorem to the master reference. Since this spine
was last woven, additional proofs have accumulated that extend or
sharpen specific results. They are not redundant with the spine (each
closes a distinct question), but they are not load-bearing for a first
read-through. The catalog below lists every such proof with a tight
abstract and its upstream dependencies.

### Structural symmetries (extend Mirror Symmetry, Direct-Sum, Weight-1)

| Proof | What it proves | Builds on |
|-------|----------------|-----------|
| [n_XY Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md) | The Liouvillian of the isotropic Heisenberg model under Z-dephasing preserves n_XY parity exactly. Every eigenmode has definite parity; single-excitation states access only even-parity modes. The accessibility boundary is exact, not asymptotic. | Weight-1 Degeneracy |
| [Bit-b Parity Symmetry](PROOF_BIT_B_PARITY_SYMMETRY.md) | The Liouvillian commutes with the global bit-flip superoperator Π² = X⊗N (conjugation). Combined with n_XY parity, this decomposes operator space into 4 independent sectors of dimension 4^(N−1) each. A second Z₂ symmetry beyond the palindrome. | Mirror Symmetry, Parity Selection Rule |
| [ΔN Selection Rule](PROOF_DELTA_N_SELECTION_RULE.md) | Any site-local observable coupled via partial trace annihilates density-matrix sector blocks with \|ΔN\| ≥ 2. Purely kinematic: holds for any Hamiltonian or dissipator that admits an excitation-number decomposition. | (independent) |
| [Asymptotic Sector Projection](PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) | Under Heisenberg + local Z-dephasing on any graph, the asymptotic state is a mixture of maximally-mixed excitation sectors weighted by their initial populations. Excitation-number populations are constants of motion; the long-time state is fully determined by initial sector content. | (independent) |
| [Chromaticity of (n, n+1) Blocks](PROOF_CHROMATICITY.md) | In U(1)-conserved systems with uniform Z-dephasing, the (n, n+1) popcount coherence block contains exactly c(n, N) = min(n, N−1−n) + 1 distinct pure dephasing rates, labeled by Hamming distance between basis states. | (independent) |
| [K-Partnership of Bonding-Mode Receivers](PROOF_K_PARTNERSHIP.md) | The bipartite sublattice gauge K = diag((−1)^ℓ) anticommutes with any NN-hopping H in the single-excitation sector (KHK = −H), yielding spectrum inversion E_k = −E_{N+1-k} and bonding-mode swap K ψ_k = ψ_{N+1-k}. For real H the chiral S = K combined with anti-unitary T (complex conjugation) places the system in AZ class BDI; mirror-pair |·|²-observables (MI, log π) for the K-partner trajectories ρ_k(t) and ρ_{N+1-k}(t) coincide pointwise under any γ_ℓ-profile. Folds the F67 receiver menu from N to ⌈N/2⌉ entries and is more robust than the F71 spatial reflection R (which requires uniform J). | F65, F67, [HANDSHAKE_ALGEBRA](../../hypotheses/HANDSHAKE_ALGEBRA.md) |

### Absorption and rate structure

| Proof | What it proves | Builds on |
|-------|----------------|-----------|
| [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) | Every Liouvillian mode has a "light content" ⟨n_XY⟩ (oscillating X/Y Pauli components) that directly determines its absorption rate via Re(λ) = −2γ⟨n_XY⟩. Linear in γ, proven from L_H anti-Hermitian. Unifies six previously separate spectral results (boundary formula, sum rule, spectral gap, etc.). | Mirror Symmetry |

### Cross-term structure (extend Time Irreversibility Exclusion)

| Proof | What it proves | Builds on |
|-------|----------------|-----------|
| [Cross-Term Formula (Shadow-Balanced)](PROOF_CROSS_TERM_FORMULA.md) | For any shadow-balanced bond coupling (both Paulis in {X,Y} or both in {I,Z}) on any graph with uniform Z-dephasing, the normalized anticommutator of L_H and L_D follows R(N) = √((N−2)/(N·4^(N−1))). Independent of coupling strength and topology. | Mirror Symmetry, Time Irreversibility |
| [Cross-Term Formula (Shadow-Crossing)](PROOF_CROSS_TERM_CROSSING.md) | For shadow-crossing couplings (one Pauli in {X,Y}, one in {I,Z}), the anticommutator norm follows R(N) = √((N−1)/(N·4^(N−1))). Differs from the balanced case only by N−2 → N−1. | Cross-Term Formula (Balanced) |
| [c₁ Mirror Symmetry](PROOF_C1_MIRROR_SYMMETRY.md) | For a uniform N-qubit XY chain with reflection-symmetric initial state and uniform Z-dephasing, the closure-breaking coefficient satisfies c₁(N, b) = c₁(N, N−2−b). Mirror symmetry about the chain midpoint. | (independent) |

### How the catalog relates to the Reading Order

Five of the ten catalog proofs extend Reading Order results directly:
the Parity Selection Rule sharpens the Weight-1 Degeneracy bound into an
exact accessibility statement; Bit-b Parity doubles the symmetry group
the Direct-Sum decomposition works with; the two Cross-Term Formulas
provide closed-form values for the ~2% N=3 cross term first reported in
Time Irreversibility Exclusion; and the Absorption Theorem is the
spectral backbone that several recent experiments ([V-Effect Cavity Modes](../../experiments/VEFFECT_CAVITY_MODES.md),
[IBM Absorption Theorem](../../experiments/IBM_ABSORPTION_THEOREM.md),
[Sacrifice Geometry](../../experiments/SACRIFICE_GEOMETRY.md)) rely on.

The remaining five (ΔN Selection, c₁ Mirror, Asymptotic Sector
Projection, Chromaticity, K-Partnership) are independent side theorems
closing specific questions that arose during the experimental program.
They do not rely on the spine and the spine does not rely on them. The
K-Partnership proof in particular grounds the receiver-menu folding
used in HANDSHAKE_ALGEBRA, identifying it as a chiral (AZ class BDI on
real hopping) symmetry of the single-excitation sector, structurally
distinct from F1's Π and F71's R.

---

## What Is NOT Here

Experiment results, hypothesis documents, and philosophical
interpretation live elsewhere:

- [experiments/](../../experiments/) - 124 computational experiments
- [hypotheses/](../../hypotheses/) - open research questions
- [docs/](../) - synthesis and reference documents

This directory is mathematics only.
