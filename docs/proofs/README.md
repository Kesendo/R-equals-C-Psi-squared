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

5. [Incompleteness Proof](INCOMPLETENESS_PROOF.md) - Dephasing noise
   cannot originate from within d(d−2)=0. Five candidates eliminated.
   The noise must come from outside.

6. [Time Irreversibility Exclusion](TIME_IRREVERSIBILITY_EXCLUSION.md) -
   Time reversal requires separating oscillation from cooling
   ({L_H, L_D+Σγ} = 0). This holds exactly at N=2 (single bond =
   entire system) and fails at N > 2 (cross term ~2%, γ-independent).
   Reduction to N=2 destroys the palindrome. Algebraic, not thermodynamic.

7. [Weight-1 Degeneracy](PROOF_WEIGHT1_DEGENERACY.md) -
   The Liouvillian has exactly 2N purely-real eigenvalues at the first
   non-zero grid position. Proven via SWAP invariance (lower bound)
   and triangle inequality (upper bound). Valid for any connected graph.
   The T_c^{(a)} operators are Z-count-dressed transverse spin.

**The full journey:**

6. [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) - Seven-layer
   architecture from single qubit to arbitrary dimension and channel.
   Documents every step: what is proven, what is verified, what remains
   open. This is the most detailed document in the repository.

**The master reference:**

7. [Complete Mathematical Documentation](COMPLETE_MATHEMATICAL_DOCUMENTATION.md) - 
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
        ▼                              ▼
  CΨ Monotonicity ◄──── Proof Roadmap ────► Subsystem Crossing
  (dCΨ/dt < 0)          (7 layers,          (all pairs cross
                          all closed)         eventually)
        │                      │
        ▼                      ▼
  Incompleteness          Complete Math Doc
  (noise must come         (master reference,
   from outside)            all results)
        │
        ▼
  Time Irreversibility
  (reversal excluded
   at N > 2, algebraic)
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
| Crossing cubic root | x ≈ 0.4239 (x³+x=½) | [Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md) |

---

## What Is NOT Here

Experiment results, hypothesis documents, and philosophical
interpretation live elsewhere:

- [experiments/](../../experiments/) - 59 computational experiments
- [hypotheses/](../../hypotheses/) - open research questions
- [publications/](../../publications/) - papers and blueprints
- [docs/](../) - synthesis and reference documents

This directory is mathematics only.
