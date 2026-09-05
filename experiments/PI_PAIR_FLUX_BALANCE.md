# F1 Pair-Rate Balance Under a Bond Perturbation

**Status:** Current-state account of a finite `N=5` perturbation reading for ten
unambiguous matches. The full orbit counts are supplied by the later
`N=2...7` census.
**Date:** 2026-04-20; orbit language corrected 2026-09-06
**Authors:** Tom, Claude Opus 4.7 (1M)
**Relates to:** [the absorption theorem proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md),
[the palindromic orbit census](FACTOR_TWO_STANDING_WAVES.md), and
[the orthogonality-selection family](ORTHOGONALITY_SELECTION_FAMILY.md)

---

## What this run establishes

For the ten unambiguously continued `N=5` matches displayed by the producer, a
bond perturbation `delta J=0.01` preserves the complementary real-part sum to
machine precision:

```text
Re(lambda_slow)+Re(lambda_fast) = -2 Sigma_gamma = -0.5.
```

Equivalently, the inferred XY-weight sum remains five. “Flux” in the original
title was an optical reading of this algebraic balance; the calculation did
not measure transport between modes.

The exact F1 operation is linear:

```text
lambda -> -lambda-2 Sigma_gamma,
mu=lambda+Sigma_gamma -> -mu.
```

Hermiticity preservation separately gives conjugate closure. Their composition
`lambda -> -conj(lambda)-2 Sigma_gamma` is a different involution. In
particular, F1 is fixed only at `lambda=-Sigma_gamma`; the composite is fixed
on the full line `Re(lambda)=-Sigma_gamma`. The earlier analysis merged these
two orbit structures, so its self-pair counts and the proposed golden-ratio/
mod-10 mechanism are withdrawn rather than reinterpreted.

## Setup

- `N=5`, bond 0, `delta J=0.01`, `gamma_0=0.05`,
  `Sigma_gamma=0.25`.
- Diagonalize the uniform-chain generator `L_A` and the bond-perturbed
  generator `L_B+`.
- Identify complementary-rate partners in `L_A`, then continue candidates into
  `L_B+` by greedy nearest-neighbour matching.

The matching is reliable only for the ten isolated examples reported below.
Degeneracies make a global eigenvalue-by-eigenvalue continuation ambiguous.

## Ten unambiguous matches

```text
pair  Re(A_s)    Re(A_f)    Sum A       Re(B_s)    Re(B_f)    Sum B       delta sum
 1   +0.00000   -0.50000   -0.500000   +0.00000   -0.50000   -0.500000   +3.3e-15
 2   +0.00000   -0.50000   -0.500000   -0.00000   -0.50000   -0.500000   +4.8e-15
 3   -0.10000   -0.40000   -0.500000   -0.10000   -0.40000   -0.500000   -1.7e-16
 ...
10   -0.10000   -0.40000   -0.500000   -0.10000   -0.40000   -0.500000   +3.6e-16
```

The corresponding inferred XY weights sum to five before and after the
perturbation; their reported changes are between roughly `3e-15` and `5e-14`.
The absorption theorem supplies the relation between real part and XY weight.

The aggregate greedy match reports a maximum real-part-sum residual of
`1.4e-3` under `L_B+`. That number cannot distinguish a physical failure from
a matching artifact. A global perturbative statement would require
overlap-based continuation of invariant subspaces, including generalized
eigenspaces where blocks are defective. F1 fixes algebraic multiplicities and
maps the generalized eigenspace at `lambda` to the one at
`-lambda-2 Sigma_gamma`; it does not canonically pair basis vectors inside a
degenerate space.

## Correct full census

The later committed census across `N=2,...,7` supersedes the old `N=3...6`
pair table:

- Linear F1: 10,903 unordered two-member orbits plus 34 fixed eigenvalues at
  `lambda=-Sigma_gamma`.
- Conjugate-composite: 9,921 unordered two-member orbits plus 1,998 fixed
  eigenvalues on `Re(lambda)=-Sigma_gamma`.

Both totals include algebraic multiplicity and account for all 21,840
eigenvalues. They are orbit combinatorics, not evidence for binary inheritance,
a spatial standing wave, or an energy current.

## Reproduction and retained evidence

- `simulations/eq018_pi_pair_flow.py`
- `simulations/results/eq018_pi_pair_flow/pi_pair_flow.json`
- `simulations/results/eq018_pi_pair_flow/run.log`
- [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)
- [`simulations/results/factor_two_standing_waves.txt`](../simulations/results/factor_two_standing_waves.txt)

The historical parity and golden-ratio scripts remain reproducibility records
of the former classification, not support for a current theorem:

- `simulations/eq018_pi_parity_scan.py`
- `simulations/eq018_golden_ratio_check.py`
- `simulations/eq018_double_involution_scan.py`

---

*The established perturbative reading is the ten-pair rate sum. Any statement
about spatial waves, transported flux, or individual eigenvector partners needs
additional evidence.*
