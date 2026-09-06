# The Sum-versus-Product Question

**Status:** Null result: the Liouvillian palindrome does not select the sum
**Date:** March 14, 2026

## Question

Two candidate expressions were compared:

```text
R_sum  = C * (Psi_A + Psi_B)^2
R_prod = C * Psi_A * Psi_B.
```

The question was whether information conservation or the palindromic
Liouvillian spectrum forces `R_sum` and excludes `R_prod`.

## Verdict

No such derivation is present. Both expressions are invariant under exchanging
the labels `A` and `B`. The affine eigenvalue pairing established by the
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) constrains the
spectrum; it does not pair observable amplitudes, identify the labels with
forward/backward spatial propagation, or preserve the cross-term in `R_sum`.

Closed-system unitarity cannot be used by itself to forbid information loss in
the reduced state of an open Lindblad system. Algebraically, a sum of amplitudes
can also vanish by cancellation unless additional positivity and phase
conditions are imposed. Therefore neither information conservation nor the
non-vanishing of `R_sum` follows from the displayed formula alone.

## What remains

The sum-squared expression may be studied as a chosen model or operational
ansatz, but the current evidence does not make it unique. Claims of a physical
handoff at the `1/4` boundary, a standing-wave mechanism, a black-hole analogy,
Born-rule deviations, or the product's incompatibility with unitarity require
separate derivations with a specified state, channel, observable, and phase
convention.
