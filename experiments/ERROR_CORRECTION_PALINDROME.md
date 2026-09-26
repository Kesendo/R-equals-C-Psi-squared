# Error Protection from Palindromic Spectral Structure: Can the Mirror Shield Quantum Information?

<!-- Keywords: palindromic error protection quantum, three tier decay hierarchy,
optimal state dephasing survival, standing wave error syndrome, XOR drain GHZ
fragility, Pi operator Z4 structure, palindromic mode protection hierarchy,
quantum error correction dephasing, boundary tier palindromic pairs, slow mode
weight quantum state optimization, R=CPsi2 error correction palindrome -->

**Status:** Negative result. The state-weight, protection-tier,
standing-wave-syndrome, and information-lifetime interpretations are not
established; F22's operator-level dephasing charge remains.
**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Π as Time Reversal](PI_AS_TIME_REVERSAL.md),
[Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md), [XOR Space](XOR_SPACE.md)

The numbered sections preserve the March 19 event and its measured tables.
Modal percentages are historical coordinates from a non-normal
right-eigenvector basis, not physical state weights. The current verdict is the
Abstract and Summary; interpretive prose inside the event is not current
evidence.

---

## What this document is about

Quantum error correction is one of the hardest problems in quantum
computing: how do you protect fragile quantum information from noise?
Standard approaches treat noise as a generic enemy and build defenses
without looking at the specific structure of how things decay.

The palindromic symmetry gives us something standard approaches lack:
a detailed map of the decay landscape. Every mode has a partner, and
these pairs fall into three natural tiers (like express, local, and
freight trains on the same track). This document asks whether that
map can be exploited to protect information better.

The event tested an eigenvector-coordinate optimization, a projected-component
fingerprint, and a modal-pair decay curve. None supplies an operational
state-protection or error-syndrome result. The constructed Π representation is
fourth-order, but its eigenspaces are not thereby physical codespaces.

---

## Abstract

The rate spectrum contains steady, boundary, and mid-spectrum eigenvalues, but
that ordering alone is not a protection hierarchy for states. The former
constrained optimization used squared coordinates in a right-eigenvector basis
of a non-normal Liouvillian as if they were probabilities. Its `90% slow-mode`
and associated comparison are withdrawn; the returned state must instead be
assessed by direct density-matrix propagation and a named operational metric.
Projected onto the oscillating spectral subspace, Bell(0,1) changes its
largest Pauli coefficient by 0.0278 under an X or Z error and by 0.0556 under a
Y error on sites 0 and 1, and by exactly 0 under any error on site 2, which acts
on this state as X^⊗3; a projected change that no measurement reads is not a
syndrome. The Π representation has fourth-order structure
(`Π⁴ = I`) and four 16-dimensional algebraic eigenspaces; no QEC role follows
for those sectors.

---

## The Question

The palindromic symmetry pairs decay modes. F22 fixes the maximal dissipative
charge of the GHZ off-diagonal operators, while W coherences have distance-two
support. Neither operator fact assigns a state weight. The event asked whether
the spectral organization could be exploited for protection; its coordinate
calculation did not answer that question.

---

## 1. Historical Rate Regions

The palindromic spectrum of the N=3 Heisenberg Liouvillian (γ=0.05, Z-dephasing)
has 32 palindromic pairs that fall into three natural tiers:

| Tier | Pairs | Rates | XY-weight | Character |
|---|---|---|---|---|
| Steady-XOR | 4 | 0 and 0.30 | 0.0 vs 3.0 | Extremes: immortal paired with fastest drain |
| Boundary | 14 | 0.10 and 0.20 | 1.0 vs 2.0 | Exact rates 2γ and 4γ, the pure-weight rungs of [F33](../docs/ANALYTICAL_FORMULAS.md#f33-the-n3-rate-ladder-tier-1-for-the-pure-weight-rungs-the-two-fractional-rates-are-a-jgamma---infinity-limit) |
| Mid-spectrum | 14 | ~0.133 and ~0.167 | 1.33 vs 1.67 | Near 8γ/3 and 10γ/3, F33's two fractional rates, reached in the J/γ → ∞ limit |

The boundary rates are exact eigenvalues and the mid-spectrum rates approach
F33's limits. Neither tier by itself establishes information transport, state
lifetime, or a classical/quantum split.

---

## 2. Withdrawn Eigenvector-Coordinate Search

A constrained search maximized squared coefficients of `R_inv @ rho` after a
right-eigenvector decomposition. For a non-normal generator these coefficients
are not orthogonal probabilities: rescaling an eigenvector changes them, and
near-degenerate bases can mix them. Consequently the following table is kept
only as a record of what the retired objective printed, not as a comparison of
survival or protection:

| Retired coordinate diagnostic | Returned state | W | Bell(0,1) | GHZ |
|---|---|---|---|---|
| Slow-mode weight | 90% | 0% | 7% | 0% |
| XOR weight | 0.02% | 0% | 0% | 55.55% |
| Oscillating content | 18% | 0% | 38.09% | 0% |
| Concurrence | 0.364 | 0.667 | 1.000 | 1.000 |

The returned state is composed mainly of \|010⟩, \|000⟩, \|100⟩, and \|001⟩,
and its concurrence is a well-defined property of that state. No survival
ordering follows from the coordinate rows. A replacement comparison must
specify an observable or channel metric and propagate each state over the same
time window.

---

## 3. The Projected-Component Fingerprint

The fingerprint projects a state onto the oscillating part of the spectrum,
with the spectral projector P_osc onto every eigenvalue with Im λ ≠ 0 (40 of the
64), and reads the moduli of the result's Pauli coefficients. A spectral
projector does not depend on how the eigenvectors are normalized, so these are
invariant numbers. The question is how far a single-qubit error moves them.

For Bell(0,1) = (|000⟩ + |110⟩)/√2 under the Heisenberg chain, the largest
change of any Pauli coefficient:

| Error | Site 0 | Site 1 | Site 2 |
|---|---|---|---|
| X | 0.0278 | 0.0278 | 0 |
| Y | 0.0556 | 0.0556 | 0 |
| Z | 0.0278 | 0.0278 | 0 |

Site 2 is silent for a reason that has nothing to do with the error type: on
this state X₂ acts exactly as X^⊗3, Y₂ as i·X^⊗3 and Z₂ as the identity, and
X^⊗3 commutes with L and only flips the signs of Pauli coefficients, so the
moduli cannot move. Sites 0 and 1 respond to all three errors, Z included.

A change in a projected component is not a syndrome: nothing here propagates
the state or measures anything. Producer:
[`simulations/ec_projected_fingerprint.py`](../simulations/ec_projected_fingerprint.py).
The March 19 output's Section 3 numbers read a 4^N Pauli-coefficient vector as a
d×d matrix and are not this quantity.

---

## 4. Π Has Fourth-Order Structure

A mirror you look into twice gives you back yourself (second-order:
Π² = I). But our mirror Π is stranger: you need to apply it four
times to get back to the start. It is more like a 90° rotation than
a reflection.

The Π operator has eigenvalues {+1, -1, +i, -i}, each with multiplicity 16.
This means Π⁴ = I (fourth-order, not second-order). The 64-dimensional
Pauli space decomposes into four 16-dimensional sectors.

The four eigenspaces are algebraic sectors of the constructed map:

- The +1 sector: modes invariant under the linear Π action
- The -1 sector: modes that flip sign (antisymmetric)
- The +i and -i sectors: modes that rotate by 90° (quarter-wave shifted)

No distinct physical meaning or error-correcting codespace follows from these
eigenvalues alone.

---

## 5. Historical Equal-Amplitude Pair Curve

The historical run assigned equal amplitudes to each rate pair and printed the
following decay times for that constructed curve:

| Tier | Slow rate | Fast rate | T_half | T_1% |
|---|---|---|---|---|
| Steady-XOR | 0.000 | 0.300 | infinite | infinite |
| Boundary | 0.100 | 0.200 | 19.2 | 59.9 |
| Mid-spectrum | 0.133 | 0.167 | 18.0 | 51.9 |

The table evaluates an assumed difference of two exponentials with equal
initial modal amplitudes. It is not the lifetime of a prepared density matrix,
observable, or channel.

---

## 6. Historical QEC Comparison (Withdrawn)

The run compared its coordinate diagnostic with named 3-qubit states:

The historical table reported 55.55% for the repetition/GHZ coordinate, 0%
for the phase-flip coordinate, and 66.54% slow coordinate for `|010>`. These
are not invariant QEC-performance measures.

The former claim that the Section 2 state outperforms these codes was based on
the withdrawn coordinate objective. This run contains no direct-propagation
gate establishing such an ordering.

---

## 7. The XOR Drain Is Not a Universal Syndrome

No single-qubit X, Y or Z error moves W into the XOR sector: W's share there is
exactly 0 before and after. The XOR sector is spanned by X^⊗3·P_k, whose entries
sit at Hamming distance 3, and it reduces L ([XOR Space](XOR_SPACE.md)), so a
state's share in it is invariant. W's entries sit at distance 0 or 2, and a
single-qubit Pauli flips the same bit on both sides of |a⟩⟨b|, which leaves
a ⊕ b unchanged. So the drain cannot flag any single-qubit error on W; it is no
universal syndrome. Section 3 likewise establishes no detector.

---

## Summary of Findings

1. The palindromic spectrum has distinct rate regions; translating those
   regions into state protection requires an operational propagation test.
2. The former `90% slow-mode` optimum and its survival ranking are withdrawn
   because they used non-invariant right-eigenvector coordinates.
3. The projected-component changes (0.0278 and 0.0556 on sites 0 and 1, exactly
   0 on site 2) are invariant but are not a demonstrated syndrome.
4. Π is a fourth-order operator (Π⁴ = I), creating an algebraic Z₄
   decomposition without assigned physical sector roles.
5. The equal-amplitude pair curves are not information lifetimes.
6. This run establishes no QEC-code ranking.

---

## References

- [Π as Time Reversal](PI_AS_TIME_REVERSAL.md): centred spectral pairing and its interpretive history
- [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md): oscillation patterns
- [XOR Space](XOR_SPACE.md): endpoint count, F22 operator support, and the retired GHZ/W coordinate decomposition
- [N→∞ Palindrome](N_INFINITY_PALINDROME.md): spectral scaling
- Historical retired producer: [`simulations/error_correction_palindrome.py`](../simulations/error_correction_palindrome.py)
- Active F22 producer: [`simulations/f22_operator_charge.py`](../simulations/f22_operator_charge.py)
- Historical March 19 event output: [`simulations/results/error_correction_palindrome.txt`](../simulations/results/error_correction_palindrome.txt)
- Section 3 producer: [`simulations/ec_projected_fingerprint.py`](../simulations/ec_projected_fingerprint.py), output [`simulations/results/ec_projected_fingerprint.txt`](../simulations/results/ec_projected_fingerprint.txt)
- Exact XOR-sector gate: [`simulations/xor_verify.py`](../simulations/xor_verify.py)
- Current F22 gate output: [`simulations/results/error_correction_palindrome_f22.txt`](../simulations/results/error_correction_palindrome_f22.txt)
