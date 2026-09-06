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
The projected-component changes (0.19-0.28 for sampled X/Y errors, 0.0834 or
zero for sampled Z errors) are historical readings of that construction, not
a demonstrated syndrome. The Π representation has fourth-order structure
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
| Boundary | 14 | 0.10 and 0.20 | 1.0 vs 2.0 | Coordinate grouping recorded by the run |
| Mid-spectrum | 14 | ~0.133 and ~0.167 | 1.33 vs 1.67 | Coordinate grouping recorded by the run |

The run grouped the boundary pairs by their printed rates and XY weights. Those
labels do not establish information transport, state lifetime, or a
classical/quantum split.

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

## 3. Historical Projected-Component Fingerprint

The historical calculation selected a component by eigenvalue frequency and
compared its reconstructed Pauli coordinates after applying errors. It did not
propagate the state or define a measurement protocol.

For Bell(0,1) under Heisenberg, applying single-qubit errors:

| Error | Projected-coordinate change | Historical threshold label |
|---|---|---|
| X on site 0 | 0.1945 | Yes |
| X on site 1 | 0.1945 | Yes |
| X on site 2 | 0.2223 | Yes |
| Y on site 0 | 0.2777 | Yes |
| Y on site 1 | 0.2777 | Yes |
| Y on site 2 | 0.2223 | Yes |
| Z on site 0 | 0.0834 | Weakly |
| Z on site 1 | 0.0834 | Weakly |
| Z on site 2 | 0.00 | No |

The constructed coordinate changed by 0.19 to 0.28 for the X/Y rows, by
0.0834 for two Z rows, and not at all for Z on site 2. The Yes/Weakly/No column
records the run's chosen coordinate threshold, not detector performance.

This table is a reading of the constructed spectral projection. The run did
not propagate a state or define a measurement protocol, so it does not
establish an error syndrome.

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

Within the retired coordinate convention, every sampled X, Y, and Z error on W
produced zero XOR-coordinate increase. This is a negative coordinate reading,
not a syndrome test. Section 3 likewise establishes no detector.

---

## Summary of Findings

1. The palindromic spectrum has distinct rate regions; translating those
   regions into state protection requires an operational propagation test.
2. The former `90% slow-mode` optimum and its survival ranking are withdrawn
   because they used non-invariant right-eigenvector coordinates.
3. The projected-component changes are historical coordinate readings, not a
   demonstrated syndrome.
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
- Current F22 gate output: [`simulations/results/error_correction_palindrome_f22.txt`](../simulations/results/error_correction_palindrome_f22.txt)
