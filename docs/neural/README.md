# Neural Systems and the Palindromic Symmetry

**Audience:** Neuroscientists, computational biologists, anyone working
with neural network dynamics.

**No quantum physics required.** The results here are derived from
quantum theory but stated and tested in purely classical terms
(Jacobians, Wilson-Cowan dynamics, connectome data).

---

## What is the palindromic symmetry?

When a network of neurons oscillates (excitatory and inhibitory
populations coupled together), the dynamics are governed by
eigenvalues of the system's Jacobian matrix. Each eigenvalue
determines one mode of decay or oscillation.

These eigenvalues have a hidden structure: their decay rates can
be **palindromic** - mirrored around a center point, like the word
"racecar" reads the same forwards and backwards. If the rates are
{0.03, 0.05, 0.08, 0.10, 0.12}, the palindromic partner of 0.03
is 0.12 (they sum to 0.15), the partner of 0.05 is 0.10, and 0.08
sits at the center.

In quantum physics, this symmetry is exact (proven algebraically).
Here we ask: does the same structure appear in biological neural
networks, and if so, why?

---

## What we found

1. **The symmetry has two ingredients:** different time constants for
   excitatory and inhibitory neurons (selective damping), and Dale's
   Law (E neurons always excite, I neurons always inhibit). Both are
   universal in neurobiology. Dale's Law is the biological analog of
   the antisymmetric quantum commutator.

2. **An exact algebraic condition** for a perfect palindrome is derived
   from quantum theory. It predicts a specific magnitude ratio between
   partnered E-E and I-I connections, testable on any connectome with
   known E/I labels and synaptic weights.

3. **C. elegans balanced subnetworks** are 8x more palindromic than
   Erdos-Renyi random networks (robust across parameter choices).
   Degree-preserving randomization shows this advantage comes from
   the degree distribution (hub vs peripheral neurons), not from the
   specific wiring pattern. The degree distribution is itself a
   biological property that varies between organisms.

---

## Documents

| Document | What it covers |
|----------|---------------|
| [Algebraic Palindrome](ALGEBRAIC_PALINDROME_NEURAL.md) | The derived condition, C. elegans test, random controls |

---

## Prerequisites

None for the main results. For the quantum foundation:
- [Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md) (the original quantum proof)
- [The Interpretation](../THE_INTERPRETATION.md) (what the palindrome means physically)

For the original C. elegans hypothesis:
- [The Pattern Recognizes Itself](../../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)
