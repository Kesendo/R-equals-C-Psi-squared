# The Polarity Layer: Inheritance, +0/−0, and Self-Description from Inside

**Date:** 2026-04-30 (updated 2026-05-01)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Status:** Tier 4 (interpretive synthesis grounded in Tier 1-2 results F77-F85). The mathematics is unchanged from the established F-chain; this document re-reads it. **Update 2026-05-01:** the original "+0/−0 polarity layer" claim is sharpened to a multi-axis Klein-Vierergruppe Z₂² (k=2) / Z₂³ (k≥3) polarity structure; operational consequences are now checkable via the `PauliHamiltonian` class and the `diagnose_hardware` workflow; the trace-back primitive `recover_H_odd_from_M_anti` listed as open in the original is closed.
**Depends on:**
[F-chain F77-F85](../docs/ANALYTICAL_FORMULAS.md),
[ON_THE_RESIDUAL](../reflections/ON_THE_RESIDUAL.md),
[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md),
[The Primordial Qubit](PRIMORDIAL_QUBIT.md),
[Resonance Not Channel](RESONANCE_NOT_CHANNEL.md),
[Gamma Is Light](GAMMA_IS_LIGHT.md),
`framework.PauliHamiltonian` (added 2026-05-01),
`framework.diagnostics.recover_H_odd_from_M_anti` (added 2026-04-30),
`framework.workflows.diagnose_hardware` (added 2026-04-30, refactored to PauliHamiltonian 2026-05-01)

---

## What this document is about

A brainstorming session on 2026-04-30 produced three reframings of the
R=CΨ² framework that are mathematically consistent with the F-chain
F77-F85 but operate in different vocabulary than the framework's
prior "memory of the mirror" reading. The reframings interlock:

1. **The system is a self-coupling channel, not a memory.** Π·L·Π⁻¹ + L + 2Σγ·I = M describes a regenerative cavity that keeps the signal alive through self-recirculation, not a memory device that stores then retrieves.
2. **The "two sides" are bra and ket of one ρ, not two subjects.** The framework describes a single self-coupling system whose internal two-indexedness Π conjugates. We are not external observers of a mirror; we are the system.
3. **The "0" axis has internal polarity structure.** d²−2d=0 splits not into "qubit (d=2)" versus "void (d=0)" but into "qubit (d=2)" versus "+0/−0 polarity layer (d=0 with internal differentiation)". The qubit is a window onto this polarity layer, accessed naturally through the X-basis.

A meta-claim ties them together: **self-description from inside is
possible because rules inherit through layers.** The F-chain F77→F85
is the operational proof.

---

## Abstract

The framework's interpretation has carried a residue of observer-dualism:
"the mirror reflects toward us; we measure what comes back; M is what
the mirror remembers." We propose this reading is Shannon-storage
vocabulary contaminating the cavity-light vocabulary established by
[Gamma Is Light](GAMMA_IS_LIGHT.md) and [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md).

Cavity vocabulary is more accurate operationally. A cavity does not
store; it actively maintains. The palindrome equation describes a
self-coupling resonator whose output reenters its input through Π
conjugation. M is the operational residue of incomplete recirculation,
not a memory contents. Truly Hamiltonians (M = 0) describe the perfect
closed cycle; Π²-odd Hamiltonians describe a self-pumping cycle where
the dynamics generator IS the recirculation operator (F81's
M_anti = L_{H_odd}).

Within this single self-coupling system, the "two sides" are not two
external subjects but the bra and ket indices of one density matrix.
The d²−2d=0 condition expresses that operator space (d² = 4^N for N
qubits) is the natural arena: bra-ket co-existence on a single state.

The d=0 axis is then not a passive midpoint between two observers;
it is the active substrate from which the system draws coherence —
the "Stromanschluss" in Tom's metaphor, analogous to the QFT vacuum's
zero-point activity. The d=0 axis itself has internal structure: a
+0/−0 polarity that constitutes (rather than passively underlies)
the layer on which the qubit lives. The qubit's natural projection
onto this layer is the X-basis, where |+⟩ = +0 and |−⟩ = −0 are the
two polarities.

This explains operationally why X-basis observables dominate the
F-toolkit's discriminators, why |+,−,+⟩ X-Néel is the canonical
hardware initial state, and why amplitude damping (σ⁻) on truly
Hamiltonians produces the cleanest F82/F84 signature: σ⁻ attacks
the polarity layer directly, breaking the ⟨Z⟩-conservation that the
truly Hamiltonian alone preserves.

---

## The inheritance correction

A common claim is that a formal system cannot describe itself from
inside (Gödel, Russell). This is true for **flat** formal systems
without internal layered structure. It is **not** a constraint on
systems where rules inherit through layers.

The F-chain F77→F85 is constructed by inheritance:

- **F77** trichotomy classifier (truly / soft / hard) on 2-body Pauli
  pairs.
- **F78/F79** structural M-decomposition for 1-body and 2-body.
- **F80** spectrum identity Spec(M) = 2i · Spec(H_non-truly).
- **F81** Π-decomposition M_anti = L_{H_odd}.
- **F82** T1 closed form ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N-1).
- **F83** anti-fraction r = ‖H_even_nontruly‖² / ‖H_odd‖².
- **F84** Pauli-channel cancellation lemma.
- **F85** k-body generalization: truly criterion = "#Y even AND #Z even" inherits to all body counts; F77-F84 generalize verbatim.

Each layer inherits from the previous and refines it. The 2-body
n_YZ-counting that worked for F49 turned out to be a coincidence at
the 2-body layer; F85 lifted it into the body-count-independent
Π²-class structure. The inheritance is the operative content, not a
metaphor.

A system whose F-chain inherits through layers can describe itself
from inside one layer at a time. F77 describes Hamiltonian
classification. F80 describes M's spectrum given F77's classification.
F82 describes T1's contribution given F81's decomposition. The whole
chain describes the system at increasing depth. The inheritance
braids the layers together; no single layer needs to describe everything,
and no closure paradox arises.

---

## The +0/−0 polarity layer

The d²−2d=0 condition has solutions d=0 and d=2. d=2 is the qubit's
Hilbert dimension. d=0 has been read as "the mirror axis" or "the
through-line" but its internal structure has not been named.

We propose: **d=0 is a polarity-birthing axis with internal structure.**
The polarity layer is constituted by polarity differentiation; a layer
exists *because* polarity differentiates. Without differentiation, there
is the pre-polarized vacuum substrate ("Stromanschluss"). With it, there
is a layer.

**Update 2026-05-01: the polarity is multi-axis, not bilateral.**

The original "+0/−0" framing identified one axis (the X-basis polarity).
This is the **bit_a axis** of the framework's operator-letter parity
structure. The full polarity at k=2 is the **Klein-Vierergruppe Z₂×Z₂**
with two independent axes: bit_a (X-eigenstate polarity, accessed via
X-basis) and bit_b (Π² parity, accessed via Z-basis). At k≥3 a third
independent axis (Y-parity) appears: the polarity is **Z₂³** with 8
sectors. The Klein-Vierergruppe is the k=2 collapse of the Z₂³ structure.

The +0/−0 reading is therefore correct for the **bit_a axis**:

  |+⟩ = (|0⟩ + |1⟩)/√2     ← X-eigenstate +1 = +0 (bit_a-positive polarity)
  |−⟩ = (|0⟩ − |1⟩)/√2     ← X-eigenstate −1 = −0 (bit_a-negative polarity)

But there is also a **bit_b axis** with its own polarity structure:

  |0⟩                        ← Z-eigenstate +1 (bit_b-trivial Mother slot under bit_a)
  |1⟩                        ← Z-eigenstate −1 (bit_b-trivial Mother slot under bit_a)

The X-basis is the natural diagonalization of bit_a. The Z-basis is
the natural diagonalization of bit_b (the framework's canonical
dephasing axis). Numerically, +0 = −0 = 0 in any of these axes; as
polarities, they are operative projections — but on the FULL polarity
layer, there are multiple axes.

Operationally, this is now visible in the `PauliHamiltonian` class:

```python
H = fw.PauliHamiltonian.from_letter_tuples(
    [('X', 'Y'), ('Y', 'X')], chain_length=3
)
H.klein_set                          # → {(0, 1)}  — bit_a=0, bit_b=1
H.is_klein_homogeneous               # → True
H.per_term_full_z2_signatures        # → [(0,1,1), (0,1,1)]  — bit_a, bit_b, Y-par
H.is_z2_homogeneous                  # → True at k=2 (redundant); independent at k≥3
```

The framework operates on this multi-axis polarity throughout:

- **Initial state |+, −, +⟩ X-Néel** projects onto the bit_a axis
  (X-eigenstates) palindromically across the chain. A Z-Néel initial
  state would project onto the bit_b axis.
- **Π reflection** under F71 chain-mirror swaps the polarities along
  the chain.
- **F-discriminating observables** cluster on X-rotated combinations
  because the F83 hardware test selects the bit_a-polarity diagonalization;
  a different choice of test could probe the bit_b axis directly.
- **Klein-homogeneity rule** (verified at k=2 full enumeration, k=3
  sample): Hamiltonians whose terms all share one Klein index are
  always F77 soft or truly. This is a structural fact, exposed via
  `PauliHamiltonian.is_klein_homogeneous`.

The X-basis is not "the" polarity layer — it is the natural
diagonalization of one of the two (k=2) or three (k≥3) Z₂ axes that
constitute the polarity layer.

---

## The two readings of polarity-layer

The +0/−0 structure can be read two ways:

**(a) Palindromic over the layer.** Π couples +0 (one site) to −0
(its reflected site). The layer exists; the polarities live on it
and palindromically pair across it.

**(b) Are the layer.** The +0/−0 differentiation IS what constitutes
a layer. Below this differentiation, there is the pre-polarized
substrate (vacuum, d=0 axis as Stromanschluss). The layer emerges
from polarity.

Both readings are mathematically consistent with the F-chain. (b) is
the more radical: layers are not pre-existing scaffolds for polarity
but emergent structures from polarity. This connects to:

- **Cavity QED:** the field IS the cavity's coupled structure; there
  is no field separate from the cavity-vacuum coupling.
- **Active matter / dissipative structures:** structure exists by
  continuous flow; static substrate has no structure.
- **Hippocampal memory:** "stored" memories are continuously
  reactivated patterns, not static traces.

---

## Operational consequences

The math is unchanged. The vocabulary shifts, and with it the way
results are read. **Update 2026-05-01: each consequence is now directly
checkable in code via `PauliHamiltonian` properties and the
`diagnose_hardware` workflow — the reframing is no longer interpretive
text but operational structure.**

| Old (memory/observer) | New (channel / multi-axis polarity) | Operational primitive |
|------------------------|--------------------------------------|------------------------|
| "The mirror remembers" | "The channel maintains" | `pi_decompose_M` returns M_anti = recirculation |
| "Two sides observe each other" | "Bra and ket of one ρ" | d² = 4^N operator-space arena |
| "0 is the midpoint" | "0 is the active substrate; multi-axis polarity layer is Z₂² (k=2) or Z₂³ (k≥3)" | `PauliHamiltonian.klein_set`, `.full_z2_signature_set` |
| "Trace back = remember" | "Trace back = extract the recirculation operator" | `recover_H_odd_from_M_anti` (built in commit 2061fa2) |
| "Memory of the dynamics" | "L_{H_odd} IS the recirculation; the dynamics drives its own maintenance" | `pi_decompose_M`'s `M_anti` field |
| "Truly = no memory" | "Truly = perfectly closed cavity; Q→∞ idealized" | `classify_pauli_pair → 'truly'`; `PauliHamiltonian.has_truly_term` |
| "Qubit is \|0⟩, \|1⟩" | "Qubit is a window onto multi-axis polarity; X-basis diagonalizes bit_a, Z-basis diagonalizes bit_b" | `PauliTerm.klein_index`, `.full_z2_signature` |

This makes prior puzzles legible:

- **Why X-basis tomography is more sensitive than Z-basis for truly
  Hamiltonians.** X projects onto the natural polarity layer where
  the dynamics live. Z projects onto the computational view (where
  XX+YY-truly's ⟨Z⟩-conservation hides the dynamics).
- **Why amplitude damping (σ⁻) breaks truly's ⟨Z⟩-conservation
  cleanly.** σ⁻ destroys |1⟩ population, which in X-basis is
  ½(|+⟩ − |−⟩). σ⁻ attacks the −0 component, breaks the polarity
  balance, and through that breaks the ⟨Z⟩-conservation that the
  truly Hamiltonian alone protects. The ~60% ⟨Z,Z⟩ damping observed
  in the 2026-04-30 F83 hardware run on Marrakesh (Job d7pol1e7g7gs73cf7j90;
  damping fraction is 60% relative to γ_Z=0.05 path-fit baseline,
  56% relative to γ_Z=0.1 framework-default baseline) is the F82/F84
  signature seen through the polarity-layer lens.
- **Why ZZ-crosstalk is not a Π-breaker.** ZZ leaves the polarity
  layer untouched (Z and X are bit_b-different in the F77 classifier).
- **Why the F-toolkit's hardware-confirmed predictions cluster on
  X-rotated observables.** They project onto the natural layer.

---

## Relation to existing hypotheses

This document does not replace any prior hypothesis; it connects them:

- [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) computes the palindrome
  at Σγ = 0 (no noise) and below. We sharpen the meaning of "0 is the
  mirror": the 0-axis is not a passive midpoint but the polarity
  layer's birthing substrate.
- [The Primordial Qubit](PRIMORDIAL_QUBIT.md) proposes "system and
  noise are two readings of a single structure split by a Z₂ mirror
  at zero." We add: the Z₂ at zero is the +0/−0 polarity, and the
  qubit's X-basis is the natural projection.
- [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md) establishes the
  cavity reading. We extend: the cavity is self-coupling through Π
  (channel-and-cavity, not channel-or-cavity).
- [Gamma Is Light](GAMMA_IS_LIGHT.md) establishes γ as light. We
  extend: light enters at the +0/−0 polarity layer; the qubit is the
  cavity's window onto this layer.

---

## What this document does NOT establish

- **No new mathematics.** F77-F85 stand unchanged. The hypothesis
  reframes how those theorems are read, not what they say.
- **No new operational primitives.** A future operational primitive
  `recover_H_odd_from_M_anti(chain, M_anti)` would test the
  trace-back direction explicitly; not built yet (would live in
  framework/diagnostics/f81_pi_decomposition.py).
- **No experimental falsification path.** The reframing predicts
  the same hardware signatures the prior reading predicted; only
  the language differs. A genuine test would require finding a
  signature that the polarity-layer reading predicts but the
  memory-reading does not.
- **No claim about consciousness, observation, or measurement.**
  The "we are the system" framing is mathematical (we = ρ, the
  bra-ket structure on d²) and structural, not phenomenological.

---

## Open directions

- **~~Operational trace-back primitive.~~** ✓ **CLOSED 2026-04-30
  (commit 2061fa2).** `recover_H_odd_from_M_anti(chain, M_anti)` lives
  in `framework.diagnostics.f81_pi_decomposition`. Round-trip verified
  bit-exact via least-squares projection onto the basis of
  commutator-action operators L_{P_α} for Π²-odd Pauli strings; pytest
  test `test_F81_recover_H_odd_from_M_anti_round_trip` covers the
  truly / pure Π²-odd / pure Π²-even non-truly / mixed cases.
- **Polarity-layer reading of higher-body F-theorems.** ◐ **PARTIALLY
  CLOSED 2026-05-01.** At k≥3, Y-parity becomes an independent third
  Z₂ axis: the polarity is Z₂³ (8 sectors), not Z₂² (4 Klein slots).
  Verified empirically — XYZ at k=3 has Klein index (0,0) but
  pi2_class 'pi2_even_nontruly' (NOT truly), distinguishable only via
  Y-parity. `PauliHamiltonian.is_z2_homogeneous` exposes this third
  axis. The k=2 Klein-Vierergruppe is the collapse of the full Z₂³
  structure under Y-parity = bit_a XOR bit_b (which holds at k=2 but
  not at k≥3). Bond-position structure within a fixed Z₂³ sector
  remains open.
- **Inheritance proof structure.** ▭ **OPEN, with new structural finding
  2026-05-01.** The inheritance argument here remains informal. The
  Klein-homogeneity rule "all terms in same Klein index → F77 soft or
  truly, never hard" is a candidate structural invariant.

  **Verified strict at k=2** (full enumeration of 36 V-Effect pairs):
  every Klein-homogeneous k=2 Hamiltonian on a chain with Z-dephasing
  is F77 truly or soft — 0/6 are hard.

  **Tendential but not strict at k=3.** Full sweep at k=3 N=4 with
  Z-dephasing across all 240 Z₂³-homogeneous pairs in {I,X,Y,Z}^3:
  90 truly + 104 soft + 46 hard. The 46 hard cases reveal a sharp
  asymmetry — they all lie in Klein **(0, 1)**, which is the same
  Klein index as the Z-dephasing dissipator's letter Z.

  Klein-(0, 0): 0 hard / 60 — Mother sector, fully clean.
  Klein-(0, 1): 46 hard / 60 — same Klein as Z-dephasing.
  Klein-(1, 0): 0 hard / 60.
  Klein-(1, 1): 0 hard / 60.

  This suggests the dissipator's Klein index identifies a "preferred"
  sector where the Klein-homogeneity rule degrades. The polarity layer's
  symmetry under Klein-homogeneity is **broken by the dissipator
  selection**: choosing Z-dephasing privileges Klein (0, 1) at k=3.

  Testable conjecture: with X-dephasing (X has Klein (1, 0)) the hard
  cases would shift to Klein (1, 0); with Y-dephasing (Klein (1, 1))
  to Klein (1, 1). The structural rule may then be: "Klein-homogeneous
  AND Klein-orthogonal-to-dissipator → F77 soft/truly always."
- **Polarity-layer reading of the d=0 axis at higher d.** ▭ **OPEN.**
  d² − 2d = 0 has solutions d=0 and d=2 because the qubit dimension is 2.
  For qudit (d-level) systems, the analog condition would be different
  (d² − cd = 0 for some c). The polarity-layer reading might generalize
  to qudits with multi-polarity structure.

---

## Closing note

Three readings developed in one session. None are physics, all are
language. The hypothesis is that *this* language is the one that
makes the F-chain operationally legible. The math sits underneath
unchanged. The vocabulary shift is what removes the observer-dualism
that contaminated the prior reading, and what makes the
self-description-from-inside legitimate via the inheritance structure
that the F-chain itself instantiates.

> "Wir SIND das Selbst-koppelnde System. 0 ist der Stromanschluss.
> +0/−0 ist der Layer."
> — Tom Wicht, 2026-04-30
