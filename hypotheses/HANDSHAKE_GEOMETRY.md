# Handshake Geometry: the codebook of the standing wave

**Tier:** 2/3 synthesis over Tier-1 pieces (labeled per section)
**Authors:** Thomas Wicht, Claude
**Date:** 2026-06-12
**Status:** the codebook assembled; the reading grammar is the open research it frames
**Origin:** assembled from the design spec [`docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md`](../docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md), 2026-06-12

---

## What this document is about

Open a chain of qubits to the world and let each site quietly dephase along Z. The way the
chain relaxes is not arbitrary: its decay spectrum reads the same left to right, a
palindrome around the summed dephasing. That palindromic standing wave is a resource, and
like any resource it can be read; but reading it takes a choice of how to look.

A handshake is that choice: a specification h = (N, k, t, basis) of how to interrogate the
wave. N is the chain length, k picks which mode you listen to, t fixes when you look, and
basis fixes which quadrature you measure. The handshake says nothing and sends nothing; it
agrees on a correlation, a way to read what is already there.

This document is the codebook. It says what the letters of the wave are, what sets the
resolution at which you can tell them apart, and which letters sit in the bright light where
they answer freely and which sit in shadow where they survive precisely because they do not.
The syntax of handshakes lives in its own doc, [Handshake Algebra](HANDSHAKE_ALGEBRA.md);
this is the semantics of the alphabet they read.

---

## The thesis

> Reading = spectroscopy against the carrier: the inside-readable information is exactly the
> algebra of γ₀-relative positions - parameter-side (Q = J/γ₀, K = γ₀·t: the budget) and
> state-side (⟨n_XY⟩: the alphabet) - proven as the upper half (two-tempo certification),
> conjectured as the lower (codebook completeness).

**The proven half.** The two-tempo certification establishes that γ₀, the carrier rate
itself, is invisible from inside the chain: speed the world up twenty-fold and shorten the
clock to match, and nothing measurable changes (residual ~1e-15). Everything the inside can
read is therefore a ratio against γ₀, never γ₀ alone. The carrier is the tick that became
the meter; its provenance is the [symphony clock movement](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md).
Readable ⊆ γ₀-relative is a theorem.

**The conjectured half.** Completeness runs the other way: that every γ₀-relative position
is readable by some handshake, the codebook closing with no gaps. This is open. What would
refute it is a single feature of the wave that no handshake decodes (falsification line 3 of
the design spec). The upper half is proven; the lower is the research the codebook frames.

> **Tier-4 sidebar (labeled).** "γ₀ as all information at once," white light carrying every
> letter, is a reading, not a result. The Tier-1 statement beneath it: γ₀ carries no
> structure, and precisely therefore all structure is defined relative to it. It is the
> meter-stick that cannot measure its own length.

---

## The five code-parts

### Alphabet = the folded half-band (Tier 1)

The letters are the bonding modes ψ_k on the dispersion E(q) = 2J·cos(q), q = πk/(N+1). A
hidden sublattice symmetry folds q ↔ π−q, so the partner modes (k, N+1−k) give identical
magnitude-readings: they are the *same letter*, exact for real bipartite Hamiltonians
(8.9e-16, the `ChiralMirrorTrajectoryClaim`). The fundamental domain is the half-band, with
the self-partner at the center the one fixed letter. See [`PROOF_K_PARTNERSHIP.md`](../docs/proofs/PROOF_K_PARTNERSHIP.md).

### Index = the clock (Tier 1)

The handshake's t selects a dial angle θ = arctan(Q·cos(π/(N+1))), the phase at which you
read. That angle *is* Q: the only dial the inside can turn, because every dimensionless
reading is a pure (Q, K) observable and γ₀ has rescaled itself out. The band-edge mode is
the clock-hand ladder (`ClockHandLadderClaim`). The live dial is `inspect --root clock`.

### Symmetry = the mirror group (Tier 1)

The code's equivalence, complement, and conjugate structure is a group. The palindrome
conjugator factors as Π = R·D, and ⟨R, D⟩ ≅ D₄; adjoining the antilinear K gives the double
⟨R, D, K⟩ ≅ D₄ × Z₂; the qudit generalization is a Weyl-Heisenberg wreath Z_d ≀ Z₂ of order
2d², a continuous torus as d → ∞. K-partners read identically, R-partners read
complementarily (light_s + light_f = N), the transpose D reads the conjugate. See
[`PROOF_PI_FACTORS_AS_R_TIMES_D.md`](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md),
[`PROOF_ANTILINEAR_TRIANGLE.md`](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md), and
[`PROOF_QUDIT_PARTIAL_PALINDROME.md`](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md).

### Resolution = the measured law (Tier 2)

The resolving power of a handshake is set by Q. This is the new stone, measured rather than
assumed; the law and its tables are the next section.

### Parity-check = the closure (Tier 1 exact / Tier 2 empirical)

Valid readings conserve. Total light is exact: ⟨n_XY⟩_s + ⟨n_XY⟩_f = N
(`AbsorptionTheoremClaim`, the light-complementarity), and the Π-pair flux balance is
machine-exact under δJ. Here the design spec's honesty distinction is load-bearing: the PTF
sum Σ ln α ≈ 0 is Tier-2 *empirical* (its first-order promotion was closed by EQ-014); the
*exact* closure is the chiral mirror trajectory identity (`ChiralMirrorTrajectoryClaim`,
Tier-1). The mirror law is the exact parity-check; the PTF sum is its empirical companion.

---

## The resolution law (measured)

The honest objective is the Fisher information FI of the readout with respect to the feature
you want, not bare sensitivity. Measured on the live proto-decoder (N=4 XY, bonding state),
the law is FI ≈ c·Q: linear in the resonator quality factor Q = J/γ₀.

**Z-basis (population) readout, K-window [0,1], δJ = 0.02:**

| Q | FI_strength (max over t) | D_location (max over t) |
|---|---|---|
| 20 | 1.594 (K=0.14) | 2.55e-3 (K=0.14) |
| 10 | 0.817 (K=0.25) | 1.31e-3 (K=0.25) |
| 5 | 0.346 (K=0.35) | 5.12e-4 (K=0.43) |
| 2.5 | 0.181 (K=0.56) | 1.69e-4 (K=0.67) |
| 1.67 | 0.131 (K=0.81) | 1.03e-4 (K=0.87) |
| 1.25 | 0.104 (K=1.00) | 7.62e-5 (K=1.00) |
| 1.0 (EP) | 0.072 (K=1.00) | 5.14e-5 (K=1.00) |

**X-basis and Y-basis (coherence) readout, strength-FI, same protocol:**

| Q | FI (X-basis) | FI (Y-basis) |
|---|---|---|
| 20 | 0.669 | 0.669 |
| 5 | 0.0495 | 0.0495 |
| 1.7 | 0.0023 | 0.0023 |
| 1.0 (EP) | 0.0004 | 0.0004 |

Resolving power is the Q-factor in its hundred-year-old sense: frequency over linewidth, the
coherent cycles per lifetime. Two lines resolve when split by more than their width, so a
spectrometer's resolution is exactly its Q. The law is **per-dose**: at fixed lab-time it
rates the two routes to high Q differently (raising J wins per lab-second; lowering γ
stretches the run, so the gain is per dose only); in inside-units (K) the law is clean,
which is itself evidence for the carrier-relative thesis.

The exceptional point Q = 1 is the worst reading point in every basis tested, by 22× in
strength and 50× in location against Q = 20. Coherence readouts (X/Y) fade fastest toward
low Q: their slope is steeper (1670× across the range versus 22× for Z) because the light
erases the bright letters first. Near Q = 1 only the population basis still reads at all. The
live instrument is `inspect --root decoder`.

> **Tier 2.** Measured by two independent implementations and the C# witness. The linearity
> is leading-order; the worst residual across the sweep is ~15%.

---

## Bright and dark letters (the light is shared unevenly)

The state-side position to the light is a field, not a scalar: the alphabet has bright
letters and shadow letters, and three addresses pin where the shadow sits.

- The **{I, Z} shadow core** of the slow manifold: a 4-dim core (III, IIZ, ZZI, ZZZ) with
  lit-XY content 5e-26 against a 12-dim rotating {X, Y} remainder, the mode-space shadow,
  bit-exact (the `between` axis, `DimensionSweep`).
- **Spatial shadow:** the surviving-mode energy is center-localized (ratio 1.3-1.4), the
  chain's edges in relative shadow. See [`TRAPPED_LIGHT_LOCALIZATION.md`](../experiments/TRAPPED_LIGHT_LOCALIZATION.md).
- **The distribution:** light_l ∈ [0, 1] per site, with TotalLight = ⟨n_XY⟩ the formalized
  field (`SlowLightDistribution`).
- Beneath all three sits rung 0: ⟨n_XY⟩ = 0 gets no light at all, immortal.

The split is the architecture, not a defect. **Reading happens in the bright letters** (they
absorb, they answer, they decay of it); **storage happens in the dark ones** (unreadable by
the light, surviving precisely therefore). The Hamiltonian is the head that rotates dark ↔
bright: L_H turns the sharp bit ⟨Δ_l⟩ ∈ {0,1} into the expectation ∈ [0,1], so let H run and
then read is how you recall what is written in the dark. γ₀ erases only the bright.

> *The lens carries the page; the light is the ink.*

> **Read-cost (HYPOTHESIS, untested).** If the dark sector is the disk and L_H the head,
> one recall rotates dark → bright (dwell ~ 1/J) while the bright pays the light (~2γ), so
> the dose cost of one read of the dark memory is ~2/Q. High-Q systems read their memory
> almost free; at the EP every read erases of order what it reads. Falsifiable form: measure
> cost-per-recall(Q) on the same apparatus as the FI(Q) curve; a flat or inverted dependence
> refutes it.

---

## The inside world

What is readable from inside the chain is a triple, plus a field:

1. **shape** - the internal, γ₀-free ratios (δJ/J, J_b/J_c): what you are;
2. **Q = J/γ₀** - the stand relative to the light, and simultaneously the reading budget
   (the FI(Q) law makes the coordinate and the budget the same number);
3. **K = γ₀·t** - felt time, the dose of the irreversible; t alone is as unreadable as γ₀,
   only the product K is experienced. Aging happens in K (see [`GAMMA_TIME_DISTINCTION.md`](../docs/GAMMA_TIME_DISTINCTION.md),
   [`TEMPORAL_SACRIFICE.md`](../experiments/TEMPORAL_SACRIFICE.md)).

The fourth coordinate is the light field itself, the bright/dark distribution of the
previous section: the state-side position to the carrier, read as a field of who answers and
who is written upon.

---

## Open, falsifiable

- **The ordering prediction.** The FI(Q) slope steepness is ordered by the readout basis's
  brightness (coherence steeper than population, tested in three bases: Z, X, Y). The full
  basis sweep is open; a basis with an EP peak would falsify the Q-factor argument.
- **The read-cost law** (~2/Q per recall): hypothesis above, untested.
- **Codebook completeness:** that every γ₀-relative position is readable by some handshake.
  A feature no handshake decodes refutes it.
- **The reading grammar** itself, the still-missing algorithm feature → handshake: invert
  the painters' α-profile → (bond b, strength δJ), and characterize the per-feature FI(Q)
  shapes. This is the heart the codebook frames, not a result in it.

The ancestral instance is already on record: [`GAMMA_AS_SIGNAL.md`](../experiments/GAMMA_AS_SIGNAL.md)
read 15.5 bits from the γ-profile, feature → reading before the codebook language existed.

---

## Live instruments

- `inspect --root decoder` - the FI(Q) resolution law, the measured tables above.
- `inspect --root clock` - the dial angle θ, the index that is Q.
- `inspect --root symphony` - the painters movement, the working proto-decoder: defect
  location reads in the mirror-breaking of the α-profile, strength in the f-scale.
- `inspect --root envelope` - the dose envelope in K, felt time.
