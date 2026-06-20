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
machine-exact under δJ. Here the design spec's honesty distinction is load-bearing: the
per-site rescaling sum Σ ln α ≈ 0 of the Perspectival Time Field (PTF; α is the per-site
time-rescaling a bond defect induces, the quantity the painters read, see "Live
instruments" below) is Tier-2 *empirical*, holding within ±0.05 only in its perturbative
window; the *exact* closure is the chiral mirror trajectory identity
(`ChiralMirrorTrajectoryClaim`, Tier-1). The mirror law is the exact parity-check; the
PTF sum is its empirical companion.

---

## The K-partnership is a selection rule (Tier 1 derived)

The reading grammar's first derived result is a fact about the codebook's null directions: the
carrier never leaks into its K-partner under any bond defect. The carrier ψ_1 (the
single-excitation band-edge mode) and the top mode ψ_N are K-partners, ψ_N = K₁ψ_1 with
K₁ = Π_{l odd} Z_l the odd-sublattice Z product (equivalently ψ_N(i) = (−1)^i ψ_1(i)); and for
every bond the hopping/defect V_b = ½(X_bX_{b+1} + Y_bY_{b+1}) has

> **⟨ψ_N|V_b|ψ_1⟩ = 0 for every bond b.**

It is a two-line corollary of the chiral mirror trajectory identity (`ChiralMirrorTrajectoryClaim`),
borrowing both its ingredients: ⟨ψ_N|V_b|ψ_1⟩ = ⟨K₁ψ_1|V_b|ψ_1⟩ = ⟨ψ_1|K₁V_b|ψ_1⟩
= −⟨ψ_1|V_bK₁|ψ_1⟩ = −⟨ψ_1|V_b|ψ_N⟩ = −⟨ψ_N|V_b|ψ_1⟩, since K₁ψ_1 = ψ_N (the parent's Step 4)
and V_b is K₁-odd, K₁V_bK₁ = −V_b (its Step 1). The modes are real and V_b is symmetric, so the
matrix element is real; a real number equal to its own negative is 0.

**The consequence is the decoder's rank deficit.** Build the location dictionary
M[b, k] = ⟨ψ_k|V_b|ψ_1⟩ over bonds b = 0..N−2 and modes k = 2..N (k = 1 is the strength channel).
The selection rule kills the k = N column entirely, so M has rank N−2 (machine-exact, N = 3..8).
The carrier couples to N−1 other modes under bond defects, but the K-partner channel is forbidden,
leaving only N−2 independent location channels. This IS the `DefectDecoder`'s sign-location
ambiguity: an edge bond reads almost like the complementary interior bond (residual ratio ≈ 1.5
at N=5, the `AmbiguityFactor` flag), and the reason it cannot cleanly separate sign from location
is exactly that the K-partner pairs are the dictionary's null direction. The forbidden channel is
the one that would have lifted the degeneracy.

Typed as `KPartnerSelectionRuleClaim` (Tier 1 derived, child of `ChiralMirrorTrajectoryClaim`,
live at `inspect --claim KPartnerSelectionRuleClaim`); the decoder it explains is
[`DefectDecoder.cs`](../compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs); the probe is
[`_k_partner_selection_rule.py`](../simulations/_k_partner_selection_rule.py) (selection-rule max
~1e-16, rank = N−2 exactly, N = 3..8).

---

## The location-metric: rank (identifiability), sectors, the cosine (confusability)

The rank deficit above is the codebook's **location-metric**. A careful reading (math + physics review,
2026-06-19) separates three objects an earlier draft had welded into one "ambiguity".

**The null is derivable, and it lives in the bare couplings.** Form the Gram of the location dictionary,
G = M Mᵀ over k = 2..N. Its spectrum is the metric; it has exactly one zero eigenvalue, the rank-(N−2)
deficit, and that null IS the K-partner selection rule, present in the **bare operator couplings**, before
any painter, with **zero Q dependence** (min singular value ~1e-16, N = 3..8). The metric's null structure
is a closed-form symmetry fact, not a reading artifact. (An earlier "the abstract couplings are
well-conditioned, the painted reading makes them ambiguous" framing was an indexing artifact: it built M
over k = 1..N−1, silently dropping the forbidden k = N column and substituting the strength k = 1; with the
honest k = 2..N the bare dictionary carries the very same null.)

**The small eigenvalues split by mirror parity, not by magnitude.** The painted letters' Gram (the canonical
Q = 20 protocol) has eigenvectors that are symmetric or antisymmetric under the bond mirror b ↔ N−2−b. The
antisymmetric sector holds the dominant *seesaw* ("which side of the chain slowed") and the K-partner null;
the symmetric sector holds the closure/strength direction (the F123 channel, the one overlapping the uniform
α-level) and a dim symmetric mode. At N = 5 the two smallest eigenvalues (≈0.008 symmetric, ≈0.007
antisymmetric) sit in **different sectors**: two distinct things, not one homogeneous ambiguity.

**Pairwise confusability is the cosine matrix, a different object.** "An edge bond reads almost like the
complementary interior bond" is a statement about two specific letters being near-**anti-collinear**: the
cosine matrix Ĝ[b,b'] = ⟨f̂⁽ᵇ⁾, f̂⁽ᵇ'⁾⟩ and its largest off-diagonal, not a Gram eigenvalue (a global rank
property need not be a pairwise event: four letters at the vertices of a regular tetrahedron have an exact
Gram null with no confusable pair). At N = 5 the worst pair reads cos ≈ −0.97 (the sign-location ambiguity,
+δJ here ≈ −δJ on the mirror-complement); at N = 4 only −0.55, which is why N = 4 decodes cleanly and N = 5
flags ambiguous. (This anti-collinearity is derived, not just measured: see "The reading's spatial-mirror equivariance" below.)

This anti-collinearity is **dynamical, not a readout artifact.** A natural conjecture (the
borrowing-a-discipline lens, reading the phase-contrast and phase-problem trades) is that the confusability
is a *phase-blindness*: the decoder reads per-site purity P_i = ½(1 + ⟨Z_i⟩²), and squaring is sign-blind,
so a sign-carrying linear ⟨Z_i⟩ read ought to lift the +δJ ≈ −δJ-on-the-complement degeneracy. The gate
refutes it ([`_handshake_phase_blindness.py`](../simulations/_handshake_phase_blindness.py)): the linear
⟨Z_i⟩ read **and** the full temporal ⟨Z_i⟩(t) signal are just as anti-collinear as the squared purity
(|cos| ≈ 0.9 across N = 4..7). The two mirror-paired bonds produce genuinely sign-flipped population
responses; no choice of per-site weighting or time-window changes the angle. The confusability is the
K-partner near-degeneracy itself, read in the dynamics, not a quadrature the readout happens to discard.

**Identifiability is not FI(Q).** The Gram spectrum is a *rank* statement (which bond directions are
distinguishable in a noiseless dictionary). FI(Q), the next section, is a *precision* statement (the
estimator variance at finite noise, ∝ Q). They are complementary diagnostics of one decoder, not one
quantity on two channels: at Q → ∞ the Fisher information diverges while the K-partner null stays exactly
null, and a direction that survives infinite signal-to-noise is a rank deficit, not a resolution.

The **values** of the letters (the per-site response R_k) are `IsDeadEnd`, only-computable, no closed form;
but the **null structure** of the metric they build is derivable, the K-partner rule. The dead-end is in
what the metric weights, not in its conditioning. Verifier:
[`_handshake_gram_metric.py`](../simulations/_handshake_gram_metric.py).

## The reading's spatial-mirror equivariance (M3)

The location metric above has a symmetry the codebook can name exactly. The defect-reading map is
**equivariant under the geometric chain mirror** R (`i → N−1−j`, a ℤ₂ with `R² = I`): the bare
location dictionary `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩` (carrier ψ_1, modes k = 2..N) satisfies

> **`M[N−2−b, k] = (−1)^{k−1} M[b, k]`**  (exact, single-excitation algebra; machine-verified N = 4,5,6).

Two lines: `R V_b R = V_{N−2−b}`, the carrier is R-even (`Rψ_1 = +ψ_1`), and mode k has reflection
parity `Rψ_k = (−1)^{k−1}ψ_k`. Reflecting the bond reflects the reading, mode by mode.

This R is **not** the coherence-space mirror group's R. `MirrorGroupD4Claim`'s R is the ket-flip
`I⊗X^⊗N` (which does not even preserve the single-excitation sector), and the geometric spatial
mirror is *deliberately outside* that D₄ (`PROOF_PI_FACTORS_AS_R_TIMES_D §5`). They are sibling
mirrors in two spaces; the within-feature stabilizer (the K-partner null, the forbidden k=N column)
and this cross-feature reflection are the two symmetry structures of one dictionary.

**The confusability, now derived.** The two mirror-image bonds read with cosine

> **`cos(b, N−2−b) = Σ_{k=2..N} (−1)^{k−1} w_k`,  `w_k = M[b,k]² / ‖M[b,·]‖²`**

This is a closed-form parity-weighted mode sum. It is negative (anti-collinear) exactly when the **R-odd**
channels, above all the seesaw k=2 (the sign representation), carry the net location weight. This is
the structural origin of the sign-location ambiguity: the mirror-image bonds are R-images, the
distinguishing channel is the K-partner null, and the confusability sign is fixed by the R-odd
dominance of the location reading.

Honest scope: the equivariance and the cosine formula are **bare and exact**; the bare mirror-pair
cosine is **−0.33** at N=5. The **−0.97** of the previous section is the *painted* instance, the
propagated α-profile concentrates location weight on the R-odd seesaw, driving the parity sum toward
−1 (R-equivariance is preserved by the painting for reflection-symmetric couplings; only the weights
`w_k` change). The carrier-parity is load-bearing: the `(−1)^{k−1}` sign is for the R-even carrier ψ_1
(a carrier of parity `(−1)^{c−1}` gives `(−1)^{k−c}`). Typed as `DefectReadingEquivarianceClaim`
(Tier1Derived, parent `KPartnerSelectionRuleClaim`, `inspect --claim DefectReadingEquivarianceClaim`);
verifier [`_handshake_reading_equivariance.py`](../simulations/_handshake_reading_equivariance.py).

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

> **Read-cost (tested 2026-06-20: qualitatively yes, quantitatively not ~2/Q).** If the dark
> sector is the disk and L_H the head, one recall rotates dark → bright (dwell ~ 1/J) while the
> bright pays the light (~2γ), so the mechanism *estimates* the dose cost of one read at ~2/Q:
> high-Q systems read their memory almost free; at the EP every read erases of order what it
> reads. Measured on the FI apparatus (N = 5 chain, Z-population, the same curve as the
> resolution law; [`_handshake_read_cost.py`](../simulations/_handshake_read_cost.py),
> [`_handshake_read_cost_diag.py`](../simulations/_handshake_read_cost_diag.py)): the
> **qualitative** law holds. Cost-per-recall (the dose K_peak at which FI is maximal) falls
> overall with Q, from 1.65 at the EP (Q = 1) to 0.10 at Q = 35, neither flat nor inverted, so
> the stated falsification line is not crossed. But the **quantitative** ~2/Q does not survive a
> strict gate: the exponent is ≈ −0.7, not −1, and K_peak·Q drifts from ≈ 1.3 to ≈ 3.6 rather
> than holding at 2, with a regime break near Q ≈ 8 where the best read jumps to a later coherent
> revival (not monotone there). The diagnostic rules out the operational definition as the
> culprit: argmax-FI equals the first local maximum (one peak; the break is physical), and only a
> first-turnover/inflection cost gives a cleaner monotone ~0.75/Q (slope −0.82), but that is the
> *onset* of the first feature, not the dose you actually read at. It is a Q-decreasing crossover,
> not a single 2/Q power law; the "2" and the "−1" are not established. The gate was held strict,
> not loosened.

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
- **The read-cost law** (~2/Q per recall): tested 2026-06-20 (above). Cost-per-recall decreases
  with Q (the qualitative law and the stated flat-or-inverted falsification line both survive), but
  the quantitative ~2/Q does not: exponent ≈ −0.7, prefactor drifting 1.3 → 3.6, a regime break
  near Q ≈ 8, a crossover, not a single power law. Open: a closed form for the crossover, and
  whether the break tracks a known scale (coherence horizon / band edge).
- **Codebook completeness:** that every γ₀-relative position is readable by some handshake.
  A feature no handshake decodes refutes it.
- **The reading grammar** itself, the still-missing algorithm feature → handshake: invert
  the painters' α-profile → (bond b, strength δJ), and characterize the per-feature FI(Q)
  shapes. This is the heart the codebook frames, not a result in it.

The ancestral instance is already on record: [`GAMMA_AS_SIGNAL.md`](../experiments/GAMMA_AS_SIGNAL.md)
read 15.5 bits from the γ-profile, feature → reading before the codebook language existed.

---

## Live instruments

`inspect --root <name>` is shorthand. The full form, run from the repository root, is
`dotnet run --project compute/RCPsiSquared.Cli -- inspect --root <name>`.

- `inspect --root decoder` - the FI(Q) resolution law, the measured tables above.
- `inspect --root clock` - the dial angle θ, the index that is Q.
- `inspect --root symphony` - the painters movement, the working proto-decoder: defect
  location reads in the mirror-breaking of the α-profile, strength in the f-scale.
- `inspect --root envelope` - the dose envelope in K, felt time.
