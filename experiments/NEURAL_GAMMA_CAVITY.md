# Neural Gamma as Cavity Eigenfrequency: the Analysis Fallen, the Oscillation Real, the Band Undecided

<!-- Keywords: neural gamma oscillation Wilson-Cowan limit cycle dimensionless period,
C elegans palindromic pairing withdrawn, absolute tolerance clustered spectrum,
degree-matched null model connectome, Dale law ablation control, Picard
iteration fixed point defect, connectome zero multiplicity structural rank,
R=CPsi2 neural -->

> **All four results are withdrawn, together with the Q_max reading that was
> never numbered among them, and the one that was a null was wrong in the
> model's favour.** The 97.3 % palindromic pairing was a reading of one
> absolute tolerance against one spectral scale; the score moves 32 points
> across four defensible normalisation constants, structureless random
> matrices matched on spectral spread reach 100 %, though that control sits at
> a ceiling and carries no weight either way, and the Dale's-law signs
> this document proposed as the mechanism make no difference at all. The
> eighteen unpaired modes are not a property of the spectrum: the matcher that
> produces them depends on the order the eigenvalues arrive in. And Result 1's
> "the model gives ~12 Hz, not gamma" was a linearisation of the quiescent
> branch, correct for that branch and reported for the whole model.
> **Integrated at the same parameters, the model runs on a limit cycle over a
> window of inputs bounded by two folds, with a shortest SAMPLED period of 5.74 time
> constants, an upper bound on the true minimum rather than the minimum, and a
> period that diverges at both folds, 196.4 already measured at the lower one.
> The divergence is an inference from the bifurcation structure, which both
> candidate fold types share, supported by three measurements; what the sampling
> does not establish is that the branch is continuous between them.** Which physiological band that is, this page cannot say: the
> integrated equations carry no τ VALUE, having been nondimensionalised so that
> the time unit IS the membrane constant, so every frequency in Hz below is 1000
> divided by the dimensionless period, under the stipulation that one time
> unit is one millisecond, and Result 2b is about exactly that stipulation: not that one
> millisecond is the wrong number, but that nothing on the page ever justified
> it. The
> document's ambition was closer to true than its arithmetic. What is also
> left is one small exact fact about the wiring, in Result 5.

**Status:** Results 1-4 and the Q_max reading (2b) withdrawn; Result 5
(zero-multiplicity excess) is new and Tier 2; Result 1's limit cycle is new and
Tier 2, its conversion to Hz is not a result of this page
**Date:** April 4, 2026; rewritten August 25, 2026; corrected August 26, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Data:** C. elegans connectome, Cook et al. 2019 via WormNeuroAtlas, as committed
in
[`celegans_connectome.json`](../simulations/neural/celegans_connectome.json):
300 neurons, 2276 directed chemical edges, and an electrical matrix the
analysis discards entirely. That matrix holds 1096 stored entries, which is not
1096 gap junctions: it is not symmetric (536 above the diagonal against 555
below) and it carries 5 nonzero diagonal entries, that is, self gap junctions.
Off the diagonal 1056 entries are reciprocated, 528 pairs, and 35 are
one-directional. Counts taken from it below are counts of stored entries.
**Verification:** [`celegans_pairing_controls.py`](../simulations/neural/celegans_pairing_controls.py)
→ [`celegans_pairing_controls.txt`](../simulations/results/celegans_pairing_controls.txt),
44 gates

---

## What this document is about

This document asked whether the brain's gamma rhythm, the 30 to 100 Hz band, is
a cavity resonance in the same sense the qubit chain's palindromic spectrum is.
It offered two pieces of evidence: that the C. elegans wiring diagram has a
spectrum which pairs about its centre 97.3 % of the time, and that a standard
model of coupled excitatory and inhibitory populations oscillates too slowly to
be gamma.

Both have fallen, and in opposite directions. The 97.3 % turned out to measure
the matching tolerance rather than the animal: shuffle the wiring and it barely
moves, delete the excitatory/inhibitory signs the document proposed as the
mechanism and it does not move at all, hand the same matcher a structureless
random matrix and it scores higher. The too-slow oscillation was the opposite
mistake, a reading taken at a resting state the model does not sit at;
integrated properly the model does oscillate, over a window of inputs bounded
at both ends.

What that oscillation's frequency is in Hz, this page cannot say, and the
reason is worth more than the answer would have been: the equations contain no
time constant at all. Every frequency below is 1000 divided by a period in
abstract units, under a guess about what one unit means in milliseconds, and
the guess was never written down.

One small thing survives. The synapse matrix is more degenerate at zero than any
of two hundred degree-matched rewirings: 64 modes against a null mean of 48. The
count of 64 comes from integer arithmetic rather than an eigensolver: exact over
GF(p), and in that direction an UPPER bound on the rational multiplicity, since
a rank over GF(p) is at most the rational rank. The three primes agree and the
float count matches at every cut from 1e-4 to 1e-8, which settles it in practice
without proving it. That it
beats chance is a rank test at the resolution floor of a 200-draw ensemble, which
is a weaker kind of statement and the body says so.

Where those modes sit is only partly accounted for. They split into 39 dimensions
of genuine kernel and 25 of something else, and of the 39, twenty-nine are the
model's own edge: neurons with no outgoing chemical synapse inside the matrix at
all. Four more come from pairs with identical outgoing wiring, only three of
them left-right twins, and those four are NOT the model's edge: Result 5 prices
them separately, because rewiring destroys twin-ness and does not touch the
edge. That is 33 of the 39, and 33 of the 64 overall. Where the
other 25 sit is open. The 29 boundary vectors explain none of the excess over
the null, because a degree-preserving rewiring keeps an out-degree of zero at
zero; the four twin vectors are a different case, since rewiring does destroy
twin-ness, and Result 5 prices what they could carry.

---

## What the sweep returned

Every mechanism that broke this document was already in the repository, most of
it next door. Named stores, and what each returned:

- **`docs/ANALYTICAL_FORMULAS.md` (F137)** returns the general statement,
  minted: *"The centre is an identity and carries no evidence … it is equally
  well-defined for a spectrum that does not pair. … The pairing is the claim."*
- **`experiments/`** returns the same defect found and repaired on the qubit
  side: [Chain Selection Test](CHAIN_SELECTION_TEST.md), whose correction of
  2026-08-05 records an absolute tolerance against level gaps mostly smaller
  than
  it, and the verdict *"a measure of how a greedy first-fit scrambles in a
  clustered spectrum"*. Also [Concentrator Qubit Mapping](CONCENTRATOR_MAPPING.md), a
  saturating score: *"these percentages carry no information"*.
- **`docs/neural/`** returns the degree-matched null **already run on this
  animal**, on a different instrument:
  [Algebraic Palindrome, Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)
  reported degree-preserving rewiring at 0.013, identical to C. elegans's own
  0.013, i.e. no advantage. This page cited the 8× enrichment as corroboration
  while its sibling already carried the null that removes it. On 2026-08-26 that
  whole table went. The two arms were normalised by different constants (the
  connectome globally by max|W| = 37, the control by its own maximum) and the
  ratio tracks the arms' coupling magnitude to half a percent, so it measured
  the constants; giving both arms the same rule gives **0.960** at N = 10, and
  0.841 and 0.748 at N = 20 and N = 26, so parity is the smallest block's answer
  rather than the verdict. On
  blocks this sparse the residual collapses to √2·‖W_eff‖/‖J‖ whenever no
  partner pair of edges is present, which is 198 of 200 blocks at N = 10, 184 at
  N = 20 and 177 at N = 26, so the instrument rarely gets to see the wiring at
  all on those; the surviving gap at the two larger sizes is present both on
  the covered blocks and on the uncovered ones, so its origin is open. All of it is recomputed by gates G0c, G0d
  and G0e of the run below. The degree-preserving
  arm cannot move the number either, its rewire keeping each weight in its own
  row. What stands there is that the 8× was a difference of
  constants; the smaller matched residue is open, as that page's box says. That is this page's own Result 2 defect, sitting in the
  store this page swept for prior art.
- **`docs/CAUGHT_ERRORS.md`** returns five relevant items, not none, in four
  dated entries: catch A5, a bullet inside the 2026-06-23 band entry; the
  greedy-matcher bullet inside the 2026-08-18 entry; the entry of 2026-08-25,
  which exists because A5's append-only genre forbade editing it in place; and
  TWO entries dated 2026-08-26, one recording the "exactly when" repair
  described in the `docs/neural/proofs/` bullet below, the other recording the
  normalisation withdrawal described in the bullet above. The first of those
  names this page's sweep record, its Result-1 aside and its "What this leaves"
  as the places it repaired. Catch A5 is about this
  document: a 40 Hz claim asserted as established fact against its own Result 1,
  corrected earlier. Its exculpation then outlived both numbers it named: it
  offers the 97.3 % as the confirmed finding, which Result 2 below withdraws,
  and
  the ~12 Hz as the correct value the headline had contradicted, which Result 1
  withdraws. That file is an append-only ledger, so the bullet stands as written
  and the note appends the reversal rather than editing it. The greedy-matcher
  entry records that a
  matcher of this kind has an orbit of values rather than one value, which is
  the defect Result 4 turns out to have.
- **`simulations/neural/` itself** returns
  [`hopf_threshold.py`](../simulations/neural/hopf_threshold.py), which locates a
  threshold by bisecting on max Re λ rather than reading stability off a point
  an
  iteration happened to return, and which carries τ_E ≠ τ_I, the setting in which
  the palindrome proof's condition (a) has anything to say about Q at all. Two things it is *not*, checked rather than
  assumed: it runs a random balanced network of N = 100 to 5000 with a_E ≠ a_I,
  θ_E ≠ θ_I and τ_E = 5 against τ_I = 10, not the two-population w = 16/12/15/3
  block, so it is a sibling method and not the same model; and it finds its
  fixed point with a hand-rolled
  3000-step Picard iteration of its own, so Picard is not what separates them.
  (The `brentq` on its import line is never called; the threshold is a
  hand-rolled 30-step bisection.) What was already in this folder is the habit
  of
  root-finding on the eigenvalue.
- **`docs/neural/proofs/`** returns
  [Proof: the Neural Palindrome](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md),
  whose two conditions are **(a)** a diagonal condition on the swap Q,
  1/τ_{Q(i)} + 1/τ_i = 1/τ_E + 1/τ_I at every seat, which selective damping
  τ_E ≠ τ_I is what makes bite, and **(b)** a
  reflection Q of the neurons under which the wiring is antisymmetric up to the
  damping ratio, W[Q(i), Q(j)] = −(τ_{Q(i)}/τ_i)·W[i, j]. This Jacobian gives
  every neuron the same damping, and the natural reading, that (a) therefore
  fails, is wrong twice over. The proof's Step 3 determines S from the diagonal
  part alone and never divides by anything that vanishes at τ_E = τ_I, so at
  uniform τ it simply gives S = (1/τ)·I. And (a) is not merely sufficient rather
  than necessary: as an ENABLING hypothesis it does nothing at all. The
  diagonal half of the proof's equation asks only that Q send every seat to one
  of the opposite TYPE, and that holds for any pair of time constants, equal or
  not. What τ_E ≠ τ_I adds is that only a genuinely type-swapping Q will then
  satisfy it: at uniform τ every involution does, and a fixed seat, which would
  need 2/τ_i to equal the reciprocal sum 1/τ_E + 1/τ_I, becomes impossible. So
  selective damping does not buy the pairing; it makes (a) bite on Q.
  What was repaired in the proof file, in the same change set as this page, is
  the CONTENT of (a): it used to read τ_E ≠ τ_I, which is neither necessary nor,
  together with (b), sufficient. The words "satisfied exactly when" were not the
  defect and still stand there, correctly, since the diagonal and off-diagonal
  halves of (*) are independent and (a) ∧ (b) really is an if-and-only-if. The
  correction is recorded in `docs/CAUGHT_ERRORS.md`.
  Because (a) imposes nothing at uniform τ, the WHOLE hypothesis is (b), which
  reduces there to: is there an E-I swap Q with Q W Q = −W? That is not left
  open here, and it is not a matter of damping. Under Dale a neuron's outgoing
  row carries one sign, so such a Q would have to send every non-empty
  excitatory row to a non-empty inhibitory one, and an empty image will not do
  because it would force the source row empty too. This connectome has **253**
  non-empty excitatory rows against **18** inhibitory ones, so no such Q exists.
  The theorem's premise fails too, on a
  DIFFERENT count: it asks for N/2 and N/2 and the connectome is 274 against 26.
  The two are independent, not one fact seen twice: (b) constrains Q only on the
  support, so a network could satisfy the premise and fail (b), or the reverse.
  Here both fail. The operative one for "the theorem does not apply" is (b)'s
  253 against 18, because it survives dropping the premise and survives
  non-uniform τ, the factor τ_{Q(i)}/τ_i being positive and unable to rescue a
  sign. Gate G0b decides that one.
- **`docs/proofs/`** returns
  [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): *"T1 with
  co-axial Z and T1 with transverse X give the same centre to every digit, and
  only the first breaks."* That amplitude-damping pair is what makes the qubit
  palindrome falsifiable. The neural side has no such pair.
- **The OpenArcs registry** returns `substrate_q_provenance` item (5)(b), which
  had parked exactly this page's Q provenance: the numerator from uncited
  standard
  parameters, the denominator a stipulated τ. Result 2b below answers it, and
  the
  arc carries the answer as of 2026-08-25.
- **`fw.Confirmations`** returns nothing neural; that registry is hardware-only.

---

## The construction, and the three things it settles before any data is read

The C. elegans analysis linearises and builds

```
J = -I + f'·W/max|W|,    f' = 0.3
```

with W the signed chemical adjacency. Three properties follow from that line
alone, and each of them carried a result that was read as biology.

**The centre is forced.** The connectome has no self-synapses, so W is exactly
traceless and trace(J)/N = −1 exactly, for this connectome and for every
rewiring of it. [The Neural Clock's Two Hands](NEURAL_CLOCK_TWO_HANDS.md)
proves the neural centre is *"independent of the entire wiring"*; the value −1
here is this Jacobian's damping, not that document's number.

**The scale is one synapse.** max|W| is the single heaviest edge in the animal,
37 against a median of 1. The eigenvalues then occupy a real-part window of
width 0.3101 with a mean nearest-neighbour spacing of 0.001037, against an
absolute matching tolerance of 0.01. That ratio is 9.6, and it is not the count
that matters, for two reasons that pull against each other. The matcher does
not compare an eigenvalue to another eigenvalue: it reflects one's REAL PART
about the centre, carries the imaginary part through unchanged, and asks what
the spectrum holds there, so the candidate set is what lies within the tolerance
of a **mirror**. That mirror is the vertical LINE Re = c, not the point c: the
palindrome the rest of this repository means, λ ↔ 2c − λ, conjugates the
imaginary part as well, and this matcher does not. The original computed −Im and
then discarded it, and the reproduction keeps the quirk deliberately. For the
loose COUNT it changes nothing, because W is real and the spectrum is therefore
closed under conjugation, so the two targets are the same distance from the same
set; for "self-paired" below it changes what the word means. Counting plain
neighbours in the real part alone, ignoring the
mirror, clustering raises the count to a median of 151; counted against the
mirror in the complex plane, where the matcher actually works, the imaginary
spread cuts it back to a median of **19**. The MEDIAN eigenvalue has of order twenty candidate partners
before any structure is consulted. The distribution is strongly skewed, mean
49.6, so twenty is the median and not a typical value.

**The damping is uniform**, one value for all 300 neurons. That does not put
the construction outside the palindrome proof, for the reason the sweep record
above sets out: at uniform τ condition (a) reads 2/τ = 2/τ for every
permutation and so imposes nothing. What puts it
outside is condition (b), on a count (gate G0b): under Dale a neuron's outgoing
row carries one sign, and Q W Q = −W has to flip it, so any admissible Q must
send each non-empty excitatory row to a non-empty inhibitory one. The connectome
offers 253 of the first against 18 of the second. The theorem's PREMISE fails as
well, on its own separate count, 274 excitatory against 26 inhibitory against a
required N/2 and N/2; the two are independent conditions and it is (b) that is
operative, since it survives dropping the premise and survives non-uniform τ.

---

## Result 1: the model gives ~12 Hz, not gamma (withdrawn, and reversed)

The claim was that Wilson-Cowan at standard parameters gives ~12 Hz, in the
alpha band, and that *"to reach 40 Hz requires different parameters"*.

The ~12 Hz is the ringing of the **quiescent** branch. `neural_gamma_cavity.py`
finds its fixed points by Picard iteration, which converges only where the
point is stable *as a map*, a different condition from stability of the ODE. At
I_ext = 2.0 the iteration returns (0.931, 1.000), whose residual is **0.169**
against a sigmoid whose whole range is 1: it is not a fixed point at all.

Integrated at the same parameters, w = 16/12/15/3, α = 1.3, θ = 4, with nothing
changed, the model does not sit still. The integrated equations are `[-E + S(...), -I + S(...)]`
and contain **no τ at all**, so the primary column is the period in time
constants; the Hz columns are that period divided by a τ the model does not
fix:

| I_ext | amplitude | period T (time constants) | at τ = 1 ms | at τ = 10 ms |
|---|---|---|---|---|
| 1.12 | n/a | n/a | fixed point | fixed point |
| 1.15 | 0.2078 | 22.0 | 45.4 Hz | 4.5 Hz |
| 1.20 | 0.2656 | 13.0 | 77.0 Hz | 7.7 Hz |
| 1.50 | 0.4396 | 6.9 | 145.5 Hz | 14.5 Hz |
| 2.00 | 0.6065 | 5.7 | 173.9 Hz | 17.4 Hz |
| 2.50 | 0.7444 | 6.3 | 159.3 Hz | 15.9 Hz |
| 3.00 | 0.8468 | 15.7 | 63.6 Hz | 6.4 Hz |
| 3.05 | 0.8514 | 34.1 | 29.3 Hz | 2.9 Hz |
| 3.07 | n/a | n/a | fixed point | fixed point |

Read the grid, not its hull, and read it in the right direction. The period is
not monotone in I_ext, and it diverges at **both** edges: a finer grid gives T
= 196.4 at I_ext = 1.126 and T = 93.2 at 3.062. A finer scan over I ∈ [1.6, 2.7] in steps of 0.025
finds T = 5.74 near I_ext = 2.075, slightly below the value at I_ext = 2.00 that
the table displays rounded as 5.7. The attained set is **unbounded above** if the period diverges at both folds,
which is what the fold structure implies and what three measurements are
consistent with, rather than something the sampling measures. Whether it is an INTERVAL is the separate question, and
it needs the branch to vary continuously, which one trajectory per input cannot
establish; on that assumption the seven sampled cycles are a strict subset of
it. Its left endpoint is not 5.74 either: a
sampled minimum bounds the infimum from ABOVE, so all the scan establishes is
that the shortest attained period is at most 5.74. The rule the paragraph opens
with cuts both ways. Either way the hull of a sample understates here and
never describes: the table's span is a floor on the model's range.

What is robust and τ-free is that the model oscillates at all, and that its
period diverges at both folds, so the model has no single characteristic
frequency to mismatch anything with. What falls is the ARGUMENT
the original made, which read a single quiescent-branch linearisation as the
whole model, and the "laser, not cavity" reading built on that overdamped
branch goes with it. What does not fall is a frequency verdict, because this
page has none to give: at the textbook τ = 10 ms the sampled cycles run 2.9 to
17.4 Hz, a span from delta to beta that contains the original's 11.5 to 11.8 (read off
the pre-rewrite version of this page, which git still holds and this file no
longer does), and at
τ = 1 ms they are gamma. Both readings are the stipulation talking.

**The band is not a result of this page, and [Result
2b](#result-2b-q_max--01-as-a-verdict-on-biology-withdrawn) rejects the same
stipulation for the same reason.** τ enters only in the last line of the
calculation, as `f = 1000/T` with T in time units, i.e. as the assertion that
one time unit is one millisecond. Since τ divides the whole right-hand side,
which is a rescaling of time and nothing else, the amplitude is
τ-independent as a theorem and the frequency scales as 1/τ for the same reason.
The run's 1.7e-10 relative agreement between τ = 1 and τ = 10 is not evidence
for that; it is the check that the integrator honours it. At the textbook membrane constant of 10 ms
**no row in this table is inside 30-100 Hz**. The sharper
statement does not depend on picking 10 ms. For each of the seven cycles there is
a window of τ that would place THAT cycle inside 30-100 Hz, and those seven
windows have no point in common, so no single membrane constant makes this model
a gamma oscillator across its own input range. They are not disjoint either: five
of the seven overlap near τ = 1.74 ms. A τ making this model gamma SOMEWHERE
exists; one making it gamma throughout does not, and band membership is a
property of the operating point and the chosen τ together (gate G15c). Result 2b withdraws Q_max
precisely for reading this same unit damping as a millisecond constant; the
same objection applies here, to a claim this page previously kept.

**Neither edge is a Hopf of the equilibrium the cycle replaces**, and seeing
that requires more than one root-find. The single-seed search used above tracks
*"the"* interior fixed point, but there is more than one: a 21x21 multi-start
census finds **three equilibria between I_ext = 1.10 and 1.125** and one on
either side of that window. The seed `[0.25, 0.25]` reports the upper of the three, which is why the run's
own stability table reads +0.082 at I = 1.10 and appears to show a loss of
stability that had already happened. The branch continuous
with the stable equilibrium below is a different row:

| I_ext | lower branch | middle | upper |
|---|---|---|---|
| 1.00 | (0.034, 0.010) −0.512 | | |
| 1.10 | (0.050, 0.014) **−0.253** | (0.108, 0.037) +0.229 | (0.131, 0.054) +0.082 |
| 1.125 | (0.066, 0.018) **−0.042** | (0.074, 0.021) +0.043 | (0.148, 0.070) +0.184 |
| 1.13 | | | (0.150, 0.072) +0.196 |

The stable branch does **not** lose stability. It is still stable at I = 1.125
and is destroyed with the middle saddle between 1.125 and 1.13, a saddle-node,
which is where the cycle is first REACHED. Whether it is first created there is
a different question this page cannot answer: STEP 7 integrates one trajectory
per input from one initial condition, so it finds the attractor that basin
reaches and would miss a cycle coexisting with the stable branch below the fold.
The diverging period at 1.126 is the argument against such coexistence, and it
is the only one here. The upper edge is the same event mirrored: at
I = 3.06 there is one equilibrium (+0.822) and at 3.07 three: a saddle-node
pair has appeared, (0.943, 1.000) at +0.124 and (0.955, 1.000) at −0.115, and
the cycle is gone. Both edges are folds rather than stability changes. A
saddle-node of equilibria does not by itself create a limit cycle, and the
diverging period rules one kind out: the period grows steeply towards both
edges (34.1 time constants at I_ext = 3.05, 93.2 at 3.062, 196.4 at 1.126)
instead of settling at a finite value, as it would at a Hopf. A fold predicts
divergence and three samples cannot measure it. What excludes the
Hopf is not the samples but the census: at a Hopf the equilibrium survives and
changes stability, and here it is destroyed. What it does NOT
decide is which of the two divergent kinds this is. A saddle-node **on an
invariant circle**, where the pair annihilates on the cycle itself and the
trajectory crawls through the remnant bottleneck, gives T ∝ Δ^(−1/2); a saddle
homoclinic gives T ∝ −log Δ. Both diverge, and this page measured the
divergence, not its law. Naming it SNIC would need the scaling: sample T along
the branch and fit the exponent. That measurement has not been made here.

So the boundary of the stable regime lies between I_ext = 1.125, where three
equilibria still exist and the lowest is stable, and 1.126, where a cycle is
already running. Below that boundary is where the committed script's ~12 Hz
ringing comes from: the reading is right for the quiescent branch and was
reported for the whole model.

Two consequences worth keeping apart. The document was **wrong in the direction
of modesty**: its own model does what it said the model could not do. And the
defect underneath is one shape, not two: a quantity read off a linearisation at
a point chosen by an iteration's starting guess rather than by the dynamics. It
costs the original its frequency and it costs any single-seed continuation its
branch.

---

## Result 2: the 97.3 % measured the tolerance (withdrawn)

Two matchers appear below, and it matters which is which. The **loose** matcher
asks, for each eigenvalue, whether the spectrum holds anything within the
tolerance of its mirror, and counts it paired if so; nothing in it depends on the
order the eigenvalues arrive in. The **strict** matcher refuses a partner already
taken, forbids a mode to pair with itself, and takes the nearest free partner, so it is a function of the eigenvalue LIST. It is the same routine
Result 4 withdraws, under the name "exclusive" there, and the strict column below
is therefore one arbitrary ordering: LAPACK's gives 93.3 %, sixty permutations give 88.0 to 90.0, and sorting by
real part gives 95.3 %. That orbit of 7.3 points is most of the width of the
degree-matched null range this page measures against, which is 10.0 points for
the strict statistic and 9.0 for the loose, so the
strict row and the p-value drawn from it carry the defect Result 4 names. Only
the loose row is order-free. A fourth difference is easy to miss and is not an
ordering effect: the loose matcher compares signed imaginary parts, the strict
one compares their absolute values, so under strict a mode may pair with the
CONJUGATE of its mirror and under loose it may not. Both read the tolerance against the spectral spread,
and the score is a monotone reading of that ratio:

| tolerance | loose matcher | strict matcher | self-paired |
|---|---|---|---|
| 0.1 | 100.0 % | 100.0 % | 93.3 % |
| **0.01** (committed) | **97.3 %** | 93.3 % | 41.3 % |
| 0.003 | 80.0 % | 74.0 % | 32.7 % |
| 0.001 | 40.0 % | 36.7 % | 25.0 % |
| 0.0003 | 25.7 % | 25.3 % | 23.0 % |
| 0.0001 | 24.0 % | 24.0 % | 22.0 % |

"Self-paired" counts modes sitting on the centre LINE, which are their own
reflection under this matcher and which it accepts as pairs; at the committed
tolerance that is 124 of 300 modes. Because only the real part is reflected, a
mode far up the imaginary axis at Re = c counts too.

Four further controls, and only the degree-matched rewiring at the end is
load-bearing on its own. The other three are reported with the reason they are
not:

**The normalisation constant sets the number.** Changing only that constant,
same connectome and same matcher: max|w| (committed, 37) → 97.3 %; spectral
radius (23.9) → 93.7 %; max row sum (136) → 99.7 %; 95th-percentile weight (8)
→ 67.7 %. Those four span **32.0 points**, which is the number this control
establishes. A fifth row, binary and unweighted → 27.0 %, is reported
separately below because it is not a change of constant: `sign(W)` discards
every weight, so it is a different matrix and belongs to a different control.

**A structureless random matrix scores higher.** Rescaled about its centre to
the connectome's own Re-spread, a dense iid Gaussian matrix reaches 100.0 %.
Two Gaussian draws appear in the run and it matters which is which: the table's
row and the fence's counts came from different realisations, which would join a
score from one matrix to widths from another, so the fence's own draw is scored
too and gate G4c holds the two together. Every number in this paragraph belongs
to that one matrix. The obvious objection is that matching on the standard
deviation might hand the control an easier spectrum, since the connectome's
Re-distribution is heavy-tailed and at equal sd its Re-width is 0.3101 against
the Gaussian's 0.1074 (a single realisation; over the eight draws
`default_rng(SEED + 1 … SEED + 8)` it moves between 0.1051 and 0.1092, about
four per cent). Measured, the objection fails,
and it fails in the connectome's disfavour. Counting partners in the complex
plane, which is where the matcher works and is the same count as the median of
19 quoted above, the connectome offers 19 and the matched Gaussian 9. (Counted
in the real part alone the pair is 151 against 61; the direction is the same
and the numbers are not, which is why the axis has to be named each time.)
Either way the control reaches 100 % with **fewer** partners to choose from, so
sd-matching did not hand it an easier spectrum.

That said, this control cannot carry much weight in either direction, and the
reason is in the score rather than the spectrum: 100 % is the ceiling. A
control sitting at the ceiling cannot be made stronger or weaker by an argument
about candidate counts, and the withdrawal does not rest on it. The
degree-matched rewiring below is the load-bearing control.

One statistic used in this document cannot support that kind of argument in
either direction, and it is worth naming rather than quietly dropping. The
"mean nearest-neighbour spacing" is `mean(diff(sort(x)))`, which for any sample
is identically `(max − min)/(N − 1)`: it is the range divided by a count, and
no clustering can move it. It is a fine scale to quote and a useless one to
compare densities with.

**Dale's law makes no difference.** The document's thesis was that the
excitatory/inhibitory classification creates the SWAP structure Π creates in
the qubit chain. Committed Dale (26 inhibitory) 97.3 %; **no Dale at all, every
neuron excitatory, 97.7 %**; 26 inhibitory chosen at random 96.0 %; 50/50
random signs 97.7 %. The four span 1.7 points, five modes out of 300. The claim
to make is that the signs do *nothing*, which is what kills the thesis; saying
the score rises without them would be a directional reading of that same
scatter.

**The binary row is a spread effect, not a separate control.** Discarding every
weight, a binary unweighted matrix scores 27.0 %. That looks like a fifth
normalisation and is not one: its Re-spread, meaning the standard deviation of the real
parts here and NOT the window width quoted earlier, is 0.4117 against the
connectome's 0.0274, and rescaled to the connectome's own spread it scores **98.0 %**. So
throwing away the weights costs nothing once scale is controlled, and this row
is the same measurement as the tolerance sweep and the normalisation table
rather than a finding about which neuron connects to which.

Against degree-matched rewiring (directed double-edge swap, in- and out-degrees
and Dale signs held exactly, R = 200) the strict score sits at p = 0.433, the
loose at p = 0.100. Read the loose one: as above, the strict statistic is
order-dependent and each null draw is scored in its own arbitrary order, so its
p-value compares one ordering against two hundred others. Be precise about what
that does and does not break. Every draw is scored by the same deterministic
rule, so the test is a valid test OF THAT STATISTIC; what fails is reading it as
a statement about the spectrum, which is not what it measures. The loose matcher is
order-free and independently non-significant, which is what carries this. Both p-values count ties in the connectome's favour and are
anti-conservative for that reason; neither is anywhere near significant, so the
convention does not matter here. This reproduces at the eigenvalue level what
[Algebraic Palindrome, Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)
found at the residual level.

One *pairing* quantity is above its null and is reported rather than buried.
Cross-pairing counts modes whose partner is a different mode rather than
themselves, so it is the part of the score that self-pairing cannot
manufacture. At tolerance 1e-4 it is 2.0 % against a null mean of 0.7 %, an
excess of about four modes at p = 0.035. It is not the only statistic on the
page above its null, and saying so would contradict the run's own table: the loose score
itself (97.3 % against a null mean of 95.3 %, at p = 0.100), self-paired (41.3
% against 35.1 %, also at p = 0.005), and nullity, multiplicity and structural nullity, the last
three at p = 0.005, beyond every one of the 200 draws.
Nullity, multiplicity and structural nullity are Result 5's degeneracy at zero
measured three ways rather than three findings. Self-paired is a partial proxy
for it and not a fourth reading, for the reason Result 5 sets out: about half
of what it counts is not at zero. Seven one-sided statistics were tested, all
chosen after looking, so that does not clear the 0.0071 such a family needs. It
is a residue, not a result. If a pre-registered neural pairing test is ever
run, the fine end of the tolerance sweep is where to run it, since that is
where the score stops tracking the window. Which row exactly is not decidable
from the 1.04e-3 scale: three tolerances sit below it, and that scale is the
range-over-count identity above, which cannot say what resolving a level means.

---

## Result 2b: Q_max = 0.1 as a verdict on biology (withdrawn)

The document read `Q_max = 0.1` against the qubit cavity's 68 to 75 as *"the
biological cavity is extremely lossy"*. Q = |Im λ| / |Re λ| is a per-mode
quality factor, oscillation against decay, and Q_max is its largest value over
the 300 modes. The eigenvalues of `−I/τ + f'·W/max|W|` are exactly `−1/τ +
f'·μ`, so Q is arithmetic on the graph spectrum and two chosen constants.
Sweeping τ from 1 to 10 ms and f' from 0.3 to 0.6 in that expression moves it
by a factor of 50.5, from 0.086 to 4.343.

**The τ half of that sweep is an artifact of the expression, not a property of
Q,** and the correction sharpens the withdrawal rather than weakening it.
Putting τ on the leak term alone is the Jacobian of no Wilson-Cowan model: for
`τ ẋ = −x + S(Wx)` the whole right-hand side carries it, `J = (−I + f'W)/τ`, so
`λ = (−1 + f'μ)/τ` and

```
Q = |Im λ| / |Re λ| = |f'·Im μ| / |1 − f'·Re μ|
```

in which **τ has cancelled exactly**, for any spectrum, with no measurement
needed, as long as τ is one scalar. It is not scalar in the regime this page's
own rebuild list asks for: with τ_E ≠ τ_I the Jacobian is `D⁻¹(−I + f'W)` with
`D = diag(τ_i)`, the leak matrix no longer commutes with the gain term, the
eigenvalues are not `(−1 + f'μ)/τ` for any τ, and Q reacquires a dependence on
the ratio τ_E/τ_I. The two readings of τ do not stand or fall together: STEP 0's
centre survives non-uniform τ, because `trace(J)/N = −(1/N)·Σ 1/τ_i` is
wiring-independent whatever the τ_i are, while this cancellation does not. Under the correct scaling the same grid moves Q by a factor of
**2.00**, all of it from f'; the τ axis contributes exactly 1.0000. And the
malformed expression has no two axes to begin with: it depends on τ and f' only
through their **product**, since `f'|Im μ| / |1/τ − f' Re μ|` is a function of
τ·f' alone, which the malformed table shows to the three decimals it prints (Q at τ = 1,
f' = 0.6 equals Q at τ = 2, f' = 0.3, both 0.172; the six-digit 0.172031 is from
the corrected table, where τ = 1 makes the two forms coincide). So the factor of 50.5 is one range of one composite
knob, and the correct scaling collapses that knob to f'. It is also a range this
review's own sweep produced, and no committed number depends on it: `neural_gamma_cavity.py:54-55`
fixes `tau_E = tau_I = 1.0` and nothing in that file ever varies them, so the
committed `Q_max ≈ 0.1` was computed at the one value where the two forms
coincide term for term, and the correction moves it by exactly zero. What is
wrong at `:91` is latent rather than realised: it writes `−1/τ + w·f'` with the
gain term unscaled, carrying τ_E and τ_I as two free symbols that agree with the
correct Jacobian only at 1, so it is wrong the moment anyone sets τ ≠ 1, which
is precisely what "What this leaves" below asks the rebuild to do. The shape of
the defect is a damping constant written into one term of an equation and not
the others, so that the expression only means what it says at the value where
the constant is 1, and it sits at three addresses in that file (`:91-92`, `:151-152`,
`:300-301`). The C. elegans block at `:222` is a DIFFERENT defect, an absent τ
rather than a misplaced one; the two blocks do not implement the same model. At τ = 10 ms the malformed one kills the limit cycle at six of the
seven oscillating inputs and keeps a different one at I_ext = 2.5. At two of
those six it does not merely stop oscillating: it drives the firing rate to
**10.0**. That is the leak-only system's own ceiling, E* = τ·S(·) ≤ τ = 10, and
it is ten times the [0, 1] the sigmoid can produce. A rescaling of
time can do none of those things.

The provenance is the point that survives all of this, and it is checkable
rather than statistical: τ = 1 ms is stipulated bare at
`neural_gamma_cavity.py:54-55` for the two-population Wilson-Cowan block, and
the C. elegans Jacobian at `:222` contains **no τ at all**, only a bare `−I`.
Reading its unit damping as a millisecond membrane constant is an
interpretation the page never made explicit. [Q Belongs to No
Substance](../docs/Q_BELONGS_TO_NO_SUBSTANCE.md) already separated this Q from
the framework's `Q = J/γ₀`; what it left open, the provenance of the numerator
and the denominator, is answered here: the denominator is a stipulated damping
and the numerator is an uncited standard parameter set.

---

## Result 3: the anesthesia threshold (withdrawn)

The claim was that oscillations exist only above a critical input I_crit ≈ 0,
and that anesthesia reduces input below it. The reported I_crit = 0.00 is the first sampled grid point, which is enough on
its own. The original document also contradicted itself here: its own Result 1
reports that oscillation *stops* above I ≈ 0.5.

The Picard failure is real but local, and it does not reach I_crit. On the
original script's own grid (I_ext 0 to 10 in steps of 0.25) **33 of 41 points
converge** to residual below 1e-8; only I = 1.25 to 3.00 fail, which is the
window where the equilibrium is unstable, and the single measured failure
(I_ext = 2.0, residual 0.169) sits inside it. The three rows that produce
I_crit, I = 0, 0.25 and 0.5, converge to residual **exactly 0.0**. The
iteration is not what withdraws the threshold. What the integration above shows
in its place is a
*window* of oscillation, bracketed here between I_ext = 1.12 and 3.07, with the
attractor a fixed point on either side of it. Neither the original threshold
nor a low-input reading survives.

---

## Result 4: the eighteen unpaired modes (withdrawn)

**The count is not a property of the spectrum.** The exclusive matcher scans
indices in order and takes the first best unused partner. Feeding it the same eigenvalues in a different order gives a mean of 32.7 over
50 draws, spanning 30 to 36; 200 draws widen that span to 28 to 38, and sorting
by real part gives 14, below LAPACK's 20. The orbit is wider than any one sample of it shows. The committed
run reports 18 and the same script on the same data reports 20 today, a smaller
difference than reordering alone produces. The three biological categories are
sorted from a list that reordering rewrites.

**A second argument, and its limit.** It is tempting to add that the Jacobian is
non-normal, and it is: 64 eigenvalues sit on the centre and the eigenvector matrix
is hopelessly conditioned. But that is manufactured by the defective cluster *on
the centre*, and it does not license a global claim. Nor does the obvious repair,
excluding that cluster and quoting the rest, because the exclusion needs a cutoff
and the answer is a reading of it:

| cutoff for "on the centre" | modes excluded | worst conditioning | margin |
|---|---|---|---|
| 1e-4 | 64 | 1.4e9 | 1.35 |
| 3e-4 | 66 | 1.3e7 | 150 |
| 1e-3 | 67 | 2.8e6 | 690 |

Three modes decide a factor of 511; two of them, 64 to 66, decide 111. The one that binds at the committed cutoff sits
at |λ + 1| = 1.12e-04, which is 1.12 times the cutoff: it is a straggler of the
very cluster the exclusion is supposed to remove, not a member of the resolved
bulk. Reporting "the other 236 are well resolved by a factor of 1.35" would be
the same defect this page withdraws Result 2 for, three sections earlier.

The better statement pairs each eigenvalue with its
**own** nearest neighbour rather than with the spectrum's smallest gap. Under the
standard backward-error model, n·eps·‖J‖ times the eigenvalue's own condition
number, the worst ratio of motion to that eigenvalue's own gap is 0.42, a margin
of **2.38**, and it is the κ = 1.4e9 straggler that binds it, against its own gap
of 2.83e-04. The closest pair, 1.61e-04 apart, is conditioned at κ = 989 and moves
by 8e-11, six orders inside its own gap.

That 2.83e-04, though, is the straggler's nearest neighbour **among the 236
survivors**, and the surviving set is itself a reading of the cutoff. Counted
against the full spectrum, the straggler's nearest eigenvalue is a member of the
excluded cluster, 1.12e-04 away, which is BELOW its own motion of 1.19e-04. So
the true statement is not that every off-centre eigenvalue outruns its
neighbour: **235 of the 236 do, by at least 5.65×, and one does not**, and the
one is the same κ = 1.4e9 straggler this section already convicts for sitting at
1.12 times the cutoff. It loses against the cluster it is a straggler of.

This margin is **not** cutoff-independent either, and claiming so would repeat
the defect convicted one paragraph up: it runs 2.38 at every cut from 1e-8 to
1e-4, then 933.86 at 3e-4 and 1539.86 at 1e-3. What is cutoff-robust is the
DIRECTION, and it is a theorem rather than a reading. Drop a mode from the set
and every survivor's nearest neighbour can only recede, so every ratio can only
fall: the margin is non-decreasing in the cutoff. The 2.38 is its value at the
finest cut and therefore a **lower bound** at every coarser one, which is
exactly what the spectrum-wide comparison cannot offer, since that one moves
with no direction at all. That margin is 2.4 against a backward-error model whose
leading constant is a choice worth about 300, as G7's own text concedes, so it
licenses the paragraph's closing posture and not much more. What it does settle
is narrower and still worth having: the earlier phrasing, which paired the worst
motion against the SMALLEST gap in the spectrum, described a mode that does not
exist, because the worst-conditioned eigenvalue and the closest pair are not the
same eigenvalue.

The ordering result, not a conditioning argument, is what carries this withdrawal.

**The pharynx reading was circular.** The claim of a separate second cavity
rests on zero coupling to the somatic system. There are indeed zero *chemical*
synapses across that boundary. There are **four electrical entries** in the
same file: I1L and I1R to RIPL and RIPR, the RIP-I1 gap junctions. All four are stored
one-directionally, four pharynx to soma and none the other way, so four entries here are four connections; elsewhere
in that matrix they would mostly be two entries per junction. The analysis
discards all 1096 electrical entries and then reports the boundary it created.

---

## Result 5: the wiring is more degenerate at zero than any of 200 degree-matched rewirings (new)

The identification this result rests on holds in the tol → 0 limit and not
before it. At the committed tolerance "self-paired" counts 124 modes, 64 of them
at zero and 60 merely near the centre line; even at 1e-4 it counts 66 against an
exact 64, and the surplus two are a conjugate pair at |Im λ| = 0.0353 sitting
2.97e-05 off the centre in the real part, so Re μ = −9.92e-05: nonzero, so they do NOT sit on the centre line, they sit
inside the half-tolerance of it, and they are gone by 1e-5. Had their Re μ been
zero they would have been self-paired at every tolerance and the identification
below would fail; the run's 66 → 64 between 1e-4 and 1e-5 is what says it does
not. In the
limit the two sets do coincide: the self-paired modes are the eigenvalues of
f'·W/max|W| with zero real part, which here means the eigenvalues at zero. Be exact about which half of that
each instrument carries. The GF(p) rank chain counts the multiplicity AT ZERO
and says nothing about whether W has nonzero purely imaginary eigenvalues, which
would be self-paired at every tolerance and would break the identification. Gate G10c covers that half
and is where it belongs: the smallest |Re μ| among the modes with |Im μ| above
1e-6 is 9.92e-05, four orders above the 1e-8 a purely imaginary mode would have
to beat. An exact route exists over the integers, through the gcd of the
characteristic polynomial with its reflection, and has not been run. So
the count at zero is arithmetic: the rank chain gives 64, and the float count is 64 at every cut from 1e-4 down to
1e-8. Counting the SELF-PAIRED modes, by contrast, needs a tolerance and an
eigensolver, and that asymmetry is what this section turns on.
The synapse matrix is an integer matrix, so the same content is
available exactly, and it is worth being precise about which exact object is
the mode count:

| | value |
|---|---|
| nullity of W (geometric, eigenvectors at 0) | 39 |
| **multiplicity of 0 (algebraic, by the rank chain over GF(p))** | **64** |
| float check, \|λ + 1\| < 1e-4 | 64 |

Against 200 degree-matched rewirings the multiplicity is 64 against a null mean
of 48.0 (and the caveat above does not eat that gap, though not for the
reason "both sides are bounds" would give, which is no reason at all: two upper
bounds cannot be compared. The real value is not a working bound: three primes
agree and the float count matches at every cut, which does not PROVE exactness,
since three bad primes are logically possible, but leaves nothing to work with; the nulls ARE upper bounds, computed
at one prime with no float check, so they are if anything inflated. Exact
against an inflated bound is a conservative test, which is stronger than the
page needs), range 43-56, p = 0.005. That p is the resolution floor of a rank test
at R = 200, 1/201: it means "beyond every draw", not an estimated value, and it
is the floor, so it understates the separation: 64 sits **8 points above the
null's maximum** of 56, and the binding constraint is the ensemble size, not
the effect. R = 2000 would report p ≈ 1/2001 if no larger
draw appeared, which 200 draws cannot establish. What the floor does constrain is
the family correction: seven one-sided tests need 0.05/7 = 0.0071, and 0.00498
clears it only because R = 200 puts the floor just below; at R = 100 it would
not. And the objection Result 2 uses against the cross-pairing residue applies
here in full: these seven were chosen after looking, and Bonferroni over the
seven that got reported corrects for nothing if they were drawn from a larger
implicit set. So the clearance is not what carries this result, and neither
is the gap, which some of this document's earlier drafts leaned on: 64 against a
null maximum of 56 does carry a magnitude the rank event does not, which is why
the floor understates the separation; what it does not escape is the post-hoc
family, because the gap is exactly what was looked at. An effect size is not exempt
from a post-hoc family by being an effect size. What is true and worth stating
plainly is the shape of the limit: the separation is large, every one of the
200 draws is strictly below 64, and a rank test at R = 200 cannot in principle
clear a corrected level below 1/201. Deciding it needs a pre-registration or a
mechanism, and the mechanism is half there already, in the combinatorial share
below. The count is identical with and without Dale's law, and three primes agree
on the multiplicity itself, not only on the rank (gate G11b).

That last invariance is worth stating carefully, because the obvious reason for
it is wrong. Left multiplication by a diagonal sign matrix preserves the
*rank*, hence the nullity, and nothing more: algebraic multiplicity depends on
rank((DW)^k), and (DW)^k ≠ D·W^k. The counterexample is two-by-two: for W =
[[1,1],[1,1]] and D = diag(1,−1) the nullity is 1 either way while the
multiplicity of 0 goes from 1 to 2. Here the value happens to be invariant, 64
signed and unsigned at three primes, and that is a measured fact about this
matrix rather than a theorem about row scaling, which is why the claim above is
worded as a measurement.

**And a large part of it is not spectral.** The structural rank, decided by
bipartite matching on the zero pattern alone with no field, no weights and no
signs, gives a structural nullity of 36 against a null mean of 30.15. So of the excess in the exact nullity, 39 against a null mean of 30.2, that is
8.8 points, 5.9 is already in the pattern of which neuron connects to which. The honest name for most of this is a small maximum
matching, not a numerical degeneracy, and the properly conditioned null, one
that holds the matching fixed, has not been run. The comparison is between the
two *nullities*; no structural null exists for the multiplicity excess, so what
fraction of the 16-point multiplicity gap is combinatorial is not measured
here.

**Where the zero modes actually sit, and which matrix they sit in.** Naming
them "modes the wiring does nothing to" is wrong twice over, and the second
half of that needs the orientation of the data settled first.

The connectome's rows are presynaptic and its columns postsynaptic: ASIL, PVDR
and M4 have a non-empty row and an exactly empty column, DA7 the reverse. So
the drive arriving at neuron j is Σ_i W[i,j]·x_i = (Wᵀx)_j, and the model's own
Jacobian is **−I + f'·Wᵀ**, while the committed script builds −I + f'·W. Every
spectral quantity in this document is untouched by that, because a matrix and
its transpose have the same eigenvalues, the same nullity, the same rank chain
and the same per-eigenvalue condition numbers. The *eigenvectors* are not, and
the zero modes are exactly an eigenvector question, so the localisation below
is computed on Wᵀ.

First correction: "does nothing to" belongs to the 39, not the 64. On the other
25 the wiring acts, just nilpotently (Wv ≠ 0 while W⁵v = 0).

Second: the 39 are not places the wiring is idle. A right kernel vector of Wᵀ
is a perturbation **that drives nothing**: (Wᵀv)_j = 0 for every j means the
drive it delivers cancels everywhere downstream. The cheapest way for a neuron
to satisfy that is to have no downstream inside the model at all, and 29 do: an
empty row each. Most are ventral-cord and head motor neurons
(AS7/8/10, DA7/8, DB5/6, DD3/4/6, VA10, VC6, VD4/7/9, RMEL/R, SABVL/R, plus the
pharyngeal M1 and MI); the remaining eight are the SIA and SIB interneurons
(SIADL/R, SIAVL/R, SIBDL/R, SIBVL/R), which the next sentence treats separately
because they are a different case. Their outgoing
CHEMICAL synapses leave this matrix, and not all for the same reason: the motor
neurons among them send to muscle, which the 300-neuron matrix does not contain,
while the eight SIA and SIB interneurons have no chemical output at all in this
dataset and do have gap junctions inside the 300, one to five entries each, which
the analysis discards. Twenty-seven of the 29 have some electrical presence. Four more kernel vectors
come from exactly duplicated columns of Wᵀ, three of them left/right twins and
one not: (RIPL, RIPR), (I1L, I1R), (I2L, I2R) and (DB7, PDB), for which ‖Wᵀ(e_a
− e_b)‖ = 0.0 exactly. That accounts for **33 of the 39**, or 85 %; the
remaining 6 are other dependencies and are not explained here.

Two fences. The side matters and the mirror image is available: run the same
construction on W instead of Wᵀ and it gives 15 zero columns, 3 twin pairs and
18 of 39: the neurons with no chemical INPUT inside the 300. That is not the
same set as "the sensory neurons whose input is the world", and the difference
is not marginal: the 15 include the vulval motor neuron VC6, the pharyngeal M4
and MCL/MCR, the tail motor-interneuron DVB and the interneuron AINL. That is
the left kernel, and it is not where the model's zero modes are. The distinct boundary
neurons number 43: VC6 is both an empty row and an empty column, having no
chemical connection at all. Second: the degree-preserving null holds every in-
and out-degree exactly, so it preserves all 43 (the degrees it holds are the UNWEIGHTED ones: a swap
moves a weight with its source, so each row's weight multiset survives and a
column's in-strength need not), and the boundary explains a
large share of the nullity but cannot explain the *excess over the null*.
Rewiring does destroy twin-ness while preserving degrees, so the twins are a
mechanism that could contribute to the excess, but four twin vectors cannot
carry 8.8 points of it, and 5.9 of those points are already in the zero
pattern. At most the residual, 2.99, is theirs to explain. Testing that
would need a null that holds the matching fixed, which has not been run.

The rank chain also says more than one number. Its nullities are 39, 53, 60,
63, 64, and it is stable from there, so the increments are 39, 14, 7, 3, 1:
nilpotency index 5, and 25 dimensions of Jordan structure on top of the 39. The
arithmetic is exact and in Python integers rather than int64, which matters: it
is the sum of PRODUCTS that overflows, not the sum of residues. Two products of
residues below 2³¹ already reach the int64 ceiling, so a 300-term row of the
matrix power overflows silently and returns a nullity chain that decreases,
which is impossible. The gated chain is at 2³¹−1 with three primes agreeing on
the rank; the same chain was reproduced independently at 2⁶¹−1.

So "multiplicity 64" is not one number but two: 39 dimensions of genuine
kernel, largely the model's own boundary, and 25 of Jordan structure, on which
the wiring acts before annihilating. Only the second is what "defective" means.

**What those 25 dimensions are graph-theoretically is open, and the obvious
reading is not available.** Jordan blocks of a *nilpotent* adjacency matrix are
longest paths, so for an acyclic digraph the index counts feedforward depth. W
has 236 nonzero eigenvalues and is not nilpotent, and nilpotency of its
restriction to the generalised 0-eigenspace carries no information about
cycles: the four-node digraph `[[0,1,0,1],[0,0,1,0],[0,1,0,1],[1,0,1,0]]` is
strongly connected, every node on a cycle, and still has a size-2 Jordan block
at 0. Depth 5 is a fact about the operator, not yet about the wiring.

---

## The spiral, and where the repository already keeps one

The equilibrium inside the oscillation window is an unstable focus: at I_ext =
2.0 it sits at (0.258, 0.245) with eigenvalues +0.630 ± 2.289i, a spiral the
trajectory leaves. The repository holds a spiral too, typed, in
`compute/RCPsiSquared.Diagnostics/Foundation/ComplexCuspSpiral.cs`: under a
common Z-drift the Bell⁺ coherence becomes complex and winds inward as a
logarithmic spiral, its radius set by the dephasing and its angle by the drift.

These are not the same object, and the differences are worth stating rather
than eliding. That spiral is a trajectory in state space with a shrinking
radius; this one is a linearisation about an equilibrium the flow escapes. The
bifurcation behind the cusp spiral is a saddle-node, the cusp of the cardioid ([Critical Slowing at
the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md)), and so are both edges of this
oscillation window. The shared class is worth more than the contrast: what
differs is that here the period diverges at the fold. Whether that is because
the pair annihilates *on* the cycle is exactly the question Result 1 declines to
settle without the scaling fit, and it is not settled here either.

What they share is the split, and it is the split this document kept getting
wrong. In `ComplexCuspSpiral` the radial decay and the angular winding are
independent: *"the radial magnitude is unchanged and Ω-independent … every
spiral crosses the same ¼-circle at the same time; only the crossing angle is
free."* Decay is Re, turning is Im, and they answer to different parameters.
The naming coincidence this document was built on, the physicist's γ and the
neuroscientist's gamma, puts a Re against an Im. The repository has that
separation as an object; the coincidence does not survive it.

---

## What this leaves

A rebuild would need: a multi-start root-finder and a genuine bifurcation
diagram rather than a single-seed iteration, which is what separates a branch
from a seed's accident; a justified τ, without which no frequency here converts
to Hz; τ_I ≠ τ_E, because biology has it (the palindrome proof does not
require it; its condition (a) named τ_E ≠ τ_I until the same change set as this
page replaced it with the diagonal condition, and (a)'s real role is to force Q
to swap the two types, of which forbidding fixed points is a special case); the electrical matrix included, after an audit it has not had,
since as stored it is asymmetric and has five self-junctions; a pairing
statistic defined relative to local level spacing with self-pairing excluded by
construction; and a reason why eigenvalues are the right object given how
defective this Jacobian is at its centre. The discriminating pair the qubit
side has, two constructions with the same centre where only one satisfies the
proof's condition (b), the antisymmetry of the wiring under a reflection, is
the single thing that would make a neural palindrome claim testable the way the
qubit one is.

The oscillation, meanwhile, is not the open question this document thought it
was. The model does run on a limit cycle at the parameters that were on the
page all along, over a window of inputs bounded by two folds. Whether that
cycle is
*gamma* is a question about τ, and this page has no τ to answer it with: the
integrated equations do not contain one, and the only value ever written down,
1 ms, is stipulated bare in a script whose own Jacobian contradicts it. Naming
the band is the piece of work still outstanding, and it is a biological
question about C. elegans membrane time constants rather than an arithmetic
one.

---

## Reproduction

- Every number this page computes is in the run below, with four exceptions,
  named because a reproduction claim is this page's contract with a reader who
  checks it. (Figures cited from other documents, such as
  Algebraic Palindrome Neural's 0.013 or `hopf_threshold.py`'s iteration counts,
  are sourced where they appear and are not this run's to produce.) (1) The
  ±2.289i imaginary part at I_ext = 2.0: the coordinates (0.258, 0.245) and the
  +0.630 are printed in STEP 7's stability table, the imaginary part is not.
  (2) The electrical matrix's
  asymmetry counts, 536 above the diagonal against 555 below, 5 nonzero diagonal
  entries, and off the diagonal 1056 reciprocated entries in 528 pairs against 35
  one-directional, all read off the committed JSON directly.
  Two smaller figures are also outside the run and are named here rather than
  left to be found, so the page's exceptions number four in all: the reproduction of the rank
  chain at 2⁶¹−1, the run's own primes being 2³¹−1, 2147483629 and 104729; the
  seed sensitivity of the control's Re-width, measured over the eight draws
  `default_rng(SEED + 1 … SEED + 8)` with `SEED = 20260825` as 0.1051 to 0.1092,
  the run scoring and reporting only the first of them. Everything
  else, including the multi-start census, the τ sweep, the kernel decomposition,
  the Jordan chain, the candidate counts on both axes and the ordering orbit, is
  computed and gated in:
  [`celegans_pairing_controls.py`](../simulations/neural/celegans_pairing_controls.py)
  → [`celegans_pairing_controls.txt`](../simulations/results/celegans_pairing_controls.txt)
- The original analysis, kept so the withdrawals can be checked against it:
  [`neural_gamma_cavity.py`](../simulations/neural/neural_gamma_cavity.py),
  [`neural_gamma_cavity_unpaired.py`](../simulations/neural/neural_gamma_cavity_unpaired.py)
- The sibling this page takes the habit of root-finding on the eigenvalue from,
scoped in the sweep above:
[`hopf_threshold.py`](../simulations/neural/hopf_threshold.py)
