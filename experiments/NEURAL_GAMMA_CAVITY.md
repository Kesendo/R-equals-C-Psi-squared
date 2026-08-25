# Neural Gamma as Cavity Eigenfrequency: the Analysis Fallen, the Gamma Band Reached

<!-- Keywords: neural gamma oscillation Wilson-Cowan Hopf limit cycle,
C elegans palindromic pairing withdrawn, absolute tolerance clustered spectrum,
degree-matched null model connectome, Dale law ablation control, Picard iteration
fixed point defect, connectome zero multiplicity structural rank, R=CPsi2 neural -->

> **All four results are withdrawn, and the one that was a null was wrong in the
> animal's favour.** The 97.3 % palindromic pairing was a reading of one absolute
> tolerance against one spectral scale; the score is tunable from 27 % to 99.7 %
> by the normalisation constant alone, structureless random matrices matched on
> spectral spread reach 100 %, and
> the Dale's-law signs this document proposed as the mechanism make no difference
> at all. The eighteen unpaired modes are not a property of the spectrum: the
> matcher that produces them depends on the order the eigenvalues arrive in. And
> Result 1's "the model gives ~12 Hz, not gamma" was a linearisation taken at a
> point the model does not sit at. **Integrated properly, at the same parameters,
> the model has a Hopf bifurcation and a limit cycle covering 64 to 174 Hz, with
> two of its own operating points inside the gamma band.** The document's
> ambition was closer to true than its arithmetic. What is also left is one small
> exact fact about the wiring, in Result 5.

**Status:** Results 1-4 withdrawn; Result 5 (zero-multiplicity excess) is new
and Tier 2; the Hopf reading of Result 1 is new and Tier 2
**Date:** April 4, 2026; rewritten August 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Data:** C. elegans connectome, Cook et al. 2019 via WormNeuroAtlas, as committed
in [`celegans_connectome.json`](../simulations/neural/celegans_connectome.json):
300 neurons, 2276 directed chemical edges, and a 1096-entry electrical matrix the
analysis discards
**Verification:** [`celegans_pairing_controls.py`](../simulations/neural/celegans_pairing_controls.py)
→ [`celegans_pairing_controls.txt`](../simulations/results/celegans_pairing_controls.txt),
18 gates

---

## What the sweep returned

Every mechanism that broke this document was already in the repository, most of
it next door. Named stores, and what each returned:

- **`docs/ANALYTICAL_FORMULAS.md` (F137)** returns the general statement,
  minted: *"The centre is an identity and carries no evidence … it is equally
  well-defined for a spectrum that does not pair. … The pairing is the claim."*
- **`experiments/`** returns the same defect found and repaired on the qubit
  side: [CHAIN_SELECTION_TEST](CHAIN_SELECTION_TEST.md), whose correction of
  2026-08-05 records an absolute tolerance against level gaps mostly smaller than
  it, and the verdict *"a measure of how a greedy first-fit scrambles in a
  clustered spectrum"*. Also [CONCENTRATOR_MAPPING](CONCENTRATOR_MAPPING.md), a
  saturating score: *"these percentages carry no information"*.
- **`docs/neural/`** returns the degree-matched null **already run on this
  animal**, on a different instrument:
  [ALGEBRAIC_PALINDROME_NEURAL](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)
  reports Erdős-Rényi 0.108 against degree-preserving rewiring 0.013, identical
  to C. elegans's own 0.013, and concludes *"the degree distribution fully
  explains the palindrome advantage"*. This page cited the 8× enrichment as
  corroboration while its sibling had already dissolved it.
- **`docs/CAUGHT_ERRORS.md`** returns **two** entries, not none. A5 is about this
  document: a 40 Hz claim asserted as established fact against its own Result 1,
  corrected earlier, and itself now carrying a second error: it names *"the
  confirmed finding is the 97.3% palindromic cavity structure"*, which is the
  claim Result 2 below withdraws. And the greedy-matcher entry records that a
  matcher of this kind has an orbit of values rather than one value, which is
  the defect Result 4 turns out to have.
- **`simulations/neural/` itself** returns
  [`hopf_threshold.py`](../simulations/neural/hopf_threshold.py), which
  integrates this same Wilson-Cowan model with a real ODE solver and a real
  root-finder. The Hopf bifurcation was understood in this folder while the
  cavity script beside it hand-rolled a Picard iteration and read the stability
  backwards.
- **`docs/neural/proofs/`** returns
  [PROOF_PALINDROME_NEURAL](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md),
  whose condition (a) is selective damping, τ_E ≠ τ_I. This Jacobian gives every
  neuron the same damping, so the proof does not reach this construction; its
  content sits in condition (b), which nothing here tests.
- **`docs/proofs/`** returns
  [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): *"T1 with
  co-axial Z and T1 with transverse X give the same centre to every digit, and
  only the first breaks."* That amplitude-damping pair is what makes the qubit
  palindrome falsifiable. The neural side has no such pair.
- **`fw.Confirmations`** returns nothing neural; that registry is hardware-only.

---

## The construction, and the three things it fixes before any data is read

The C. elegans analysis linearises and builds

```
J = -I + f'·W/max|W|,    f' = 0.3
```

with W the signed chemical adjacency. Three properties follow from that line
alone, and each of them carried a result that was read as biology.

**The centre is forced.** The connectome has no self-synapses, so W is exactly
traceless and trace(J)/N = −1 exactly, for this connectome and for every
rewiring of it. [NEURAL_CLOCK_TWO_HANDS](NEURAL_CLOCK_TWO_HANDS.md) proves the
neural centre is *"independent of the entire wiring"*; the value −1 here is this
Jacobian's damping, not that document's number.

**The scale is one synapse.** max|W| is the single heaviest edge in the animal,
37 against a median of 1. The eigenvalues then occupy a real-part window of width
0.3101 with a mean nearest-neighbour spacing of 0.001037, against an absolute
matching tolerance of 0.01. Every eigenvalue has of order ten candidate partners
inside the window before any structure is consulted.

**The damping is uniform**, one value for all 300 neurons, which is the case the
palindrome proof's condition (a) excludes.

---

## Result 1: the model reaches gamma (withdrawn, and reversed)

The claim was that Wilson-Cowan at standard parameters gives ~12 Hz, in the alpha
band, and that *"to reach 40 Hz requires different parameters"*.

The ~12 Hz is the ringing of the **quiescent** branch. `neural_gamma_cavity.py`
finds its fixed points by Picard iteration, which converges only where the point
is stable *as a map*, a different condition from stability of the ODE. At
I_ext = 2.0 the iteration returns (0.931, 1.000), whose residual is **0.169**
against a sigmoid whose whole range is 1: it is not a fixed point at all.

Integrating the same model with the same parameters, w = 16/12/15/3, α = 1.3,
θ = 4, τ_E = τ_I = 1 ms, nothing changed:

| I_ext | limit-cycle amplitude | frequency | band |
|---|---|---|---|
| 1.0 | n/a | n/a | fixed point |
| **1.2** | 0.266 | **77.0 Hz** | **gamma** |
| 1.5 | 0.440 | 145.5 Hz | above |
| 2.0 | 0.607 | 173.9 Hz | above |
| 2.5 | 0.744 | 159.3 Hz | above |
| **3.0** | 0.847 | **63.6 Hz** | **gamma** |
| 3.5 | n/a | n/a | fixed point |

There is a Hopf bifurcation between I_ext = 1.0 and 1.2, a limit cycle from there
to about 3.0, and two of the sampled operating points sit inside 30-100 Hz. The
frequency mismatch this document reported does not exist at its own parameters,
and the "laser, not cavity" reading built on the overdamped branch goes with it.

Two consequences worth keeping apart. The document was **wrong in the direction
of modesty**: its own model does what it said the model could not do. And the
defect is the same one that sinks Results 2 to 4: a quantity read off a
linearisation at a point that is not the system's state.

---

## Result 2: the 97.3 % measured the tolerance (withdrawn)

The score is a monotone reading of the tolerance against the spectral spread:

| tolerance | loose matcher | strict matcher | self-paired |
|---|---|---|---|
| 0.1 | 100.0 % | 100.0 % | 93.3 % |
| **0.01** (committed) | **97.3 %** | 93.3 % | 41.3 % |
| 0.003 | 80.0 % | 74.0 % | 32.7 % |
| 0.001 | 40.0 % | 36.7 % | 25.0 % |
| 0.0003 | 25.7 % | 25.3 % | 23.0 % |
| 0.0001 | 24.0 % | 24.0 % | 22.0 % |

"Self-paired" counts modes sitting on the centre, which are their own reflection
and which the matcher accepts as pairs; at the committed tolerance that is 124 of
300 modes.

Three further controls, each sufficient on its own:

**The normalisation constant sets the number.** Changing only that constant, same connectome
and same matcher: max|w| (committed, 37) → 97.3 %; spectral radius (23.9) →
93.7 %; max row sum (136) → 99.7 %; 95th-percentile weight (8) → 67.7 %; binary
and unweighted → 27.0 %.

**A structureless random matrix scores higher.** Rescaled about its centre to
the connectome's own Re-spread, which is the only quantity the matcher responds
to, a dense iid Gaussian matrix reaches 100.0 %. (This is a random-matrix
control, not a statement about γ, which in this repository is light.)

**Dale's law makes no difference.** The document's thesis was that the
excitatory/inhibitory classification creates the SWAP structure Π creates in the
qubit chain. Committed Dale (26 inhibitory) 97.3 %; **no Dale at all, every
neuron excitatory, 97.7 %**; 26 inhibitory chosen at random 96.0 %; 50/50 random
signs 97.7 %. The four span 1.7 points, five modes out of 300. The claim to make
is that the signs do *nothing*, which is what kills the thesis; saying the score
rises without them would be a directional reading of that same scatter.

Against degree-matched rewiring (directed double-edge swap, in- and out-degrees
and Dale signs held exactly, R = 200) the strict score sits at p = 0.433, the
loose at p = 0.100. This reproduces at the eigenvalue level what
ALGEBRAIC_PALINDROME_NEURAL found at the residual level.

One quantity is above its null and is reported rather than buried: cross-pairing
at tolerance 1e-4 is 2.0 % against a null mean of 0.7 %, an excess of about four
modes at p = 0.035. Seven one-sided statistics were tested, all chosen after
looking, so that does not clear the 0.007 such a family needs. It is a residue,
not a result.

---

## Result 2b: Q_max = 0.1 as a verdict on biology (withdrawn)

The document read `Q_max = 0.1` against the qubit cavity's 68 to 75 as *"the
biological cavity is extremely lossy"*. Q = |Im λ| / |Re λ|, and the eigenvalues
of `−I/τ + f'·W/max|W|` are exactly `−1/τ + f'·μ`, so Q is arithmetic on the
graph spectrum and two chosen constants. Across τ from 1 to 10 ms and f' from 0.3
to 0.6 it moves by a factor of 51, from 0.086 to 4.343.

The provenance is the sharper point, and it is checkable rather than statistical:
τ = 1 ms is stipulated bare at `neural_gamma_cavity.py:48-49` for the
two-population Wilson-Cowan block, and the C. elegans Jacobian at `:216` contains
**no τ at all**, only a bare `−I`. Reading its unit damping as a millisecond
membrane constant is an interpretation the page never made explicit.
[Q Belongs to No Substance](../docs/Q_BELONGS_TO_NO_SUBSTANCE.md) already
separated this Q from the framework's `Q = J/γ₀`; what it left open, the
provenance of the numerator and the denominator, is answered here: the
denominator is a stipulated damping and the numerator is an uncited standard
parameter set.

---

## Result 3: the anesthesia threshold (withdrawn)

The claim was that oscillations exist only above a critical input I_crit ≈ 0, and
that anesthesia reduces input below it. The reported I_crit = 0.00 is the first
sampled grid point. More decisively, every stability reading in that section is
taken at points the Picard iteration never found, so the threshold is not wrong
so much as undetermined. What the integration above shows in its place is a
*band*: oscillation between roughly I_ext = 1.15 and 3.0, with fixed points on
both sides. Neither the original threshold nor a low-input window survives.

---

## Result 4: the eighteen unpaired modes (withdrawn)

**The count is not a property of the spectrum.** The exclusive matcher scans
indices in order and takes the first best unused partner. Feeding it the same
eigenvalues in a different order gives 30 to 36 unpaired modes, mean 32.7;
LAPACK's ordering gives 20. The committed run reports 18 and the same script on
the same data reports 20 today, a smaller difference than reordering alone
produces. The three biological categories are sorted from a list that reordering
rewrites.

**A second argument, and its limit.** It is tempting to add that the Jacobian is
non-normal, and it is: 64 eigenvalues sit on the centre and the eigenvector
matrix is hopelessly conditioned. But that is manufactured by the defective
cluster *on the centre*, and it does not license a global claim. The
per-eigenvalue condition numbers of the other 236 have a median of 128, which
puts their induced motion at 3.58e-14 against a mean spacing of 1.04e-03. Those
eigenvalues are well resolved. The ordering result, not a conditioning
argument, is what carries this withdrawal.

**The pharynx reading was circular.** The claim of a separate second cavity rests
on zero coupling to the somatic system. There are indeed zero *chemical* synapses
across that boundary. There are **four electrical ones** in the same file: I1L
and I1R to RIPL and RIPR, the RIP-I1 gap junctions and the known route by which
the somatic system modulates pumping. The analysis discards all 1096 electrical
entries and then reports the boundary it created.

---

## Result 5: the wiring is more degenerate at zero than chance (new)

The self-paired modes are eigenvalues of f'·W/max|W| at zero: modes the wiring
does nothing to. Counting them needs a tolerance and an eigensolver. The synapse
matrix is an integer matrix, so the same content is available exactly, and it is
worth being precise about which exact object is the mode count:

| | value |
|---|---|
| nullity of W (geometric, eigenvectors at 0) | 39 |
| **multiplicity of 0 (algebraic, by the rank chain over GF(p))** | **64** |
| float check, \|λ + 1\| < 1e-4 | 64 |

Against 200 degree-matched rewirings the multiplicity is 64 against a null mean
of 48.0, range 43-56, p = 0.005. Three primes agree on the rank, and the count is
identical with and without Dale's law, as row scaling must leave it: the word
"signed" does not belong in this claim.

**And most of it is not spectral.** The structural rank, decided by bipartite
matching on the zero pattern alone with no field, no weights and no signs, gives
a structural nullity of 36 against a null mean of 30.15. So of the 8.8-point
excess in the exact nullity, 5.9 is already in the pattern of which neuron
connects to which. The honest name for most of this is a small maximum matching,
not a numerical degeneracy, and the properly conditioned null, one that holds the
matching fixed, has not been run.

---

## The spiral, and where the repository already keeps one

The Hopf above is an unstable focus: at I_ext = 2.0 the true fixed point is
(0.258, 0.245) with eigenvalues +0.630 ± 2.289i, a spiral the trajectory leaves.
The repository holds a spiral too, typed, in
`compute/RCPsiSquared.Diagnostics/Foundation/ComplexCuspSpiral.cs`: under a
common Z-drift the Bell⁺ coherence becomes complex and winds inward as a
logarithmic spiral, its radius set by the dephasing and its angle by the drift.

These are not the same object, and the differences are worth stating rather than
eliding. That spiral is a trajectory in state space with a shrinking radius; this
one is a linearisation about an equilibrium the flow escapes. That bifurcation is
a saddle-node, the cusp of the cardioid
([CRITICAL_SLOWING_AT_THE_CUSP](CRITICAL_SLOWING_AT_THE_CUSP.md)); this one is a
Hopf.

What they share is the split, and it is the split this document kept getting
wrong. In `ComplexCuspSpiral` the radial decay and the angular winding are
independent: *"the radial magnitude is unchanged and Ω-independent … every spiral
crosses the same ¼-circle at the same time; only the crossing angle is free."*
Decay is Re, turning is Im, and they answer to different parameters. The naming
coincidence this document was built on, the physicist's γ and the
neuroscientist's gamma, puts a Re against an Im. The repository has that
separation as an object; the coincidence does not survive it.

---

## What this leaves

A rebuild would need: a root-finder and a bifurcation diagram rather than an
iteration; τ_I ≠ τ_E, both because biology has it and because the palindrome
proof requires it; the electrical matrix included; a pairing statistic defined
relative to local level spacing with self-pairing excluded by construction; and a
reason why eigenvalues are the right object given how defective this Jacobian is
at its centre. The discriminating pair the qubit side has, two constructions with
the same centre where only one satisfies the proof's condition (b), is the single
thing that would make a neural palindrome claim testable the way the qubit one
is.

The gamma band, meanwhile, is not the open question this document thought it was.
It is reached, at the parameters that were on the page all along.

---

## Reproduction

- Every number above:
  [`celegans_pairing_controls.py`](../simulations/neural/celegans_pairing_controls.py)
  → [`celegans_pairing_controls.txt`](../simulations/results/celegans_pairing_controls.txt)
- The original analysis, kept so the withdrawals can be checked against it:
  [`neural_gamma_cavity.py`](../simulations/neural/neural_gamma_cavity.py),
  [`neural_gamma_cavity_unpaired.py`](../simulations/neural/neural_gamma_cavity_unpaired.py)
- The sibling that had the method right:
  [`hopf_threshold.py`](../simulations/neural/hopf_threshold.py)
