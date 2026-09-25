# The Fragile Bridge: How Much Amplification Before the System Explodes?

*How much amplification can a resonator tolerate before it explodes?*

**Status:** Computed (Tier 2). Four scripts compute the quantum gain-loss bridge, all four at two qubits per chain; `fragile_bridge_bifurcation.py` adds three qubits per chain and `fragile_bridge_n4.py` three and four (Section 4). One of the four certifies the two-qubit threshold as an exceptional point at J_bridge = 1.0 and 1.9 and locates the five couplings where it vanishes; a fifth script takes a first look at a neural node (Section 5).
**Date:** March 29, 2026; the threshold's certificate, zeros and sweep September 2026
**Scripts:**
- [fragile_bridge_bifurcation.py](../simulations/fragile_bridge_bifurcation.py)
- [fragile_bridge_anomaly.py](../simulations/fragile_bridge_anomaly.py)
- [fragile_bridge_ep_signature.py](../simulations/fragile_bridge_ep_signature.py) (the threshold's certificate, zeros and sweep, Sections 2 and 3)
- [fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)
- [fragile_bridge_n4.py](../simulations/fragile_bridge_n4.py) (sparse, N=4)

---

## What this document is about

[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) showed that the palindrome
exists on both sides of zero: decay (positive noise) and amplification
(negative noise, gain). This raises an obvious question: what happens
when you connect a decaying system to an amplifying one?

The answer: it works, but only within limits. Too little coupling and
the two sides cannot interact. Too much coupling and the system
explodes, like a microphone placed too close to a loudspeaker. There
is a sweet spot in between, at roughly twice the internal coupling
strength, where the balance holds and the system is maximally stable. And
at five exact couplings, four below the sweet spot and one above it (two
qubits per chain), there is no balance to hold: any gain, however small,
breaks it.

This document maps where the quantum bridge holds and where it breaks.
Then we asked the same question of the brain's standard two-population
model, a Wilson-Cowan E/I node, whose sigmoid keeps its activities
bounded. What that first scan can and cannot tell us is Section 5.
Before writing it we went through what the repository already holds on
the neural side: the mirror condition [F36/F37](../docs/ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra)
and [its proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md), the
connectome support null in [Neural Gamma Cavity](../experiments/NEURAL_GAMMA_CAVITY.md),
the [quantum proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) the bridge
stands on, and [fw.Confirmations](../simulations/framework/confirmations.py),
which holds no neural hardware result. That side of the story is told in
[the neural account](../docs/neural/README.md) and its
[mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md).

---

## Abstract

Two identical quantum systems, one decaying and one amplifying, are connected by
a single coupling (the "bridge"). The bridge has three stability regimes: weak
coupling stabilizes close to linearly as a trend (γ_crit ≈ 0.189 ×
J_bridge^1.035), the maximum is bracketed near twice the internal strength, and
strong coupling destabilizes as 1/J_bridge, with γ_crit × J_bridge still
drifting downward (0.578 to 0.508 over J_bridge = 10 to 100) where the sweep
ends. For two qubits per chain the threshold falls to exactly zero at five
couplings, J_bridge = 3/4, √5/2, 4/3 and the real roots of two cubics (≈ 1.0294
and ≈ 3.3289): there any gain, however small, makes the bridge unstable (Section
2.1). Away from these five, wherever it was read, the instability is a
second-order exceptional point on the real γ axis (Section 3): certified at
J_bridge = 1.0 and 1.9, and read as one at every coupling of a 4000-point sweep
from 0.005 to 20 except the grid point 3/4, one of the five, all at two qubits
per chain. In the first block of the Liouvillian to go
unstable, two eigenvalues on the imaginary axis meet in a Jordan block at
γ_crit and leave the axis as a mirror pair. The exact Π mirror pairing survives
through the transition. Unlike the fold threshold (~0.5% of J for the product
state, flat in N over the measured N = 2-5, geometric), bridge stability depends
strongly on system size (at J_bridge = 0.1, N=3 is 33.6× less stable than N=2),
making it a topological rather than geometric property.

The Wilson-Cowan comparison is an exploratory scan of a different generator.
Its cross-coupling multiplier is not the quantum bridge/internal-coupling
ratio. A shared optimum, Hopf mechanism or finite window is not established.

---

## 1. Setup

Two identical Heisenberg chains (N qubits each), one with dephasing
(+γ, decay) and one with gain (-γ, amplification), connected by a
single Heisenberg bridge coupling J_bridge between the boundary qubits.

This is the quantum analog of a microphone-loudspeaker system:
decay side = loudspeaker (loses energy), gain side = microphone
(amplifies), bridge = the door between them.

Question: At what gain rate γ_crit does the coupled system become
unstable?

## 2. Results

### 2.1 Linear stability window (weak bridge)

For J_bridge < 2J (bridge weaker than internal coupling):

**γ_crit = 0.1891 × J_bridge^1.0346** (coefficient of determination R² = 0.99984, eight points at
J_bridge ≤ 2 in [fragile_bridge_bifurcation.txt](../simulations/results/fragile_bridge_bifurcation.txt),
AUFGABE B)

Across the regime the critical gain grows roughly in proportion to the
bridge strength, as a trend, with a constant of ~0.19 for N=2 per chain.

The line is clean on that grid and much less clean on a finer one. The 14
points below the peak in [fragile_bridge_anomaly.txt](../simulations/results/fragile_bridge_anomaly.txt)
give a through-origin slope of 0.1845 with R² = 0.6149: same quantity, closer
spacing, and R² falls from 0.9998 to 0.61.

The scatter has a cause, and at its lowest points it is exact. Resolved into its
blocks (Section 3 introduces them), the Liouvillian at zero gain holds
oscillation frequencies that are even or odd under the chain reflection that
swaps the two chains, and the gain-loss, odd under that reflection, couples an
even frequency only to an odd one. Near a coupling where one block holds an even
and an odd frequency that nearly coincide, the threshold is, to first order, the
gain it takes to overcome their mismatch, |ω_even − ω_odd|/(2s), with s the
strength of the gain-loss coupling between the two per unit γ. Where they
coincide exactly there is nothing to overcome: the gain-loss splits the pair at
first order, Re λ = ±s·γ + O(γ²), and any gain, however small, makes the bridge
unstable; the producer finds max Re λ > 0 at every gain it samples at these
couplings, from 10⁻⁶ to 2. For two qubits per chain that happens at five
couplings (section 3 of the [producer's
output](../simulations/results/fragile_bridge_ep_signature.txt)), four of them
inside this regime: J_bridge = 3/4, the real root of 12x³ − 37x² + 72x − 48
(x = J_bridge/J; ≈ 1.0294), √5/2 and 4/3. The fifth, the real root of
4x³ − 21x² + 40x − 48 (≈ 3.3289), is the reversal of Section 2.3. An exact search over every such pair of
levels in every block finds these five and no others at any J_bridge > 0. Around
each one γ_crit rises linearly from zero on both sides, a V whose slope is the
rate at which the mismatch opens divided by 2s, from 0.18 per unit of J_bridge
at 3.3289 to 9.1 at 1.0294. The finer grid's low points sit on these flanks:
recomputed block by block they are 0.0374 at 0.7 and 0.0406 at 0.8 on either
side of 3/4, 0.0381 at 1.1 beside √5/2, and 0.1105 at 1.3 beside 4/3.

So the fitted law is a trend across the regime, not a law each individual
J_bridge obeys. Its eight points (0.01 to 2.0) miss all four zeros inside the
range, and between them γ_crit falls to zero four times.

### 2.2 Optimal bridge near J_bridge ≈ 1.95J

Maximum stability (γ_crit = 0.405849 by the bisection in
[fragile_bridge_anomaly.txt](../simulations/results/fragile_bridge_anomaly.txt);
the collision itself sits at 0.405853, Section 3) is sampled at J_bridge = 1.9, between
neighbours at 0.403820 (J_bridge = 1.8) and 0.383839 (J_bridge = 2.0). The grid
is spaced 0.1 there, so the sweep brackets the maximum and resolves it no
further. What the sweep does say is that J_bridge = 2J is already 5.4% below the
peak, so reading the optimum as "the bridge equals the sum of the internal
couplings" is an interpretation this data does not single out.

Resolved into blocks, the maximum is a kink (sections 2 and 6 of the
[producer's output](../simulations/results/fragile_bridge_ep_signature.txt)).
From J_bridge ≈ 1.4907, where the (0,1) orbit hands over (Section 3), up to
≈ 1.95 the threshold is set by the spin-flip-odd half of the (2,2) block, which
rises across that stretch of the sweep; beyond it the flip-even half takes over
and falls, down to the zero near 3.3289. The two halves cross at
J_bridge = 1.952892, γ = 0.407019, above the sampled 0.4058 and short of 2J,
with every other block higher there; at 2J the threshold is 5.70% below that
crossing value. On the sweep of section 6, 4000 couplings
from 0.005 to 20, no grid value exceeds that crossing: the largest is 0.406954,
at 1.950. Why the two halves cross where they do is open.

### 2.3 Strong bridge destabilizes (J_bridge > 2J)

Above the optimum γ_crit falls, but not by one law the whole way. Between
J_bridge = 3.2 and 4.4 it turns around and rises again (0.0237, 0.0447, 0.1002,
0.1462) before resuming its descent, and that reversal is the fifth zero of
Section 2.1. At the real root of 4x³ − 21x² + 40x − 48, J_bridge ≈ 3.3289, the
threshold of the flip-even half of the (2,2) block touches zero. From there it
leaves zero at 0.18 per unit of J_bridge and keeps rising more slowly, 0.135 per
unit on average, until it meets the threshold of the (1,2) orbit at
J_bridge = 4.4594, γ = 0.1524 (section 3 of the producer's output). On the
sweep of section 6 the (1,2) orbit then goes first from 4.460 to 20, its
threshold falling from 0.1523 to 0.0269. The reversal is what
[fragile_bridge_anomaly.py](../simulations/fragile_bridge_anomaly.py) is named
for, and it is not a 1/J_bridge tail.

Far above the optimum the 1/J_bridge scaling does hold, and the product is
still drifting where the sweep stops:

**γ_crit × J_bridge: 0.578 at J_bridge = 10, falling monotonically to 0.508
at J_bridge = 100** (not converged; a single asymptotic constant is not
measured)

The sampled threshold becomes smaller at the far end of the strong-coupling
sweep. This scan measures γ_crit(J_bridge) only; it does not establish chain
merging, a missing spatial buffer, dimer formation, or resonance as the cause.
The block resolution names both features instead: the turnover is the crossing
of the two (2,2) halves near 1.953 (Section 2.2) and the reversal is the zero
at 3.3289.

The sampled large-J_bridge trend is compatible with γ_crit tending to zero as
J_bridge grows, but the finite sweep does not establish the infinite-coupling
limit or immediate instability at every nonzero gain.

### 2.4 Three regimes summary

| Regime | Condition | γ_crit | Physics |
|--------|-----------|--------|---------|
| Weak sampled range | J_bridge ≤ 2J | trend 0.189 × J_bridge^1.035 | a trend between exact zeros, not a per-point law |
| Exact zeros, across the regimes (two qubits per chain) | J_bridge = 3/4, ≈ 1.0294, √5/2, 4/3, ≈ 3.3289 | 0, rising linearly on both sides | an even and an odd frequency of one block coincide at zero gain |
| Sampled maximum | J_bridge ∈ [1.8, 2.0] | 0.4058; the kink at J_bridge = 1.9529 reaches 0.4070 | a crossing of the two (2,2) halves; its mechanism open |
| Strong sampled tail | J_bridge ≫ 2J | γ_crit·J_bridge = 0.508 at J_bridge = 100 | decreasing threshold; limit and mechanism open |

The 0.1 grid brackets the maximum in [1.8, 2.0] and no finer; resolved into
blocks, the turnover sits at the crossing near J_bridge = 1.953 rather than where
the bridge equals the total internal coupling.

## 3. An exceptional point on the real γ axis, with the mirror intact

The instability is an **oscillating** one: a pair that was already
oscillating at Re = 0 acquires a positive real part, so the system does not
drift away quietly, it screeches like microphone feedback. That much is solid,
and it is why a Hopf bifurcation is the first name that comes to mind.

The symmetry leaves a Hopf no room. At Σγ = 0 the conjugation operator Π
forces exact inversion pairing λ ↔ −λ through the transition. Hermiticity
preservation separately gives λ ↔ λ*. Together the two relations produce a
generic off-axis quartet {λ, λ*, −λ, −λ*}; the across-axis partner −λ* is not
supplied by Π alone. This P-type relation becomes involutive after Π²-parity
resolution; the full SRP sector class remains OPEN.

What holds an eigenvalue on the axis works one block of L at a time. The
Heisenberg Hamiltonian and the Z-dephasing both conserve how many spins point up
in the ket and in the bra, so L splits into 25 blocks labelled by those two
counts (p, q). Write G for the chain reflection that swaps the gain chain and
the loss chain, applied to both sides of ρ. G composed with complex conjugation
(call it GK) maps every block onto itself and sends λ to −λ*. For an eigenvalue
sitting on the imaginary axis that is the eigenvalue itself, so a simple
eigenvalue of a block cannot leave the axis on its own: leaving takes two
eigenvalues meeting, and by Krein's rule two whose signatures under the
reflection are opposite, since levels of equal signature pass through each other
on the axis. Below γ_crit every eigenvalue lies on the imaginary axis; above
γ_crit at least one off-axis quartet appears. What kind of meeting it is decides
what the threshold is. The same split of the popcount blocks by the reflection
organizes [The Reflection Keeps Past and Future
Apart](../experiments/THE_REFLECTION_KEEPS_PAST_AND_FUTURE_APART.md) on the
uniformly dephased chain, where the reflection commutes with all of L; here the
gain-loss is odd under it, which makes G a Krein metric instead.

The blocks come in families with one spectrum. Flipping every spin maps block
(p, q) onto (4 − p, 4 − q) with the same spectrum. So does S = G∘R, G composed
with the one-sided flip R(ρ) = ρ·X^⊗N: it commutes with L exactly at Σγ = 0 and
maps (p, q) onto (p, 4 − q) (section 2 of the output checks both to 0.0). R is
an edge mirror of [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md)'s
mirror group, the R in Π = R·D, and the spin flip is that group's centre, Π².
Hermitian conjugation maps (p, q) onto (q, p) with the conjugate spectrum. The
middle block (2,2) is its own image under all three, and the spin flip splits it
into an even and an odd half. The threshold γ_crit is where, in the first of
these blocks to go unstable, a pair of opposite signature meets. Away from the
five couplings of Section 2.1 the two must first approach, and they meet at
γ* > 0; from here on γ* names that collision point. At those five they coincide
already at zero gain, and γ_crit = 0.

The producer locates the meeting inside the first sector and reads it there, on
two two-qubit chains at J_bridge = 1.0 and 1.9, where that sector is a (2,2) half,
flip-even at 1.0 and flip-odd at 1.9
([`fragile_bridge_ep_signature.py`](../simulations/fragile_bridge_ep_signature.py),
section 1 of [its output](../simulations/results/fragile_bridge_ep_signature.txt)).
The squared gap f = (λ_a − λ_b)² of the colliding pair is real for real γ: minus
the squared distance between the two imaginary eigenvalues below the collision,
(2 Re λ)² above it. So f changes sign at the collision and its zero lies on the
real axis; a root finder on that sign change gives γ*, and there the pair is one
eigenvalue with one eigenvector:

| | J_bridge = 1.0 | J_bridge = 1.9 |
|---|---|---|
| collision γ* | 0.187310108345 | 0.405853184723 |
| eigenvalue at the collision | −2.651477232 i | −6.983542344 i |
| pair gap there | 3.3·10⁻⁸ | 1.0·10⁻⁷ |
| eigenvalues of the full L within 10⁻⁵ of it | 2 | 2 |
| in the sector: smallest singular values of L − λ·I | 6.5·10⁻¹⁶, then 0.176 | 2.7·10⁻¹⁵, then 1.55 |
| half-split / √\|δ\| at δ = ±10⁻⁶ | 0.20293 | 1.1035 |

Two eigenvalues sit within 10⁻⁵ of each other (their gap is of the order of the
square root of machine precision, which is what a located EP2 leaves), yet inside their sector
L − λ·I has a single null direction: one singular value at the rounding level,
the next far above it. Algebraic multiplicity 2 with geometric multiplicity 1 is
a 2×2 Jordan block. The full L agrees, with two eigenvalues there and one
singular value at the rounding level; its next two (both 2.9·10⁻³ at
J_bridge = 1.0) belong to two eigenvalues 3.07·10⁻³ away, one in the (1,2)
block and its flip copy in (3,2). The threshold is a
**second-order exceptional point**, and it sits on the real γ axis at γ* itself.

The approach and the departure are one branch. Below γ* the two are distinct
eigenvalues on the imaginary axis moving toward each other (−2.651274 i and
−2.651680 i at δ = −10⁻⁶, writing δ = γ/γ* − 1); above it they share one
imaginary part and move apart along the real direction (±2.03·10⁻⁴ at
δ = +10⁻⁶). Half the split over √|δ| is the same number on both sides: the
Puiseux law λ − λ_EP = ±c·√(γ/γ* − 1), imaginary below the threshold and real
above it. The coefficient is c = 0.20293 at J_bridge = 1.0 and 1.1035 at 1.9, on
both sides at δ = ±10⁻⁶ and above the collision at every sampled δ from 10⁻⁴
down; at 1.9 it reads 1.1038 at δ = 10⁻³ and 1.1070 at 10⁻², where the next
Puiseux order shows. The meeting is a Krein collision. The chain reflection that
swaps the gain chain and the loss chain commutes with the Hamiltonian part of L
and flips the sign of the gain-loss part, and L is complex symmetric, so iL is
pseudo-Hermitian under that reflection at every γ, and the eigenvalues on the
axis carry Krein signatures. The two approaching eigenvectors carry opposite ones
(at J_bridge = 1.0, ±0.0158 at δ = −10⁻⁴ and ±0.0016 at −10⁻⁶; at 1.9, ±0.0126
and ±0.0013), shrinking to zero at the collision.

At these two couplings the colliding pair lives in a (2,2) half, which the spin
flip and S map to themselves, so the full L carries a single Jordan pair at the
collision and one quartet leaves the axis. Where the first block to go unstable
has partners under the spin flip and S, the full L carries the collision once in
each of them (section 2 of the output): twice for the (0,2) orbit at
J_bridge = 0.1 and for the (1,2) orbit at 0.5, 1.2, 1.4, 5 and 10, and four
times for the (0,1) orbit, which goes first for J_bridge between 1.4344 and
1.4907. On the sweep of section 6 that window is the only stretch with four
copies; the (1,1) orbit, which would have four as well, never goes first there.
At 1.46 the full L holds eight eigenvalues within 10⁻⁵ of λ* and four
singular values at the rounding level, and four quartets leave together. The
(0,1) block holds the coherences between the empty state and one excitation: a
four-site hopping chain with loss 2γ on the decaying chain and gain 2γ on the
amplifying one. Inside the first block the collision is an EP2 at each of the
thirteen sampled couplings (section 2), and inside the symmetry sector of it that
goes first at every coupling of the sweep but the grid point 3/4 (section 6).

A Hopf bifurcation moves a simple pair across Re = 0 with a finite slope. Here
GK forbids a simple eigenvalue of any block from leaving the axis alone, at
every coupling, so the threshold is never a Hopf in that sense. Where the pair
that leaves is a Jordan pair, the departure has infinite slope and the
threshold is an exceptional point. In a family with one parameter, two levels
of opposite signature that meet generically form a Jordan pair: a semisimple
meeting, or a third level joining, takes a further condition, generically met
only at isolated couplings. Five such couplings are known, the zeros of Section
2.1, all at zero gain: there the pair is degenerate but semisimple and leaves
with the finite slope s, not a Hopf either, since the pair is not simple, and
not an exceptional point. Section 6 of the output reads the collision in the
symmetry sector that goes first at every coupling of its sweep, 4000 couplings
from 0.005 to 20, and finds an EP2 at all of them except the grid point 3/4,
one of the five zeros: inside that sector one singular value of L − λ*·I at
the rounding level (at most 1.8·10⁻¹⁴), the next at least 7.9·10⁻⁴, the third
eigenvalue at least 1.2·10⁻² away, and one Puiseux branch. (In Hamiltonian
dynamics a Krein collision that releases a quartet is called a Hamiltonian-Hopf
bifurcation; the Hopf ruled out here is the dissipative one this section opened
with.)

The finite-offset readings above γ_crit fit that picture (section 4 of the
output, with δ measured from the bisected γ_crit = 0.1873101 and 0.4058524):

| δ | max Re λ (J_b=1.0) | Re/√δ | gap to across-axis partner −λ* | Petermann K |
|---|---|---|---|---|
| 10⁻² | 0.02033551 | 0.2034 | 4.07·10⁻² | 40.9 |
| 10⁻³ | 0.00641845 | 0.2030 | 1.28·10⁻² | 403.2 |
| 10⁻⁴ | 0.00202899 | 0.2029 | 4.06·10⁻³ | 4027 |
| 10⁻⁵ | 0.00064059 | 0.2026 | 1.28·10⁻³ | 4.04·10⁴ |

- **Re λ/√δ is nearly constant** over the three decades (four samples) of this
  table, where a transversal crossing would move it by a factor of 32. At
  J_b = 1.9 the same column reads 1.107, 1.103, 1.093, 0.991; its last two rows
  carry the bisection offset described next.
- **The across-axis partner gap is 2|Re λ|.** Its decrease is algebraically
  equivalent to the max-Re onset: the derived partner gap adds no independent
  evidence.
- **The Petermann readings grow one decade per decade**: 40.9, 403, 4027,
  4.04·10⁴, the K ∝ 1/δ of an EP2. These readings suggest an exceptional point;
  the Jordan rank at the collision establishes it.

The bisected γ_crit sits a hair below the collision: 3.5·10⁻⁸ below γ* (relative)
at J_b = 1.0 and 1.9·10⁻⁶ at 1.9. A bisection settles where its own sequence of
max Re λ readings crosses 10⁻¹², and near the collision the two approaching
eigenvalues are ill-conditioned, so rounding noise in their real parts crosses
that level while they are still split on the axis (by 7.6·10⁻⁵ and 3.1·10⁻³ at
the bisected values). At J_b = 1.9 the offset is 1.9 % of δ at 10⁻⁴ and a fifth
of δ at 10⁻⁵, which is what bends the last two rows there; measured from γ* the
coefficient is 1.1035 for δ ≤ 10⁻⁴.

Below threshold the axis is clean to machine precision: max Re λ is
2.4·10⁻¹⁴ at δ = −10⁻³ and 1.6·10⁻¹⁴ at δ = −10⁻² at J_bridge = 1.0, and
4.5·10⁻¹⁴ and 2.0·10⁻¹⁴ at 1.9.

So the geometry is Hamiltonian PT breaking rotated 90°: there two real
eigenvalues meet in a Jordan block and become a complex pair, here, generically,
two imaginary eigenvalues meet in a Jordan block and become a mirror pair
across the imaginary axis. Π itself is linear rather than anti-linear, so
it supplies a P-type relation rather than PT; after Π²-parity resolution it is
involutive, and the full SRP sector class is OPEN. The anti-linear symmetry at
work inside every block is a different one: GK, under which iL is PT-symmetric
in the textbook sense (the Heisenberg Hamiltonian is real). Below γ_crit it is
unbroken; the threshold is where it breaks. A second anti-linear symmetry sends
λ to −λ* too: Π composed with Hermitian conjugation, Π's antiunitary twin in
F119's antilinear double. With Π = R·D and D∘† = conj it is R∘conj, the
one-sided flip with entrywise conjugation (section 2 of the output checks it to
0.0). It maps block (p, q) onto (p, 4 − q), so it acts inside a block only when
q = 2, the (2,2) block of both certified collisions among them, and its product
with GK is S = G∘R.

Each piece of that symmetry is one the repo already owns. The identity G·L(γ)·G
= L(−γ) is the conjugation identity
[F131](../docs/proofs/PROOF_MIRROR_ORDER_SORTING.md) starts from (F131 and F71
write the site reversal R; this page keeps R for F118's one-sided flip): its
unitary column (§2), stated there for the XX chain, which needs only a
reflection-even Hamiltonian, and the bridged Heisenberg chain with bonds (1,
J_bridge, 1) is one. The gain-loss profile (+γ, +γ, −γ, −γ) sums to zero over
every mirror pair of sites, so it lies in the anti-palindromic class of
[F91](../docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md) about its zero mean, a
fixed point of F91's parameter reshuffle R₉₀. The complex-conjugation half is the transport law of
[F119](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md): conjugation flips the sign
of a real Hamiltonian's part of L and fixes the dephasing. The γ ↦ −γ half is
the gain turn of MirrorWorld's `GammaFold`, which at Σγ = 0 coincides with its
anti-watch turn. The one-sided flip R inside S is one of F118's edge mirrors: it
fixes the Hamiltonian part and reflects the dissipator, R·L_diss·R = −L_diss −
2Σγ·I, so that R·L(γ)·R = L(−γ) at Σγ = 0; MirrorWorld's `Lattice` runs the same
map as its bra-side reading. S is F131's reflection composed with it. And since
the reflection anticommutes with the gain-loss part while it commutes with the
Hamiltonian part, the gain-loss couples only reflection-even to reflection-odd
states. Krein's rule makes the opposite signatures of a meeting pair necessary;
this coupling is why such a pair attracts as the gain rises, and where the two
already share a frequency at zero gain it splits them at once: the zeros of
Section 2.1.
See [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md).

The oscillation frequency at threshold, |Im λ*|, follows the first block to go
unstable and is not monotone in the bridge (section 2 of the output):

| J_bridge | first block | \|Im λ*\| |
|----------|-------------|-----------|
| 1.46 | (0,1) orbit | 0.722 |
| 1.5 | (2,2) flip-odd | 5.956 |
| 1.9 (sampled optimum) | (2,2) flip-odd | 6.984 |
| 2.0 | (2,2) flip-even | 7.214 |
| 5.0 | (1,2) orbit | 0.260 |
| 10.0 | (1,2) orbit | 0.115 |

It is 5.956 to 7.214 where a (2,2) half goes first (1.5 to 2.0 in the table)
and below 1 where the (0,1) or (1,2) orbit does (0.722 at 1.46; 0.260 and 0.115
at 5 and 10). At 5 and 10 the first block is the (1,2) orbit, so there the
collision happens in two blocks at once. Whether the threshold frequency tends
to zero as the bridge grows is open.

### 3.1 Relation to the local-EP question

This file's own exceptional point is certified at J_bridge = 1.0 and 1.9 and
read at every coupling of the section-6 sweep but the grid point 3/4, all with
two qubits per chain (Section 3): it sits on the real γ axis at the collision γ*
itself, in whichever block of L goes unstable first. At the five zeros of
Section 2.1 there is none. The Petermann
readings follow from it. K diverges as 1/δ at the EP, so a scan that samples at
δ ≈ 10⁻³ reads about 400 and one at 10⁻⁴ reads about 4000: the K = 403.2 at
γ/γ_crit = 1.001 is K at the scan's own distance from the threshold, a
finite-offset reading of a divergence, not another locus and not a property of
the threshold. Coarse-scan rows whose selected eigenvalue is degenerate omit
single-vector K because it is basis-dependent within that eigenspace.

What is shared is the algebra read at two residuals of the F1 palindrome
`Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ vs Σγ = 0), and three different
objects carry exceptional points. F86's toy 2×2 rate-channel reduction has its
EP at Q_EP = 2/g_eff, a double root at the positive decay rate 4γ₀·k (4γ₀ for
the slowest pair, k = 1). This file's
SEPARATE Σγ = 0 gain-loss system has the EP2 of Section 3 (two qubits per
chain), a double root at zero decay on the imaginary axis. F89 separately certifies narrow real-q defective
EP2s of the full Σγ = N·γ₀ block at N=5,7,9; its all-odd endpoint-nullity
surplus does not extend that character verdict. Whether the full Σγ = N·γ₀
block shares a defective-EP structure off the real axis is open (the nearest
complex-Q coalescences found 2026-06-21 are themselves diabolic, ‖P‖ = 1).
Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs`
(**OpenQuestion**; the four PetermannSpikeWitness rows are kept only as a
cautionary non-normality record, not EP evidence). The claim carries this system's collision point, and
`FragileBridgeThresholdTests` recomputes it from the C# builders.

## 4. N-dependence: topological, not geometric

How does the stability limit change with system size? At J_bridge = 0.1,
the one coupling computed at three qubits per chain, it depends strongly on
how many qubits are in each chain: three qubits per chain are 33.6× less
stable than two. Whether even-length chains are more stable than odd-length
ones is a question the N=4 point below cannot yet answer.

γ_crit depends strongly on chain length at the one coupling computed:

| N (per chain) | γ_crit (J_bridge=0.10) | Ratio to N=2 | instrument |
|---------------|----------------------|-------------|------------|
| 2 | 0.017292 | 1.000 | dense eigenvalues |
| 3 | 0.000515 | 0.030 | dense eigenvalues |
| 4 | 0.001186 | 0.069 | norm-fit estimator, error uncharacterised |

The first two rows are exact eigenvalues of the coupled generator: at 4 and 6
qubits the Liouvillian is 256² and 4096², and `max_re_sparse` solves it densely.
N=3 is 33.6× less stable than N=2 on those numbers. The
[bifurcation producer](../simulations/fragile_bridge_bifurcation.py) reaches
0.000500 for the same point on a tighter bisection, i.e. 34.6×; the two agree to
3%, which is the resolution of the claim.

The third row is a different instrument and carries a different weight. At 8
qubits the generator is 65536², beyond a dense solve here, so γ_crit(N=4) is
bisected on a norm fit: propagate a random vector, take the slope of log‖v‖ over
five samples. Where that fit CAN be checked against exact eigenvalues, at 4 and 6
qubits, it under-reports max Re(λ) by 8.8% to 100%, always in the same direction,
which biases a bisection's γ_crit high. At 8 qubits there is no reference, so the
error is not measured and 0.001186 is a reading of the estimator rather than of
the spectrum.

**So the even/odd reading is not supported by this table.** "N=4 is 2.3× more
stable than N=3" compares a dense number with an estimator number whose one-sided
bias points the same way as the claimed effect. A parity law in the chain length
would be an attractive thing to have, and it needs the N=4 point measured by an
instrument with a known error before it can be asserted; a sector-resolved solve
or a Krylov eigensolver with a residual bound would supply one. No simple power
law or exponential fits the three points either.

This is fundamentally different from the fold threshold
(Σγ_crit/J ≈ 0.5% for the product state, flat in N over the measured N = 2-5). The fold is a **geometric**
property of the palindrome. The bridge stability is a **topological**
property: it depends on chain length, at the one coupling where chain lengths
were compared. Whether it also depends on the number of bridge connections is
open question 2; one bridge is all that has been computed.

| Property | Fold threshold | Bridge stability |
|----------|---------------|-----------------|
| N-dependence | Independent | Strongly dependent at J_bridge = 0.1 |
| Type | Geometric constant | Topological |
| Bifurcation | Fold (saddle-node) | two qubits per chain: generically an exceptional point on the real γ axis (an EP2 in the first block to go unstable, oscillating); a linear onset from zero gain at five exact couplings |
| What it measures | When irreversibility begins | When coupled gain-loss explodes |
| Determined by | Palindrome geometry | Chain length and bridge coupling (one bridge computed) |

## 5. Neural comparison: a first scan, and what it would need

Does the same stability problem appear in the brain? That was our next
question on March 29, and we went straight to the Wilson-Cowan model, a
standard model of how a population of excitatory neurons (which amplify,
like the gain side) meets a population of inhibitory neurons (which damp,
like the decay side) through synaptic connections. It looked like the
biological twin of the quantum setup, so we asked whether the three
properties of the bridge would show up there too.

Script: [fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)

### 5.1 Test design

A single Wilson-Cowan E-I node: two activities, a 2×2 Jacobian (the matrix
of partial derivatives that decides local stability), τ_E = 8 and τ_I = 18.
The internal couplings w_EE = 16 and w_II = 3 stay fixed; the two
cross-couplings are scaled together by a factor s, as w_EI = 12s and
w_IE = 15s. P is the external input, and P_crit(s) is the input the bisection
returns between a stable P = 0 and an unstable P = 10.

One difference from the quantum bridge is built into this design. The
multiplier s scales two cross-weights with different fixed bases, against
internal weights of 16 and 3; it is not the ratio J_bridge/J, and s = 2
does not mean "twice the internal coupling" in the quantum sense.

### 5.2 What the scan printed

| s | P_crit | Im λ at 1.01·P_crit | Freq (Hz) |
|-----|--------|-----------|-----------|
| 0.1–1.5 | none in 0 ≤ P ≤ 10 | - | - |
| 2.0 | 2.28 | 0.229 | 36.5 |
| 2.3 | 4.78 | 0.419 | 66.8 |
| 2.5 | 8.33 | 0.593 | 94.4 |
| 3.0–10.0 | none in 0 ≤ P ≤ 10 | - | - |

([fragile_bridge_neural.txt](../simulations/results/fragile_bridge_neural.txt),
with the fine sweep s = 2.0–2.5 in steps of 0.1; the imaginary part is
read just above the threshold, at 1.01·P_crit.) Where the scan found an
instability, the leading eigenvalue pair was complex. Read at face value the
table shows a threshold that climbs from s = 2.0 to 2.5 and is gone by
3.0, its last value 8.33 already near the scan's P = 10 ceiling, with an
oscillatory onset like an EEG rhythm. It invites the three things the
quantum bridge had shown: a special coupling near 2×, an oscillating
instability, and a finite window.

The instrument cannot carry that reading yet, and the reasons are in the
code. `find_fixed_point` runs a damped iteration for a fixed number of steps
and returns wherever it stands, without checking the equation residual, so
a returned point need not be an equilibrium at all. `find_P_crit` tests the
two endpoints P = 0 and P = 10 and bisects only when they differ; when both
are stable it reports "none", which cannot exclude an unstable interval
between them. Nothing continues an equilibrium branch, and a complex
eigenvalue at a nearby sampled point supplies no crossing or nondegeneracy
test. The same trap caught a larger network in the
[mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md#a-bounded-iteration-cannot-locate-a-hopf-bifurcation),
where fresh residuals at the returned points came out between 0.5 and 1.

So the three properties do not transfer as measured. The "2×" coincides in
name only, since s is not the quantum ratio and the quantum optimum sits
at J_bridge ≈ 1.953 with 2J already 5.70% below it (Section 2.2). The oscillatory
onset is what the scan saw, but a Hopf threshold is a verdict this scan
cannot give. The finite window may be real, or may be the two endpoints
agreeing, with the last threshold sitting close to the P = 10 ceiling the
scan never looks past.

### 5.3 The sigmoid bounds the model

One property of the model holds by a short argument. For either activity x
the equation is τ dx/dt = −x + S(input), with τ > 0 and 0 < S < 1. At x = 0
the derivative points inward, and at x = 1 it points inward too, so the
continuous model keeps both activities inside [0,1]² once they start there.
The density matrix of the quantum bridge has no such wall and can run away
under gain. It is tempting to call the sigmoid a biological safety mechanism, the
thing that keeps the neural equivalent of a laser from exploding. What the
argument gives is smaller: boundedness permits fixed points, transients and
oscillations alike, and it proves neither that an equilibrium is stable nor
that anything is protected against a pathology.

### 5.4 What a neural threshold would take

A defensible threshold needs a converged equilibrium branch over the
declared ranges of P and s, a fresh residual at every point, a nonzero
imaginary eigenvalue pair tracked through Re λ = 0, a transverse crossing
and the nonlinear nondegeneracy checks. The oscillation and its onset then
need timestep and duration checks in the time domain. These are the next
gates, and none of them has run on this probe yet.

A comparison through the neural mirror condition F36 asks for more again: a
specified involution and a scalar centre satisfying both the diagonal and
the effective-coupling conditions at the operating point. The
[neural translation gate](../simulations/neural/neural_translation_gate.py)
holds exact palindromes with complex spectra and constructed unstable
instances side by side; pairing by itself gives no silence and no
stability.

## 6. Open questions

1. **N-scaling law (open, and the N=4 point needs a better instrument first):**
   N=4 was computed at 65536×65536 through `expm_multiply`, a SciPy routine that
   applies the matrix exponential to a vector without forming the operator, with
   γ_crit bisected on the growth rate that fit returns. Measured against exact
   eigenvalues at the two smaller sizes, that fit under-reports max Re(λ) by 8.8%
   to 100%, one-sided, which pushes a bisected γ_crit up. The suspected even/odd
   parity effect is exactly what such a bias would manufacture, so N=5 is not the
   next step: re-measuring N=4 with a bounded-error instrument is. A Krylov
   eigensolver carrying a residual bound, or a solve restricted to the sector the
   leading mode lives in, would give one.

2. **Multiple bridges:** What if the two chains are connected by
   more than one qubit pair? Does γ_crit recover N-independence
   when bridges scale with N? Real neural networks have distributed
   rather than point-to-point E-I coupling, which makes the question
   tempting on that side too; a network there would need its own
   generator and its own stability analysis rather than this curve.

3. **Cascade stability:** If each level of a frequency hierarchy
   were a coupled gain-loss pair, each would have its
   own bridge window. The weak-regime trend γ_crit ≈ 0.19 × J_bridge
   would suggest one condition per level, but the trend does not hold per
   coupling: a level whose bridge sits at one of the exact zeros of
   Section 2.1 is unstable at every gain sampled there, from 10⁻⁶ to 2. Whether
   windows compose is a question for the combined generator: build the
   chain of bridges and compare its spectrum with the isolated-bridge
   prediction.

4. **Large-bridge asymptotics:** Does γ_crit × J_bridge converge as
   J_bridge grows? The executed values decrease from 0.578 at J_bridge=10 to
   0.508 at J_bridge=100. Section 7 of the producer's output reads it at six
   couplings, with the (1,2) orbit first at each: 0.5384, 0.5254, 0.5151, 0.5075,
   0.5025 and 0.5008 at J_bridge = 20, 30, 50, 100, 300 and 1000. The samples
   still decrease and establish no limit. If a limit exists, is it 1/2? That
   value is a hypothesis to derive or falsify, not a measured constant.

5. **Saturation as design principle:** The sigmoid keeps Wilson-Cowan
   activities inside [0,1]². Is there a quantum analog? A state-dependent
   γ could act as a quantum sigmoid and might bend the bell curve into a
   window. To find out, write that gain model down with its state domain
   and run a stability test on it; the bounded neural model predicts
   neither the quantum curve's shape nor a biological safety mechanism.

---

*Computed March 29, 2026. Three scripts on the quantum side map the
bridge: a rising weak regime, an optimum near 1.9J, and a falling tail that
is not one law the whole way. A fourth, in September 2026, certifies the
threshold as an exceptional point at J_bridge = 1.0 and 1.9, finds the five couplings where it
is zero and resolves the optimum into a kink near 1.95J, all at two qubits per
chain. The neural
scan asked whether the brain meets
the same limit; it found an oscillatory window whose equilibria it never
checked, so that answer waits for the gates in Section 5.4.*
