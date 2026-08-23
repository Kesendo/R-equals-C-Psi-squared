# The Blind Site: counting what dephasing at one seat cannot reach

**Date:** August 23, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [`simulations/blind_site.py`](../simulations/blind_site.py), which imports its propagator from [`simulations/bridge_sector.py`](../simulations/bridge_sector.py)
**Data:** [`simulations/results/blind_site/blind_site_run.txt`](../simulations/results/blind_site/blind_site_run.txt)

Z-dephasing at a single site of an open Heisenberg chain has a subspace of
**single-excitation** states it cannot touch: not weakly, not slowly, exactly.
Which seat you choose decides whether such states exist and how many, and the
count is a divisor question, in the ordinary arithmetic sense and not F140's
frozen `Divisor`, which §5 keeps separate,

  **dim (blind single-excitation subspace at site j) = (gcd(2j+1, N) − 1) / 2,**

sites indexed j = 0..N−1, zero unless 2j+1 and N share a common factor above 1.
On the eleven-site mediator bridge the mediator is site 5, gcd(11, 11) = 11, and
the subspace has dimension 5. Dephasing there at γ = 50, fifty times the bond
coupling J = 1, moves such a state's trajectory by 6.2e-16; the same rate one
seat over moves it by 0.40.

**The mechanism is not new and the count is.** That a mode with a node at the
dephasing site has rate zero is **F64**, exact at every γ_B and stated for
"Z-dephasing on any single site B", with `review/EMERGING_QUESTIONS.md:542`
spelling it out: "F64 captures protection, not just dissipation:
rate = 0 ⟺ |v(B)|² = 0 (mode has a node at the dephasing site)." The eigenbasis
is `PROOF_UNIFORM_LAW.md`'s, gated there entry-exactly. The value ⌊N/2⌋ at the
centre of an odd chain is `MEDIATOR_NOISE_GATE_LEVEL_THREE.md`'s, five days old,
in this same arc. Even the word is borrowed twice over. What is this page's own
is the **count as a function of the site**, its even-N and composite-N half, the
intersection law for several dephased seats, and the measured mediator table. §10
is the ledger.

---

## What the repo already held, before any of this was run

The sweep, store by store. It moved the result three times and every time
downward, which is why this section is long.

**`docs/ANALYTICAL_FORMULAS.md`.** **F64** (line 1852) is the nearest prior
result and owns the mechanism: −Re(λ_k) = 2γ_B·|v_k(B)|² with v_k the
coherence-sector Liouvillian eigenvector, "mode by mode, max relative error
6.2·10⁻¹³", scoped (line 1869) to "any graph topology … XX+YY or Heisenberg
single-excitation, Z-dephasing on any single site B". Its closure at
`review/EMERGING_QUESTIONS.md:542` names the zero case as protection. **One
thing, and only one, separates this page from F64: F64 says when a rate is zero,
and this page says how many such modes a given seat has.** F64 lives in the
(0,1) coherence block while §§1-4 here live in the (1,1) density block, but that
difference is no defence and is not offered as one: §5 shows the same count
holds in F64's own block, by a two-line bijection off F64's own identity.
`EMERGING_QUESTIONS` reaches its protected modes through **degenerate**
topologies (ring, star, Y, K₅) where "B-decoupled rotations" can be built inside
a degenerate subspace; the open chain's spectrum is nondegenerate, so the modes
counted here come from the mode profile itself and not from any freedom of
basis.

**F66 (lines 1906-1922) owns the genre**, though only at one seat: it counts
what a single-site Z-dephasing channel cannot touch, "completely shielded from
Z-dephasing at B", multiplicity N+1, **verified for B at the endpoint only**.
Its scope paragraph leaves open exactly this direction: "Whether the same
structure … persists … for interior B-positions is open", noting the α = 0
multiplicity at the centre of an N = 5 chain is 64 rather than 6. **This page
does not answer that open item.** F66 is the XY chain counting full-Liouville
modes in the extreme XY-weight sectors, which the registry says the
single-excitation sector never reaches for N ≥ 3; this is the Heisenberg chain
counting states in one N-dimensional sector. It is a neighbouring question, to
be read next to F66 and not as closing it.

F126 (header line 5182, the Green's function at line 5193) solves this sector
exactly by a renewal representation. So this page does not solve the sector; it
asks which states the solution never moves. F2 (line 111) names the mechanism
that fixes the eigenbasis, "the ZZ term in the Heisenberg Hamiltonian supplies
the diagonal shift that turns the adjacency matrix into that Laplacian, and it
is absent in the XY case". F152 (line 7058) gives the (0,1) block generator used
in §5. F11 and F65 give mode profiles and rates but not a per-seat count.

**`docs/proofs/`.** `PROOF_UNIFORM_LAW.md` B0 (lines 182-186) owns the starting
point, in this page's own 0-based indexing, and supplies both facts §5 needs:
"The sector-1 restriction of the chain H is (N−1)·J·Id − 2J·L with L the
path-graph Laplacian (the gate compares the matrices entry-exactly). Its
eigenvectors are the Neumann cosine modes of D10, ψ_m(j) ∝ cos(πm(2j+1)/(2N))
for m = 0..N−1, with distinct eigenvalues … so the sector-1 spectrum is
nondegenerate." The first run here re-derived that identity independently; it is
**not** a finding of this page and now stands in the script only as a visible
check, where the residual is exactly 0.00e+00 at every N tested.
`PROOF_R90_FROZEN_DIVISOR.md` Lemma 5 (line 176, Heisenberg bullet at line 178)
carries the same basis with eigenvalues written out, λ_k = 4cos(kπ/N) + N − 5.
D10 (line 180) states the contrast that makes the count a Neumann count: "The
Dirichlet family sin(πkj/N) is **NOT** the eigenbasis here."

`PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md` proves that the decoherence-free
subalgebra of the dephasing dissipator is the diagonal operators, a
d_w-dimensional algebra. §6 both uses it and reports an internal inconsistency
in it that this page's own state exhibits.

**`hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md` has the word and the observation.**
Line 133: "when the dephased site sits on a node of a Hamiltonian eigenmode,
that mode is **blind to dephasing**, its partners degenerate together, and the
fine structure collapses", from an N = 5 B-position scan; line 144 adds "for
interior B the slowest S-coherence mode can be exactly dark, |a_B|² = 0". So the
title word is borrowed, and so is the qualitative interior-B observation. What
that scan reports is the collapse of the α-value count, and those numbers must
be quoted in their corrected form: the familiar 57 at an endpoint, 7 at
positions 1 and 3, 8 at the centre are counts of distinct *first-order* values
binned at tolerance γ₀·10⁻³, and the same line's 2026-08-09 correction gives
110, 25 and 12 distinct eigenvalues at finite γ₀ = 0.1. Either way it counts
α-values, not a subspace dimension per seat.

**`experiments/`.** `MEDIATOR_NOISE_GATE_LEVEL_THREE.md` (lines 257-260, commit
`751ec96`, five days old and in this arc) already states the centre case: "every
antisymmetric one vanishes at the centre site exactly … ⌊N/2⌋ of the N modes
have an exact node at the centre, 5 of 11 here and 2 of 5 at N=5", and its own
per-site table at line 267 carries a column headed "exact nodes there", gated by
the committed `simulations/bridge_node_weights.py`. **So the N = 11 value 5 is
that sentence, not this page's.**

`ORTHOGONALITY_SELECTION_FAMILY.md:37` defines "orthogonal complement H_M^⊥ = the
**blind subspace**", and line 279 gives the Meta-Theorem as "conservation law +
summed measurement → blind **channel**", so this page's title is borrowed and
its result is an instance, with conservation = H-invariance and
measurement = n_k. Its line 349 counts "modes with ⟨n_XY⟩ = 0 (pure lens, fully
γ-immune, decoherence-free subspace) … exactly N+1 at every N"; those are
Liouvillian **modes**, not Pauli strings, and one sector apiece, a different
object from a subspace inside one sector. `J_BLIND_RECEIVER_CLASSES.md` (lines
32-40) already carries the two-condition logic this page re-derives, "(i)
L_D[ρ_0] = 0 … (ii) ρ_0 is a simultaneous eigenstate of every bond Hamiltonian",
with the explicit warning that "(i) alone is NOT sufficient".

`GAMMA_CONTROL.md:31` points the other way for its own observable and is
confronted in §7. `SYMMETRY_CENSUS.md:101` reports one steady state per sector
and "no … dark states", reconciled in §6. `RESONANT_RETURN.md` falsifies
spatially symmetric **γ profiles**, a different lever from the symmetric
**state** prepared here: its 43 % (line 578) is N = 7, and the N = 9 figure of
47 % is at line 582.

`MATHEMATICAL_FINDINGS.md:400-402` states "noise distribution matters when the
Hamiltonian breaks the symmetry between subsystems, not when the noise does."
§3 and §4 look like a counterexample and are not one: that sentence closes a
bipartition question, C_int against C_ext, and speaks of symmetry between
**subsystems**. Here the Hamiltonian is symmetric, where γ is placed decides
everything, and what is broken is the chain's own reflection rather than a split
between two parties. The two do not overlap, but the resemblance is close enough
to say so out loud.

**The OpenArcs registry.** `the_gate_that_does_not_gate` (line 7896) owns the
question and records at lines 7963-7970 the site-by-site response profile at
N = 11 and the conclusion that the mediator's node structure "is measurably
inert". §7 says how that stands. `site_resolved_vacuum_block` (lines 4842-4850)
settles the labels below and warns of "the trap inside the trap".

**`docs/GLOSSARY.md`: nothing.** One occurrence of "standing wave", line 190,
under **Resonator**, whose boundaries are CΨ values rather than positions. Zero
for node, antinode, group velocity, wavepacket.

**`fw.Confirmations`: nothing.** The nearest entry, `ibm_ep_onset_may2026`,
measures per-site populations of a single-excitation walk on three sites, with
no dephasing-support variation.

**`docs/CAUGHT_ERRORS.md`: four shapes, all of which bit.** Line 619, a linear
law shipped four months under a γ² heading because a script divided by an
assumed power and eyeballed the residual, answered here by exact ranks rather
than fitted ones. Line 877, "a profile comparison is a measurement of Σγ unless
Σγ is matched by construction", which caught a defect in this page's own second
draft: see §4. Line 424, integer counts belong to the bare coherence and not to
dressed eigenmodes at J > 0, which also caught this page calling
`ORTHOGONALITY_SELECTION_FAMILY`'s dressed modes "Pauli strings" two paragraphs
after citing the entry. And line 758, this page's risk written down in advance:
F143 recorded as "the narrower XY re-derivation of a lemma the repo already
held, minted five days afterwards without citing it, which is exactly the prey
class the hunt was looking for". The lemma re-derived here is the same one and
the five days are the same five days.

**Which objects, because there are three and the arc insists on the split.** The
**propagation** of §§1-4 runs on ρ in the **(1,1) joint-popcount block**, the
N × N Haken-Strobl density block. The **counting** of §5 runs on **h_SE, the
N-dimensional single-excitation Hamiltonian**, and the blind dimension is a
dimension of that Hilbert space: a count of **states**, not of operators. The
blind operator space inside the (1,1) block is larger and is not computed here.
The third is the **(0,1) vacuum-coherence block** of F152 and F64, where §5
shows the same count holds. None of the three is the Pauli w = 1 sector, which
is not L-invariant.

Two more words below are the repo's and are used in the ordinary sense, so they
are pinned here once. **Mirror** always means the site reversal j ↦ N−1−j,
object 1 of the ⌊N/2⌋ map, and never Π or `MirrorWorld`'s `Mirror`, the
block-lattice group of eight. **Block** does mean `MirrorWorld`'s `Block`, the
joint-popcount block (p, q), and is used for nothing else.

The Hamiltonian is **Heisenberg**, Σ_bonds J(XX + YY + ZZ), carrying the ZZ
degree diagonal, so h_SE = (N−1)J·Id − 2J·𝓛 and the modes are Neumann at
modulus N. It is *not* the XY h_SE of F2b at modulus N+1, which is what
`compute/MirrorWorld/Cone.cs` builds (`hop[a,b] = j`, nothing ever written to
`hop[a,a]`) and what `simulations/cone_defect_arrival.py` builds; the latter
also damps with a single uniform scalar and so cannot express single-site
dephasing at all.

**The lever this page pulls is the assignment, not the budget.** Every arm
below holds one γ and moves it from seat to seat, so Σγ is matched by
construction and the comparison is of a spatial pattern rather than of a total.
That is the distinction `docs/CAUGHT_ERRORS.md:877` was written about, and it is
also why "total noise" is the wrong register for this page: a scalar cannot say
which seat.

**The γ book depends on the support.** Two single-excitation configurations
differ in exactly two bits, so a coherence between them decays at 2(γ_a + γ_b).
With γ on every site that is 4γ, which is F126's Γ = 4γ. With γ on one site f,
this page's case, the largest rate present is **2γ_f and never 4γ_f**, so
F126's book does not describe §§3-4 even though the engine convention is
identical. That convention is `bridge_sector.py`'s
mask = −2·Σ_{differing bits} γ_q, matching γ_k(Z_k ρ Z_k − ρ).

---

## 1. The preparation that had not been run

Two earlier tests, both headed "STANDING WAVE ACROSS THE BRIDGE", looked for
counter-propagating structure without ever preparing two equal counter-running
waves. `RunTest3` (`Program.cs:897`, N = 11) starts from `BellPair(0,1,n)` with
nothing on sites 2 to 10, genuinely one-sided. `run_test_4`
(`mediator_bridge.py:642`, N = 5) is **not** one-sided: `bell_A_0M_pp_B` is
Bell(0,1) ⊗ |0⟩₂ ⊗ |+⟩₃ ⊗ |+⟩₄, so sites 3 and 4 carry weight from t = 0. It is
asymmetric and unentangled on the right and spreads across every popcount, which
makes it a different kind of unfit rather than the same one. Neither prepares a
mirror-symmetric pair, and `make_initial_state('bell_A_0M_bell_B')`, which
would, is defined in the same file at line 233 and referenced exactly once, in
Test 2.

Four preparations run here at N = 11, uniform chain, J = 1, t ≤ 20, RK4 at
dt = 0.01:

| | state |
|---|---|
| ONE | \|0⟩, the earlier control |
| SYM | (\|0⟩ + \|10⟩)/√2 |
| ANTI | (\|0⟩ − \|10⟩)/√2 |
| MIX | ½(\|0⟩⟨0\| + \|10⟩⟨10\|) |

MIX carries the same two wavepackets, the same energy and the same Σγ as SYM
and ANTI with no phase relation between them, so it cannot interfere. By
linearity MIX = ½(SYM + ANTI) exactly; the run checks that rather than assuming
it, and gets 6.3e-16.

## 2. At γ = 0 the middle seat is exactly empty, or exactly twice as full

Site 5 occupancy, γ = 0:

| t | ONE | MIX | SYM | ANTI | SYM/MIX |
|---|---|---|---|---|---|
| 2 | 0.148424 | 0.148424 | 0.296849 | −0.000000 | 2.000000 |
| 4 | 0.003418 | 0.003418 | 0.006837 | −0.000000 | 2.000000 |
| 6 | 0.154835 | 0.154835 | 0.309669 | 0.000000 | 2.000000 |
| 16 | 0.230590 | 0.230590 | 0.461180 | −0.000000 | 2.000000 |

ANTI never puts anything on site 5, to 1.9e-16 over the whole trajectory. MIX's
own occupancy there ranges over [0.000000, 0.301154] and the ratio 2.000000
holds across that entire range rather than on a plateau. No max/min ratio is
quoted: an occupancy that passes through zero makes such a ratio arbitrarily
large, which is why an earlier draft quoted 67 from four sampled rows.

**The MIX column is degenerate at this seat, and that is worth admitting.** At
the reflection-fixed seat MIX equals ONE to all six digits at every time, by
mirror symmetry, so here the incoherent control and the one-sided control
coincide. MIX earns its place at the other ten seats, where they differ, and as
the exact linearity check.

The other ten sites do not show "nothing"; they show the complement. Over the
time-averaged profile SYM exceeds MIX by +0.08093 at the middle and falls short
by between 0.00576 and 0.01007 at each of the others, and those ten deficits sum
to 0.08093. That is trace conservation, not a second finding.

The factor is 2 and not 4, and this is a clarification rather than a correction.
`BORN_RULE_MIRROR.md` §4.3 is labelled **Tier 3** and already says the 4 "cancels
in normalization"; it never claimed a measurable 4. What the run adds is the
number that survives normalization: against an incoherent two-wave control at
equal budget, the constructive factor at the antinode is exactly 2. Since
MIX = ½(SYM + ANTI), the antinode's excess and the node's deficit are the same
number, so the two columns are not independent evidence.

## 3. Blindness is about support, not about purity

An earlier draft of this section said purity protects. It does not.

What protects is that the state's **support** lies in a subspace that is both
inside the silent part of the channel and invariant under H. The
single-excitation action of the jump is `Z_k = I − 2|k⟩⟨k|`, and at the
reflection-fixed seat f every reflection-odd state has ⟨f|ψ⟩ = 0, so
`Z_f ψ = ψ` and the dissipator is not weak on it but absent. Purity is
irrelevant: a mixture of two reflection-odd states at **purity 0.5200** is
exactly as blind, moving by 4.1e-16 at γ_5 = 0.05 and 3.7e-16 at γ_5 = 50.

Mirror covariance is not the protecting thing either, and the clean evidence is
one row below. The support {0, 10} is a reflection-symmetric placement, so
the Lindbladian is reflection-covariant there, and the node fills anyway. Covariance
is not even available for the asymmetric placements: rows three and five are not
reflection-covariant Lindbladians at all.

Max occupancy on site 5 and minimum purity, ANTI preparation, γ = 0.05 on the
listed support (purities rounded to six decimals; the run file carries nine, and
the first two rows agree there to all nine):

| dephasing support | max occ at f | min purity |
|---|---|---|
| none (γ = 0) | 1.9e-16 | 0.999999822 |
| the reflection-fixed seat {5} only | 4.0e-16 | 0.999999822 |
| the left end {0} only | 3.7e-02 | 0.475369 |
| {0, 10}, symmetric | 5.8e-02 | 0.258160 |
| {1, 9}, asymmetric | 4.6e-02 | 0.371813 |
| every site | 8.9e-02 | 0.091525 |

At this N every seat but the centre is non-blind, and §5 shows that is a fact
about 11 being prime rather than a general one: at N = 15, dephasing at site 2
alone leaves a two-dimensional blind subspace.

## 4. On the mediator bridge, with the observable that decides it

`Topology.MediatorBridge(3)` at N = 11 is mirror-symmetric about site 5, its
meta-mediator, and at uniform J = 1 that topology *is* the uniform chain, so the
first rows below re-run §3's system under a γ ladder.

**The purity column cannot establish a no-op.** Its truncation floor here is the
γ = 0 row's own deficit, 1 − 0.999999822 = 1.8e-7, so two identical nine-digit
purities bound the dissipative effect only at that level. (The mirror-kept arm
below sits on its own floor, 4.8e-7, 2.7 times larger.) The observable that
settles it is the same state run twice, once at the given γ and once at γ = 0,
compared over the whole trajectory:

| max\|ρ_γ − ρ_0\| | γ = 0.05 | γ = 0.5 | γ = 5 | γ = 50 |
|---|---|---|---|---|
| mediator (site 5) | 6.7e-16 | 5.8e-16 | 4.7e-16 | 6.2e-16 |
| one seat off (site 4) | 9.4e-02 | 3.2e-01 | 3.7e-01 | 4.0e-01 |

Flat at machine level across four decades of γ against a monotone rise to 0.40:
fifteen orders of magnitude, and no tolerance is asked to carry anything.

**A control the second draft got wrong.** That draft detuned bond (1,2) to
J = 1.6, called it "breaks the mirror and changes nothing else", and concluded
from it. That is the `CAUGHT_ERRORS.md:877` shape on the J axis: the detune
changes the couplings too, so alone it cannot say which did the work. The fourth
arm is the control that can. It detunes (1,2) **and** (8,9), a *larger*
perturbation (ΣJ = 11.2 against the broken arm's 10.6 and the intact 10.0) that
**preserves** the mirror, since under j ↦ 10 − j those two bonds exchange.

| max occupancy at site 5 | γ = 0.05 | γ = 0.5 | γ = 5 | γ = 50 |
|---|---|---|---|---|
| mediator, mirror intact | 4.0e-16 | 5.3e-17 | 6.5e-17 | 6.3e-18 |
| one seat off (site 4), mirror intact | 2.5e-02 | 6.9e-02 | 2.0e-01 | 2.8e-01 |
| mediator, mirror **broken** | 2.6e-02 | 2.2e-02 | 6.1e-02 | 6.5e-02 |
| mediator, mirror **kept**, bigger detune | 1.2e-16 | 1.1e-16 | 6.3e-17 | 3.1e-18 |

The larger, mirror-preserving perturbation leaves the mediator blind; the
smaller, mirror-breaking one hands it its bite back. **At the reflection-fixed
seat it is the mirror, not the mediator and not the coupling budget.** The scope
matters: §5's blind seats at N = 12, sites 1, 4, 7 and 10, have nothing to do
with reflection, so the mirror is the odd-N-centre special case of the divisor
condition and not the general mechanism.

Each row's minimum purity is constant *down the γ ladder* (0.999999822 for the
intact arm, 0.999999520 for the mirror-kept arm); the two arms differ from each
other, which is why the floor above is called grid-and-Hamiltonian dependent.
The comparison across γ is like-for-like: `stable_dt` would shrink dt once
max\|mask\| exceeded 150, and with γ on a single seat max\|mask\| = 2γ_f reaches
only 100 at γ = 50, so every row ran at dt = 0.01 and `stable_dt` never fires
anywhere on this page.

Two readings of "harm" part company at large γ and both are reported. The
occupancy column rises monotonically; the "one seat off" **minimum purity** does
not, running 0.574, 0.133, 0.144, 0.632 across the four rates, as strong local
dephasing begins to freeze the state it is destroying.

## 5. How many such states there are, exactly

**The condition, stated correctly.** `Z_k = I − 2n_k` has eigenvalues +1 on the
configurations avoiding k and −1 on those occupying k, and the dissipator is
silent on any ρ with no coherence **across** those two eigenspaces. It is *not*
true, as a second draft of this page had it, that silence requires lying in
ker(n_k); that is the +1 branch only, and §6 shows what the wrong version cost.

For a subspace of **pure states** the count is unaffected, which is why the
numbers below stand. A pure state blind at k must lie wholly inside one
eigenspace; in the single-excitation sector the −1 eigenspace is the single
ray |k⟩, and H hops the excitation off k whenever site k carries an incident bond
with J ≠ 0, so that branch is not H-invariant and contributes 0. That hypothesis
holds for every topology on this page but is not automatic: an isolated seat, or
a bond detuned to zero, would make the −1 branch one-dimensional. The remaining
branch is the largest H-invariant subspace inside ker(n_k), which for Hermitian
H is the orthogonal complement of the block Krylov space generated by the
configurations occupying k. That is a **rank**, taken over GF(p) on the integer
sector matrix at J = 1, with no eigensolver. One prime gives only
rank_p ≤ rank_ℚ, which would make the reported dimension an upper bound, so
every value below is computed at two primes, 2³¹ − 1 and 2⁶¹ − 1, and they agree
everywhere.

| N | dimensions, sites 0 … N−1 |
|---|---|
| 8 | 0 0 0 0 0 0 0 0 |
| 9 | 0 1 0 0 **4** 0 0 1 0 |
| 11 | 0 0 0 0 0 **5** 0 0 0 0 0 |
| 12 | 0 1 0 0 1 0 0 1 0 0 1 0 |
| 15 | 0 1 2 0 1 0 0 **7** 0 0 1 0 2 1 0 |
| 21 | 0 1 0 3 1 0 0 1 0 0 **10** 0 0 1 0 0 1 3 0 1 0 |

**dim = (gcd(2j+1, N) − 1)/2**, at every site for N = 3..21.

**Derivation, from `PROOF_UNIFORM_LAW.md` B0.** The spectrum there is
nondegenerate, so every H-invariant subspace is spanned by eigenvectors and the
blind subspace is spanned by exactly the modes with a node at site j. With
ψ_m(j) ∝ cos(πm(2j+1)/(2N)) and m = 0..N−1 that node condition is
(2j+1)m ≡ N (mod 2N). Write g := gcd(2j+1, 2N); since 2j+1 is odd, g is odd and
divides N, so g = gcd(2j+1, N) and no parity case split on N is needed. The
congruence has exactly g solutions in [0, 2N). Of those, m = 0 is never one (it
would need N ≡ 0 mod 2N), m = N always is ((2j+1)N − N = 2jN ≡ 0 mod 2N), and
m ↦ 2N − m maps solutions to solutions while fixing only 0 and N. So the g
solutions are m = N together with (g−1)/2 mirror pairs, exactly one member of
each lying in [1, N−1]. Counting those gives (g−1)/2.

**The same count holds in F64's own block, and that is where a reader of F64
would look for it.** The exact identity is EQ-015's, −Re(λ_k) = 2γ_B·|v_k(B)|²,
mode by mode at every γ_B. Its generator L_coh = i·H₁ − 2γ_B|B⟩⟨B|
(`EMERGING_QUESTIONS.md:508`) is derived for an **XY** Hamiltonian; on the
Heisenberg chain the vacuum bra carries E_vac = (N−1)J, so the block generator
is M = −i(h_SE − E_vac·Id) − 2γ_B|B⟩⟨B|, which is exactly **F152**'s
M = −2i·𝓛_J − 2·diag(γ) once h_SE = (N−1)J·Id − 2J·𝓛 is substituted. The run
uses F152's form; the XY form differs by a sign on the Hamiltonian part and an
imaginary shift, neither of which touches Re λ, so the count is the same either
way. The bijection is two lines and needs no eigensolver: if v is an H₁
eigenvector with v(B) = 0 then Mv is i times a real multiple of v, so Re λ = 0;
conversely Re λ = 0 forces |v(B)|² = 0 by the identity, and v is then an H₁
eigenvector with a node at B. So the undamped modes of the coherence block are
exactly the node modes counted above, and **(gcd(2j+1, N) − 1)/2 counts them
there too**. The `coherence` run checks this at every site for N = 9, 11, 12 and
15 at γ_B = 0.01 and 0.7, with zero mismatches. It uses a general eigensolver on
a non-normal matrix, so it is a check and the bijection carries the claim; its
zero-cut is reported beside the gap it sits in, the largest |Re| counted as zero
(at most 3.9e-15) against the smallest excluded (at least 2.6e-04). Twelve
decades, not a threshold chosen to pass.

**Why the composite half was invisible.** The existing data live at N = 5 and
N = 11, both **prime**, where gcd(2j+1, N) = 1 at every non-centre seat. On a
prime chain the reflection-parity reading and the divisor law give identical
answers at every seat, so no run at N = 5 or N = 11 could distinguish them. The
first discriminating case is composite, and the smallest is N = 9 at site 1.

**Several dephased seats.** A blind subspace is a span of modes, so for a support S
it is the intersection, and

  **dim blind(S) = #{ m ∈ 0..N−1 : (2j+1)m ≡ N (mod 2N) for every j ∈ S }.**

Checked against the rank for every 1-, 2- and 3-site support at N = 9, 11, 12
and 15, with no exceptions. At N = 15 site 7 carries modes {1,3,5,7,9,11,13},
sites 2 and 12 carry the *same* set {3, 9}, so {2, 12} still gives 2 and {7, 2}
gives 2, while {1, 4, 13} gives {5}, dimension 1. A second draft of this page
called two seats "no closed form yet" and, before that, claimed any second seat
collapses the subspace to zero, which is true at N = 11 and false at N = 9 and
N = 15.

**The odd-N centre.** There 2j+1 = N, the condition collapses to m odd, and the
dimension is (N−1)/2 = ⌊N/2⌋, which is
`MEDIATOR_NOISE_GATE_LEVEL_THREE.md:260`'s number. That subspace **is** the
reflection-odd site space as a set and not merely as an integer, verified in
exact integer arithmetic for the odd N from 5 to 21, and derivable in one line
from `PROOF_R90_FROZEN_DIVISOR.md:193`, where site reversal acts on these modes
as the diagonal sign (−1)^k, so the odd-k modes are exactly the reflection-odd
site space of line 195. This is object 1 of the ⌊N/2⌋ map, and its unpaired seat
is the dephased seat itself, so the fixed-seat parity agrees as well as the value.
At even N there is no such seat, blind seats still exist (N = 12 at sites 1, 4,
7, 10, one dimension each), and the reflection-odd space is the wrong answer
there, six dimensions against one. F140's ⌊N/2⌋ is a third thing again, which
`PROOF_R90_FROZEN_DIVISOR.md:335` disclaims as "not a decoherence-free
structure".

## 6. The all-sites theorem, and an inconsistency this state exhibits

A second draft of this page wrote that at full support "the blind dimension is 0
in every sector, and trivially so", and reconciled the committed results with
that. It is wrong, and the counterexample is the simplest state there is:
ρ = I/N on the single-excitation sector, under uniform γ = 0.5 on all eleven
sites, gives max |ρ(t) − I/N| = **0.000e+00**. Exactly blind, full support.

That is not a contradiction of the committed work; it *is* the committed work.
`SYMMETRY_CENSUS.md:101` reports exactly one steady state per sector, the
maximally mixed one, and `PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md:70` gives the
full-support decoherence-free subalgebra as the diagonal operators, a
d_w-dimensional algebra. The two numbers, d_w and 1, are reconciled by the
second condition rather than in conflict: the diagonal algebra is where the
dissipator vanishes *instantaneously*, and requiring the trajectory to stay
there adds H-invariance, which on a connected graph cuts d_w down to the single
ray P_w/d_w. So the correct statement is narrow and worth stating narrowly: **at
full support the only blind trajectories are the sector steady states; no
non-stationary trajectory survives.** What this page counts is a blind space of
non-stationary states, and it needs the support to be a strict subset.

**What that proof requires, and what this page's state says about it.**
`PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md` holds for a **connected** graph with
**γ_k > 0 at every site**, which is what Step 2 part (a) (line 70) and the
separate Step 2b spend, what Step 2 part (b) (line 72) needs the connectivity
for, and what its closing sentence (line 92) and scope block (line 112) state.
This page's state is one reason the γ half cannot simply be widened to
γ_k ≥ 0: γ = 0.5 on site 5 alone, run to t = 200, leaves max |ρ − I/N| = 0.268
at purity 0.999998, so that diagonal sector block does not reach P_w/d_w.

The hypothesis is sufficient rather than necessary, and where it fails the
attractor **refines rather than dissolves**. The mechanism is one commutator:
with M the site reversal, Z_k commutes with M only at the reflection-fixed seat
(measured 0 there against 2 at every other seat of an N = 5 chain), so a support
inside that fixed set leaves the two parity projectors conserved. Those are
constants of motion beyond the N+1 sector populations, and the limit is then
maximally mixed per **(popcount, parity)** block instead of per popcount sector.
At N = 5 with γ on seat 2 alone: distance to the per-sector limit 2.604·10⁻²,
distance to the per-parity-block limit 4.0·10⁻¹⁰ (`blind_site.py scope`). The
deviation is exact rather than approximate, which is what says it is structure:
the centre seat gains exactly 1/48 of population, the four others lose exactly
1/192 each, and 5/192 survives as coherence between the reflection-partner
sites. The conserved parity-odd block is exactly this page's blind subspace.

## 7. The two mediator-inertness findings, and how they stand

`MEDIATOR_NOISE_GATE_LEVEL_THREE.md` measures site 5 at 15.57 % response and
concludes the node structure is "measurably inert". `GAMMA_CONTROL.md:31` says
"the mediator's own noise always harms; it merely harms least". Both stand and
this page corrects neither, because both measure a **different preparation**: a
Bell-on-vacuum state living in popcount {0, 2}.

The rank alone would not settle that. In the popcount-2 sector the blind
dimension at N = 11 is 0 at every one of the eleven seats, but that rank covers
the (2,2) block while such a state also occupies (0,0) and (0,2). The direct
check does settle it: propagated whole at γ = 0.5, the Bell-on-vacuum state
moves by 2.84e-01 with γ on the mediator and 2.82e-01 with it on site 0.
There is no blind state anywhere in that family, which is why the mediator
looked ordinary to it.

What changes is the reading, not the measurement. The mediator's node structure
is not inert as a structure; it is invisible to a state prepared outside it. A
cheap prediction follows: repeat the site-by-site γ sweep of NextStep (1) with a
reflection-odd single-excitation preparation, and the response at site 5 should
be exactly zero rather than fourth of nine.

## 8. What this is not

It is not a standing wave. Restricted to git-tracked `.md` files, that phrase
appears in 125 documents in this repo and means a Π-paired eigenvalue pair, a
Pauli-string oscillation, or a two-observer cross term, none of them spatial,
and "node" there already means a Pauli string that does not oscillate; this page
uses the phrase only when quoting the two test titles and in this sentence. It
is not "dark" either, which here means the ⟨n_XY⟩ = 0 of the Absorption Theorem
(**AT**, `docs/ANALYTICAL_FORMULAS.md:316`, not F123, which is the closure
functional) and the F135/F136 record classes. It is not a new solution of the
sector, which F126 solved. It is not a statement about `Cone.cs` or
`cone_defect_arrival.py`, both of which build the XY chain.

It is not a protection scheme. The blind subspace is blind to a **named
support** that must be known and strict by construction. In any setting where
γ sits on every seat it is gone, by §6, so it is not hardware-relevant as it
stands.

## 9. Seven errors this run made, and how each was caught

Three in the arithmetic. The first rank routine used an SVD and returned blind
dimensions that were not mirror-symmetric on a mirror-symmetric chain, which is
impossible, and is why the exact GF(p) route replaced it. The replacement
iterated the block Krylov to a fixed point of the *row count* rather than the
rank; an unreduced block has more rows than rank, so it exited early and
invented a blind dimension of 21 at the chain ends. A naive float power basis
[e, He, H²e, …] used for one cross-check reported an overlap of 0.73 where exact
integer arithmetic gives 0, and dimension 3 where the answer is 7. None was
found by inspection; each was found by a symmetry or an exactness the answer had
to satisfy.

Two in the physics, and these are the ones worth remembering. The blindness
condition was stated as an iff on ker(n_k), which is only one of the two
branches of Z_k, and §6's reconciliation was then built on the wrong version and
came out false. And "purity protects" was a mislabel for "the support does",
which a mixture at purity 0.52 refutes in one line. Both were found by an
outside reader asking for a counterexample, not by re-reading.

Two in the prose: the §5 generalisation from two measured cases, and §4's
unmatched control, the same shape as the `CAUGHT_ERRORS.md:877` entry this page
cites in its own sweep.

## 10. What is borrowed and what is not

| | source |
|---|---|
| rate = 0 ⟺ node at the dephasing site | F64 (line 1852), `EMERGING_QUESTIONS.md:542` |
| "blind to dephasing" as the name, and interior-B nodes | `PRIMORDIAL_GAMMA_CONSTANT.md:133, 144` |
| "blind subspace" as a typed term with a Meta-Theorem | `ORTHOGONALITY_SELECTION_FAMILY.md:37, 279` |
| the two-condition DFS logic, and that (i) alone fails | `J_BLIND_RECEIVER_CLASSES.md:32-40` |
| h_SE = (N−1)J·Id − 2J·𝓛, Neumann modes, nondegeneracy | `PROOF_UNIFORM_LAW.md:182-186`, gated entry-exactly |
| λ_k = 4cos(kπ/N) + N − 5 | `PROOF_R90_FROZEN_DIVISOR.md:176-178` |
| the ZZ term as the mechanism, Neumann against Dirichlet | F2 line 111, D10 line 180 |
| the (0,1) block generator M = −2i·𝓛_J − 2·diag(γ) | F152 line 7058 |
| ⌊N/2⌋ nodes at the centre of an odd chain | `MEDIATOR_NOISE_GATE_LEVEL_THREE.md:257-260` |
| the centre subspace = the reflection-odd site space | `PROOF_R90_FROZEN_DIVISOR.md:193, 195` |
| the genre, an immune count at one named seat | F66, verified at the endpoint, interior open |
| the sector's exact solution | F126, `PROOF_DEPHASING_FRONT_RENEWAL.md` |
| **the count as a function of the site, (gcd(2j+1,N)−1)/2** | this page |
| **blind seats at even and composite N** | this page |
| **the intersection law for several dephased seats** | this page |
| **the same count inside F64's own (0,1) block** | this page |
| **the mediator table to γ = 50 with its matched control** | this page |
| **the state that bounds the asymptotic theorem's γ hypothesis** | this page |

## 11. Open

- Whether the law survives on the **XY** chain at modulus N+1, where the node
  condition becomes m(j+1) ≡ 0 (mod N+1) and a different divisor answer is
  expected. Neither `Cone.cs` nor `cone_defect_arrival.py` can decide it: both
  damp with one uniform scalar and cannot place γ on a single seat.
  `PRIMORDIAL_GAMMA_CONSTANT.md`'s N = 5 B-position scan is already data on that
  lattice.
- The blind **operator** space inside the (1,1) block, which is larger than the
  state count here and is the number F66 would compare against.
- Whether anything is blind in **popcount ≥ 2** for any topology, or whether the
  N = 11 zeros of §7 are a theorem. Those zeros are certified at one prime
  already, since a Krylov space of full rank mod p has full rank over ℚ.
- The §7 prediction: the site-by-site sweep repeated on the reflection-odd
  preparation.
- Whether every failure of the asymptotic sector projection's all-sites
  hypothesis **refines** the way §6's does, under some graph automorphism whose
  fixed set contains the dephasing support. At N = 5 the one failing support is
  exactly the reflection-fixed seat and the refined limit is exact, but that is
  one instance on one graph. Settling it would turn an all-sites hypothesis into
  a statement about which symmetries the support leaves alive.
- `review/OBC_SINE_BASIS_FINDINGS.md` asserted for four months that Heisenberg
  single-excitation eigenvalues "do not follow any simple cos formula" while its
  own tabulated N = 3 row, (−4, 0, 2), is λ_k = 4cos(kπ/N) + N − 5 exactly. It
  now carries the closed form, and `simulations/eq021_obc_sine_basis.py` Phase 1
  computes and gates it. Recorded here because the false clause outlived three
  commits of that file and nothing else in the repo contradicted it.
- `experiments/CONCENTRATOR_OPTICS.md:106-109` reports that at N = 5 the same γ
  budget on the middle site gives 2618 against the edge's 352, and N = 5 centre
  is j = 2 with gcd(5, 5) = 5 and blind dimension 2. That file's next line says
  "Both numbers above are real; what they measure is not settled", so this is a
  lead and not an explanation. Not checked here.
