# The Blind Site: counting what dephasing at one seat cannot reach

**Date:** August 23, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [`simulations/blind_site.py`](../simulations/blind_site.py), which imports its propagator from [`simulations/bridge_sector.py`](../simulations/bridge_sector.py)
**Data:** [`simulations/results/blind_site/blind_site_run.txt`](../simulations/results/blind_site/blind_site_run.txt)
**Second script (§7 only):** [`simulations/blind_seat_mi_sweep.py`](../simulations/blind_seat_mi_sweep.py), parts `gate | sweep | algebra | support | zeno | converge`
**Second data file:** [`simulations/results/blind_site/blind_seat_mi_sweep.txt`](../simulations/results/blind_site/blind_seat_mi_sweep.txt)

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
only 100 at γ = 50, so every row of THAT ladder ran at dt = 0.01 with `stable_dt`
inert. It is not inert everywhere on this page: the Zeno scan below reaches
γ_centre = 100 over a baseline, max|mask| passes 150, and the step is shrunk
there.

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
a bond detuned to zero, would leave that ray H-invariant, and it would then
contribute 1 rather than 0. The remaining
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
first discriminating case is composite, and the smallest is **N = 6**, whose
blind seats 1 and 4 carry one dimension each. An even chain has no
reflection-fixed seat at all, so a reflection-parity reading predicts no blind
seat there and the divisor law predicts two; the dimension table above prints
that row. N = 9 is the smallest ODD discriminating case.

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
Bell-on-vacuum state living in popcount {0, 2}. The measurements stand
untouched; one reading does not, and `GAMMA_CONTROL.md:31`'s "always" is
corrected in scope below by a counterexample on a preparation that page did not
run.

The rank alone would not settle that. In the popcount-2 sector the blind
dimension at N = 11 is 0 at every one of the eleven seats, but that rank covers
the (2,2) block while such a state also occupies (0,0) and (0,2). The direct
check does settle it: propagated whole at γ = 0.5, the Bell-on-vacuum state
moves by 2.84e-01 with γ on the mediator and 2.82e-01 with it on site 0.
There is no blind state anywhere in that family, which is why the mediator
looked ordinary to it.

What changes is the reading, not the measurement. The mediator's node structure
is not inert as a structure; it is invisible to a state prepared outside it.

### The prediction, run

This section used to end with a cheap prediction: repeat the site-by-site γ
sweep of `MEDIATOR_NOISE_GATE_LEVEL_THREE.md` on a reflection-odd
single-excitation preparation, and the response at seat 5 should be exactly zero
rather than fourth-smallest of the nine interior sites, which is where that
sweep had put it. `simulations/blind_seat_mi_sweep.py` runs it on
(|1₀⟩ − |1₁₀⟩)/√2, the excitation on the first site minus the excitation on the
last. The runner is gated before it is pointed anywhere new: on Bell-on-vacuum
it returns all eleven published rows at that run's own step dt = 0.05, worst
deviation 0.005 against a table printed to two decimals (part `gate`).
Everything after the gate runs at dt = 0.02, for the reason under convergence
below.

**Run exactly as predicted, it does not give zero: it gives 0.0809 %.** Zero is
what comes out in a configuration the prediction did not name, one with no γ
anywhere but the swept seat, and there it is exact in the strongest available
sense. So the prediction holds where it can hold exactly and misses where it was
aimed, and the six points below are what the run put in place of it.

Two definitions first, because every percentage in this section depends on them.
**A** = {0, 1, 2, 3, 4} and **B** = {6, 7, 8, 9, 10} are the two five-site halves
the centre separates, and the centre is in neither, which is why it is the
bottleneck the sixth point turns on. A **span** is what the committed sweep
reports and what every percentage here means: the seat named in a row is moved
from γ = 0 to γ = 0.5 while every other site holds the baseline, and the span is
100·(value at γ = 0 − value at γ = 0.5)/(value at γ = 0), so a positive span
means the dephasing lowered the correlation.

| baseline γ elsewhere | seat 5 | the other ten seats |
|---|---|---|
| none (only the swept seat carries γ) | −2.2·10⁻¹⁴ % | not swept |
| 0.05, the sweep's own | 0.0809 % | 12.19 % to 17.51 % |

Both rows are spans of the window mean of I(A:B) over t ∈ [0, 20] at dt = 0.02
(parts `sweep` and `converge`). The
comparison against Bell-on-vacuum is made under **that same functional**, not
against the committed peak: like for like the centre reads 7.635 % on
Bell-on-vacuum and 0.0809 % here, a factor of 94, and the Bell column then runs
5.56 % to 15.77 %. Setting the committed 15.57 % beside 0.0809 % would compare
two different functionals and overstate the contrast by about twice.

**The mechanism is the state's own support, not a parity of ρ, and this page
already said so.** The
single-excitation action of the jump is Z_k = I − 2|k⟩⟨k|; a reflection-odd
state has a node at the centre, so Z₅ acts there as the identity, measured
max|Z₅ P_odd − P_odd| = 0.0 exactly (part `algebra`), while at seats 4 and 0 the
same quantity is
1.0. The state is picked out BY its Hilbert-space reflection parity, and that is
what puts the node there (§5). What buys nothing is the parity of ρ as an
operator, and the centre-swept configurations show it cleanly. There the profile
[b, …, b, γ₅, b, …, b] is itself reflection-symmetric, so R permutes the jump set
and the Liouvillian is **exactly reflection-covariant**: measured
max|R H R − H| = 0.0 and max|R·mask·R − mask| = 0.0 at baselines 0.05 and 0.2.
And ρ₀ = |ψ⟩⟨ψ| with Rψ = −ψ is Ad_R-**even**, max|Rρ₀R − ρ₀| = 0.0. So in
exactly the configuration where the blindness is exact, the Ad_R grading is
conserved and the state sits in its trivial half; it cannot be what protects.
The other ten rows of the sweep are not reflection-covariant at all, since a
single off-centre seat's profile is not reflection-symmetric, which §3 above
already states and which is a second reason covariance cannot be the
explanation. It is not quite the distinction
[PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)
draws, and the difference is worth keeping: that proof works at the commutator,
Z_k commuting with the site reversal only at the mirror-fixed seat, and reads the
consequence as two conserved mirror-parity projectors whose odd block it names as
this page's blind subspace. Both are true of the same object. Its projectors are
Hilbert-space parity and are the right reading; the operator grading of ρ is the
one that is inert here.

**The peak functional is saturated, and the repo knew that before this run.**
The reflection-odd state begins at I(A:B) = 2 bits and never rises above it, so
its peak is the t = 0 value at every seat and every γ. With no dephasing anywhere
it does not even fall: the odd subspace is H-invariant and every state in it has
p₅ = 0, so p_A = p_B = ½ for all t and I(A:B) ≡ 2 exactly. That 2 is not a measurement here:
it is **F75** ([ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md), the
mirror-pair MI closed form) evaluated at its own saturation point p₀ = 1/2, and
F75 is independent of the mirror sign, so reflection-oddness buys the node and
not the 2 bits. That the peak of a single-excitation mirror-symmetric preparation sits at t = 0
is the observation of
[RECEIVER_VS_GAMMA_SACRIFICE](RECEIVER_VS_GAMMA_SACRIFICE.md), 2026-04-23:
"Peak MM sits at or near t = 0 … F75 predicts the peak to within 7 % without any
propagation". That is a resemblance and not the same statement: it is measured
on the **bonding modes**, which this preparation is not, and it is approximate
to 7 %, where here the peak sits at t = 0 **exactly**, forced by the 2-bit
ceiling the state already occupies. A second thing the page cannot check for
itself: F75's MM is a mirror-PAIR mutual information and the sweep's I(A:B) is a
halves mutual information, and they agree at t = 0 only because the one occupied
mirror pair straddles the A/B cut while every other site is empty and the state
is pure.
The window mean used instead is a partial answer to a request
`MEDIATOR_NOISE_GATE_LEVEL_THREE.md` already had open, for "a window-stable
functional (steady-state MI, or MI integrated over time)"; the answer is that
the integrated one is **not** window-stable, see below.

**The 0.0809 % is a cancellation residue, not a leak size.** It is tempting, and
this page's first draft did it, to read the residual as the size of the state's
leak out of the odd subspace. Restricting the baseline to named subsets refuses
that. The number in each row below is the **centre's own span**, unchanged in
definition, with only the set of other sites carrying the 0.05 baseline
changed (part `support`):

| sites carrying the baseline | the centre's span |
|---|---|
| all ten, which is the sweep itself | +0.0809 % |
| the two ends {0, 10} | −0.2595 % |
| the interior, both ends and the centre excluded | +0.3127 % |
| one end {0} | +0.0239 % |
| the centre's two neighbours {4, 6} | +2.0954 % |

The contributions carry both signs and the neighbours alone are an order of
magnitude larger than the total. What the run does establish is the exact statement: Tr(P_odd ρ) is
1.000000000000 across the whole window when nothing else dephases, with or
without γ at the centre, so the centre is a node for the whole trajectory and not
only at t = 0; with the baseline it is 0.469 by t = 20, on its way to the
unbiased share 5/11. Why the cancellation lands at 0.08 % is open.

**The sign flip has a name, and it is a counterexample.** The sweep's span is a
two-point difference, γ₅ = 0 against 0.5, of a function that is not monotone in
γ₅. Scanning it (part `zeno`): at baseline 0.05 the span rises to +0.189 % near
γ₅ = 2, then
crosses to −1.65 % at 20 and −10.76 % at 100; at baseline 0.2 it is already
negative at 0.5. Across the same scan the window-mean occupation of the centre
falls from 0.0663 to 0.0310. That is Zeno at the bottleneck: strong dephasing on
the seat every A-to-B path crosses blocks the transport, and the preparation
keeps more of the maximal I(A:B) it started with. The repo names this class
already (`docs/GLOSSARY.md` on ENAQT, `D06_SPECTRAL_GAP`,
`PROOF_ABSORPTION_THEOREM`). Stated rather than hidden: this is a
counterexample, **on this preparation**, to `GAMMA_CONTROL.md:31`'s "the
mediator's own noise always harms; it merely harms least", which was measured on
another.

**The percentages are a window, the ratio is the law.** Unlike the committed
peak, which that page reports as window-converged at 15.566 % for t_max = 20, 40
and 80, the window mean cannot be: past the transient every γ profile decays
toward the same limit, so a longer window dilutes every span alike. The centre
reads 0.0809 % at t_max = 20, 0.0497 % at 40 and 0.0272 % at 80, and the other
seats fall in step. What holds still is the ratio of the smallest other seat to
the centre: 150.7, 148.5, 149.3, stable to about one percent over a factor of
four in the window. That ratio is the quotable number. In the integrator the
spans are converged (0.0808, 0.0809, 0.0809, 0.0809 at dt = 0.05, 0.02, 0.01,
0.005) even though the curve itself is not: with no dephasing anywhere the run
is pure, so I(A:B) = h(p_A) + h(p_B) − h(p₅) ≤ 2 bits exactly, and at dt = 0.05
the **peak** of that pure arm exceeds the bound by 1.3·10⁻³. That excess is
common to both arms and cancels in the difference, which is why the spans are
converged anyway; the sweep is run at dt = 0.02 regardless, where it is
3.5·10⁻⁵ (part `converge`).

The −2.2·10⁻¹⁴ in the first row is **not** an integrator floor, and calling it
one would be a mislabel: with nothing on the other sites the two arms are the
same computation, since Z₅ is the identity on the subspace the state never
leaves. The two arms are not bit-identical, and the printed column says so: −2.220·10⁻¹⁴,
+2.220·10⁻¹⁴, −2.220·10⁻¹³ at t_max 20, 40 and 80, a sign flip and a tenfold
growth with the window. What is measured is a difference at the rounding floor,
consistent with exact blindness; the exactness itself rests on the argument
above, that Z₅ is the identity there, and not on this number.

The odd column is exactly reflection-symmetric, seat j against seat 10 − j. That
is forced rather than observed, and not by the covariance of a single
run: R carries the seat-j profile to the seat-(10 − j) profile, ρ₀ is Ad_R-even,
and R swaps A with B while I(A:B) is symmetric, so the two runs are
Ad_R-conjugate to each other. The Bell column is not symmetric because its preparation
sits on sites 0 and 1.

**What the exact-blindness argument BELOW uses**, beyond "an odd chain":
couplings that are reflection-symmetric about the seat, excitation number
conserved so the single-excitation sector is invariant, and γ = 0 at every
non-centre seat. The three are used together, and the argument does not run
without all of them; what §11 now adds is that the reflection is not NECESSARY
for the blindness, being neither necessary nor sufficient for the divisor law's
value and buying only the one seat it is symmetric about. Excitation-number
conservation is a different matter: the criterion is stated inside the
single-excitation sector, so that hypothesis is not one §11 lifts. At the reflection-fixed seat it does **not** need Heisenberg
over XY, but not for any reason about the ZZ term, which is exactly what moves
the nodes elsewhere by turning the adjacency matrix into the Laplacian and the
modulus from N+1 to N (F2 line 111, D10 line 180, and the XY question §11 used to carry, now answered there).
The reason is weaker and sufficient: any reflection-symmetric single-excitation
Hamiltonian leaves the reflection-odd space invariant, and every state in that
space has amplitude zero at the fixed seat by definition, ψ(f) = −ψ(f). That
argument names no eigenbasis, so it carries to XY unchanged. Away from the fixed
seat the count does depend on the modulus, and §11 now carries both closed
forms. Nor should it be specific to Z-dephasing: any jump L that
ANNIHILATES the odd subspace leaves it alone, and so does any L satisfying BOTH
Lψ = λψ and L†Lψ = |λ|²ψ on it. Amplitude damping at the centre is the
annihilating branch. Both halves of the second condition are needed and neither
implies the other, which is what Z_k supplies by fixing the subspace pointwise
AND being unitary. Drop the isometry: L = P_odd + |v⟩⟨e_c|, with v odd and c the
reflection-fixed seat, has Lψ = ψ for every odd ψ, yet L†Lψ = ψ + ⟨v|ψ⟩·e_c, so
for ρ = |v⟩⟨v| the dissipator is −½(|e_c⟩⟨v| + |v⟩⟨e_c|)·⟨v|v⟩ and does not
vanish. Drop the eigenvector half: where the odd space has two directions, at
N = 5 and above, L = |w⟩⟨v| with v and w distinct unit odd vectors has
L†Lv = v while sending ρ = |v⟩⟨v| to |w⟩⟨w| − |v⟩⟨v|. So "acts as an isometry
there" is not the condition on its own either. All of that is an argument, not
a measurement, and is the one claim in this section the run does not check. The
open boundary is needed for the cosine derivation.

Run: `python simulations/blind_seat_mi_sweep.py`; output in
`simulations/results/blind_site/blind_seat_mi_sweep.txt`.

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

## 9. Eight errors this run made, and how each was caught

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

**An eighth, from the §7 run, and it is the worst of the eight because the page
had already made it once.** The first draft of the §7 result explained the
residual as "Z_k does not preserve reflection parity, so the baseline leaks the
state out". The wrong object was the parity of ρ: where the blindness is
exact the Liouvillian is exactly reflection-covariant and ρ₀ is Ad_R-even, both
measured at 0.0, so the grading the sentence invoked is conserved and the state
sits in its trivial half. Hilbert-space reflection parity is a different object
and is not conserved in the swept configuration at all: Tr(P_odd ρ) falls to
0.469 by t = 20. The correct statement is about the state's support, which is §3
of this same page, three sections above where the error was written.
The same conflation, superoperator grading against Hilbert-space grading, had
been caught and reverted the day before in a different set of documents; it came
back the moment a new result needed a sentence. The lesson is not "check
parity claims" but that a mechanism sentence written in a *grading* word is worth
suspecting on sight when the object at hand is a subspace.

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
| MI = 2 bits for a mirror pair at p = 1/2, and its η-independence | **F75**, `docs/ANALYTICAL_FORMULAS.md` |
| the peak of a mirror-symmetric single excitation sitting at t = 0 | `RECEIVER_VS_GAMMA_SACRIFICE.md:281`, as a resemblance: measured on bonding modes, approximate to 7 % |
| mirror-parity projectors conserved for a fixed-seat support, and their odd block named as this page's | `PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md` |
| the Zeno / ENAQT class the §7 sign flip belongs to | `docs/GLOSSARY.md`, `D06_SPECTRAL_GAP`, `PROOF_ABSORPTION_THEOREM.md` |
| **the count as a function of the site, (gcd(2j+1,N)−1)/2** | this page |
| **blind seats at even and composite N** | this page |
| **the intersection law for several dephased seats** | this page |
| **the same count inside F64's own (0,1) block** | this page |
| **the mediator table to γ = 50 with its matched control** | this page |
| **the state that bounds the asymptotic theorem's γ hypothesis** | this page |
| **the §7 sweep on a blind preparation, gated against the committed eleven rows** | this page |
| **the like-for-like contrast under one functional, 7.635 % against 0.0809 %** | this page |
| **the residual as a cancellation of both signs rather than a leak size** | this page |
| **the sign flip as a counterexample, on this preparation, to `GAMMA_CONTROL.md:31`'s "always"** | this page |

## 11. Open

**The XY chain is no longer open, and this page's own law turns out to be one
evaluation of a wider one.** [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) settles
both. Its §4 gives the XY answer predicted here, on the **uniform** chain: the
node condition is m(j+1) ≡ 0 (mod N+1) and the count is **gcd(j+1, N+1) − 1**,
with no halving. The proof is two lines. With d = gcd(j+1, N+1) the solutions are
the multiples of (N+1)/d, and d−1 of them lie in 1..N. Verified three ways, the
third being an exact integer enumeration to N = 200 which certifies the count,
not the kernel. The two laws disagree loudly: at N = 6 and N = 12 the Heisenberg
chain has blind seats and the XY chain has none, while at N = 11 Heisenberg has
one and XY has seven. At the reflection-fixed seat of an odd chain they agree,
both giving (N−1)/2, which is §7's claim reached a second way. The reflection
argument stays the better one, because it names no eigenbasis; what the closed
forms add is independent agreement with it.

**And the divisor law is a uniform-chain law**, which this page did not say and
should have. That page's §2 gives the general form:
blind(j) = deg gcd(χ(H_left), χ(H_right)), exactly, for any bond profile with no
zero bond, where χ is the characteristic polynomial and H_left, H_right are the
two principal submatrices of the single-excitation Hamiltonian that the seat
leaves behind. Its reason: H_SE on such an open chain is a Jacobi matrix, so an
eigenvector vanishes at a site exactly when its eigenvalue is shared by those two
blocks. Two words of that page's own correction belong here, because this
paragraph first carried both errors: the dephasing does **not** cut the chain, it
clears the watched seat's row and column and leaves coherences that cross the
seat untouched; and the two blocks are principal submatrices, not the
free-standing subchains; the two readings agree on the XY chain, where there is
no diagonal to differ in, and disagree on every Heisenberg row that page
measures. The no-zero-bond hypothesis is F143's fence for a neighbouring
object, and here a zero bond makes the criterion wrong rather than inapplicable.
Reflection symmetry is neither necessary nor sufficient: a reflection-symmetric
N = 7 chain has more blindness than the divisor law and a reflection-symmetric
N = 9 chain has less, while an asymmetric N = 5 chain reaches the law's own
centre value. What a reflection **about seat j** buys is that one seat, which is
what §7 uses it for and all it may be used for.

What is still open:

- The blind **operator** space inside the (1,1) block, which is larger than the
  state count here and is the number F66 would compare against.
  [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) §3 measures a neighbouring
  object, the **stationary** operators of the whole Liouvillian on the
  single-excitation sector, and finds their count is exactly 1 + the blind state
  dimension, for any bond profile with no zero bond, with the same Jacobi
  simplicity as the reason. The two objects are still not shown to coincide: one is killed by the
  dissipator, the other by the dissipator and the commutator together.
- Whether anything is blind in **popcount ≥ 2** for any topology, or whether the
  N = 11 zeros of §7 are a theorem. Those zeros are certified at one prime
  already, since a Krylov space of full rank mod p has full rank over ℚ.
- The §7 prediction is **run** and reported there. What it leaves open is the
  size of the residual at the centre, 0.0809 % at baseline γ = 0.05. §7 shows it
  is a cancellation between contributions of both signs (ends −0.2595 %,
  interior +0.3127 %, the centre's neighbours alone +2.0954 %), so there is no
  single quantity it is the size of, and nothing here predicts where the
  cancellation lands. The spans at the centre are 7.5·10⁻³, 5.1·10⁻², 8.1·10⁻²
  percent at baselines 0.001, 0.01, 0.05 and −2.4·10⁻² at 0.2, all at
  t_max = 20, dt = 0.02. The 0.001 column is the one that is not
  step-converged, moving by 6 % between dt = 0.05 and dt = 0.005; the other
  three are stable in the fourth digit.
- Whether a **discriminating** window-stable functional exists for this
  comparison at all. Both candidates `MEDIATOR_NOISE_GATE_LEVEL_THREE.md` named
  are now spent. §7 rules out MI integrated over time: the integral's mean
  dilutes with the window and every span falls with it, so only the ratio
  survives. [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) §1 rules out
  steady-state MI for the arm this page's sweep compares: wherever nothing is
  blind the stationary state of the **single-excitation** sector is I/N at every
  seat and every rate, its mutual information on halves of size (N−1)/2 with the
  centre in neither is log₂N − ((N+1)/N)log₂((N+1)/2) bits, and the span is zero
  everywhere. It is window-stable in the way a closed eye is steady. It does not
  cover §7's first row, "none (only the swept seat carries γ)", where the kernel
  is degenerate and the limit depends on the preparation. What that page found
  instead is that the stationary manifold's **dimension** carries the whole
  blindness count, which is window-free and rate-free but is an integer and not a
  mutual information. A third candidate, the spectral gap, is named in an
  untracked design spec under `docs/superpowers/` and untried on this
  comparison.
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
