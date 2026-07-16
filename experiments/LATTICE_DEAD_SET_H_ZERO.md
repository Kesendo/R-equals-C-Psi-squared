# The muteness axis: the dead-set law is generic in h, and at h = 0 it refines

**Date:** 2026-07-16 (the same day's third act: the h thread found the zeros, the dead-set rule closed them, this maps how the rule moves along the axis it was found on)
**Status:** verified from below. Gate: [`simulations/lattice_dead_set_h_zero.py`](../simulations/lattice_dead_set_h_zero.py) (29 checks, all PASS, exit 0). Refinement of [F132](../docs/ANALYTICAL_FORMULAS.md); no new F number, the h = 0 rule is F132's own structure at the special point.
**Where it came from:** [LATTICE_DEAD_SET_RULE](LATTICE_DEAD_SET_RULE.md) left one honest asymmetry: necessity (forbidden ⟹ identically zero) is derived, sufficiency (allowed ⟹ alive) is the gated observation, and every census behind that observation scanned exactly ONE fixed generic h direction per config. The handover named the follow-up "the muteness-axis scout": treat the mute readouts as an object and walk the h axis. Two questions fell out. Is sufficiency generic, or did the single sample get lucky? And what happens at the non-generic points of the axis?

## The two answers

**Sufficiency is generic.** Across five configs (N = 3..5, population and real
coherence, seeds s = 1 and s = 3) and six seeded random h directions and scales
each, the F132 law predicts the alive set exactly at every sample. The
separation is not marginal: every allowed readout sits at ≥ 5·10⁻³, every
forbidden one at machine zero (~2·10⁻¹⁶). The dead set is not "small", it is
exactly zero, at every generic h we threw at it. Scope: this gate runs the
one-sided watching profile throughout; the parent gate swept full non-uniform
watching, and the necessity arguments below are watching-profile-independent
(the dephasing mask is diagonal on Pauli strings for ANY profile), but the
h-axis genericity observation itself is gated at one profile per config.

**At h = 0 the law refines, and the refinement collapses again.** Staggered and
single-site fields are already generic (the law unchanged). But at h = 0 and at
uniform h the alive set shrinks by exact, reproducible counts, and the h = 0
shrinkage has a closed one-line rule of its own (below). The counts:

| Config | law alive | extra dead at h = 0 | extra dead at uniform h |
|---|---|---|---|
| N=3 population (s=1) | 15/63 | 6 | 6 |
| N=3 real coh (s=1) | 35/63 | 16 | 6 |
| N=4 population (s=1) | 55/255 | 24 | 24 |
| N=4 real coh (s=1) | 71/255 | 32 | 0 |
| N=5 real coh (s=3) | 507/1023 | 256 | 130 |

## The h = 0 mechanism: the Majorana graph disconnects

Under the left Jordan-Wigner map (the same one F132's degree d lives on), the
XX+YY hopping is quadratic with Majorana edges (2l+1, 2l+2) and (2l, 2l+3),
and the field Z_l is quadratic with the on-site edge (2l, 2l+1). Walk the
hopping edges alone and the 2N Majorana indices fall into two components that
never meet:

    E = {k : k mod 4 ∈ {0, 3}},   O = {k : k mod 4 ∈ {1, 2}},

and the field edges are the ONLY bridge between them. So at h = 0 the flow
conserves not just the total degree d but the BI-degree (d_E, d_O), the
Majorana count per component (the adjoint dissipator is diagonal on Pauli
strings and rides along, as always in this thread). The gate pins the
conservation at the operator level: [H₀, O] expands into strings of O's own
(d_E, d_O) sector with leak exactly 0.0, for all 4^N strings at N = 3 and
sampled strings at N = 4, 5.

**The collapsed h = 0 rule.** The prep populates few bi-sectors: the population
pair expands in even-n_Z diagonal strings, and every Z contributes exactly one
Majorana to each component, so the population sectors are (k, k) with k even;
the cross coherence expands in all-XY strings of total degree N, spread over
mixed sectors (d_E, N − d_E). A readout can only ever show what its own
bi-sector holds, and everything else about F132 becomes automatic per populated
sector (in a population sector (k, k), d = 2k ≡ 0 mod 4 and the sign conditions
below hold identically). One kinematic line remains:

    alive at h = 0  ⟺  (d_E, d_O)(O) is populated in ρ(0)  ∧  O is K-readable.

Exact on all six gated configs (the five above plus N = 5 s = 1 coherence),
including the ten strings the sign rule alone misses (next section). As
everywhere in this thread: the necessity direction is the conservation
argument; sufficiency is the gated observation.

## The tower repeats one level down

F132's discovery shape was: an antiunitary symmetry V (the composed mirror
W_g = U_g·X^N with conjugation) whose kill sign turned out to be the mod-4
shadow of the exact fermion degree. At h = 0 the SAME tower appears one level
finer, and it appeared in the same order:

**The bare chiral-conj symmetry.** Without the field, the sublattice gauge U_g
flips the hopping alone (every bond has exactly one site in g), so no X^N is
needed:

    V'_g(ρ) = U_g · conj(ρ) · U_g

is an exact flow symmetry at h = 0 (pinned dynamically: V' commutes with the
propagated flow, residual 0.0). Its kill sign on a Pauli string is
ε'_g = (−1)^(n_Y + xy_g), and the per-channel sign rule (population needs
ε'_g = +1 for both gauges, the coherence channel needs ε'_g to match the sign
U_g puts on the coherence pair) reproduces the h = 0 alive set exactly at
N = 3 and N = 4, both preps.

**The crack, in the same place as last time.** At N = 5 with the popcount-2
seed the sign rule over-predicts by exactly ten strings (IXXXX, IYYYY, XXIYY,
XXXXI, XXXZY, XZYYY, YYIXX, YYYYI, YYYZX, YZXXX; it never under-predicts).
F132's F layer first bit at N = 5 with the popcount-2 seed; the bi-degree
first bites at N = 5 with the popcount-2 seed. All ten are strings whose
(d_E, d_O) sector is unpopulated, e.g. IXXXX carries {3, 4, 7, 8}, all four
Majoranas in component E, so (4, 0), while the population prep only ever fills
(k, k).

**The shadow identity, one level down.** The per-index sign of V'_g is uniform
on each component (for the even-sites gauge: −1 on E, +1 on O; mirrored for
the odd gauge), which gives

    ε'_even = (−1)^(d(d−1)/2 + d_E),   ε'_odd = (−1)^(d(d−1)/2 + d_O),

checked against the letter formula for all 4^N strings × both gauges,
N = 2..6: 10912 pairs, zero mismatches. Compare F132's
ε_even = (−1)^(d(d+1)/2), ε_odd = (−1)^(d(d−1)/2): the composed mirror's sign
is a function of the total degree, the bare mirror's sign is the same phase
plus one component's count. V' is the mod-4 shadow of the bi-degree exactly as
V is the mod-4 shadow of the degree, and the two shadows differ by one factor:
ε_g/ε'_g = (−1)^(d_O) for BOTH gauges, because X^N itself acts per Majorana
index with the component sign (+1 on E, −1 on O; from a_2l ↦ (−1)^l a_2l,
a_2l+1 ↦ (−1)^(l+1) a_2l+1). Composing with X^N is what fuses the fine shadow
into the total-degree one.

## Uniform h: half-way back to generic

The uniform field is a function of total popcount, so it commutes with the
hopping and acts as a pure phase per popcount block. On the population channel
(diagonal blocks) the phases cancel: the dynamics are IDENTICAL to h = 0, and
the h = 0 deaths persist. On the coherence blocks the phase rotates
(e^(2i·Δp·h·t)), turning the real coherence complex and reviving every
sign-killed coherence readout. The gate pins the union form per config:

    alive(uniform h) = alive(h = 0) ∪ {coherence on ∧ K_coh ∧ d = N},

exact on all five configs (e.g. N = 5 s = 3: 251 ∪ 252 with overlap 126 gives
377). This is why the table above shows the population configs unchanged
between h = 0 and uniform h, and the N = 4 coherence config fully generic at
uniform h: there the whole h = 0 shrinkage sat in the coherence channel.

## What this adds to F132

- The sufficiency face, F132's honest Tier-2 residue, now rests on thirty
  random h samples with machine-zero dead sets instead of one fixed direction,
  for the five configs gated here (the scope note above). Still an
  observation, but a generic one.
- Within the swept family of field profiles (random, staggered, single-site,
  uniform, zero), exactly two members are non-generic (h = 0 and uniform h),
  each with a closed, gated rule. This is not an exhaustive scan of the
  continuous axis; it is the swept family's verdict.
- The mod-4-shadow-of-a-conserved-degree pattern is not a one-off: it is a
  tower that repeats at the finer grading the moment the parameter point turns
  the finer grading on. Leaving h = 0 fuses the two Majorana components, the
  bi-degree coarsens to the degree, the bare mirror dies and only the composed
  one survives: the muteness axis is the fusion knob.

## Files

- Gate: [`simulations/lattice_dead_set_h_zero.py`](../simulations/lattice_dead_set_h_zero.py)
  (T1 generic sufficiency, T2 the structured points, T3 the bare chiral-conj
  symmetry incl. the ten-string over-prediction pin, T4 bi-degree conservation
  and the collapsed h = 0 rule, T5 the shadow identity, T6 the uniform-h union
  rule).
- The law this refines: [LATTICE_DEAD_SET_RULE](LATTICE_DEAD_SET_RULE.md)
  (F132; the typed claim `DeadSetLawClaim` carries the law itself, this
  refinement is doc + gate level).
- The setup and the axis: [LATTICE_H_THREAD](LATTICE_H_THREAD.md) (the h scan,
  X^N as the third F131 mirror).
