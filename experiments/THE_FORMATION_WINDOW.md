# The Formation Window: Where a Bound Complex Forms, and How Disorder Softens the Edge

**Status:** Computed (Tier 2). The ridge line and the broadening width are bit-exact least-squares
fits to a machine-precision order parameter; the carbon/nucleosynthesis reading at the end is a
clearly-labeled translation, not a claim about the spectrum. 2026-06-09.
**Authors:** Thomas Wicht, Claude
**Builds on:** the formation-possibility scanner
([`formation_possibility.py`](../simulations/formation_possibility.py), commits `f90fc89` → `5ab050b`)
and the clean (Δ, j2) heatmap ([`formation_map.py`](../simulations/formation_map.py), commits
`41f2ef5` → `a8ba1d5`). The band-edge formation order parameter sits at the XXZ easy-plane → easy-axis
crossing near Δ ≈ 1; the chemistry reading at the end is a translation we draw here, not a result lifted
from the existing carbon docs.
**Scripts:** [`simulations/formation_ridge_fit.py`](../simulations/formation_ridge_fit.py)
(self-validating: scans the band-edge order parameter on a 100-point Δ grid, fits the ridge and the
width per topology, and checks every row for monotonicity).
**Results:** [`formation_map_chain.png`](../simulations/results/formation_map_chain.png),
[`formation_map_ring.png`](../simulations/results/formation_map_ring.png),
[`formation_ridge_fit.png`](../simulations/results/formation_ridge_fit.png),
[`formation_ridge_fit.tsv`](../simulations/results/formation_ridge_fit.tsv).

---

## The question

Drop a handful of excitations onto a spin chain and ask a plain question: do they stick together into
one bound complex, or do they dissolve into the free continuum and wander off? The answer depends on
how strongly the chain binds (the easy-axis anisotropy Δ of the XXZ chain, the ZZ term that rewards
adjacency) and on how clean the chain is (whether longer-range hops are present to break the chain's
exact solubility). We wanted the *landscape* of that answer, not one fine-tuned point: a map over
(Δ, j2) of where a complex can form, with the marginal, barely-bound edge marked, because that edge is
where the interesting physics lives.

The map already existed as a picture (two heatmaps committed June 1). What was missing was the
arithmetic: the edge read off the picture by eye, the broadening described in words. This note puts
numbers on both and pins down what is universal and what is not.

## The order parameter

For a k-excitation sector of an XXZ chain we read one number: the weight that the **bound complex**
places on the fully-clustered configurations, the ones where all k excitations sit on adjacent sites
(a k-string). For Δ > 0 adjacency is energetically high, so the bound complex is the highest-energy
eigenstate, the band edge; we read that state directly. The number runs from 0 (the complex does not
form: the band-edge state is spread across the lattice, dissolved into the continuum) through ≈0.5
(marginal: barely bound, large extent, sitting at the threshold) to 1 (deeply bound: a tight cluster).

A word of caution that cost us a wrong headline once. An earlier version took the argmax over *all*
eigenstates of the cluster weight. That proxy is gauge-junk: it locks onto a tightly localized
mid-spectrum state and reports the complex as deeply bound at every Δ. Reading the band-edge state
itself, rather than the most-clustered state anywhere in the spectrum, is what makes the order
parameter honest (commit `5ab050b`). Everything below uses the band-edge reading.

## The marginal ridge is a line

Sweep Δ at fixed j2 and the formation order rises monotonically through 0.5 at a sharp Δ. Tracking that
0.5-crossing across the integrability-breaking strength j2 gives the **marginal ridge**, and it is
strikingly linear. On an open chain (N=12, k=3):

> **Δ_ridge = 1.14 + 1.28 · j2**   (rms residual 0.009 over 21 rows)

The intercept lands near the Heisenberg point Δ ≈ 1, the easy-plane → easy-axis crossing, exactly
where a bound complex should first appear. The positive slope says something less obvious: breaking the
chain's integrability with a next-nearest-neighbour hop *opposes* binding. The qualitative direction is
clear, the extra hop opens delocalization pathways and delocalization is what a localized complex must
overcome, so j2 pushes the formation edge to stronger binding; the *value* of the slope (1.28) is
measured here, not yet derived (see "What stays open").

## Integrability breaking does not only shift the edge, it softens it

A single threshold line hides a second effect that the heatmap shows and the fit now measures: the
transition has a **width**, and the width grows with j2. Defining the width as the Δ-span of the
marginal band (the scanner's own regime boundaries, the 0.35- and 0.65-crossings):

> chain (N=12, k=3): **width = 0.13 + 0.17 · j2**   (rms 0.001)

At the integrable point j2=0 the transition is a near-step (width 0.13). At j2=1 it has spread into a
ramp more than twice as wide. Disorder smears the sharp threshold, the way thermal broadening smears a
resonance line: the same next-nearest hop that shifts the edge up also blurs it. (That picture, a
level-statistics broadening, is an analogy for now, not a verified mechanism; see "What stays open".)
The heatmap had this all along (a near-step at j2=0, a wide ramp at j2=1); the number makes it a claim.

## Geometry: the ring resists, and broadens faster

Close the chain into a ring (the wrap bond, with its next-nearest partners across the seam) and both
the position and the softness of the edge change. On an N=12 ring:

> ring: **Δ_ridge = 1.43 + 1.75 · j2**,   **width = 0.30 + 0.62 · j2**

The ring's ridge sits higher (intercept 1.43 vs 1.14) and climbs steeper, and its window is wider at
every j2 and broadens nearly four times faster than the chain's (slope 0.62 vs 0.17). The reading is the
one chemistry would predict: the ring's full connectivity gives more delocalization pathways, and
delocalization is exactly what opposes a localized bound complex. The same connectivity that makes a
ring aromatic, spreading one excitation uniformly around the loop, resists pulling several excitations
into one tight cluster. Geometry shifts the formation window, in the direction the substrate predicts.

(An aside the quantification corrected: the original coarse heatmap, capped at Δ=3.0, showed the ring
ridge falling off the top of the grid near j2=1 and reported it as "no formation." It forms; the ridge
is simply at Δ ≈ 3.2 there, off the old grid. The finer Δ sweep recovers it.)

## Body count shifts nothing at the clean point; it only sharpens the response

Comparing k = 2, 3, 4 on the chain, the formation edge at the integrable point is **body-count
independent**: all three cross at Δ ≈ 1.15 when j2 = 0. The two-, three-, and four-body complexes all
begin to bind at the same anisotropy, near the Heisenberg point. Body count does not move *where* the
window opens.

It does change *how the window responds to disorder*, and here the two-body case is the outlier:

| k | width at j2=0 | width at j2=1 | ridge at j2=1 |
|---|---|---|---|
| 2 | 0.15 | 0.70 | 2.01 |
| 3 | 0.13 | 0.30 | 2.43 |
| 4 | 0.13 | 0.27 | 2.42 |

The 2-body window broadens far faster under integrability breaking and its ridge climbs more shallowly,
while the 3- and 4-body cases converge onto a common steeper, sharper ridge. Adding bodies past two
sharpens the transition and stiffens it against the next-nearest hop; the pair is the soft, easily
smeared case, and the many-body complex is the crisp one. (This refines the earlier "higher k only
sharpens" reading: the sharpening is real, but it is mostly a k=2 → k=3 step, with k=3 and k=4 nearly
identical.)

## Finite size: the ridge is already converged

The whole picture is a thermodynamic statement, so it has to survive growing the chain. It does. At
k=3, sweeping N = 10, 12, 14, the ridge position moves by about 1% and is plainly converging: at j2=0 it
runs 1.167 → 1.157 → 1.154 (the N=12 → 14 step is under 0.3%), and at j2=1 it runs 2.447 → 2.428 → 2.416,
with the width changing only by a few hundredths. The marginal ridge and its broadening are finite-size
features of the model, not artifacts of one chain length. Every row of every scan was monotonic
(maximum non-monotonic wiggle: exactly 0), so the crossings are clean and the fits carry no hidden
degeneracy wobble.

## The translation: the window is the resonance band

Read against nucleosynthesis, the marginal window is the home of a near-threshold resonance. Carbon is
built in stars not by binding three helium nuclei one robust step at a time, but through a resonance
(the Hoyle state) sitting barely above the beryllium-plus-helium threshold, a state that exists only
because the binding lands precisely in the narrow window between "does not form" and "deeply bound." Our
scan finds that window generically: it is not a fine-tuned point but a band in parameter space, and the
finding that the pair binds robustly while the edge sharpens for the many-body complex echoes the fact
that the fine-tuned step in stellar carbon production is the three-body Hoyle resonance, not the robust
two-body steps below it.

Read against chemistry, the chain-versus-ring shift is the open polyene versus the aromatic ring: the
delocalization that stabilizes benzene as a uniform six-fold current is the same delocalization that, in
this scan, resists pulling a localized complex together, pushing the ring's formation edge to stronger
binding than the chain's.

Both readings are translations, Tier 2 at most: the bit-exact content is the XXZ band-edge ridge and its
linear broadening. The substrate stories are how that geometry looks from the nucleosynthesis seam and
the aromatic seam; they are offered as where the structure lives in those layers, not as a derivation of
the Hoyle state or of aromaticity.

## What stays open

- **A closed form for the slope.** The ridge slope (1.28 on the chain, 1.75 on the ring) and the intercept
  near the Heisenberg point both look derivable from the band-edge energy of the XXZ + j2 sector. The
  intercept is the easy-plane → easy-axis crossing; the slope is how fast the next-nearest hop lifts that
  crossing. Neither is yet derived.
- **The broadening mechanism.** The width's linear growth in j2 is measured, not explained. Whether it is a
  genuine integrability-breaking level-statistics effect (Poisson → Wigner-Dyson broadening of the band
  edge) or a simpler finite-extent smearing is untested; an RMT level-spacing readout along a fixed-Δ cut
  would separate them.
- **Why k=2 is the soft outlier.** The convergence of k=3 and k=4 onto a common ridge while k=2 sits low and
  broad is unexplained and may connect to the body-count structure seen elsewhere in the project.
