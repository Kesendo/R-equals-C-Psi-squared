# The N=6 Route-B A₂ Locus Atlas

**Status:** Executed finite-N atlas. The exact schema-3 carrier supplies 266 isolated direct-*t* A₂ loci; the C# full-sector classifier consumes every locus once and reports 266 diabolic crossings, 0 defective points and 0 unresolved loci. This is an N=6 result, not an all-*N* theorem.

**Date:** 2026-09-08

**Authors:** Thomas Wicht and Codex

**Physical scope and coordinates:** Open uniform nearest-neighbour XY chain
(Δ=0, zero field) with uniform local Z-dephasing γ=1, restricted to the
SE-ket/DE-bra coherence block. The carrier uses the same `qCSharp=J/γ`
normalization as the earlier N=4 work, but stores the N=6 search plane as the
exact quarter-turn *t*=i·`qCSharp`; its spectral coordinate is Λ=2λ. Thus a
reported isolation margin in physical λ units is the carrier's Λ margin
divided by two. Complex *t* (equivalently complex *q*) is analytic
continuation, not a physical coupling.

## A map after the search

The Route-B calculation did not merely return a count. It left 266 small
isolation boxes in the complex *t* plane: 133 pairwise disjoint boxes in each
R-parity, each one carrying an algebraic double root and two exact involutions:
complex conjugation within one R-parity and
parity transport *t* → −*t*. The atlas is the moment those boxes stop being a
list and become a place we can revisit.

![The 266 N=6 Route-B A2 loci in the complex t plane.](../visualizations/route_b_a2_n6_constellation.png)

*Every point is drawn at the midpoint of its exact rational isolation box,
converted to a floating-point display coordinate. The marker is therefore a
representative of a certified box, not the algebraic root itself. Cyan is
R-even, violet R-odd; diamonds are the 118 real-*t* `HermitianAxis` readings,
circles the 148 off-axis `EpCharacterStable` readings. The gold outline records
the executed verdict Diabolic. Magenta is reserved for Defective, of which the
atlas contains none.*

The colour and shape carry different questions. Parity says which residual
sector owns the locus. Shape says how its local character was established.
Neither is allowed to stand in for the other, and the absence of magenta is a
measured outcome rather than an empty design choice.

## Ninety-six rooms, not 266 unrelated points

The two partner maps close the inventory into 96 orbits. The 59 real-axis
orbits have size two: conjugation fixes each point and parity transport joins
the even and odd loci at opposite *t*. The remaining 37 orbits have size four:
conjugation supplies the upper/lower pair and parity transport supplies its
opposite-parity image.

![The conjugation and parity-partner orbits of the N=6 A2 loci.](../visualizations/route_b_a2_n6_orbit_map.png)

*Solid grey joins the carrier's conjugation partners; dotted gold joins its
parity partners. These are the only edges present in the schema-3 artifact.
The figure deliberately draws no “crossfold” relation: adding one from visual
symmetry alone would turn the atlas into a second, unsupported dataset.*

This is the useful compression: 266 certified local events become 96 rooms
without discarding a locus. It is also a practical search index. A later
perturbation can be asked whether it moves an entire room coherently, splits a
conjugate pair, or breaks the parity transport. The map does not answer those
questions in advance; it makes their controls visible.

## How we know the crossings are silent

Algebraic multiplicity two comes from the exact A₂ carrier. Geometric
multiplicity is a separate question and passes through the full
45-dimensional R-parity block at the exported seed.

![Executed isolation and three-contour character evidence.](../visualizations/route_b_a2_n6_evidence_profile.png)

On the real-*t* axis the full block is Hermitian to an executed residual of
exactly 0.0 in this run, so every algebraic double root there is semisimple.
Off the axis, `EpCharacter` reads the isolated two-dimensional restriction at
three radii. Every radius must enclose exactly the same pair and all three must
agree on character. The diabolic gate requires relative departure below
10⁻⁶; the defective gate would require geometric multiplicity one and
departure above 5·10⁻². The smallest middle-contour isolation margin is
0.06589058580248003 in physical λ units at γ=1. No exact-rank fallback was
executed, and the source artifact's `exactRankCertificates` array is empty.
For 22 off-axis loci the recorded maximum relative departure is exactly 0.0;
the logarithmic evidence plot places those markers at a 10⁻¹⁸ display floor,
which is not a measured departure.

| Inventory reading | R-even | R-odd | total |
|---|---:|---:|---:|
| certified direct-*t* A₂ loci | 133 | 133 | 266 |
| `HermitianAxis` | 59 | 59 | 118 |
| `EpCharacterStable` | 74 | 74 | 148 |
| `ExactRankExecuted` | 0 | 0 | 0 |
| Diabolic, alg=geo=2 | 133 | 133 | 266 |
| Defective / unresolved | 0 | 0 | 0 |

## Provenance and reproduction

The source of geometry is
[`simulations/results/route_b_a2_n6.json`](../simulations/results/route_b_a2_n6.json),
SHA-256
`c4b19015a750e56a55b0c8fdba3d05534c0c133bc5355a74eb89b77020532d47`.
The broader double-root interpretation and the boundary to other N live in
[`THE_DOUBLE_ROOT.md`](../docs/THE_DOUBLE_ROOT.md).

The atlas manifest binds that digest to every accepted character reading and
orbit:

```powershell
dotnet run --project compute/RCPsiSquared.Cli -c Release -- route-b-n6-atlas --out simulations/results/route_b_a2_n6_atlas.json
python simulations/route_b_a2_n6_atlas.py --manifest simulations/results/route_b_a2_n6_atlas.json --out-dir visualizations
```

The manifest is not an independent measurement. It is a deterministic view of
the exact inventory after the existing reconciliation gate has accepted all
266 executed readings. The renderer refuses changed source bytes, count drift,
an incomplete orbit partition, a non-diabolic verdict, or fewer than three
contour readings.

## The earlier flashlight

The visual ancestor is the N=4
[`f89_octic_branch_locus.png`](../visualizations/f89_octic_branch_locus.png).
That image searched the complex *q* plane with a min-gap heatmap and showed a
mixed branch locus: 20 exceptional points and four diabolic points. This atlas
does something narrower and later. It maps only the N=6 direct-*t* A₂ layer,
with no heatmap and no A₁ exceptional points. The two figures share a visual
grammar: magenta for defective seams, gold for silent crossings, and the same
`qCSharp` normalization; the displayed N=6 *t* plane is exactly the N=4 *q*
plane quarter-turned by multiplication with i. They do not share N, sector,
dataset, or claim.

## Boundary of the atlas

The map establishes the finite N=6 inventory, its two exact partner
involutions, and the executed local character at the committed seeds. It does
not turn box midpoints into exact roots, prove stability under Δ, disorder or
topology changes, infer a crossfold edge, or settle the F₅₃ N=7 doubled layer.

What it gives us is a better place from which to ask those questions. Instead
of hunting “another point,” we can now choose a room, its mirror, its parity
partner and the evidence route that made it visible, and design a perturbation
that is capable of breaking exactly one of those relations.
