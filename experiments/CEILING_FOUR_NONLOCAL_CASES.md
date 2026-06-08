# The Four Non-Local Cases: the palindrome ceiling, narrowed to two

**Date:** 2026-06-06 (updated 2026-06-07: the 4 → 2 step, two more turned out local)
**Status:** the routing is constructively verified; the genuine non-locality narrowed from four cases to
two, and that pair is the open analytic frontier.
**Scripts:** [`simulations/ceiling_6to4_verification.py`](../simulations/ceiling_6to4_verification.py)
(the 6 → 4 step), [`simulations/ceiling_4to2_iheavy_local.py`](../simulations/ceiling_4to2_iheavy_local.py)
(the 4 → 2 step)

> **Update, the arc is 6 → 4 → 2.** This writeup keeps its title (the historical four), but the count has
> moved on twice. The first correction found two of an original six were local (XIX+XIY+YIX, YIY+XIY+YIX,
> routed by a continuous-uniform per-site map). The second correction, recorded here, finds two of the four
> are local as well, in exactly the same spirit: the I-heavy pair IXI+IIY+YII and IYI+IIX+XII palindromize
> through a constructive single-site-field route. That leaves the two genuine non-local cases, the Z-middle
> pair XZX+XZY+YZX and YZY+XZY+YZX. So the real ceiling today is two, and the sections below carry both
> corrections in order.

## The correction: six were really four

We had a ceiling, and the ceiling was too high. For a while the certifier's story ended with six k-body
term-sets that no per-site product Q could palindromize, the cases where the mirror is only ever an
entangled, non-local Π. Six "non-local" sets, neatly named XZX+XZY+YZX and its siblings. The trouble is
that the count was read off the wrong test. Two of the six were never given the chance the others got.

So we gave all six the same chance, under one reliable test: the complete uniform-continuous per-site map
family. Concretely, a 16-real-parameter palindrome objective, Q = M^⊗N built from a single 4×4 site map M,
scored by the residual ‖Q L Q⁻¹ −(−L−2σ)‖ relative to its target, and minimized over many restarts. Run
against all six candidate sets, the verdict splits 2 and 4. Two route: their residual falls to ≈8×10⁻¹⁴
with an invertible M, and it stays there at N=3, 4, 5. Four resist: the best the optimizer can manage is
bounded well away from zero (0.576 and 0.814, depending on the case).

The reason the split is trustworthy is that the optimizer is demonstrably capable of finding the routers
when they exist. It found the two. The same search, run on the four, simply cannot drive the residual
down; the floor is the structure pushing back, not the optimizer giving up. At this step the ceiling stood
at four, not six. (That four would later drop again to two, once a second, non-uniform route was tried on
the I-heavy pair; the uniform search here cannot see it. See "The four, on a second look" below.)

## The two that are local: XIX+XIY+YIX, YIY+XIY+YIX

What hid these two from the earlier count is that they route in a way the per-term test could not see. The
per-term routing strategy (Stufe B) asks for a single per-site Q that palindromizes every term on its own.
These two do not satisfy that; no fixed Klein router works term by term. They satisfy something weaker and
quieter: the continuous-sum cancellation {M^⊗k, Σ_terms [T, ·]_k} = 0. The router only has to balance the
whole sum of commutators, not each one, and a continuous map can thread that needle where a per-term map
cannot. This is why they also fail the discrete {P1, P4, M2} routers outright (residual ≈1.15): the
discrete Klein elements act per term, and per term there is nothing clean to grab.

So these are local, with no caveat: a per-site product Q palindromizes them, verified at N=3, 4, and 5.
But they are local in an awkward register. The routers the optimizer lands on are arbitrary continuous
maps; they are not order-4, and they differ from case to case. That makes them poor citizens for a clean
"addable" certifier, which wants a small fixed alphabet of routers, not a bespoke continuous map per
Hamiltonian. The honest reading is that their absence from the earlier list was a coverage gap, not a
discovery about them: the per-term lens was blind to continuous-sum routing. Fill the gap and they fall
out of the non-local class. They were never non-local; they were unseen.

## The four, on a second look: two local, two genuinely not

The original list of four was XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII, IYI+IIX+XII. A closer look splits it
again. Two of the four, the I-heavy pair IXI+IIY+YII and IYI+IIX+XII, turned out LOCAL, by a constructive
route the first pass missed; the next subsection makes that explicit and corrects the obstruction reading
it had picked up. The other two, the Z-middle pair XZX+XZY+YZX and YZY+XZY+YZX, are the genuine non-local
cases: they admit no per-site product Q at all, not in any family we can reliably close. The
uniform-continuous search leaves them bounded away from a palindrome, and the discrete-periodic Klein
routers miss them too. What is interesting is not only that they resist but how they resist; the
obstruction has a clean mechanism, described after the I-heavy correction.

### The two that are local after all: IXI+IIY+YII, IYI+IIX+XII (the 4 → 2 step)

The first pass read the I-heavy pair as non-local, with a "per-site sector overload" obstruction and a
joint floor near 0.47, "a single product map cannot turn both at once." That is true of the UNIFORM map,
the one the 16-parameter search above optimizes, where a single 4×4 site map M is repeated on every site.
But it is not the last word, because the I-heavy pair does not need a uniform map. It needs a SITE-VARYING
per-site product, and that exists in closed form.

Here is why. Summed over the windows of an N-chain, each I-heavy term-set is a sum of weight-1 transverse
fields: the chain Hamiltonian is H = Σ_i (a_i X_i + b_i Y_i), a little transverse field a_i X_i + b_i Y_i
sitting on each site. Single-site Pauli operators on different sites commute, so the Liouvillian splits as
a sum L = Σ_i L_i over commuting single-site pieces, one per site. Now each single-site transverse field is
something we already know how to palindromize: rotate it about the Z axis until it points along X (an
R_z-rotation, which leaves the Z-dephasing untouched because R_z commutes with Z), and it becomes a plain
single-site X-field, whose Liouvillian spectrum {0, −2γ, −γ ± 2i} is a clean palindrome about the centre
−γ. Call the resulting per-site map M_i (it carries its own rotation angle θ_i = atan2(b_i, a_i), so it is
genuinely site-varying, different on each site). The product Q = ⊗_i M_i then palindromizes the whole
chain at once, because the chain is just the commuting sum of pieces each M_i handles. The construction is
explicit, it is N-independent, and it checks to machine precision (residual ~1e-14) at N=4, 5, and 6
([`simulations/ceiling_4to2_iheavy_local.py`](../simulations/ceiling_4to2_iheavy_local.py)). In C# the
`SingleSiteField` strategy carries exactly this certificate.

The 0.47 floor was therefore the floor of the wrong family. It is real for the uniform map (one M repeated
everywhere genuinely cannot turn two sectors at once), but the I-heavy pair was never asking for a uniform
map; a product of distinct single-site crossover maps clears it cleanly. This is the same shape of
correction as the 6 → 4 step: a coverage gap, not a fact about the Hamiltonians. The per-term and uniform
lenses were blind to the site-varying single-site route, and once it is supplied the I-heavy pair falls out
of the non-local class. They were local all along.

One guardrail is worth stating, because it is what makes the certificate sound rather than greedy. The
locality is specific to TRANSVERSE fields. A single-site X or Y field is soft, its spectrum
{0, −2γ, −γ ± 2i} palindromic about −γ; but a single-site Z (longitudinal) field is HARD, its spectrum
{0, 0, −2γ ± 2i} leaving the 0 eigenvalue with no partner about −γ. So the route certifies a sum of
transverse single-site fields and pointedly excludes Z. The I-heavy pair is all-transverse, which is
exactly why it qualifies.

### The two that are not: XZX+XZY+YZX, YZY+XZY+YZX

The Z-middle pair is the genuine non-local ceiling. It admits no per-site product Q at all, not in any
family we can reliably close, and the way it resists has a clean mechanism.

The cleanest way to see the obstruction is the F1 / F1² lens. The canonical residual M = Π L Π⁻¹ + L + 2σ
is anti-Hermitian everywhere it does not vanish. Anti-Hermitian means a rotation: to undo it, a per-site
map has to supply an angle, not a reflection, and the question becomes whether one product map carries
enough of the right angle to cancel it. Splitting M by its behaviour under Π² (the Π²-even sector is real,
the Π²-odd sector complex) is what exposes where the angle runs out.

For the Z-middle pair, XZX+XZY+YZX and YZY+XZY+YZX, the residual M carries an extra rotation that the
routable cases never have. On top of the Π²-odd (complex) rotation that the local cases route, the
Z-middle cases pick up an additional Π²-even (real) rotation. The best per-site map can clear the Π²-even
sector cleanly, but once it is spent on that, the Π²-odd sector still holds an intrinsic floor of about
0.30 that no per-site map can reach. The Z in the middle is not a spectator; it loads a second rotation
into the residual, and a single product map cannot turn both at once.

(The I-heavy pair shows a different, milder pattern under the same UNIFORM lens: each Π²-sector is
individually clearable but not both at once for a single repeated map, a per-site sector overload with a
joint floor near 0.47. That floor is real for the uniform map, but it is not the ceiling it once looked
like; the site-varying single-site-field route above clears the I-heavy pair outright, which is why they
are local and only the Z-middle pair remains here.)

So for the Z-middle pair the conclusion stands: a product of single-site maps does not have the freedom the
loaded second rotation needs. Their palindrome exists only as an entangled, non-local Π.

## The explicit frontier (no hidden bottom)

A word on what "non-local" means here, because it is a bounded claim and we want the boundary visible. For
the two Z-middle cases that remain, non-local means: admits no per-site product Q in the reliably-testable
families. That is two families, both of which we can actually close by optimization, the uniform-continuous
family (the complete 16-parameter one above) and the discrete-periodic Klein routers. Across both, the
Z-middle pair is bounded away from a palindrome. (The I-heavy pair, by contrast, is now closed on the
constructive side: a site-varying single-site-field product palindromizes it outright, so it carries no
frontier asterisk at all.)

There is one family we deliberately do not claim to have closed for the Z-middle pair: the
continuous-periodic maps, period-2 and higher, where M itself varies from site to site. That objective is
32-dimensional and globally non-convex, and an optimizer that fails to find a router there has not proven
one absent; it has only failed to find one. So we name it rather than absorb it. The continuous-periodic
family is the explicit open frontier of the Z-middle non-locality. It is a visible edge, not a hidden trap.
By contrast the local cases carry no such asterisk: the XIX/YIY pair routes in the uniform-continuous
family outright, and the I-heavy pair routes by the constructive single-site-field product, both
definitively local.

Beyond even that frontier sits the harder question we are not attempting here: the analytic all-Q
obstruction, a proof of why the Z-middle floor holds for every Q, of any family, periodic or not. That is
the banked hard problem. We reference it; we do not open it. What this experiment establishes is the
constructive half, told honestly: across the 6 → 4 → 2 arc, four route and are local (two by continuous-sum
cancellation, two by the single-site-field product), and the two Z-middle cases resist over every test we
can trust; the line between "verified" and "open" is drawn exactly where the optimization stops being
reliable.

## Links
- Mechanism + classification: [PROOF_F103 §7.12](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (the single-site transverse-field lemma)
- The certifier + the 2-case Z-middle ceiling: `compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs` (the `SingleSiteField` strategy), `PalindromeSoftCertifierClaim.cs`, `KBodyPalindromeRouting.cs`
- The 4 → 2 verification: [`simulations/ceiling_4to2_iheavy_local.py`](../simulations/ceiling_4to2_iheavy_local.py)
- Related: [SOFTNESS_IS_N_DEPENDENT](SOFTNESS_IS_N_DEPENDENT.md)
