# The Four Non-Local Cases: the genuine palindrome ceiling

**Date:** 2026-06-06
**Status:** the routing is constructively verified; the full non-locality is the open analytic frontier.
**Script:** [`simulations/ceiling_6to4_verification.py`](../simulations/ceiling_6to4_verification.py)

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
down; the floor is the structure pushing back, not the optimizer giving up. The ceiling is real, but its
height is four, not six.

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

## The four that are not: XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII, IYI+IIX+XII

The remaining four admit no per-site product Q at all, not in any family we can reliably close. The
uniform-continuous search leaves them bounded away from a palindrome, and the discrete-periodic Klein
routers miss them too. What is interesting is not only that they resist but how they resist, and there
are two different mechanisms cornering them.

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

For the I-heavy pair, IXI+IIY+YII and IYI+IIX+XII, the mechanism is different: a per-site sector overload.
Here each Π²-sector is individually clearable, the map can zero either one alone, but not both
simultaneously. Aim the map at one sector and the other lights up; the two demands are incompatible for a
product map, and the joint floor sits at about 0.47. Where Z-middle has one sector with an irreducible
residue, I-heavy has two sectors that each go to zero but never together.

Both mechanisms say the same thing in the end: a product of single-site maps does not have the freedom
these four need. The palindrome exists only as an entangled, non-local Π.

## The explicit frontier (no hidden bottom)

A word on what "non-local" means here, because it is a bounded claim and we want the boundary visible. For
the four, non-local means: admits no per-site product Q in the reliably-testable families. That is two
families, both of which we can actually close by optimization, the uniform-continuous family (the complete
16-parameter one above) and the discrete-periodic Klein routers. Across both, the four are bounded away
from a palindrome.

There is one family we deliberately do not claim to have closed: the continuous-periodic maps, period-2
and higher, where M itself varies from site to site. That objective is 32-dimensional and globally
non-convex, and an optimizer that fails to find a router there has not proven one absent; it has only
failed to find one. So we name it rather than absorb it. The continuous-periodic family is the explicit
open frontier of the four's non-locality. It is a visible edge, not a hidden trap. By contrast the two
local cases carry no such asterisk: they route in the uniform-continuous family outright, definitively
local.

Beyond even that frontier sits the harder question we are not attempting here: the analytic all-Q
obstruction, a proof of why the Z-middle floor and the I-heavy overload hold for every Q, of any family,
periodic or not. That is the banked hard problem. We reference it; we do not open it. What this experiment
establishes is the constructive half, told honestly: two route and are local, four resist over every test
we can trust, and the line between "verified" and "open" is drawn exactly where the optimization stops
being reliable.

## Links
- Mechanism + classification: [PROOF_F103 §7.12](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
- The certifier + the 4-case ceiling: `compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs`, `PalindromeSoftCertifierClaim.cs`, `KBodyPalindromeRouting.cs`
- Related: [SOFTNESS_IS_N_DEPENDENT](SOFTNESS_IS_N_DEPENDENT.md)
