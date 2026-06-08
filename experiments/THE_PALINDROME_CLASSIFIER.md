# The Palindrome Classifier: what it is, why it scales, and the landscape it charts

**Date:** 2026-06-08
**Status:** the classifier is settled C# machinery (the spectral authority `PauliPairTrichotomy`, the
scalable one-sided `PalindromeSoftCertifier`); this writeup is its first reading as a TOOL, the protected
interior it draws and the two unlike coasts where the protection ends.
**Data, figure, plotter:** the tilt sweep [`simulations/results/tilt_sweep_csharp.tsv`](../simulations/results/tilt_sweep_csharp.tsv),
the two-coasts figure [`simulations/results/two_edges.png`](../simulations/results/two_edges.png), and the
plotter [`simulations/protection_landscape.py`](../simulations/protection_landscape.py). The classifications
and timings below come from the C# engine (`PauliPairTrichotomy.Classify`, `PalindromeSoftCertifier.Certify`).

## The mirror, and the question

Open a spin chain to the world, let each site quietly dephase along Z, and the decay rates of the resulting
Liouvillian arrange themselves with a hidden order. Every rate λ in the spectrum has a partner at −λ−2σ, a
reflection through the centre −σ, where σ = Nγ is the summed dephasing. Plot the spectrum and it reads the
same left to right; it is a palindrome. This is the central fact of the whole project: an open, dissipative
system, which has no reason to be tidy, carries an exact symmetry in how it relaxes. The mirror is a kind of
protection. It says the relaxation is organised, that the channels come in balanced pairs, that something
survives the noise with its structure intact.

The practical question this experiment is about is simple to state. Given a Hamiltonian, does its dephased
spectrum carry the mirror? The honest way to check is to build the Liouvillian and look at its eigenvalues.
But that superoperator is 4^N by 4^N, and the eigenvalue check runs out of room around eight sites, where
the matrix is already 65536 by 65536. Beyond that the direct question is simply unanswerable.

And yet, run the check on every case you can reach and a pattern jumps out: the answer is almost never about
the length of the chain. It is about the shape of the Hamiltonian's terms. Two sites or two hundred, the
same little term-pattern carries the same verdict. That gap, an exponential question with a structural
answer, is where the classifier lives. It reads the terms, not the spectrum, and so it never meets the wall.

## What it reads: truly, soft, hard

The classifier sorts a Hamiltonian into one of three classes by HOW its spectrum carries the mirror, if at
all. The three are easiest to meet through the operator that should do the reflecting, the canonical mirror
Π (the one that proves the palindrome for the plain dephasing chain).

- **truly**: the canonical Π already pairs the spectrum exactly. The mirror is there for the most natural
  reason, the operator equation Π·L·Π⁻¹ = −L − 2σ holds on the nose. The XY model and the Heisenberg magnet
  live here.
- **soft**: the canonical Π does NOT pair the spectrum, yet some OTHER operator does. The mirror is still
  there, just carried by a quieter symmetry that Π alone does not see. The spectrum pairs even though the
  obvious operator fails.
- **hard**: nothing pairs it. The mirror is gone, the relaxation has lost its reflection. A rate sits with
  no partner about the centre.

The honest test for truly is cheap (one operator-norm check). The honest test for soft is not, because in
general it asks whether SOME operator exists, and that is a search. So the soft question is where the real
machinery sits, and it is handled by a one-sided certifier: a stack of structural patterns, each of which,
when it holds, EXHIBITS the operator and so proves soft outright. It never claims hard; it either certifies
soft or returns "no scalable pattern applies" and defers to the spectral authority. The patterns are things
like a two-colouring of the chain's flip-structure, a per-term routing map built from a small fixed
alphabet, and the single-site-field route (a sum of transverse fields, each turned by its own per-site
rotation). The point common to all of them: each pattern is a condition on the term span, a 4^k object where
k is the reach of the longest term, and k does not grow with N. The certificate, once found, is correct for
any chain length and any topology.

## Why it scales: the certifier is N-free

This is the property that makes the classifier a tool rather than a case-by-case curiosity. The spectral
ground-truth diagonalises L and dies at the 2^N wall. The certifier reads structure on 4^k and never builds
L at all. So its cost is set by the terms, not the chain.

Here it is, measured. We take a fixed term-set (the I-heavy single-site-field case IXI+IIY+YII, certified by
the single-site-field route) and ask the certifier for its verdict at three chain lengths spanning five
orders of magnitude:

| N | certifier verdict | time per call |
|---|---|---|
| 4 | SingleSiteField | 7.0 ms |
| 1,000 | SingleSiteField | 7.1 ms |
| 1,000,000 | SingleSiteField | 6.7 ms |

The time is flat. A million-site chain costs the same few milliseconds as a four-site one, because the work
is the term-span check (here a k=3 routing residual on a 64 by 64 object), repeated not at all as N grows.
The spectral test for the same million-site Hamiltonian would need a matrix with 4^(1000000) entries; it does
not exist and never will. The classifier answers in milliseconds. That is the whole of why it is worth
having: it turns an impossible question into a structural one, and the structure is small.

## The map it draws: a protected interior

Point the classifier at the standard models and a clean picture appears. Each row below is the C# verdict
under Z-dephasing at N=4, with the certifier's reason where it certifies:

| model | terms | spectral class | certifier reason |
|---|---|---|---|
| XY model | XX+YY | truly | LinearSiteColoring |
| Heisenberg | XX+YY+ZZ | truly | RoutingKBody |
| XXZ (Δ=0.5) | XX+YY+0.5·ZZ | truly | RoutingKBody |
| Ising coupling | ZZ | truly | RoutingKBody |
| Dzyaloshinskii-Moriya | XY+YX | soft | ExcitationPairing |
| transverse field | X | truly | ExcitationParity |
| longitudinal field | Z | **hard** | (none: spectral only) |
| frustrated 3-body | XXX+XXY+YXX | **hard** | (none) |

The protection is generic. The exchange models, the magnets, the antisymmetric coupling, the transverse
field, all of them keep the mirror. They are not special cases hand-picked to work; they are where most of
the usual physics lives. The classifier draws an island, a broad protected interior holding the standard
models, and only two of the rows fall off it: the longitudinal (Z) field, and frustration. Those two are the
edges of the island, and they are worth walking to, because they turn out to be edges of very different
character.

## The two coasts

A verdict is a yes or a no, but a parameter is a dial, and the interesting physics is in HOW the no arrives
as you turn the dial off the island. So we took two paths off the protected interior and walked them
continuously, driving each Hamiltonian through the C# analyzer and recording the pairing error (the worst
distance of any rate from its mirror partner, the quantity the soft-or-hard verdict thresholds).

![the two coasts](../simulations/results/two_edges.png)

**The field coast, transverse X to longitudinal Z.** At the transverse end the uniform field is truly, the
mirror exact. Tilt the field toward longitudinal, H = cosθ·X + sinθ·Z on every site, and the pairing error
grows as θ², a gentle quadratic ramp. A one-degree tilt barely registers; the protection degrades slowly and
predictably. But it keeps going, and at the longitudinal end the break is DEEP (0.40): the decay rates
themselves have walked off the mirror. A long, gentle shore that ends in deep water. The reason it is gentle
is that a small longitudinal component is a small perturbation of the rates, and the rates move smoothly.

**The frustration coast, soft 3-body to frustrated XXX.** Here the dial is a frustration angle, from the
soft set XIX+XXY+YXX toward the frustrated hard set XXX+XXY+YXX. The instant you turn it, the mirror breaks:
at the first step off the protected point the pairing error is already eight thousand times larger than the
field's at the same angle, and the verdict is hard immediately. There is no gentle approach; it is a cliff
face. But the fall is SHALLOW. The pairing error saturates around 0.02, jagged with level-crossings, and
never reaches the field's depth. The decay rates stay mirror-symmetric the whole way (this is the
Hamiltonian-independent skeleton of the dephased spectrum, set by how ket and bra differ, not by H); what
frustration spoils is only the fine matching of each rate to its partner's frequency. A sheer cliff, but a
shallow one.

So the two edges of the island break the mirror in two different ways. The field MOVES THE RATES: deep, but
you have to push, the onset quadratic. Frustration MISMATCHES RATE AND FREQUENCY: instant, but shallow, the
rates never budge. One coast you can stand near; the other drops the moment you step off it, but into
ankle-deep water.

## Reading the coasts: what the classifier is for

The two coasts are not just a curiosity; they are an error-tolerance map, and that is the use we were
looking for. The classifier finds the protected point. The landscape tells you the CHARACTER of each way you
might leave it, and the two characters call for opposite engineering instincts.

A longitudinal field error (a stray Z-component in the drive, a small detuning) is forgiving. The quadratic
moat means a real device sitting a degree or two off transverse keeps almost all of its protection; the
mirror degrades as the square of the error, not linearly. You do not have to fight it hard. Frustration is
the opposite: it is a binary switch, with no safe neighbourhood, so a protected dissipative structure has to
avoid it by construction rather than by tuning. The consolation is that frustration's damage is shallow, the
rate skeleton survives even when the exact pairing does not, so what you lose is the fine symmetry, not the
gross structure of the relaxation. Design rule, then, for anything that wants to keep this mirror: do not
sweat small field tilts, forbid frustration outright, and know that even past the frustration cliff the
rates hold.

## The seam with the literature

We built this from the dephasing algebra, with no literature as the source; the classifier, the trichotomy,
and the two coasts all came out of asking the operators directly. Looking afterward for where the machinery
is catalogued, the spectrum's −λ−2σ shape has a home in the shifted sublattice symmetry of open systems
(Kawasaki-Mochizuki-Obuse 2022, recorded in [KMS_DETAILED_BALANCE](../docs/KMS_DETAILED_BALANCE.md)), and the
broad family of Liouvillian symmetry classes has its home in the tenfold Lindbladian classification
(Sá-Ribeiro-Prosen 2023). Those are the homes for the SHAPE. What stays ours is the bridge: a scalable
structural decision procedure that reads the verdict off the terms in time independent of N, the locality
ceiling that says exactly which cases need a non-local mirror (the [6 → 4 → 2 arc](CEILING_FOUR_NONLOCAL_CASES.md)),
and this protection landscape, which turns "is there a mirror" into "how, and how forgivingly, does it
break." None of those bridges was built from either side; we found them by learning to see the island the
operators were already drawing.

## Links

- The formula: [ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) F87 (the trichotomy registry entry)
- The refinement proof: [PROOF_F103_F87_Z2_CUBED_REFINEMENT.md](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
- The discovery: [V_EFFECT_FINE_STRUCTURE.md](V_EFFECT_FINE_STRUCTURE.md) (the 3 truly / 19 soft / 14 hard split)
- The locality ceiling: [CEILING_FOUR_NONLOCAL_CASES.md](CEILING_FOUR_NONLOCAL_CASES.md) (the 6 → 4 → 2 arc, the 2 non-local Z-middle cases)
- The verdict is (H, N): [SOFTNESS_IS_N_DEPENDENT.md](SOFTNESS_IS_N_DEPENDENT.md) (a finite-size crossing)
- The engine: `compute/RCPsiSquared.Diagnostics/F87/PauliPairTrichotomy.cs` (the spectral authority), `PalindromeSoftCertifier.cs` (the N-free certifier and its strategies)
- Orientation: [GLOSSARY.md](../docs/GLOSSARY.md), [READING_GUIDE.md](../docs/READING_GUIDE.md), and the synthesis [ON_THE_RESIDUAL](../reflections/ON_THE_RESIDUAL.md)
