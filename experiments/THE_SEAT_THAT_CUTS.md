# The Seat That Cuts: a dephased site is blind to exactly what the two blocks it leaves behind have in common

**Date:** August 23-24, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [`simulations/seat_cut_blindness.py`](../simulations/seat_cut_blindness.py), parts `steady | kernel | scope | xy | full | criterion | graphs | sector | deleted`
**Data:** [`simulations/results/seat_cut_blindness/seat_cut_blindness_run.txt`](../simulations/results/seat_cut_blindness/seat_cut_blindness_run.txt)

[MEDIATOR_NOISE_GATE_LEVEL_THREE](MEDIATOR_NOISE_GATE_LEVEL_THREE.md) asked for a
window-stable functional and named two candidates, mutual information integrated
over time and steady-state mutual information. [The Blind Site](THE_BLIND_SITE.md)
§7 ruled out the first: the window mean dilutes as the window grows. This page
tries the second, finds it silent, and finds the law that was standing behind it.

§1 to §4 are read in the **single-excitation sector**; §5 alone leaves it and says
so where it starts. Every mutual information here is in **bits**.

**The steady-state functional is window-stable because it has stopped looking.**
Wherever dephasing at a SINGLE seat leaves no blind subspace, the single-excitation
Liouvillian has a one-dimensional kernel and its steady state is I/N at every seat
and every rate. §1's sweep carries γ on ten or eleven seats rather than one, and
this sentence does not reach that far: what licenses §1's own one-dimensional
kernels is the exact GF(p) rank taken there profile by profile, and the
several-seat statement below is a connectivity CLAUSE with a counterexample, not
a form that could license them. For one seat
that holds on any graph and needs no simplicity: blind(j) = 0 says e_j is cyclic for
H, so anything commuting with H and with the seat's projector sends e_j to a multiple
of itself and is therefore a multiple of the identity. For a support of several seats
it needs connectivity as well, and the run prints the pair: at N = 4 with bonds
[1, 0, 1] and γ on seats 0 and 2 no state is blind to both seats, yet the kernel is 2,
the two disjoint pairs each keeping their own scale, while the uniform [1, 1, 1] at
the same two seats gives 1.

**What the seat can be blind to is decided by a shared factor.** The statement is
about a DIMENSION and not about "a state", because plenty of states survive with
amplitude at the seat. On an open chain **with no zero bond** the stationary space
of the sector is spanned by Q, the sum of the projectors onto the modes that do not
vanish at the seat, together with one projector per node-mode. Q is unnormalised
and carries Q[j, j] = 1 exactly at the seat, so the state it becomes puts
1/(N−blind) there. What is bounded is the part of that space carrying **no**
amplitude at the seat, and its dimension is

    **blind(j) = deg gcd( χ(H_left), χ(H_right) ).**

Here χ(M) is the characteristic polynomial of M; H_SE is the single-excitation
Hamiltonian, the N × N matrix the chain's one-excitation states live on; and
H_left, H_right are its two principal submatrices on the sites strictly left and
strictly right of the seat. So: take those two blocks and count the roots their
characteristic polynomials share. The seat is blind to the eigenvalues those two
blocks have in common, and to nothing else.

**The no-zero-bond hypothesis is load-bearing twice on the chain**, for the
criterion and for the span, and §7 measures what each half is really about: for the
span it is not the zero bond at all. The run sweeps the criterion's half rather than
exhibiting it, over every profile in {0,1,2}^(N−1) for N = 3..6 and every seat, **on
both books**. On HEISENBERG the failure is total: wrong on **all 1682** (profile,
seat) pairs that carry a zero bond, right on **all 316** that do not. **On XY it is
not total**: the same zero-bond set leaves **60** pairs the criterion still gets
right, while the 316 zero-free pairs are right on both books.

**So the fence is not really about the zero bond; it is about the degeneracy a zero
bond forces.** On Heisenberg it always forces one, since each component contributes
the one-magnon descendant of its own ferromagnetic vacuum at the same eigenvalue and
two components therefore repeat a level; on XY there is no such term and a cut chain
can keep a simple spectrum. Simplicity is necessary there and not sufficient, and over
the swept range, N = 3..6, that is counted rather than asserted: of the 1682 XY zero-
bond pairs the four cells (simple, degenerate) × (right, wrong) hold **60, 342, 0,
1280**. The empty cell is the necessity, since no pair with a degenerate spectrum is
ever right; the 342 is why simplicity is not sufficient. A zero bond disconnects the chain, so H_SE stops being
unreduced, and the kernel then outgrows the span too: 7 against 5 at N = 6 with bonds
[1, 1, 0, 1, 1], **at seats 1 and 4**. That profile splits H_SE into two identical
halves, so every eigenvalue is doubled; which seat is watched then decides whether
the excess appears, and the chain's other four seats give 4 against 4. Off that
fence: verified against the exact kernel at every seat of 330 (profile, seat) pairs
over N = 3..8, on uniform, ramp, palindromic and pseudo-random integer profiles, on
EACH book, zero mismatches on both, in exact arithmetic with no eigensolver.

**On the uniform chain that shared factor is the divisor law**, which is the
committed reading of [The Blind Site](THE_BLIND_SITE.md):

    **uniform Heisenberg: blind(j) = (gcd(2j+1, N) − 1)/2**

    **uniform XY:          blind(j) = gcd(j+1, N+1) − 1**

the second closing what stood as §11's first open item of that page at commit
`c9c16d9`. And under the same hypothesis the stationary manifold of the sector, the
kernel of the Liouvillian L_SE restricted to it, has

    **dim ker L_SE(j) = 1 + blind(j)**,

the +1 being the Q direction, the one that does sit on the seat.

## What this is about

Put one site of a chain under a light. That is all dephasing at a single seat is:
one place where the world is looked at, continuously, while everywhere else stays
unwatched.

The light does not cut the chain, and it does not cut the state either. The
excitation still hops straight through the watched site, and the Hamiltonian never
notices. Nor does the watching separate the two sides from each other: it erases
only what touches the watched place, its own row and its own column, and a
surviving state is free to stay coherent straight across it. The **blind** ones all
are, and not by accident. On a Jacobi matrix no eigenvector can vanish at an END of
the chain, because the three-term recursion would then drag the whole vector to
zero, so every blind mode of every zero-free profile carries a non-zero coherence
between the two ends. That is a theorem, not a measurement; the run prints the
three instances at N = 7 with the light on the middle seat, ρ[0, 6] = −0.0538,
−0.1746 and −0.2716. It does not extend to every surviving state: the sector's
maximally mixed state survives too, and being diagonal it has no end-to-end entry
at all.

What the light cuts is the room the surviving states have. Plenty of states survive
with weight on the watched site; what the light takes away is the freedom to be
anywhere else *independently* of it. The part of the surviving space that owes the
watched place nothing at all is what gets counted here, and its price is to be
**nowhere at that place**: one wave with a node there, not two waves that have
stopped talking.

The node is what forces the arithmetic. A wave with a node at the seat is pinned to
end at zero on the left of it and to begin at zero on the right of it, so each side
has to be able to carry that pitch under exactly those end conditions, and the
pitch has to be the same one on both sides. **So the watcher is blind to exactly
the pitches its two sides can hold in common**, and that is §2's shared factor in
words.

One qualifier, because the obvious reading of that sentence is wrong for half of
this page. "The two sides" are not the two free-standing chains you would get by
snipping the seat out. They are the principal submatrices the seat leaves behind,
and with the ZZ term those carry a boundary term at the cut and a different shift
on each side. For the **XY** chain of §4 the two readings coincide at every seat;
for the **Heisenberg** chain of §1 to §3 they do not, and the run prints both. At
N = 11 the free-standing reading invents two blind seats where the chain has none.
The picture is right about the node and about the sharing; it is the halves that
are not free-standing.

**The mirror is not the protection, the mirror is the agreement.** A chain mirrored
about the watched seat has two identical halves, so they agree about everything, so
the watching is as blind as it can be. That is why the middle seat of an odd chain
is the blindest seat. Symmetry is the crudest way to force agreement, not the
reason agreement matters, and the arithmetic of the divisor law is what that
forcing looks like when the chain is regular enough to count.

**That splits the dark states in two.** There are forced ones, where a mirror
guarantees the agreement, and met ones, where two unequal halves happen to find the
same pitch. Happening to agree is a coincidence a generic chain does not have, so a
genuinely irregular chain will usually have no blind seat at all. That is a
plausible reading and not a theorem here. Whether the two kinds behave differently
under disorder, the forced ones holding while the met ones vanish at the first
detuned bond, is not measured here and is the sharpest question this page leaves.

**What it lets someone do.** Two doors, both untried, both directions and not
results. A blind subspace can be *built* at a chosen seat of a chain with no
symmetry at all, by tuning the couplings until the two blocks share as many pitches
as one wants protected states; the only route known before was to make the chain
symmetric. And read backwards, how blind a seat is reports the spectral overlap of
the two blocks it leaves behind, without either block ever being addressed on its
own. What the instrument cannot do: the criterion is a yes or no about an exact
coincidence, so it is structurally blind to a chain that *nearly* agrees, however
nearly.

**Three words here are on loan.** **Light**: γ as illumination is a Tier-4 reading
in [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md) and nothing above leans on it.
**Mirror**: `compute/MirrorWorld/Mirror.cs` and F1 own it as the block-lattice
group and the palindromizer; here it means a chain read backwards. **Divisor**:
F140 and `Divisor.cs` own the frozen divisor on the R₉₀ locus; the "divisor law"
inherited from [The Blind Site](THE_BLIND_SITE.md) is a greatest common divisor of
integers.

## What the repo already held

The sweep the look-back gate asks for. Each store, and what it returned.

**`docs/proofs/`.** `PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md`, item 2 under its
`## Consequences`: "Asymptotic state is a function of (p_0, ..., p_N) alone … The
vector (p_0, ..., p_N) is a complete invariant for the purpose of predicting ρ(∞)."
No γ in that conclusion, so §1's rate-blindness is a consequence of a proven theorem
rather than a measurement. Its Step 2 states §3's and §5's structure in one line:
"The fixed-point algebra of the restricted Lindblad generator L_w = −i[H_w, ·] + D_w
is the intersection of (a) the decoherence-free subalgebra of D_w and (b) the
commutant of H_w." Its one exception is narrower than it looks: what it states and
measures is the SITE-REVERSAL at N = 5; the general form, "under some automorphism
whose fixed set contains the support", it declares open.
`PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md` (Tier 1 derived) proves
`dim ker L = Π_c (|c|+1)` under **uniform** Z-dephasing, resting on
[DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md) Result 2: the only operators with
[H, Q] = 0 and D(Q) = 0 are the identity plus the N popcount projectors.
`PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md` §(b) **Lemma A** gives §2's simplicity
step, and its §(g) gives the reason the non-vanishing bonds matter, "at a vanishing
bond the chain is cut and the pieces can repeat each other". `PROOF_UNIFORM_LAW.md`
B0 holds the sector-1 nondegeneracy and the Neumann cosine modes of the Heisenberg
chain.

**`docs/ANALYTICAL_FORMULAS.md`.** **F4** returned the open question §5 walks into
and its two numbers. **F64** owns the criterion §2 counts, in the words of its
[EQ-015 closure](../review/EMERGING_QUESTIONS.md#eq-015): "**F64 captures
protection,** not just dissipation: rate = 0 ⟺ |v(B)|² = 0 (mode has a node at the
dephasing site)", for Z-dephasing on a single site, exact at every γ, across chain,
ring, star, Y-junction and K₅; the entry is Tier 1-2, analytical. **F65** (Tier 1,
proven, N = 3..30; titled there the uniform open **XX** chain, the object §4 calls
XY) holds the eigenbasis ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)) that makes §4's
proof two lines. **F66**'s Scope clause had already measured that the seat decides:
"at B = center of N=5 chain, α = 0 has multiplicity 64 (not 6), so the N+1 count is
endpoint-specific". **F143** returned §2's fence, same hypothesis and same reason
for a neighbouring object; the only difference is that a zero bond makes F143's
question ill-posed and makes this criterion outright wrong. **F76** gives a
resembling shape on a two-site pair rather than on blocks, and its own `Valid for`
line fences the t → ∞ endpoint, so it is a resemblance to §1 rather than a prior
form of it. **F122** returned the star's zero-sum leaf space, which it calls the
"(N−2)-fold 0-eigenvalue leaf manifold"; the name *dark leaf manifold* is
[THE_HUB_KILLS_THE_HORIZON](THE_HUB_KILLS_THE_HORIZON.md)'s. What the registry did
**not** hold at `c9c16d9`: any steady-state mutual information of a general state,
or any formula for the seat dependence of a kernel dimension.

**`experiments/`.** [CUSP_LENS_CONNECTION](CUSP_LENS_CONNECTION.md) already holds
§1's limit: "the density matrix becomes the maximally mixed state within the SE
subspace", "Purity: 1/N". So ρ_∞ = I/N is not this page's.
[SYMMETRY_CENSUS](SYMMETRY_CENSUS.md) holds the attractor census and the flat
sentence "Cross-sector coherences are always destroyed asymptotically"; its section
is headed `(N=5, uniform gamma)` and it surveyed ring and star too, so the axis §5
pushes on is uniform γ against a SINGLE-SEAT support.
[WEIGHT2_KERNEL](WEIGHT2_KERNEL.md) already reports that the weight-2 kernel
dimension is topology-dependent with no closed form, the same wall §5 hits one
sector up. [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) is
the original owner of "blind subspace", defined there generally as the orthogonal
complement H_M^⊥ of a measurement's detector subspace; this page borrows it in
[The Blind Site](THE_BLIND_SITE.md)'s narrower sense. **That page also holds the
fence-free form of §2's middle step, and this page reached the criterion without
using it**: its §5 defines the blind space as "the largest H-invariant subspace
inside ker(n_k)", the orthogonal complement of the Krylov space the seat generates,
taken as a GF(p) rank at two primes with no eigensolver. That statement asks for no
chain, no simple spectrum and no zero-bond fence. §2's gcd of two halves is what it
becomes once a chain lets the cut matrix fall into two pieces, and §7 measures the
two against each other. The connection was available from the start: this page cites
The Blind Site throughout and never once for this. **No prior computation of an
asymptotic mutual information**: every MI in that directory is read at a finite
time, and [SCALING_CURVE](SCALING_CURVE.md)'s column "MI_steady (t=20)" is a fixed
finite time under a misleading label.

**`docs/CAUGHT_ERRORS.md`** returned the genre §1 belongs to, in its own words:
"when a fixture is invariant under the very thing a test is meant to detect, the
test is decoration"; "A witness that cannot disagree is worse than no witness"; "A
green gate is not evidence until you can say what would make it fail." It also
returned a live disagreement, which §5 settles.

**`docs/GLOSSARY.md`** returned the concept: "When the sum is zero there is nothing
to compare and the reading is empty, which is a third answer alongside balanced and
broken." The word **SILENT** is not in the glossary; it is a verdict string in
`compute/RCPsiSquared.Diagnostics/Foundation/PhysicalGeneratorPolarityBreakWitness.cs`
and a doc comment in
`compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs`, both object-level. This page borrows it for an
**instrument** whose reading is empty, which is an extension.

**The OpenArcs registry** returned `the_gate_that_does_not_gate`, whose 2026-08-23
annotation said "the integrated one is now ruled out, the steady-state one is
untried" (quoted from `c9c16d9`; this change set replaces it). That is the CANDIDATE
this page spends. The arc's item stays open, and §7 says so: two candidates are now
ruled out and no instrument has replaced them.

**Code.** `compute/RCPsiSquared.Diagnostics/DZero/StationaryModes.cs` (reached through
`DZeroDecomposition.cs`), its
Python original `simulations/framework/diagnostics/d_zero.py`, and
`compute/RCPsiSquared.Diagnostics/Ptf/StationaryManifold.cs` all compute a stationary
manifold by a bare tolerance
that nothing at run time compares against a gap; the third documents a measured
separation in its doc comment, so its threshold is justified even though no code
reads it. None of the three can express this page's question, because `ChainSystem`
builds the Liouvillian with a **uniform** γ₀ and cannot place dephasing on a single
seat. Two `simulations/` helpers, `reading_the_30_percent.py` and `mixed_bridge.py`,
take the single argmin|λ| eigenvector as the steady state, which is wrong whenever
the kernel is degenerate.

## 1. The steady state says nothing

Setup as in The Blind Site §7, whose J = 1 is that page's standing convention:
N = 11, uniform chain, halves A = {0,1,2,3,4} and B = {6,7,8,9,10} with the centre
in neither. The sweep moves one seat from γ = 0 to γ = 0.5 against a baseline of
0.05 on the rest. Rather than propagating, the single-excitation Liouvillian is
built as a 121 × 121 matrix and its kernel taken.

Over the 23 profiles the sweep compares, the kernel is one-dimensional and the
steady state is the maximally mixed single-excitation state to machine precision,
worst |ρ_∞ − I/N| = 2.901·10⁻¹⁶. One-dimensionality is certified profile by profile
with the exact GF(p) rank of §3 before the float solver is allowed to pick a single
direction; without that step this section would be doing what the sweep record
faults its two `simulations/` helpers for. The three artifacts faulted above are
faulted for a different thing, a bare tolerance.

The mutual information at ρ_∞ is therefore the same number for all of them, spread
2.220·10⁻¹⁵ bits, and the span is zero at every seat (largest magnitude
3.472·10⁻¹³ %). The value is a closed form in N alone. For halves of size (N−1)/2
with the centre excluded, ρ_∞ = I/N gives S_AB = log₂N (the N−1 occupied sites plus
the one vacuum outcome are all distinguishable) and S_A = S_B = ((N−1)/2N)log₂N +
((N+1)/2N)log₂(2N/(N+1)), so

    I(A:B)|_∞ = log₂N − ((N+1)/N)·log₂((N+1)/2)

which matches the value computed on I/N at N = 5, 7, 9, 11, 13, 15, the worst
difference being 1.78·10⁻¹⁵ bits at N = 15 and the rest below 10⁻¹⁵. At N = 11 the
closed form is 0.639472526941491 and the computed value 0.639472526941490, which is
the double's own resolution and not a disagreement.

Nothing in that expression is a rate, a coupling, a seat, or a preparation. **The
steady-state functional is window-stable in the way a closed eye is steady.** In the
repo's vocabulary it is SILENT at the level of the instrument: not balanced, not
broken, empty.

**One arm this does not cover**, and it is the arc's own first row. The Blind Site
§7's table begins with the profile "none (only the swept seat carries γ)". Read as a
profile that is the ZERO support, every seat at γ = 0, where the kernel is the whole
commutant of H (dimension 11 at N = 11) and the steady state is undefined without the
initial state, so the guard refuses to compute it. It is the γ = 0 endpoint that
blocks the arm and not the single-seat support: a single seat gives kernel 1 at ten
of the eleven seats and 6 at the centre. The
negative verdict is therefore about the arm where at least ten of the eleven seats
carry γ, which is the arm the sweep compares.

## 2. What the seat cuts is the set of survivors

Dephasing at seat j kills exactly the entries of ρ whose seat bit differs, which in
this sector means the seat's own row and column OFF THE DIAGONAL. The qualifier is
load-bearing: ρ[j, j] has the seat's bit on both sides, so it survives, and it is the
entry that supplies the +1 in §3's dim = 1 + blind. It does **not** kill the entries
that cross the seat: at N = 7 with γ on seat 3 the mask leaves ρ[0, 6] untouched and
kills ρ[2, 3]. Nor is it only the states with no amplitude at the seat that survive:
Q and every combination containing it survive while sitting on the seat. What the
criterion counts is the part of the stationary space that owes the seat nothing.
Made exact:

    blind(j) = deg gcd( χ(H_left), χ(H_right) )

**Why, and only the middle step is this page's.** That a mode with a node at the
dephased site is exactly a mode the dephasing cannot touch is **F64**'s. That a path
with every bond non-zero is non-derogatory, which the proof states for an unreduced
tridiagonal and which covers H_SE, is
[PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE](../docs/proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md)
§(b)'s **Lemma A**, stated there "in either Δ book". §(g) is where that document
says it of a Hermitian matrix outright, "Lemma A applies to the Hermitian Jacobi
matrix A_J as well, so a path's hopping spectrum is simple whenever A_J is Hermitian
and every J_b ≠ 0", which is exactly what H_SE is. §(g)'s own setting is the Δ = 0
adjacency block, and the argument it DISCARDS is a Δ = 1 Laplacian one; what it keeps
is Lemma A again, which is stated for either Δ book. Non-derogatory gives geometric
multiplicity one; H_SE is real
symmetric, hence diagonalisable, so geometric one forces algebraic one and the
spectrum is simple. What is added here is the step between: for a Jacobi matrix an
eigenvector vanishes at a site exactly when that eigenvalue is an eigenvalue of both
principal submatrices the site leaves behind, so counting the shared eigenvalues
counts the node-modes, and F64 turns that count into the blind dimension.

Both directions of the middle step are elementary. If λ is shared, the last
component of an eigenvector of a Jacobi H_left is non-zero, so the two blocks'
eigenvectors can be scaled to satisfy the equation at the seat with amplitude zero
there. Conversely, if an eigenvector vanishes at the seat, the three-term recursion
forces each restriction to satisfy its own block's equation, and neither restriction
can itself vanish, because the recursion would then drag the whole vector to zero.

**Verified against the kernel rather than against the eigenvectors.** The run
computes `blind_by_gcd` in exact `Fraction` arithmetic (Faddeev-LeVerrier for the
characteristic polynomials, Euclid for the gcd, no eigensolver anywhere) and
compares it to the exact kernel dimension over 330 (profile, seat) pairs across
N = 3..8, on uniform, ramp, palindromic and pseudo-random integer bond profiles, on
EACH book. **Zero mismatches on both.**

**Reflection symmetry is neither necessary nor sufficient for the uniform-chain
law's value.** The qualifier is the whole of the claim: reflection about seat j IS
sufficient for blindness AT that seat, derived below under "What reflection
symmetry buys", and what
fails is the step from the symmetry to the law's profile. [The Blind
Site](THE_BLIND_SITE.md) already scopes the mirror as "the odd-N-centre special case
of the divisor condition and not the general mechanism"; what is added here is that a
palindrome departs from the law in BOTH directions. The run says so in one table:

| profile | measured blind | uniform-chain law |
|---|---|---|
| N = 7 palindrome [1,2,1,1,2,1] | 0, **1**, 0, 3, 0, **1**, 0 | 0, 0, 0, 3, 0, 0, 0 |
| N = 9 palindrome [1,2,3,4,4,3,2,1] | 0, **0**, 0, 0, 4, 0, 0, **0**, 0 | 0, 1, 0, 0, 4, 0, 0, 1, 0 |
| N = 5 asymmetric [1,4,2,2] | 0, 0, 2, **1**, 0 | 0, 0, 2, 0, 0 |
| N = 4 asymmetric [1,3,2] | 0, **1**, 0, 0 | 0, 0, 0, 0 |

Bold marks every cell that departs from the law: a reflection-symmetric chain with
more blindness than the law, one with less, and two asymmetric chains gaining
blindness at a seat the law calls sighted. The N = 5 row carries a second point that
is not bolded because it is an agreement: that chain reaches the law's own centre
value 2 with no reflection anywhere, so the reflection is not necessary even at the
centre. The criterion reproduces all four rows exactly.

**What reflection symmetry buys**, read off the criterion: a chain
reflection-symmetric **about seat j** makes the two blocks equal up to the reversal
that relabels one onto the other, so they share a characteristic polynomial, the gcd
is the whole of it and blind(j) = j. Such a reflection needs 2j = N−1, so the only
seat it can ever speak about is the centre of an odd chain. The general ceiling is
blind(j) ≤ min(j, N−1−j), the smaller block's size.

## 3. The uniform chain, the divisor law, and the kernel

On the uniform chain the shared factor is the cosine coincidence that
[The Blind Site](THE_BLIND_SITE.md) counts, and the criterion evaluates to the
committed divisor law:

    blind(j) = (gcd(2j+1, N) − 1)/2

Verified at every seat of every N from 3 to 13 against the exact kernel, zero
mismatches, and the modular identity behind it,
#{m ∈ 0..N−1 : m(2j+1) ≡ N (mod 2N)} = (gcd(2j+1, N) − 1)/2, checked by integer
enumeration for N = 2..200.

**And the stationary manifold is one larger:** dim ker L_SE(j) = 1 + blind(j), the
+1 being the direction spanned by Q, the sum of the projectors onto the modes that
do not vanish at the seat. The sector's maximally mixed state is in the kernel too,
being Q plus every node projector, but it is not the extra direction: it has a
component in the blind part. The same Jacobi simplicity gives the identity: H_SE
simple means its commutant is the diagonal in the eigenbasis, and the seat's
constraint forces one common coefficient across every mode that does not vanish at
the seat while leaving the rest free.

**There is no tolerance in this.** The kernel dimension is an exact GF(p) rank on
integer inputs. For Hermitian ρ the commutator term of Tr(ρ·L(ρ)) is identically
zero by cyclicity, so Tr(ρ·L(ρ)) = 0 leaves a sum of non-positive terms that must
each vanish, giving ker L_SE = {ρ : [H, ρ] = 0 and ρ_ij = 0 wherever a dephased
seat's bit differs between i and j}, which for integer J is a plain integer
commutant. Two preconditions the argument needs and the code does not check:
**γ_k ≥ 0 at every site** (a gain term breaks the non-positivity the argument rests
on) and, for the "J does not enter" reading, a **uniform** rescaling of J rather
than a change of profile, since a profile is exactly what §2 shows does move the
answer.

**The rank is one-sided.** A GF(p) rank can only be smaller than the rational rank,
so a reported kernel dimension can only be too large, never too small, and only if p
divides a pivot minor. Two things weaken that without removing it: §5's `sector`
part takes every **diagonal** block at a second, unrelated prime (2³¹−1) by a
different implementation and the two agree block by block, which is shipped and
rerunnable; and the direction matters, since the one committed formula this page
refutes is refuted by a number above ours rather than below it. An independent
reviewer also rebuilt the same kernel over ℚ(i) in exact fractions over N = 3..13
with zero disagreements; that one is a review report and is not reproducible from
the artifacts this page ships.

**A wide singular-value gap is not evidence.** At N = 11 with γ on the centre, where
the answer is 6, an SVD rank returns 6 at J = 1, **21 at J = 10⁻⁵ while reporting a
gap of 5.95·10³**, and 101 at γ = 10⁹. MirrorWorld's `Divisor` records the trap: "a
floating-point rank silently miscounts it once the coupling is small and the chain
long, where the other eigenvalues crowd the root at spacing J^(2d)"
(`compute/MirrorWorld/Divisor.cs`).

## 4. The XY chain

On the **uniform** XY chain the node condition is m(j+1) ≡ 0 (mod N+1) and the count
is

    blind(j) = gcd(j+1, N+1) − 1

with no halving, against the Heisenberg (gcd(2j+1, N) − 1)/2.

**The proof is two lines and the first is already committed.** The uniform XY
single-excitation modes on an open chain are ψ_m(j) ∝ sin(πm(j+1)/(N+1)) for
m = 1..N, which is **F65**'s eigenbasis, Tier 1 and proven there for N = 3..30; the
second line is that a node at seat j means (N+1) divides m(j+1). Writing
d = gcd(j+1, N+1), those m are exactly the multiples of (N+1)/d, and d−1 of them lie
in 1..N. The argument is uniform-chain only, because it names that eigenbasis; off
the uniform chain §2's criterion is the statement, with the ZZ term dropped.

Three independent routes agree. The exact kernel dimension of §3 with the ZZ term
dropped, at N = 6, 7, 9, 11, 12, 13. The node count read off the eigenvectors,
N = 3..20. And an exact integer enumeration with no floating point anywhere,
N = 2..200: zero mismatches on both laws. The third certifies the **count**, that is
the arithmetic of the node condition, and not the kernel identity.

The two laws disagree loudly, which is what makes the pair worth having:

| N | Heisenberg blind seats | XY blind seats |
|---|---|---|
| 6 | 1:1, 4:1 | none |
| 7 | 3:3 | 1:1, 3:3, 5:1 |
| 11 | 5:5 | 1:1, 2:2, 3:3, 5:5, 7:3, 8:2, 9:1 |
| 12 | 1:1, 4:1, 7:1, 10:1 | none |
| 16 | none | none |

At the reflection-fixed seat of an odd chain the two agree, both giving (N−1)/2,
since gcd(N, N) = N on one side and gcd((N+1)/2, N+1) = (N+1)/2 on the other. That
is The Blind Site §7's claim that the exact blindness at the centre does not need
Heisenberg over XY, reached a second way. The reflection argument stays the better
one at that seat, because it names no eigenbasis; it is not the general statement,
though, since a reflection about seat j speaks about that seat and no other.

## 5. F4's open question, and a reading that failed

This section alone leaves the single-excitation sector. It works two ways: on the
whole 4^N space by an SVD, which is where the graph table and F4's own numbers come
from, and per popcount block by an exact GF(p) rank, which needs no tolerance and
reaches further. Where both apply they agree to the last digit.

`docs/ANALYTICAL_FORMULAS.md`, F4, in its own words: "on the N = 3 open chain at
J = 1, γ = 0.5 on the **end** seat alone already gives kernel 4, while the same γ on
the **middle** seat gives 6. Which seat carries the γ decides, and why is an open
question rather than a formula here." Both numbers reproduce. Two things are added,
and neither closes it.

**On the open Heisenberg chain the full-space kernel is block-diagonal in popcount,
and the qualifier that does the work is the ZZ term, not the topology on its own.**
The popcount of a computational basis state is its number of excitations, so the
operator space splits into (p, q) blocks by the popcount of the ket and of the bra,
and "block-diagonal" means the kernel lives entirely in the p = q blocks. §4's own
second object is the counterexample: drop the ZZ term and keep the same uniform open
chain at N = 5, and the kernel runs 6, 20, 24, 20, 6 by seat against diagonal sums of
6, 10, 16, 10, 6, that is cross-sector weight 0, 10, 8, 10, 0. The XY open chain
carries cross-sector stationary coherence.

Put the ZZ term back and the same chain behaves the other way. The `sector` part
ranks every one of the (N+1)² popcount blocks exactly, the cross ones included, so
the split is measured rather than inferred from diagonal sums: cross-sector weight is
**exactly zero at every seat** of the uniform chain at N = 5 and N = 7, of
[1, 4, 2, 2] at N = 5, and of [2, 1, 1, 2, 1] at N = 6, and the exact route
reproduces the SVD wherever both are run. At N = 5 with γ on the centre the twelve
Heisenberg kernel directions distribute as 1, 3, 2, 2, 3, 1 across the six sectors.
Not claimed: that this holds for every N, or off the chain, or with the ZZ term
dropped, where it is measured to fail.

Each sector contributes the commutant of that sector's Hamiltonian, intersected with
the operators vanishing wherever the dephased seat's bit differs, which is the
structure `PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md` states and
`PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md` proves the **uniform**-γ face of. In the
single-excitation sector H is nondegenerate and that dimension is 1 + blind, which is
§3. For popcount ≥ 2 the count is **larger**: sectors 2 and 3 give 2 each where
1 + blind would give 1.

**Which breaks the obvious sum rule.** The natural guess is (N+1) + Σ_w blind_w(j),
with blind_w(j) this page's count carried over to the popcount-w sector. Read
strictly it is not well posed, since blind_w is defined nowhere here and computed
nowhere; treated as the arithmetic statement that each sector contributes 1 plus
whatever is blind in it, it fails visibly. The per-sector lists at the
reflection-fixed centre are [1, 3, 2, 2, 3, 1] at N = 5 and [1, 4, 2, 2, 2, 2, 4, 1]
at N = 7. Sectors w = 1 and w = N−1 carry 1 + blind exactly (3 and 3 at N = 5, 4 and
4 at N = 7), w = 0 and w = N carry one diagonal state each, and every sector strictly
in between carries **2** where the rule allows only 1, which is where the measured 12
and 18 exceed the rule's 10 and 14.

The rule survives at N = 3, 4 and 6 for two different reasons. At N = 3 there is no
sector strictly between w = 1 and w = N−1 at all, so nothing can overflow. At N = 4
and N = 6 there is (w = 2, and w = 2, 3, 4) and each carries **1**: what those lack is
not the sector but the blindness, since the middle seat of an even chain has b = 0,
and where nothing is blind every sector carries 1 and any rule of this shape is
trivially right. The two readings are one statement from two sides: the rule fails at
exactly the N where a blind middle seat and an in-between sector exist together, and
over the N this page reaches, N = 3, 4 and 6 each miss one of the two.

**The reason is not degeneracy.** On the uniform chain the popcount-2 and popcount-3
sector Hamiltonians are simple, checked exactly rather than by a smallest gap: a
spectrum is simple exactly when the characteristic polynomial is squarefree, that is
when it shares no factor with its derivative, and that degree is 0 at N = 5, 6 and 7
in both sectors. Heisenberg degeneracy is across popcount sectors, in the SU(2)
multiplets, not inside one. That last is a statement about the UNIFORM chain: the
reflection-symmetric profiles below carry a doubled eigenvalue INSIDE popcount 2 and
popcount 3, so it is the profile that puts a degeneracy inside a sector.

**What sets a sector's count, and the hypothesis that travels with it.** On a
**simple** sector the count is the number of connected components of the graph on
that sector's eigenmodes, two modes joined when the dephased seat's occupation number
has a non-zero matrix element between them: ρ commutes with H_w, so it is diagonal in
the eigenbasis, and the seat's condition forces one common coefficient along every
edge. In the single-excitation sector the seat's occupation is a rank-one projector,
so every mode that does not vanish there fuses into one component and each node-mode
stands alone, which is §3's 1 + blind. Sectors above the first can carry more than one
component, and do not have to carry more than the first: at N = 5, centre seat, the
per-sector kernel runs 1, 3, 2, 2, 3, 1, so popcount 2 carries fewer than popcount 1.

**Off a simple sector the count depends on the basis, and this page first read that
as the mechanism failing.** It is not. On a degenerate sector ρ need not be diagonal
in the solver's eigenbasis, and the basis the mechanism's own argument entitles it to
is the one that diagonalises the seat's occupation INSIDE each degenerate eigenspace.
Take that basis and the count matches the exact kernel on both profiles where the
solver's basis undercounted: at N = 5 with bonds [4, 3, 3, 4] the popcount-2 sector
has a doubled eigenvalue, the exact kernel is 2, the solver's basis says 1 and the
adapted basis says **2**; N = 7 with bonds [4, 3, 4, 4, 3, 4] does the same in
popcount 3.

**The test this page offered against that reading could not have failed.** The
component count is re-run on fifty random orthogonal re-bases of each degenerate
eigenspace and returns one value throughout, and that was read as the count being a
property of the sector rather than of the basis. The adapted basis is measure-zero
inside a degenerate eigenspace, so a random draw misses it with probability one and
returns the generic, maximally connected value every time. A one-value set is
therefore evidence of nothing here, which is this page's own quotation from
`docs/CAUGHT_ERRORS.md` turned on itself: a witness that cannot disagree is worse
than no witness. Both columns are printed side by side rather than argued.

**The adapted basis is not the general recipe either**, and the general statement is
algebraic rather than a basis at all. The kernel is the commutant of the algebra A the
sector's H and the seat's occupation generate. Writing A's Wedderburn blocks as
(n_i, m_i) with Σ n_i m_i the sector dimension, the kernel is **Σ m_i²** while the best
component count any eigenbasis can reach is **Σ m_i**; the two agree exactly when every
m_i = 1, that is when the commutant is abelian. Fitting the basis to the seat's
occupation inside each degenerate eigenspace happens to reach the best count on these
two profiles, but it can fall short elsewhere and is undefined when the occupation is
itself degenerate there.

That settles both rows with no solver at all, and the arithmetic is the whole of it: a
multiplicity 2 costs 4, so **a kernel of 1, 2 or 3 forces every m_i = 1**, and the best
basis then realises the kernel as its component count. Both exhibited kernels are 2.
**Four is the first kernel dimension at which a multiplicity 2 fits**, so it is the
first at which the mechanism can fail at all; the smallest failure this page exhibits
sits one above it. The star at N = 4 with the HUB dephased has blocks (2,1) and (1,2)
in the single-excitation sector, so its kernel is 5 while no basis reaches more than
3 components. The `sector` part measures the four invariants that force those blocks
rather than asserting them. So the mechanism CAN genuinely undercount, the shortfall is
Σ m_i² − Σ m_i, and the two profiles this page exhibited were never the witnesses.
What stays open is the narrower question, whether a zero-free open chain can do it;
nothing here exhibits one either way.

Every kernel and every simplicity degree in that table is exact; the component counts
are the only numbers in it taken off a float eigenbasis, and after the above they
exhibit a basis effect rather than certifying anything. The script calls `eigh` in
three further places, one of them §4's node count, and `svd` in three; the tolerance
paragraph at the end of this section prices them.

Note what this costs the neighbouring question: at N = 5, popcount 2, **nothing is
blind at all** and the count is still 2, so a component is not a blind state and the
two must not be conflated.

**Off the chain the block-diagonality need not hold, and the reading offered for
which graphs keep it is a null result.** The `graphs` part measures the cross-sector weight at every seat of ten
graphs, as the trace of the kernel projector per block, which is basis-independent.
The N = 4 ring carries weight exactly 4 at every seat. The reading under test was that
the carrier is a non-adjacent pair of vertices with identical neighbourhoods, which
graph theory calls *false twins*. The table refuses it on both sides:

- **not sufficient**: the star has three such pairs at N = 4 and six at N = 5 and
  carries no cross-sector weight at any seat;
- **not necessary**: the graph {01, 02, 03, 12, 13, 24} on five vertices has no such
  pair at all and carries 6, 12, 6 at seats 2, 3, 4;
- **and no seat rule survives**: on C₄ with a pendant the weight is zero at exactly
  the pair's two seats, on K₄ minus an edge it is non-zero at exactly the pair's two
  seats and zero elsewhere, and on C₄ every seat lies in a pair and every seat
  carries 4.

What survives is the open Heisenberg chain's statement and nothing more: **which
systems carry cross-sector stationary coherence is open, and it is not a question
about topology alone, since the same chain with the ZZ term dropped already carries
it.**

The repo holds three pieces near this and none of the joins.
[PROOF_RING_N4_DIHEDRAL_LOCK](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md) has a
section "Why this is N=4-specific" holding one half, that C₄ **is** K₂,₂ and that this
is unique to N = 4; the other half is in its opening, where the Casimir spectrum has a
seven-fold zero level built from rows of different total sublattice spin, an H
degeneracy spanning popcount sectors. **F89d** is where the N = 4 self-fold is
recovered, as the degenerate partner-equals-self case of its cross-block mirror; the
sentence "at N=4 the DE sector is its own Hamming complement" is
[F89_PATH_K_DIABOLIC](F89_PATH_K_DIABOLIC.md)'s. And the star's zero-sum leaf space is
F122's, named the *dark leaf manifold* by
[THE_HUB_KILLS_THE_HORIZON](THE_HUB_KILLS_THE_HORIZON.md), the nearest thing in the
repo to an explanation of the star row. Also near: that ring proof says of its own
mechanism, "Bipartite-completeness therefore cannot be the cause of the number; it is
one route to computing it."

**A committed line is false, and this page settles it.** `docs/CAUGHT_ERRORS.md`
carries the parenthetical "kernel excess is b²+b, not gcd". Reading b as the blind
count at the seat, the `sector` part evaluates the full kernel at the reflection-fixed
centre of the odd open chain by an exact rank on every popcount block, cross blocks
included, with no eigensolver and no tolerance: **6 at N = 3, 12 at N = 5, 18 at
N = 7**, against the formula's 6, 12 and **20**. It agrees at the two smallest N and
fails at the smallest N that can tell them apart, which is also the first N the SVD
route could not reach.

Three things guard that number, and they are separate arguments.

First, it is a **full** kernel and not a diagonal sum: every cross block is ranked and
every one comes out exactly zero. That half needs no second prime and is not a
certificate but a proof. A GF(p) rank can only be smaller than the rational one, so the
nullity it reports can only be too LARGE; a measured nullity of zero therefore forces
the rational nullity to be zero. One prime settles an empty block.

Second, the **diagonal** blocks are computed twice, by the pure-Python elimination at
2⁶¹−1 and by a vectorised one at 2³¹−1, and the two lists agree block by block. Only
the diagonal: the cross blocks are ranked at 2³¹−1 alone. Two primes do not discharge
the one-sidedness §3 concedes; they make it a coincidence two unrelated primes would
have to share.

Third, and not probabilistic, the remaining error can only run one way. If 18 is wrong
it can only be too **large**. The formula says 20, above 18, so every error this method
is capable of moves our number *towards* the formula rather than away. The refutation
survives its own worst case.

3(N−1) fits all three, and three points are not a law, so that is a fit and not a
replacement. Repairing `docs/CAUGHT_ERRORS.md` is a separate act on an append-only
file, whose own rule is a note appended and never an edit in place; the note is
appended there under 2026-08-24 and anchors back to this section.

**Six counts on this page rest on a tolerance**, and it has just faulted three
artifacts for the same thing, so each is named rather than assumed. First, this
section's full-space kernel dimensions are taken by a bare TOL = 10⁻⁸ on the singular
values, with the narrowest deciding ratio over the N = 3..6 table equal to
7.35·10¹²; the exact block route recomputes the same numbers with no tolerance
wherever the exact route is also run, which is not every row the SVD covers. Second,
**the ten-graph table earlier in this section takes the same SVD route and prints no
gap at all**, and it is named here because it carries the false-twin null result. Third, §4's node count off the eigenvectors uses a bare
10⁻⁹ and also prints no gap; it is corroborated by two exact routes, which is why it
is kept, not because the threshold is defended. Fourth, the component count above is
a float eigensolver by necessity, being a claim about an eigenbasis, and it certifies
nothing. Fifth and sixth, two further `eigh` calls carry the same bare 10⁻⁹, and
each DECIDES A PARTITION rather than merely reading a vector: the one behind the
ρ[0, 6] entries at the top of this page splits the modes into blind and sighted, and
the run prints the count of blind ones from that split; the one behind the Q rows
makes the node against non-node partition those four rows display. Neither prints a
gap either.

## 6. What is borrowed and what is ours

| | |
|---|---|
| ρ(∞) is a function of the sector populations alone, no γ in the conclusion | `PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md`, item 2 under `## Consequences` |
| the fixed-point algebra as DFS ∩ commutant; and that the limit refines under the SITE-REVERSAL's fixed set at N = 5 (the general "some automorphism" form is open THERE, not established) | the same proof |
| the commutant characterisation the N+1 rests on | `PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md` + [DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md) Result 2 |
| ρ_∞ = I/N on the single-excitation sector, stated in `experiments/` already | [CUSP_LENS_CONNECTION](CUSP_LENS_CONNECTION.md) |
| the uniform divisor law (gcd(2j+1, N) − 1)/2 itself | [The Blind Site](THE_BLIND_SITE.md) |
| the node condition m(j+1) ≡ 0 (mod N+1) as the XY form to expect | The Blind Site §11, which predicted it |
| the sector-1 nondegeneracy §3's commutant argument leans on, and the Neumann cosine modes of the HEISENBERG chain | `PROOF_UNIFORM_LAW.md` B0 (it does NOT hold §4's XY sine basis, a different Hamiltonian) |
| **rate = 0 ⟺ the mode has a node at the dephasing site**, which is the test §2 counts | **F64**, via its [EQ-015 closure](../review/EMERGING_QUESTIONS.md#eq-015) |
| the uniform XY sine eigenbasis §4's two-line proof names | **F65** |
| that the seat decides the α = 0 multiplicity at all (multiplicity 64, not 6, at the N = 5 centre) | **F66** Scope |
| that a path with every bond non-zero is non-derogatory "in either Δ book", hypothesis included (simplicity then follows for our real symmetric H_SE) | [PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE](../docs/proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md) §(b) **Lemma A** |
| the reason the non-vanishing bonds matter, "at a vanishing bond the chain is cut and the pieces can repeat each other" | the same proof's §(g), where it is one of the two load-bearing hypotheses of §(g)'s own simplicity statement about the hopping spectrum, the other being Hermiticity of A_J. Of that PAIR Lemma A carries only the non-vanishing one, "needing no Hermiticity" in §(g)'s own words (its other hypothesis is pathness, which row 1 above carries); and it is Lemma A, not the simplicity statement, that holds in either Δ book. Not to be read as the proof's labelled **Corollary**, which is §(c)'s and says something else |
| the phrase "blind subspace" | [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md), where it is a measurement's H_M^⊥ generally; borrowed here in The Blind Site's sense |
| the N = 3 end-against-middle kernel numbers, 4 and 6 | **F4** |
| that the weight-2 kernel is topology-dependent with no closed form | [WEIGHT2_KERNEL](WEIGHT2_KERNEL.md) |
| cross-sector coherences destroyed asymptotically, stated unconditionally | [SYMMETRY_CENSUS](SYMMETRY_CENSUS.md) |
| the empty-reading concept, and SILENT as its name | `docs/GLOSSARY.md` for the concept, `compute/RCPsiSquared.Diagnostics/Foundation/PhysicalGeneratorPolarityBreakWitness.cs` for the word |
| "a green gate is not evidence until you can say what would make it fail" | `docs/CAUGHT_ERRORS.md` |
| the float-rank trap at small J and long chains | MirrorWorld `Divisor.cs`, whose own object (the frozen divisor on the R₉₀ locus) is unrelated |
| **blind(j) = deg gcd(χ(H_left), χ(H_right)), for any bond profile with no zero bond**: the step from F64's node test to a COUNT | this page |
| **that reflection symmetry is neither necessary nor sufficient for the uniform-chain law's value, while buying blindness outright at the seat it is symmetric about** | this page |
| **dim ker L_SE(j) = 1 + blind(j)**, argued from the simplicity Lemma A gives | this page |
| **I(A:B)\|_∞ = log₂N − ((N+1)/N)log₂((N+1)/2)** | this page |
| **the uniform XY law gcd(j+1, N+1) − 1**, the N = 200 enumeration certifying the COUNT and not the kernel | this page |
| **that the SYMMETRY_CENSUS sentence holds under a single-seat support on the open chain and FAILS on the N = 4 ring** | this page |
| **that the false-twin reading of the cross-sector weight is false both ways** | this page, after a reviewer's five-vertex sweep |
| **that the XY open chain carries cross-sector weight, so the qualifier is the ZZ term** | this page |
| **that `CAUGHT_ERRORS.md`'s "b²+b" is false at N = 7**, on a FULL kernel: cross blocks ranked and zero (one prime suffices for an empty block), diagonal blocks agreeing at two | this page, `sector` |
| **that the two palindromic profiles are NOT counterexamples to the component mechanism**, and the law that decides it: dim ker = Σ m_i² against a best component count of Σ m_i over the Wedderburn blocks of ⟨H_w, n_seat⟩, so a kernel below 4 can never be a witness and the star's hub is one, its four block invariants measured | this page, `sector` (a3) and its star rows |
| **that the two blocks are principal submatrices and NOT free-standing subchains, which coincide on XY and differ on Heisenberg** | this page, `criterion` |
| **that the fence itself is BOOK-SPECIFIC, and that what it is about is the DEGENERACY a zero bond forces rather than the zero bond**: total on Heisenberg, where each component contributes the one-magnon descendant of its own ferromagnetic vacuum at the same eigenvalue, and not total on XY, which has no such term and where a cut chain can stay simple; simplicity is there NECESSARY and not sufficient over the swept N = 3..6, counted rather than asserted (60, 342, 0, 1280 in the four cells) | this page, `criterion` |
| **that the criterion is a chain evaluation of a fence-free law**: blind(j) = deg gcd(χ(H), χ(H with row and column j struck)), matching the definition at every seat of twenty graphs on both ZZ books, degenerate spectra and zero bonds included, so neither the fence nor the simplicity hypothesis is the phenomenon's | this page, `deleted`, on the fence-free form committed in The Blind Site §5 |

## 7. Open

- **F4's question**, still. §5 gives its shape on the chain and its first sector; the
  popcount ≥ 2 contributions are a commutant dimension with no closed form here, and
  [WEIGHT2_KERNEL](WEIGHT2_KERNEL.md) reached the same wall from the other side.
- **Which systems carry cross-sector stationary coherence.** The false-twin reading
  is dead, and so is the chain-versus-graph framing: the XY open chain carries it, so
  whatever decides this is not the topology. One direction comes out of §2 rather
  than out of staring at examples, and it is a direction and not a conjecture. On a
  chain the criterion is really a statement about the PRINCIPAL SUBMATRIX left by
  striking row and column j, and on a general graph that submatrix is usually
  connected, so there are no two spectra to compare and the chain's phrasing does not
  typecheck. What carries over is that submatrix's spectrum as a whole, against the
  full one, and §2's shape is what it becomes where the cut falls into two pieces.
  The `deleted` part measures it: **blind(j) = deg gcd(χ(H), χ(H with row and column
  j struck))** reproduces the definition at every seat of all twenty graphs, on the
  ZZ book and off it, where the definition is N minus the rank of the Krylov matrix
  the seat generates, taken exactly at two primes with no eigensolver. Fourteen of
  those twenty have a degenerate spectrum on the Heisenberg side and thirteen on the
  XY side, so **the COUNT carries no simplicity hypothesis and no zero-bond fence**:
  the halves form is defined on the five paths and wrong on the three that carry a
  zero bond, while the general form is right on all three. For the count the fence
  belongs to the phrasing and not to the phenomenon. **The SPAN is a separate question
  and does not generalise**: dim ker L_SE(j) = 1 + blind(j) holds on twelve of the
  twenty graphs on the ZZ book and eleven off it, [1,1,0,1,1] at seats 1 and 4 giving
  a count of 4 against a kernel of 7. **And the zero bond is not what decides it**,
  which is worth saying because the opening summary calls the fence load-bearing twice.
  Counting the failures that carry "no zero bond" would turn on whether a missing edge is
  WRITTEN as J = 0 or simply left out, which is why the count given here is by
  CONNECTIVITY instead: of the eight failures on the ZZ book and the nine off it, five
  are on connected graphs, and it is the same five both times (both stars, K₄, K₅, the
  bridged triangles), while the zero-bond path [1, 0, 1], itself disconnected, holds. What every failure does carry is
  a degenerate spectrum, and some degenerate spectra hold anyway (the N = 4 ring
  does), so on this table **simplicity is SUFFICIENT for the span and not necessary**. Carrying the fence-free
  criterion across to the kernel would carry it past a hypothesis the kernel needs and
  the count does not. Say submatrix and not "the deleted graph":
  rebuilding the rest as a free-standing graph is the reading §2 fences off and the
  `criterion` part measures wrong at every Heisenberg row, and the distinction survives the
  generalisation, since a ZZ term leaves on the cut site's NEIGHBOURS a diagonal
  shift that an induced subgraph does not carry. What none of this decides is the open item itself,
  which asks WHICH systems carry the coherence; it removes the typecheck failure that
  stood in the way of asking on a graph.
- **The two kinds of dark state.** §2 makes blindness a coincidence between two
  spectra, which a mirror can force and an irregular chain can only meet by chance.
  Forced blindness should survive disorder and met blindness should not. Nothing here
  measures that, and it is the cheapest open item on this list: the criterion is
  exact, so a disorder sweep needs no propagation and no eigensolver.
- **Whether §2's criterion needs a proof in `docs/proofs/`.** Half of what such a
  proof would carry is already there: §(b)'s Lemma A owns the non-derogatory step in
  either Δ book, hypothesis included, §(g) owns the reason the hypothesis matters,
  and F64 owns the node criterion. What has no proof file is the step between them
  and the kernel identity dim = 1 + blind, both written here as prose and neither
  gated.
- **The blind operator space** inside the (1,1) block, which The Blind Site §11 lists
  as open. §3 counts operators killed by the whole Liouvillian, not operators killed
  by the dissipator; whether the two objects coincide is not settled.
- **Whether a discriminating window-stable functional exists at all.** Both
  candidates the arc named are spent. A third, the **spectral gap**, is window-free
  and rate-sensitive and untried on this comparison; it is named in an untracked
  design spec under `docs/superpowers/`, which is why no quotation from it appears
  here, and it is offered there for a different question.
