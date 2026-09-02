# PROOF: The Crack and the End-Pair Anisotropy Carry One Blindness Locus, Sector by Sector

**Registry:** none. The result is harvested from [F157](../ANALYTICAL_FORMULAS.md) and
[F160](../ANALYTICAL_FORMULAS.md) and adds no closed form of its own that those two do not
already carry; whether the harvest earns its own number is left open rather than decided here.
**Status:** Tier 1 derived. Lemma 1 is a matrix identity at every N and every complex u. Lemma 2a
is the fence-free Cramer argument of [The Seat That Cuts](../../experiments/THE_SEAT_THAT_CUTS.md)
§7 read per eigenvalue; Lemma 2b, the sector split, is this file's. Lemma 3 is Lemma J of the
node-lemma proof. Lemma 4 is this file's own boundary-system argument and is deliberately
independent of F160, whose
simplicity clause is fenced to u ≥ 0 while half of the u axis here is negative. The Theorem is
exact at every N ≥ 3 and both parities. Nothing here touches MULTIPLICITY: every locus below is a
SET, F157's blind COUNT at a given Δ is a multiplicity, and this file claims nothing about it.
**Date:** 2026-09-02.
**Authors:** Thomas Wicht, Claude (Opus 5).
**Script:** [`simulations/blind_seat_two_axes_proof.py`](../../simulations/blind_seat_two_axes_proof.py),
25 checks under 22 labels in five blocks (L2a fires twice, L2b three times), exact in sympy,
about 10 seconds measured quiet under sympy 1.14.0, which is not in the dependency line in
`CLAUDE.md`. Run committed at
[`blind_seat_two_axes_proof_run.txt`](../../simulations/results/blind_seat_two_axes/blind_seat_two_axes_proof_run.txt).
**Closes:** open item 1 of [The Blind Seat on the Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md),
whose §(c) was an observation with an identified mechanism at odd N = 5..17 and at 44 of 66
even-N seats. That page states the relation; this file states why it holds, why it fails at the
other 22, and why it holds at every N rather than at the ones that were reached.

## Statement

Let A be the adjacency matrix of the open N-site path, the single-excitation XY block, and put

    P = E₀,ₙ₋₁ + Eₙ₋₁,₀        the wrap bond, F160's road, u = 0 open chain, u = 1 ring
    D = E₀₀ + Eₙ₋₁,ₙ₋₁         the end-pair diagonal, F157's anisotropy axis

For a seat j write **locus(j)** for the set of real knob values at which seat j is blind, that is
the real root set of Res_x(χ(H), χ(H_j)) with H_j the matrix H with row and column j struck.
Let E(j) and O(j) be the corresponding sets computed inside the reflection-even and
reflection-odd sector of A + tD, and let −S denote {−s : s ∈ S}. Then

> **Theorem.** u-locus(j) = E(j) ∪ (−O(j)) ∪ {+1, −1} and Δ-locus(j) = E(j) ∪ O(j), at every
> N ≥ 3 and every seat. Hence
>
> - at **odd N**: u-locus(j) = Δ-locus(j) ∪ {+1, −1}, at every seat, every odd N;
> - at **even N**: u-locus(j) = E(j) ∪ {+1, −1} and Δ-locus(j) = E(j) ∪ (−E(j)), so the same
>   relation holds **exactly when** E(j) ∪ {+1, −1} is negation-closed.

The odd-N line is the companion page's §(c) headline, there measured at N = 5..17. The even-N
line is a criterion rather than the relation, and it is that page's gate C8, which measured
*every break is a non-closure* over N = 6..16 and called it a coincidence with no mechanism. It
also names the breaking seats: over N = 6..16 the criterion predicts 22 seats, and they are the
22 the companion gate found.

## What the repo already held, store by store

Swept 2026-09-02 by three agents, one per primitive (the sector reduction, the blindness
criterion, level crossings between symmetry sectors), each over the named stores rather than over
the arc's vocabulary, then corrected by three empty review rounds which moved four claims here
against the author. What was drafted as a new lemma turned out to be owned in folded form; what
was drafted as the hard step turned out to be a corollary; a blanket "nothing" over three stores
turned out to repeat the exact defect `docs/CAUGHT_ERRORS.md` logged for the companion page the
same day; and Lemma 4's original citation turned out not to reach the negative half of its own
axis. All four are recorded rather than quietly absorbed.

- **`docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md`** (F160), Corollary B: **owns the folded
  half of Lemma 1.** *"R commutes with H(u) and H splits into an R-even and an R-odd block.
  Folding the ring along R, each block is a tridiagonal matrix with every off-diagonal entry
  nonzero ... and a diagonal carrying ±u at the crack end"*, with its helper `folded_blocks`
  (whose body is literally `e[0, 0] = u; o[0, 0] = -u`) gated symbolically in u at N = 3..9 as
  gate P9 of
  [`cracked_ring_exact_curve_proof.py`](../../simulations/cracked_ring_exact_curve_proof.py). So
  Lemma 1 is not a new lemma. It is one comparison line: the same fold applied to the
  anisotropic OPEN chain puts the SAME +t in both folded blocks, which is the crack's even block
  at u = +t and its odd block at u = −t. That file's simplicity clause also gives a route to
  Lemma 4, and §(d) explains why this file does not take it.
- **`experiments/THE_SEAT_THAT_CUTS.md` §7**: owns Lemma 2a in all but the phrasing, fence-free
  and for arbitrary real symmetric H. *"By Cramer's rule the (j, j) entry of the adjugate of
  xI − H IS the characteristic polynomial of the principal submatrix, so
  χ(H_cut)/χ(H) = [(xI − H)⁻¹]_{jj} = Σ_λ ‖P_λ e_j‖²/(x − λ) ... An eigenvalue the seat projects
  onto keeps one factor fewer in the numerator than in χ(H); one the seat misses keeps at least
  as many, and the gcd takes the minimum."* That sentence is the per-eigenvalue bookkeeping
  Lemma 2a needs; the lemma is its restatement, not a derivation.
- **`docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md`**: owns the same counting as
  Corollary C, `blind(j) = Σ_λ (m_λ − s_λ)` with `s_λ ∈ {0, 1}`, and owns **Lemma J**, the
  simplicity of an unreduced Jacobi spectrum, which is the whole of Lemma 3. Its own per-λ
  eigenvector-vanishing statement (J3) is fenced to the two-halves Jacobi case, which is why
  Lemma 2a is taken from §7 above and not from here.
- **`docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md`** (F161) §(b), the lemma the companion page
  cites as *"Lemma B"*: the wrap bond in the mode basis, `⟨ψ_k|V|ψ_l⟩ = a_k·a_l·(η_k + η_l)`,
  rank one per sector. That is Lemma 1 seen in the mode basis; it is what let the companion page
  write *"two end-pair objects differing only in that factor"* without taking the last step.
- **`docs/ANALYTICAL_FORMULAS.md`**: F157 owns the Δ axis in closed form (Δ_k, the resultant
  packaging P_j, the live witness `inspect --root blindlocus`). F160's entry owns the u axis
  curve and, correcting a draft of this bullet that said otherwise, **does carry the reflection
  sectors**: *"its simplicity clause (G = 2AB, the two reflection sectors, each a nonvanishing
  prefactor times an unreduced Jacobi block's characteristic polynomial, and a Bézout identity
  forbidding a common zero unless u² = 1)"*, and it carries the fence this file crosses, *"Not
  the blind seat: THE_SEAT_THAT_CUTS's open item asks for a detuned bond under a seat cut, a
  different object"*. F161's entry names the sectors too. What no entry carries is a statement
  relating the two PERTURBATIONS.
- **`compute/MirrorWorld/BlindSeat.cs`**: owns the Δ-locus as a resultant at a named non-centre
  seat, *"at N = 9 seat 1 that resultant factors as `128*D*(D-1)*(D+1)*(D^2-3)`, so the seat is
  blind at Delta = 0, +-1 and +-sqrt(3) and nowhere else"*, and it owns §(a)'s third convention
  too, the centre seat's two principal submatrices being *"conjugate by the chain reflection, so
  they carry the same characteristic polynomial and share every root, at any Delta and any
  palindromic profile (their resultant is identically zero)"*. That is §(a)'s object in the
  small. `Crack.cs` names F160's factorization in passing and does no fold.
- **`compute/RCPsiSquared.Diagnostics/Foundation/DiabolicReflectionParityWitness.cs`**: owns the
  SHAPE of Lemma 5 on a different object, the (SE,DE) coherence block: *"At EVEN nBlock the two
  sectors have equal dimension ... At ODD nBlock the reflection-fixed central site makes
  dim(R-even) − dim(R-odd) = (nBlock−1)/2 ≠ 0; the dimension mismatch forbids the cross-pairing,
  forcing the antiunitarity to act within each sector."* Same two operators in kind, an
  antiunitary and a site reflection, and the same preserve-versus-swap dichotomy by the parity
  of N. Lemma 5 is that dichotomy on the single-excitation chain; it is a second instance, not a
  new shape.
- **`docs/GLOSSARY.md`** holds the blind seat itself, under the headword *The blind seat (F157,
  August 2026)*, with `blind(j) = N − dim Krylov(e_j)` and both gcd laws: §(a)'s vocabulary is
  already written there. What it has no headword for is reflection sector, end pair, crack or
  comb.
- **`reflections/ON_LEAVING_THE_CIRCLE.md`** holds Lemma 5's dichotomy in prose, *"in an even
  chain the two halves are each other's reflection ... An odd chain has a middle seat that
  reflection cannot move"*, about the unmirrorable seat rather than about a locus.
- **`hypotheses/DIABOLIC_BY_INTEGRABILITY.md`** turns on the same Δ axis and records the
  opposite case to §(d)'s, *"No symmetry puts the two modes in different sectors so they cannot
  couple: the site-reflection R commutes but both modes are R = +1"*: a degeneracy INSIDE one
  sector, which is exactly what Lemma 3 excludes here.
- **Nothing** from `recovered/` and from `fw.Confirmations`, both checked entry by entry; the
  nearest confirmation is 24, on the n = 9 comb, a neighbour and not a hit.
- **A near-miss worth naming, because it looks like a hit:** `experiments/CROSSING_TAXONOMY.md`,
  `SUBSYSTEM_CROSSING.md` and every "crossing" row in the glossary are the CΨ boundary crossing,
  not a level crossing. The repo owns no general statement of when levels in different
  reflection sectors cross.

## (a) The object, and three conventions the numbers depend on

Two one-parameter families on the same N sites, both real symmetric, both commuting with the
chain reflection R: l ↦ N − 1 − l:

    H_u(t) = A + t·P        the cracked ring
    H_Δ(t) = A + t·D        the end-pair-anisotropic open chain

The scale is the companion page's, hop 1 and a bare t, which differs from F157's committed
convention (hop 2, 2Δ per end, a Δ(N−5)·I shift) by a common scaling and a common shift. Neither
moves a coincidence between χ(H) and χ(H_j), and the companion gate's C1 checks that against
F157's committed table rather than asserting it.

A **seat** j ∈ {0, …, N−1} is where the watching sits and blind(j) is F157's count, defined in
`docs/GLOSSARY.md`. The **locus** is the real root set of Res_x(χ(H), χ(H_j)) in t, a SET; F157's
count at a given point of that set is a multiplicity and is out of scope here. Three conventions,
because without them not one number below can be reproduced:

1. **How E(j) and O(j) are computed.** Both are the same strike-and-resultant construction
   applied INSIDE a sector: the sector's matrix in the pair basis {e_p + η·e_{N−1−p}}, with the
   row and column at the fold coordinate min(j, N−1−j) struck. A sector vector's amplitude at site
   j is
   ± that pair coordinate, the sign being −1 only in the odd sector at j > N−1−j, and vanishing
   is unaffected by a sign, so the criterion transports.
2. **The sector block is not symmetric.** The pair basis is orthogonal but not orthonormal, so
   the block is the conjugate of a real symmetric matrix by the positive diagonal
   diag(√‖·‖). That conjugation changes neither χ nor which coordinate of an eigenvector
   vanishes, which is all Lemma 2a reads, so Lemma 2a applies to it; but the block itself does
   not satisfy Lemma 2a's stated hypothesis, and saying so is not optional.
3. **The reflection-fixed centre seat of an odd chain is a degenerate case of the definition.**
   There every R-odd eigenvector vanishes (ψ_j = −ψ_j), so O(j) is all of ℝ and
   Res_x(χ(H), χ(H_j)) is identically zero: the locus is ℝ on both axes and is not the root set
   of anything. The Theorem holds there because both sides are ℝ, not because a polynomial was
   compared. The companion page counts 7 such seats in its odd range, the script carries them as an
   explicit `'ALL'` branch, and `BlindSeat.cs` already records the identically vanishing
   resultant.

## (b) Lemma 1, the sector reduction

> **Lemma 1.** On the R-even subspace P acts as D; on the R-odd subspace P acts as −D. Hence
> H_u(t)|even = H_Δ(t)|even and H_u(t)|odd = H_Δ(−t)|odd, as matrices, at every N and every t.

*Proof.* Both P and D commute with R, so both preserve the two sectors and the restrictions are
defined. Let ψ be R-even, so ψ_{N−1} = ψ_0. Then (Pψ)_0 = ψ_{N−1} = ψ_0 = (Dψ)_0 and
(Pψ)_{N−1} = ψ_0 = ψ_{N−1} = (Dψ)_{N−1}, and both vanish at every other site. So Pψ = Dψ. Let ψ
be R-odd, so ψ_{N−1} = −ψ_0; the same two lines give (Pψ)_0 = −ψ_0 = −(Dψ)_0 and
(Pψ)_{N−1} = ψ_0 = −ψ_{N−1} = −(Dψ)_{N−1}, so Pψ = −Dψ. ∎

The content is one substitution, ψ_{N−1} = ±ψ_0, and it is why the two axes can meet at all:
**the cracked ring is, sector by sector, the anisotropic open chain, at +u and at −u.** At odd N
the reflection-fixed site is untouched, neither P nor D having any support there.

Gate L1a checks the matrix identity symbolically in u at N = 3..14 in both sectors, exact zero
matrix. L1b is its anti-vacuity partner on the same code path: the WRONG sign in the odd sector
disagrees at every N, so L1a is not comparing a matrix with itself. L1c checks the mechanism
rather than a difference of values: R commutes with both P and D and does NOT commute with a
one-end diagonal, both directions asserted at every N.

## (c) Lemma 2, the locus decomposition, and where the striking goes

The obstacle the open item ran into is real and must be met head on: F157's count is defined by
striking row and column j, and **striking a non-central seat destroys the reflection**, so
Lemma 1 does not transfer to the struck matrix. It does not have to, because the striking can be
removed from the argument before the sectors are used.

> **Lemma 2a.** For real symmetric H and any seat j, λ is a common root of χ(H) and χ(H_j) if
> and only if some nonzero eigenvector of H at λ vanishes at j.

*Proof.* This is `THE_SEAT_THAT_CUTS` §7 read one eigenvalue at a time, and is repeated here
only because the two directions are used separately below. Write m_λ = dim E_λ and
s_λ = dim span{P_λ e_j} ∈ {0, 1}; the eigenvectors at λ vanishing at j are E_λ ∩ e_j^⊥, of
dimension m_λ − s_λ, because P_λ e_j spans the only direction of E_λ not orthogonal to e_j.

Should such an eigenvector ψ exist, delete its j-th entry to get ψ′ ≠ 0; for i ≠ j,
(H_j ψ′)_i = (Hψ)_i − H_{ij}ψ_j = λψ_i, so λ ∈ spec(H_j) and λ is a common root.

Conversely let λ be a common root and suppose no eigenvector at λ vanishes at j, so m_λ = s_λ,
which with s_λ ≤ 1 and λ ∈ spec(H) forces m_λ = s_λ = 1. Then the residue of the resolvent entry
at λ is ‖P_λ e_j‖² ≠ 0, so by §7's Cramer identity χ(H_j)/χ(H) has a genuine pole at λ: λ is a
root of χ(H) to one order higher than of χ(H_j), and with m_λ = 1 that means λ is not a root of
χ(H_j) at all, contradicting the assumption. ∎

Hermiticity is what the argument uses; real symmetry is what makes spec(H) ⊂ ℝ, so that a real
root set is the whole of the common-root set, and both families here have it.

> **Lemma 2b.** If H commutes with R, write E(j) for the knob values at which some nonzero
> R-even eigenvector vanishes at j, O(j) the same for R-odd, and **C** for the knob values at
> which H has any eigenvalue of multiplicity at least two. Then locus(j) = E(j) ∪ O(j) ∪ C, for
> every seat, and C does not depend on the seat.

*Proof.* For ⊇: E(j) and O(j) lie in locus(j) by Lemma 2a; and if dim ker(H − λ) ≥ 2 then that
kernel, a subspace of ℝ^N of dimension at least two, meets the hyperplane {ψ_j = 0} in a nonzero
vector, so by Lemma 2a again the value lies in locus(j), whatever j is.

For ⊆: let ψ ∈ E_λ be nonzero with ψ_j = 0 and write ψ = ψ⁺ + ψ⁻ along the R-invariant splitting
E_λ = E_λ⁺ ⊕ E_λ⁻. Should one part vanish, the other is a sector eigenvector vanishing at j and
the value lies in E(j) or O(j). Should both be nonzero, dim E_λ ≥ 2 and the value lies in C. ∎

**C carries no seat.** That is the force of the third term: a degeneracy blinds every seat at
once, so it appears in every seat's locus and in none of them for a reason having to do with that
seat. It is also the missing entry in the companion page's taxonomy, which splits a vanishing
mode into *forced* (a reflection fixes the seat) and *accidental* (a sine lands on a multiple of
π) and has no room for this. Stated on its own:

> **Remark.** For any Hermitian H, any λ with m_λ ≥ 2 is blind at EVERY seat.

It is one line from `THE_SEAT_THAT_CUTS` §7 (an eigenvalue of multiplicity m keeps at least
m − 1 factors in the gcd whether the seat projects onto it or not) and the sweep found it drawn
nowhere. It is not the same statement as C: C is the set of knob values where such a λ exists,
and the Remark is why that set enters every seat's locus.

Gate L2a checks the decomposition on both axes, 60 seats each over N = 4..11. Its three
mutations L2b drop C, O and E from the union in turn and are gated against the literal counts
54, 10 and 14 of the 60 crack seats.

## (d) The two degeneracy sets, which is where the axes differ

> **Lemma 3.** C_Δ = ∅ at every N and every real t.

*Proof.* H_Δ(t) is tridiagonal with every off-diagonal entry equal to 1, an unreduced Jacobi
matrix, whose spectrum is simple (Lemma J of `PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA`).
Unreducedness is a fact about the zero pattern and t
never touches an off-diagonal, so this holds at every t, including the |t| at which levels leave
the band. ∎

> **Lemma 4.** C_u = {+1, −1} at every N ≥ 3.

F160's Corollary B offers a route here, through its Bézout identity, and this file deliberately
does not take it: that clause is stated *"for every u ≥ 0 with u ≠ 1"*, and half of this file's u
axis is negative, so the citation would not reach the point (u = −1) at which the statement has
its content. The argument below is self-contained, needs no fold, no band restriction and no
prefactor relation, and holds at every real u.

*Proof.* Fix a value u of the knob and let H_u(u)ψ = λψ. The interior rows are the recursion
ψ_{l−1} + ψ_{l+1} = λψ_l for
l = 1, …, N−2, whose solution is ψ_l = ψ_1·U_{l−1} − ψ_0·U_{l−2} with U_n = U_n(λ/2) the
Chebyshev polynomials of the second kind and the conventional U_{−1} = 0, U_{−2} = −1.
Substituting into the two boundary rows leaves a 2×2 homogeneous system in (ψ_0, ψ_1),

    ⎡ −(λ + u·U_{N−3})     1 + u·U_{N−2} ⎤
    ⎣    U_{N−2} + u         −U_{N−1}    ⎦

and since an eigenvector is determined by (ψ_0, ψ_1), dim ker(H_u − λ) ≥ 2 holds exactly when all
four entries vanish. Now

    u·(1 + u·U_{N−2}) − (U_{N−2} + u)  =  U_{N−2}·(u² − 1)

so a common zero of the two anti-diagonal entries, (0,1) and (1,0), forces u² = 1 or
U_{N−2} = 0. On the branch
U_{N−2} = 0 the same two entries read 1 and u, so the first is never zero and the branch is
empty. Hence C_u ⊆ {+1, −1}.

For the converse, at λ = 2cos(mπ/N) with 1 ≤ m ≤ N−1 one has U_{N−1} = 0 and
U_{N−2} = sin((N−1)mπ/N)/sin(mπ/N) = (−1)^{m+1}, and U_{N−3} = (−1)^{m+1}·λ. All four entries
then vanish for u = −U_{N−2}, that is at u = +1 for even m and at u = −1 for odd m. Both
occur for every N ≥ 3, in counts ⌊(N−1)/2⌋ and ⌈(N−1)/2⌉. ∎

**This is the entire difference between the two axes.** The open chain is Jacobi; the ring is
not, because the wrap entry lifts it out of that form. What the companion page carries as the
unexplained appendix *"∪ {+1, −1}"* is the ring degeneracy, and it appears at every seat because
C carries no seat.

**Which multiplicity, and why the gate may certify it.** Lemma 4 bounds the GEOMETRIC
multiplicity, dim ker, and its forward direction is a polynomial identity that holds over ℂ. The
gates read the ALGEBRAIC multiplicity, off disc_x(χ). For real t both families are real symmetric
and the two coincide, which is what entitles a discriminant to certify a kernel statement; over ℂ
they part, and at N = 3, u = 2√2·i the discriminant vanishes with dim ker = 1, an exceptional
point. So the reality of the knob is named here rather than assumed. A free by-product of the
same 2×2: dim ker(H_u − λ) ≤ 2 at every N and every u, since an eigenvector is determined by
(ψ_0, ψ_1).

Gates C1 and C2 read both degeneracy sets off the same solved object, the real root set of
disc_x(χ(H(t))), at N = 3..14, so the empty answer and the nonempty one are not different code
paths; each is the other's control. B0 is the tie the rest of the block needs: the 2×2 is built
from the Chebyshev ansatz and never touches H_u, so a slip in either boundary row would go
unseen; det(the 2×2) = χ(H_u), exactly and symbolically, pins both rows, the ansatz and the
kernel isomorphism against the real matrix at once. B1 checks the displayed identity, B1b the
empty branch, B3 the converse in the boundary system itself, by exact reduction modulo U_{N−1},
which is also the only gate that reads entry (0,0); and B2 checks the same converse from the
outside as an exact multiplicity count in the discriminant, 2⌊(N−1)/2⌋ at u = +1 and
2⌈(N−1)/2⌉ at u = −1, the expectations combinatorial in N rather than read off the measurement.
B2 gates attainment only; the *exactly* in C_u = {±1} is C2's. C3 and C3b read the Remark on the
ring at N = 5..11 against a generic u = 1/3 on the same code path.

## (e) Parity, which is the only place N's parity enters

Let Σ = diag((−1)^l). Then ΣAΣ = −A, ΣDΣ = D and ΣPΣ = (−1)^{N−1}P, so

    Σ H_Δ(t) Σ = −H_Δ(−t)     at every N
    Σ H_u(t) Σ = −H_u(−t)     at odd N only

which is the companion gate's C5 and C5b, and the second identity is what that page calls the
staggering identity. If Hψ = λψ then H(−t)(Σψ) = −λ(Σψ), and Σ preserves |ψ_j|, so blindness at
seat j at t implies blindness at seat j at −t: **every Δ-locus is negation-closed, at every N.**
The arc `the_forced_and_the_met` already records this as the reason Δ = −1 was not worth probing.

The sector-level statement is the one this file adds, and it is a single commutation:

    Σ R = (−1)^{N−1} R Σ

> **Lemma 5.** At odd N, Σ preserves each reflection sector, so E(j) and O(j) are EACH
> negation-closed. At even N, Σ swaps the sectors, so O(j) = −E(j).

*Proof.* Both from the commutation, applied to the map ψ ↦ Σψ of Lemma 2b's sector
eigenvectors. ∎

The same preserve-versus-swap dichotomy, by the same parity of N and between the same kinds of
operator, is already live in the repo on a different object: see
`DiabolicReflectionParityWitness` in the sweep above.

Gates P1 (odd N ∈ {5, 7, 9, 11, 13}), P2 (even N ∈ {6, …, 16}) and P3, which is P1's fence: at
even N the sectors are NOT separately closed, at 30 of the seats tested, so P1 is a statement
about odd N and not a statement that happens to be checked there.

## (f) The theorem

By Lemma 1 the crack's even sector at u IS the anisotropy's even sector at u, and the crack's
odd sector at u IS the anisotropy's odd sector at −u. Eigenvector vanishing at a site is a
property of the sector alone, so

    E_u(j) = E_Δ(j)        O_u(j) = −O_Δ(j)

and with Lemmas 2b, 3 and 4,

    u-locus(j) = E_Δ(j) ∪ (−O_Δ(j)) ∪ {+1, −1}
    Δ-locus(j) = E_Δ(j) ∪ O_Δ(j)

**Odd N.** By Lemma 5, −O_Δ = O_Δ, so u-locus = E_Δ ∪ O_Δ ∪ {±1} = Δ-locus ∪ {±1}, at every
seat and every odd N. The companion page's §(c) headline is proved, and proved past the range
N = 5..17 in which it was measured.

**Even N.** By Lemma 5, O_Δ = −E_Δ, so u-locus = E_Δ ∪ {±1} and Δ-locus = E_Δ ∪ (−E_Δ). The
relation u-locus = Δ-locus ∪ {±1} is therefore equivalent to −E_Δ ⊆ E_Δ ∪ {±1}, and since {±1}
is itself negation-closed that is equivalent to E_Δ ∪ {±1} being negation-closed, that is to the
**u-locus** being negation-closed. ∎

Two things follow that the companion page could not state.

1. **Its gate C8 is a consequence, not a coincidence.** That page proves one direction (a
   non-closed u-locus cannot equal a negation-closed right-hand side) and has no mechanism for
   the converse. The converse is Lemma 5 at even N.
2. **The breaking seats are named.** Applying the criterion over N = 6..16 predicts 22 seats,
   and they are exactly the 22 the companion gate reports. Gate T4 checks this against that
   gate's list taken as a literal; nothing in this file derives the expectation from a
   measurement.

The ±1 in the criterion is load-bearing rather than decorative. Dropping it predicts 30 seats
instead of 22, the eight extra being (6,1), (6,4), (10,2), (10,7), (12,4), (12,7), (14,3),
(14,10), where E_Δ is not negation-closed on its own and the ring ends restore the closure. Gate
T5 asserts exactly that difference, against the eight seats as a literal. It is the same clause
under which the companion page's gate C9 finds its two exceptions at N = 12 seats 1 and 10, where
the Δ-locus contains ±1 itself.

Gates T1 (a consistency check, see the table), T2 (odd N, 45 seats), T3, T4 and T5.

## (g) The gates, in one table

| gate | what would make it red |
|---|---|
| L1a | the sector reduction failing at one N or one sector; exact zero matrix, symbolic in u |
| L1b | the identity holding with the WRONG odd-sector sign, which would make L1a vacuous |
| L1c | R failing to commute with P or D, or commuting with a one-end diagonal |
| L2a | one seat where locus(j) ≠ E ∪ O ∪ C, either axis |
| L2b ×3 | a dropped summand noticed by a number of seats other than 54, 10, 14 |
| C1 | one N or one Δ at which the anisotropy's spectrum degenerates |
| C2 | the crack's degeneracy set differing from {+1, −1} at one N |
| C3 | one ring seat blind by less than the number of degenerate pairs |
| C3b | every seat blind at a generic u, which would make C3 a property of the code path |
| B0 | det(the 2×2) ≠ χ(H_u) at one N, i.e. the Chebyshev reduction not describing the real matrix |
| B1 | the boundary identity failing at one N, which would break Lemma 4's ⊆ direction |
| B1b | the U_{N−2} = 0 branch turning out to be nonempty |
| B2 | a ring end not attained, or attained with a multiplicity other than the predicted one. It gates ATTAINMENT only; the *exactly* in C_u = {±1} is C2's |
| B3 | one of the four entries failing to vanish at u = −U_{N−2}; the only gate that reads entry (0,0) |
| P1 | one odd-N seat whose sector locus is not negation-closed |
| P2 | one even-N seat where O ≠ −E |
| P3 | the even-N sectors turning out separately closed after all |
| T1 | **a consistency check throughout, not independent evidence.** Its even leg cannot fail at all: E_u = E_Δ is L1a, the two blocks being the same matrix. Its odd leg tests that the locus of a family at −t is the negation of its locus, i.e. that `real_factors` and `negate` commute with t ↦ −t |
| T2 | one odd-N seat where u-locus ≠ Δ-locus ∪ {±1} |
| T3 | one even-N seat where the two-term forms fail |
| T4 | a predicted break set differing from the committed 22 |
| T5 | the ring-ends clause turning out not to change the criterion, or changing it at seats other than the eight named |

Three mutations of the objects were run by hand against an earlier build of this gate and all
three went red: the wrong odd-sector sign, the crack fed through C1's door, and the criterion
without the ring-ends clause. None of the three is a hand mutation any longer. The first is gate
L1b, the second is what C2 asserts through C1's own door, and the third is gate T5.

## (h) Scope and fences

- **Multiplicity is out of scope.** Every locus here is a SET. F157's count at a point of the
  Δ-locus is a multiplicity, load-bearing in F157 and untouched here. The theorem says the two
  axes carry the same points; it does not say a seat is blind by the same AMOUNT at a shared
  point, and no gate below tests that. The companion page records why its construction cannot
  see it: the resultant it forms is a rational multiple of F157's P_j squared.
- **F157's standing fence is not lifted here.** *"A Δ is NOT the detuned bond that The Seat That
  Cuts leaves open; do not report one as the other."* Lemma 1 identifies no operators: P is a
  bond and D a diagonal, they sit in different positions on the full space, and what coincides
  is their RESTRICTIONS to each reflection sector, at opposite signs of the knob in the odd one.
  Everything here reports loci, never one perturbation as the other.
- **The XY book only.** Every H here is the single-excitation XY block, F160's book. F157's
  Heisenberg law is on modulus N, a different comb, and is not on this road.
- **The theorem is about the two named families**, the wrap bond and the end-pair diagonal. The
  companion page measured four interior perturbations and found none that carries the whole
  Δ-locus, the best reaching 44 of 48 seats by containment. Which interior perturbations do, and
  at which seats, is untouched here. Lemma 1 says why the end pair is the one the crack meets,
  and says nothing about why a symmetric interior diagonal pair comes close.
- **"The end pair is special" is still refuted** for the locus, as the companion page found, and
  nothing here reinstates it. Lemma 1 is about the crack and the anisotropy sharing an end pair,
  not about the end pair being the only site pair that can carry a locus.
- **What is this file's own**: the comparison line in §(b), Lemma 2b, Lemma 4, Lemma 5's sector
  clause, and the Remark in §(c). Lemma 1's folded half is F160's, Lemma 2a is
  `THE_SEAT_THAT_CUTS` §7's, Lemma 3 is Lemma J of the node-lemma proof, and the Σ identities are
  the companion page's.
- **Two ordinary words here are typed objects elsewhere, and neither is meant.** A **pair** on
  this page is the site pair {p, N−1−p} the reflection joins, never
  [`Pair.cs`](../../compute/MirrorWorld/Pair.cs)'s bare coherence |i⟩⟨j| with its rate −2γk; the
  **end pair** is F157's {0, N−1}, which the companion page already fences against F140's
  "corner block". A **block** here is the matrix of H restricted to one reflection sector, never
  [`Block.cs`](../../compute/MirrorWorld/Block.cs)'s joint-popcount block (p, q). No γ and no
  Liouvillian appears anywhere in this file.
- **Bare-letter lemma names are ambiguous in this neighbourhood** and this file numbers its own
  for that reason. Both cited proofs have a Corollary B and a Corollary C, of different objects,
  and `PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md` uses "Lemma B" for a third. Every citation here
  names the file.
- **No hardware claim.** Nothing here is proposed for a flight.

## Where it came from

The companion page landed on 2026-09-02 with §(c) as an observation and open item 1 asking for a
proof of the shape u-locus = (even-sector Δ-locus) ∪ −(odd-sector Δ-locus). That shape is correct
and is Lemma 1 plus Lemma 2b. What the item did not anticipate is that the same lemmas settle the
even-N side as well, converting its own gate C8 from a measured biconditional into a consequence,
and that the appendix {+1, −1} is not an appendix but the third term of the decomposition.

Related: [The Blind Seat on the Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md) ·
[PROOF_CRACKED_RING_EXACT_CURVE](PROOF_CRACKED_RING_EXACT_CURVE.md) ·
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) ·
[PROOF_COLLISION_GAP_ODD_ORDERS](PROOF_COLLISION_GAP_ODD_ORDERS.md) ·
[The Seat That Cuts](../../experiments/THE_SEAT_THAT_CUTS.md) ·
[The Cracked Bell](../../experiments/THE_CRACKED_BELL.md)
