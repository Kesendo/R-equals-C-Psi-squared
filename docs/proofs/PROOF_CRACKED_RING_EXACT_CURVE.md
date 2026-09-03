# PROOF: The Cracked Ring Is Exactly Solvable, and Its Curve Is the Chain's Amplitude Closed

**Registry:** [F160](../ANALYTICAL_FORMULAS.md).
**Status:** Tier 1 derived. Theorem A is a polynomial identity in x and u and holds for
every N ≥ 3 and every complex u; Corollaries B and C and Theorems D, E, F are read off it,
with three matrix arguments added (Perron-Frobenius, the reflection sectors and Weyl's
inequality), and
Theorem G has a second, eigenvector route beside the polynomial one. Corollary B (the
curve) needs 0 ≤ u < 1 for the sentence *the whole spectrum lies on the open interval
(0, π)*, which is Perron-Frobenius on a nonnegative matrix. Its simplicity clause is the
factorization G = 2AB into the two reflection sectors and holds for every u ≥ 0 with
u ≠ 1. Past u = 1 the same curve continues to complex k. Corollary C (the two ends) is
exact at u = 1 and u = 0. Theorem D (the join) is stated for u > 0, where the amplitude t
is defined. Theorem E (the split's next order) is a δ → 0 statement and says nothing at
a finite δ beyond its own O(δ²). Theorem F's three factor forms are polynomial identities
in u; its departure COUNT is stated for u ≥ 0, the object's range (at negative u the count
reads differently and is not this file's), and reports the edge, where a level sits ON the
band edge, separately rather than through a tolerance. Theorem G (the velocity at the
chain end) is a first-order statement at u = 0. The gates are of three kinds: exact symbolic
identities (sympy: the multiple-angle identities at concrete N by reduction modulo
s² + c² − 1, the symbolic-N identities and the series elimination by the symbolic
simplifier, exact either way, no floating point anywhere in the script), exact integer
routes in the two C# test suites
(Bareiss, Descartes, GF(p)), and the experiment's own stage-E gates, most of them on the
curve alone and some through an eigensolver, which this file cites and does not repeat.
**Date:** 2026-09-02 (the law derived and gated 2026-08-31 on the experiment page; this
file written the day after the F160 mint)
**Authors:** Thomas Wicht, Claude (Fable 5.1)
**Script:** [`simulations/cracked_ring_exact_curve_proof.py`](../../simulations/cracked_ring_exact_curve_proof.py)
(nine symbolic gates P1 to P9, exact in sympy; run committed at
[`cracked_ring_exact_curve_proof_run.txt`](../../simulations/results/cracked_ring_exact_curve/cracked_ring_exact_curve_proof_run.txt)).
The eigensolver gates are stage E of
[`simulations/cracked_bell_gate.py`](../../simulations/cracked_bell_gate.py), the
exact-integer meeting of Theorem A is `CrackTests` on
[`compute/MirrorWorld/Crack.cs`](../../compute/MirrorWorld/Crack.cs), and the typed
claim's own from-below tests are `CrackedRingExactCurveClaimTests` in
`compute/RCPsiSquared.Core.Tests/`.
**Builds on:**
- [The Cracked Bell](../../experiments/THE_CRACKED_BELL.md) §The crack is exactly
  solvable, where the curve was found from the matching condition and gated through an
  eigensolver, and [Coupling Defect Walk-Time Step](../../experiments/COUPLING_DEFECT_WALK_TIME_STEP.md),
  which has carried the transmission amplitude of the same bond since 2026-07-12 and
  since 2026-08-31 states the join. This file supplies the determinant route, which the
  pages did not have, and the derivations the registry entry could only name.
- [PROOF_F139_SEAM_IDENTITY](PROOF_F139_SEAM_IDENTITY.md) for the Chebyshev
  polynomial of the second kind in the 2cos normalization, S_m(2cos θ) = sin((m+1)θ)/sin θ,
  which is the U_n(x/2) of this file under another name. Cited, not re-derived.
- [PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md)
  Lemma J for the unreduced Jacobi matrix (every off-diagonal entry nonzero, hence a
  simple spectrum), which is the u = 0 end of the road, and
  [PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE](PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md)
  §(b) Lemma A (the non-derogatory tridiagonal, and the house precedent of writing a
  classical tridiagonal fact out anyway, with the reason given) and §(g) (Hermitian
  simplicity).
- [PROOF_RING_GAP_DOMINANCE](PROOF_RING_GAP_DOMINANCE.md) for the uniform ring, whose
  §Scope fence (the wrap bond detuned by 1e-4) is this curve's Perron root, and
  [PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md) for the rate of the coherences
  this spectrum labels, −2γ, which this file does not touch.
- [F2b](../ANALYTICAL_FORMULAS.md) for the open chain's comb and its two virtual
  sites, [F65](../ANALYTICAL_FORMULAS.md) for the endpoint rate comb that Theorem G
  meets, [F129](../ANALYTICAL_FORMULAS.md) for the comb law that Theorem G is then
  applied to in [The Comb on the Road](../../experiments/THE_COMB_ON_THE_ROAD.md).

## What the repo already held, store by store

The sweep was run on 2026-09-02 by three agents over the markdown layer, the typed
layer and the comb laws, starting from the record the registry entry carried the
evening before and adding what a record auditor then found missing from it
(F123, PROOF_CHAIN_GAP_DOMINANCE §4.3, PROOF_DEPHASING_FRONT_RENEWAL).

- **[`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).** F160 holds the law as
  an INDEX of the two experiment pages and says so in its header; its Proof field read
  *"none in docs/proofs/ … Open: a PROOF_ home"* and its Typed field *"not yet"*, and
  both are the reason this file exists. F2b is the u = 0 end in words (the two virtual
  sites i = −1 and i = N); F65 is the rate comb α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)) that Theorem G
  finds inside the polynomial's derivative; F122 names the ring's k ↔ −k degeneracy and
  never lifts it; F157 owns the exact-integer genre (a locus polynomial over ℤ by
  Bareiss, counts by GF(p) rank) that the adopted `Crack` runs on this object; F129 is
  the comb law whose collisions Theorem G moves; F123 is Theorem G's Re-side sibling in
  the MOVE and not in the form, the same Hellmann-Feynman derivative in one bond's knob,
  there giving the squared difference of the mode's site occupations across that bond,
  (n(j) − n(j+1))², which vanishes at the ends where Theorem G's signed endpoint product
  lives; F86c's per-bond Hellmann-Feynman response K_b (PROOF_F100, its Q_peak witness) is the same derivative
  on a Liouvillian observable.
- **[`docs/proofs/`](.).** No proof holds the cyclic tridiagonal determinant; at HEAD `1a65a9c` *cyclic
  tridiagonal* appeared in no proof file, *wrap bond* only in PROOF_RING_GAP_DOMINANCE (its
  two Jordan-Wigner parity remarks and the §Scope detuning this file cites), and *transfer
  matrix* only in two unrelated senses (INCOMPLETENESS_PROOF, PROOF_SUBSYSTEM_CROSSING).
  PROOF_CHAIN_GAP_DOMINANCE §4.3 holds the chain's own single-bond fence with the same
  O(δ^1.999) drift power on the other graph. [PROOF_DEPHASING_FRONT_RENEWAL](PROOF_DEPHASING_FRONT_RENEWAL.md)
  stands ON this object since 2026-08-31 (its Status line: verified *"on the cracked
  ring with an off-diagonal observable"*), a proof about the renewal ladder and not about
  the curve.
  What the proofs do hold, and this file cites: PROOF_F139_SEAM_IDENTITY's S_m (the
  same polynomial, the same normalization, Lemma 5 there is a Chebyshev fold);
  PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA's Lemma J on unreduced Jacobi matrices;
  PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE's Lemma A (the non-derogatory tridiagonal) and
  its §(g) (Hermitian simplicity);
  PROOF_F89_PATH_D_CLOSED_FORM's Chebyshev-U sums (a different use of the same family);
  PROOF_HANDSHAKE_TRANSITION_INVARIANT reading one bond's matrix element as a scattering
  amplitude, which is the frame Theorem D's join sits in. Faddeev-LeVerrier, the route
  `Crack` takes to the left side of Theorem A, appears in no proof; it lives in F160's own
  identity paragraph and on experiment pages (THE_BLIND_SITE gate G1, THE_SEAT_THAT_CUTS,
  THE_SPREAD_IS_A_RESONANCE, two F89 pages); this file's route is the Laplace expansion
  instead, with the matrix determinant lemma and the Cassini identity as the cross-check,
  and the Cassini step is the quietest new thing in this file: it appears in no other
  proof, and until this change nowhere in the repo. [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) Lemma 1 owns
  KHK = −H for bipartite nearest-neighbour hopping, non-uniform J included, which is the
  identity the comb page's mirror-pair argument stands on (K·H(u)·K = −H(−u) at odd N is
  that lemma applied to the crack, whose wrap bond K leaves alone).
- **[`experiments/`](../../experiments/).** The two owners in the header. THE_BLIND_SITE
  §5 holds the rank-one member of the same determinant family, det(λI − H + s·E_jj) =
  χ_H + s·χ_cut, the crack being the rank-two member. XY_FROZEN_BAND names the wrap bond
  as a pricing deformation. No null result touches the object.
- **[`hypotheses/`](../../hypotheses/) and `simulations/`.** PERSPECTIVAL_TIME_FIELD §3.2
  (Tier 2) carries Theorem G's eigenvector route one bond over, *"A = 2 ψ_k(b) ψ_k(b+1) is
  the (real) bond overlap of the single-excitation mode at the defect bond (b, b+1)."*, and
  two committed scripts state the same Hellmann-Feynman element on the chain's sine modes in
  the repo's own notation, `simulations/handshake_rk_first_principles.py` (*"d eps_1/dJ_b =
  <psi_1|V_b|psi_1> = 2 psi_1(b)psi_1(b+1)"*) and `simulations/handshake_rk_block.py`. What
  Theorem G adds is the closed form at the wrap bond, the product of the two END amplitudes
  being the signed F65 rate. A first version of this record listed `hypotheses/` under
  *nothing found* and did not open `simulations/` for the move at all.
- **The OpenArcs registry** (`compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs`).
  No cracked or warble-shaped arc. `the_forced_and_the_met` NextStep (1) asks for the
  comb test on exactly this curve's u = 0 end and recorded, until this change, that none
  of the five comb laws had been run on the road; that item is taken up in
  [The Comb on the Road](../../experiments/THE_COMB_ON_THE_ROAD.md), not here, and the
  same change rewrites the arc's item.
- **[`docs/GLOSSARY.md`](../GLOSSARY.md).** No entry for the object, none for *comb*,
  none for *crack*; *Band* is spent in a pre-registration sense that this file's *band*
  (the cosine band [−2J, 2J]) is not.
- **[`docs/CAUGHT_ERRORS.md`](../CAUGHT_ERRORS.md).** Four of the eight entries of
  2026-08-31 and 2026-09-01 are this arc's: the THE_CRACKED_BELL landing round and the
  exactness landing round (the departure count's first draft forgot the parity of N: *a
  count is a claim*) on the first day; *five review rounds on one bundle* (rounds four to
  eight of the same page) and the F160 mint with the Crack adoption (one entry carrying a
  continued curve that cancelled at κ = ln u, a test row that certified the row u = 1e5,
  which returned 1e5.0107, as fine, and a limit stated in the wrong variable) on the second.
  Of the other four, the PROOF_RING_GAP_DOMINANCE header round of 2026-08-31 is this
  neighbourhood's (the proof whose §Scope fence this curve's Perron root sits in), and the
  caveat-outlived-its-gate entry of the same day, the structural-ceiling verifier and the
  F157 letter round of 2026-09-01 belong to other arcs. Theorem F below is written so that the parity is in
  the statement and not in a remark. The entry of 2026-09-02 is this file's own round.
- **`reflections/`.** ON_THE_REFUND reflects on the sibling page and the renewal proof
  (*every observation is a rebirth*), the outward reading of the object this curve sits on;
  nothing there about the curve.
- **Nothing found:** `fw.Confirmations` and the C# `ConfirmationsRegistry` (nothing
  ring-shaped, nothing defect-shaped; the F129 flight on the comb belongs to the comb
  page's sweep), `recovered/`.

**What is new here, stated narrowly.** The algebra is standard tight-binding
scattering, and the experiment page fences it so. What this file adds to the repo is
the determinant route to a curve the pages already had (Theorem A; the multiplicity
clause was on the page as gate E1c, the polynomial form in F160's entry and `Crack.cs`),
with the Cassini identity as its cross-check, the SIMPLICITY of the spectrum for every
u ≥ 0 except u = 1 (Corollary B's clause, a theorem where the page had a sweep, and the
one this file's neighbours lean on), the written-out derivations of the split's
next order and of the departure count, and Theorem G, which was not on any page: the
first-order motion of the chain comb along the road is F65's rate comb with a sign.

## (a) The object

N ≥ 3 sites on a ring, the single-excitation block of the XY Hamiltonian in the
adjacency book of PROOF_RING_GAP_DOMINANCE (H_se = J·A, NOT the Heisenberg Laplacian of
D10 and not the (1,1) Haken-Strobl block `Cone` runs). Every bond carries J except the
wrap bond N−1 ↔ 0, which carries J′ = u·J, u ≥ 0. Units of J throughout. Two
conventions for the knob are committed and both are kept: the crack WEAKENS, u = 1 − δ
(`Warble.cs`, the bell); the walk-time step STRENGTHENS, u = 1 + δ (`WalkTime.cs`, the
sibling). Every statement below is in u.

Notation: P(x) = det(x·I − H) is the characteristic polynomial of the cracked ring, monic
of degree N. D_n(x) is the determinant of the n × n path pencil (x on the diagonal, −1 on
the two off-diagonals), with D_0 = 1, D_1 = x and the recursion D_n = x·D_{n−1} − D_{n−2}.
It is U_n(x/2), the Chebyshev polynomial of the second kind, monic of degree n with
integer coefficients, and at x = 2cos k it is sin((n+1)k)/sin k (PROOF_F139's S_m, with
m for n). T_n is the first kind, 2T_n(x/2) = 2cos(nk) at x = 2cos k, used once, at the
ring end.

## (b) Theorem A, the polynomial

**Theorem A.** det(x·I − H) = U_N(x/2) − u²·U_{N−2}(x/2) − 2u.

*Proof.* Write M = x·I − H. Row 0 of M is (x, −1, 0, …, 0, −u). Expand along it. Three
minors appear.

The (0,0) minor deletes row 0 and column 0 and with them both corner entries; what is
left is the path pencil on sites 1..N−1: D_{N−1}.

The (0,1) minor deletes row 0 and column 1. Its first column is (−1, 0, …, 0, −u)ᵀ
(the entries M_{1,0} and M_{N−1,0}); expand along it. The −1 in the first position has
as cofactor the path pencil on sites 2..N−1, D_{N−2}, with sign +. The −u in the last
position has as cofactor the matrix of rows 1..N−2 against columns 2..N−1, which is
lower triangular with every diagonal entry −1 (row i meets column i+1 at M_{i,i+1} = −1,
and nothing above it survives), so its determinant is (−1)^{N−2}, and the sign of the
position is (−1)^{N−2} as well. Hence the (0,1) minor is −D_{N−2} − u.

The (0,N−1) minor deletes row 0 and column N−1. Its last row is (−u, 0, …, 0, −1) (the
entries M_{N−1,0} and M_{N−1,N−2}); expand along it. The −1 has as cofactor the matrix
of rows 1..N−2 against columns 0..N−3, upper triangular with diagonal −1, determinant
(−1)^{N−2}, sign +. The −u has as cofactor rows 1..N−2 against columns 1..N−2, the path
pencil D_{N−2}, with sign (−1)^{N−2}. Hence the (0,N−1) minor is (−1)^{N−1}·(1 + u·D_{N−2}).

Assemble with the row-0 cofactor signs (−1)^{0+j}: det M = x·D_{N−1} + (−1)¹·(−1)·(−D_{N−2} − u) +
(−1)^{N−1}·(−u)·(−1)^{N−1}·(1 + u·D_{N−2}) = x·D_{N−1} + (−D_{N−2} − u) − u·(1 + u·D_{N−2})
= x·D_{N−1} − D_{N−2} − 2u − u²·D_{N−2}, and x·D_{N−1} − D_{N−2} = D_N. ∎

*Cross-check by a second route.* The crack is the rank-two update M = M_path − u·(e_0e_{N−1}ᵀ
+ e_{N−1}e_0ᵀ). By the matrix determinant lemma, det M = D_N·det(I₂ − u·W), where W is the
2 × 2 block of M_path⁻¹ on the sites 0 and N−1 read crosswise: its diagonal entries are
(M_path⁻¹)_{N−1,0} = (M_path⁻¹)_{0,N−1} = 1/D_N (the corner of the inverse of a tridiagonal
matrix is the product of the off-diagonal entries over the determinant, here
(−1)^{N−1}·(−1)^{N−1}/D_N), its off-diagonal entries are (M_path⁻¹)_{00} =
(M_path⁻¹)_{N−1,N−1} = D_{N−1}/D_N. So det M = [(D_N − u)² − u²·D_{N−1}²]/D_N, and the
Cassini identity of the second-kind polynomials, D_{N−1}² − D_N·D_{N−2} = 1 (gate P2), turns
the bracket into D_N·(D_N − 2u − u²·D_{N−2}). The same polynomial, by a route that never
touches a minor. ∎

Gate P1 meets Theorem A symbolically for N = 3..9, as polynomials in (x, u). `Crack.cs`
meets it over the integers (the left side by Faddeev-LeVerrier on the matrix scaled by
the denominator of u, coefficient for coefficient, residual compared as 0) at twelve
(N, u) and mod two primes to N = 1001. The typed claim's own test meets it against a
fraction-free Bareiss determinant at eight (N, u) and nine integer x each.

## (c) Corollary B, the curve, and why the zero set is the whole spectrum

**Corollary B.** At x = 2cos k,

    det(2J·cos k·I − H) = J^N·G(k)/sin k,
    G(k) = sin((N+1)k) − u²·sin((N−1)k) − 2u·sin k                      (the SUM form)
         = (1 − u²)·sin(Nk)·cos k + [(1 + u²)·cos(Nk) − 2u]·sin k          (the sin(Nk) form)
         = 2·A(k)·B(k),    A(k) = cos(k(N+1)/2) − u·cos(k(N−1)/2),    B(k) = sin(k(N+1)/2) + u·sin(k(N−1)/2)   (the SECTOR form).

The three forms are named, not numbered, because the gate script and the typed claim print
them in the other order.

For 0 ≤ u < 1 the spectrum of H is the zero set of G on the OPEN interval (0, π), with
multiplicities; and for every u ≥ 0 with u ≠ 1 the spectrum is SIMPLE, so for 0 ≤ u < 1
the zero set consists of exactly N simple zeros.

*Proof.* Substitute D_n(2cos k) = sin((n+1)k)/sin k into Theorem A, which gives the sum
form; the sin(Nk) form follows from sin((N±1)k) = sin(Nk)·cos k ± cos(Nk)·sin k. For 0 ≤ u < 1 the weighted
adjacency is nonnegative and irreducible with row sums 2 except at the two crack sites,
where they are 1 + u < 2, and an irreducible nonnegative matrix reaches its maximal row
sum as spectral radius only when all row sums are equal (Perron-Frobenius), so every
eigenvalue lies strictly inside (−2, 2), that is at x = 2cos k for a unique k ∈ (0, π).
On that interval sin k is analytic and nonzero and dx/dk = −2·sin k is nonzero, so a root
of P of multiplicity r at x₀ = 2cos k₀ is a zero of G of multiplicity r at k₀ and
conversely. G vanishes at k = 0 and k = π for free (every term of the sum form and of the
sin(Nk) form does), which is why the interval is open and why the two edges need Theorem F.

*The simplicity clause.* The sector form is the product-to-sum identity read
backwards: with a = k(N+1)/2 and b = k(N−1)/2, 2AB = sin 2a − u²·sin 2b + 2u·sin(b − a) =
sin((N+1)k) − u²·sin((N−1)k) − 2u·sin k (gate P8, symbolic in N). The two factors are the
two reflection sectors, and each is a nonvanishing prefactor times a characteristic
polynomial. The reflection R: j ↦ N−1−j through the crack bond maps every chain bond to
a chain bond and the wrap bond to itself, so R commutes with H(u) and H splits into an
R-even and an R-odd block. Folding the ring along R, each block is a tridiagonal matrix
with every off-diagonal entry nonzero (the chain bonds, 1; at odd N the middle site
enters the even block through a √2) and a diagonal carrying ±u at the crack end and ±1
at the fold end (even N) or nothing at the fold end (odd N): an unreduced Jacobi matrix,
whose spectrum is simple (the node lemma of the blind-seat proof, Lemma J, or the
edge-block proof, §(b) for non-derogatory and §(g) for Hermitian simplicity). Writing
χ_e and χ_o for the two blocks' monic characteristic polynomials,

    even N:   A(k) = cos(k/2)·χ_e(2cos k),        B(k) = sin(k/2)·χ_o(2cos k),
    odd N:    A(k) = ½·χ_e(2cos k),               B(k) = sin k·χ_o(2cos k),

consistent with 2AB = sin k·χ_e·χ_o = sin k·P in both parities. Each prefactor is nonzero
on (0, π), so the zeros of A there are the roots of χ_e WITH MULTIPLICITY and are simple
by Lemma J, and likewise for B. A double zero of G on (0, π) therefore needs
A(k₀) = B(k₀) = 0, and the Bézout identity

    (cos a + u·cos b)·A(k) + (sin a − u·sin b)·B(k) = 1 − u²

(gate P8, symbolic in N) makes that impossible unless u² = 1; the prefactor relations
themselves are gate P9, the folded blocks built symbolically for N = 3..9 with u free. So
for u ≥ 0, u ≠ 1, every
zero of G on (0, π) is simple, and for 0 ≤ u < 1, where Perron-Frobenius puts all N
levels inside, G has exactly N simple zeros there. Past u = 1 the levels outside the band
are simple as well, by Theorem F's Weyl bound (at most one level beyond each edge, and
the same bound keeps a level sitting exactly ON an edge, the odd-N threshold case, apart
from the level next to it), and the in-band ones by the argument just given, so the
spectrum is simple for every u ≥ 0 with u ≠ 1. At u = 1 the pairs m ↔ N−m coincide (the
ring). ∎

Until this file the simplicity was a measurement: the experiment's gate E2d(i) sweeps it
over 1050 (N, δ) points and said so (*"measured rather than proved"*, a clause the same change
now reads as history), and the page now
reads *"nothing here proved it for every N and δ"*. The factorization proves it, and it
is what makes the labelling of the road's levels by descending energy the same as their
continuation from u = 0, which [The Comb on the Road](../../experiments/THE_COMB_ON_THE_ROAD.md)
rests on. The factorization was found by a review lens on 2026-09-02, not by the writer,
and the prefactor step that makes it a proof by a second lens the same day: the first
written version inferred simplicity of A's zeros from their coinciding with a reflection
sector's eigenvalues, which says nothing about their order.

At u = 1 the top level sits at k = 0 (and at even N the bottom one at k = π), where G
vanishes for free: the ring's edge levels are counted by Theorem F, not by the interval.
Past u = 1 the departed levels ride the same polynomial at |x| > 2, that is at
k = iκ and k = π + iκ, and `Crack.cs` reads them off the curve continued in a grouped
form with no cancellation (a first draft cancelled at κ = ln u; CAUGHT_ERRORS 2026-09-01).

Gate P3 meets Corollary B symbolically in the sum form and the sin(Nk) form for N = 3..9,
gate P8 the sector form and the Bézout identity, gate P9 the prefactor relations. The
experiment's gates meet the curve by three routes: E1a through an eigensolver (one Newton
step along G moves an eigenvalue by at most 1.8e-14 over 48 (N, δ) points), E1b on the
curve alone (exactly N sign changes of G on (0, π): for a monic real-rooted polynomial of
degree N with all roots inside, that count alone forces N simple roots at each sampled
point, since a double root contributes no sign change), E1c by a determinant (the identity
to 4.1e-13 relative, the curve pinned to the matrix). What the simplicity clause
adds is every N and every u ≥ 0 with u ≠ 1, not a blind spot closed.

## (d) Corollary C, the two ends, and what u is

**Corollary C.** P(x; 1) = 2T_N(x/2) − 2, whose roots are the ring comb 2cos(2πm/N),
m = 0..N−1 (the pairs m ↔ N−m coinciding), and P(x; 0) = U_N(x/2), whose roots are the open
chain comb 2cos(πm/(N+1)), m = 1..N, F2b.

At u = 1 the polynomial is 2T_N(x/2) − 2 and the curve reads cos(Nk) = 1, the perfect
ring's comb k = 2πm/N with its m ↔ N−m pairs. At u = 0 the polynomial is U_N(x/2) and the
curve reads sin((N+1)k) = 0, the OPEN N-site chain's comb k = πm/(N+1), which is F2b read
as the matching condition ψ_{−1} = u·ψ_{N−1}, ψ_N = u·ψ_0 with u set to zero (the two
virtual sites of F2b). So u interpolates the ring's modulus N and the chain's modulus
N+1, and it is a BOUNDARY CONDITION: for every u > 0 the graph is still a ring, only the
endpoint u = 0 is a chain. This file does not use the word *topology* for u; the repo
spends it twice already (a discrete choice of graph, and the band invariant of
TOPOLOGICAL_EDGE_MODES, a committed negative result).

Both ends are met in `CrackTests` as two different recursions (2T_N − 2 and U_N) against
the matrix; the typed claim's tests meet the ring end against the first-kind recursion and
the chain end against a Bareiss determinant of the path pencil, a route that shares no code
with the polynomial.

## (e) Theorem D, the join

**Theorem D.** Let t(k) = −2iu·sin k/(e^{−ik} − u²·e^{ik}) be the transmission amplitude of
the same bond on the infinite chain (the sibling's τ(q), carried since 2026-07-12). Then

    1/t(k) = (1 + u²)/(2u) + i·(1 − u²)·cot k/(2u),        Re[e^{−iNk}/t(k)] − 1 = G(k)/(2u·sin k),

so the quantization condition of the ring is Re[e^{−iNk}/t(k)] = 1: the chain's
scatterer with one round trip of phase.

*Proof.* e^{−ik} − u²e^{ik} = (1 − u²)·cos k − i(1 + u²)·sin k; dividing by −2iu·sin k
gives the first line. Then Re[e^{−iNk}/t] = cos(Nk)·(1 + u²)/(2u) + sin(Nk)·(1 − u²)·cot k/(2u),
and subtracting 1 and putting everything over 2u·sin k gives G's sin(Nk) form. ∎

Two different chains meet in that sentence and stay apart: the infinite one whose
leads define t, and the open N-site one the same equation reaches at u = 0, where t
itself vanishes and the condition is to be read after multiplying through by 2u·sin k,
as G(k) = 0. What the
sibling sets aside as *"the upstream O(δ) reflected wave … is not read as signal"* is on
the closed ring the whole signal. Gate P4 meets both lines (the first symbolically, the second by the same exact reduction, with the Re step as an identity in real A, B, N = 3..9); the experiment's gate E4
meets them through the eigensolver (1.4e-14 and 7.9e-10).

## (f) Theorem E, the split and its next order

**Theorem E.** Write u = 1 − δ, k_m = 2πm/N with 1 ≤ m < N/2, s = sin k_m, c = cos k_m.
The two levels the crack makes of the ring's degenerate pair m ↔ N−m are separated by

    ΔE_m = (4δJ/N)·[1 + δ·c_m + O(δ²)],        c_m = ½ − 1/(N·s²).

*Proof.* Put Nk = 2πm + x, so that sin(Nk) = sin x and cos(Nk) = cos x, and use the
sum form: sin((N+1)k) = sin(x + k), sin((N−1)k) = sin(x − k), so
G = sin(x + k) − u²·sin(x − k) − 2u·sin k. With sin(x + k) − sin(x − k) = 2·cos x·sin k
this is G = 2·sin k·(cos x − u) + (1 − u²)·sin(x − k), and G = 0 reads

    cos x − u = (1 − u²)·sin(k − x)/(2·sin k),        k = k_m + x/N.

Left side: δ − x²/2 + O(x⁴). Right side: δ(1 − δ/2)·R(x) with R(x) = sin(k_m − x(1 − 1/N))/
sin(k_m + x/N). Expanding R to second order in x, with a = 1 − 1/N and b = 1/N,

    R(x) = 1 − (c/s)·x + [(b² − a²)/2 + (c/s)²·b(a + b)]·x² + O(x³)
         = 1 − (c/s)·x + [1/N − ½ + c²/(N·s²)]·x² + O(x³),

the N-dependence cancelling at first order because a + b = 1. At order δ² the equation is
x² − 2δ·(c/s)·x − δ² = 0, whose roots x_± = δ·(c ± 1)/s are the committed branch shifts
(the flat split, since E = 2cos(k_m + x/N) = 2c − 2s·x/N + …, gives E_− − E_+ =
−(2s/N)(x_− − x_+) = 4δ/N at this order). At order δ³, writing x = a₁δ + a₂δ² and
K = 1/N − ½ + c²/(N·s²),

    a₂·(2a₁ − 2c/s) = −(c/s)·a₁ − 2K·a₁²,        so        a₂,± = ∓(s/2)·[(c/s)·a₁,± + 2K·a₁,±²].

Then, with E = 2c − 2s·(x/N) − c·(x/N)² + O(x³),

    ΔE = E_− − E_+ = −(2s/N)·(x_− − x_+) − (c/N²)·(x_−² − x_+²)
       = 4δ/N − (2s/N)·δ²·(a₂,− − a₂,+) + 4c²δ²/(N²s²) + O(δ³),

and a₂,− − a₂,+ = (s/2)·[(c/s)(a₁,− + a₁,+) + 2K(a₁,−² + a₁,+²)] = [c² + 2K(c² + 1)]/s, using
a₁,− + a₁,+ = 2c/s and a₁,−² + a₁,+² = 2(c² + 1)/s². Dividing by 4δ/N,

    ΔE/(4δ/N) = 1 − (δ/2)·[c² + 2K(c² + 1)] + δ·c²/(N·s²),

and c² + 2K(c² + 1) = c² + (c² + 1)(2/N − 1) + 2c²(c² + 1)/(N·s²) = −1 + (2(c² + 1)/N)·(1 + c²/s²)
= −1 + 2(c² + 1)/(N·s²), since 1 + c²/s² = 1/s². So the coefficient of δ is
½ − (c² + 1)/(N·s²) + c²/(N·s²) = ½ − 1/(N·s²). ∎

Gate P5 does this elimination symbolically in (N, s, c) with c² = 1 − s² and finds the
difference to the displayed c_m identically zero. The experiment's gate E5a gates it as
a LAW (the residual to c_m falls one decade per decade of δ over N = 6…300, ratios 9.58
to 10.71; a fixed bound would be the wrong shape because the residual's coefficient grows
with N), and E5b/E5c the value and the sign at finite N (both signs of cos k_m occur in
that grid; the symbolic gate P5 takes c real). The sign of c_m changes at
N·s² = 2, positive at the band centre and negative at the edge; for m = 1 the crossing is
between N = 19 and 20 (the root of N·sin²(2π/N) = 2 is 19.03; the large-N equation
4π²/N = 2 gives 2π² = 19.74, which the finite-N root undershoots by 3.6%).

**Scope of c_m.** It is the δ → 0 form and governs the split. At a finite δ the split's
own next order has already moved that zero (measured off the spectrum: between N = 14 and
15 at δ = 0.1), and the (1,1) block's zero-crossing reading of the same pair sits elsewhere
again (between N = 16 and 17 at δ = 0.1, the admixture into other pairs displacing it by
two units of N). The experiment page carries both numbers; this file carries the law.

## (g) Theorem F, the departures

**Theorem F.** Let P(x) = det(x·I − H) as in Theorem A. Then

    P(+2) = −((N−1)u + (N+1))·(u − 1)                 at every N,
    P(−2) = ((N−1)u − (N+1))·(u + 1)                  at odd N,
    P(−2) = −((N−1)u + (N+1))·(u − 1)                 at even N,

and the number of levels outside the band [−2, 2] is: none at u ≤ 1; past u = 1 exactly
one above the band at every N; below the band exactly one at every u > 1 for even N, and
for odd N none until u = (N+1)/(N−1), where the level sits exactly on the edge, and
exactly one beyond it.

*Proof.* U_n(1) = n + 1 and U_n(−1) = (−1)^n·(n + 1). So P(2) = (N+1) − u²(N−1) − 2u =
−((N−1)u² + 2u − (N+1)) = −((N−1)u + (N+1))(u − 1), and P(−2) = (−1)^N[(N+1) − u²(N−1)] − 2u,
which at even N is the same expression as P(2) and at odd N is (N−1)u² − 2u − (N+1) =
((N−1)u − (N+1))(u + 1). The parity of N enters through the sign of U_n(−1) and nowhere else.

For the count, with u ≥ 0: H is real symmetric, so P is real-rooted and monic. P(2) < 0
means an odd number of roots above 2; the sign of P at −∞ is (−1)^N, so (−1)^N·P(−2) < 0
means an odd number of roots below −2. The crack is the rank-two update
u·V, V = e_0e_{N−1}ᵀ + e_{N−1}e_0ᵀ, of the path matrix H_path (the cross-check paragraph of
Theorem A carried the same update with the opposite sign because it perturbs x·I − H;
Theorem G below uses the same V at first order), with eigenvalues (u, 0, …, 0, −u), so by
Weyl's inequality λ₂(H) ≤ λ₁(H_path) + λ₂(uV) = 2cos(π/(N+1)) < 2, and symmetrically
λ_{N−1}(H) ≥ λ_N(H_path) + λ_{N−1}(uV) = −2cos(π/(N+1)) > −2; hence at most one level
passes each edge, and an odd number at most one is exactly one. Now read the factors:
the first factor of P(2) is positive for every u ≥ 0, so P(2) < 0 iff u > 1. At even N,
P(−2) has the same form, so the bottom level leaves iff u > 1 too. At odd N,
(−1)^N·P(−2) = −((N−1)u − (N+1))(u + 1) is negative iff (N−1)u > N+1, that is u > (N+1)/(N−1),
and at u = (N+1)/(N−1) exactly P(−2) = 0: a level sits on the edge, counted here as not
departed and reported as an edge level. At u < 1 Perron-Frobenius (Corollary B) puts every
level strictly inside, consistent with P(2) and (−1)^N·P(−2) both being positive there; at
u = 1 the ring's edge levels sit ON the edges (Corollary C), P(2) = 0 and, at even N,
P(−2) = 0, and none is beyond them. ∎

Gate P6 meets the three factor forms symbolically for N = 3..12 and finds (N+1)/(N−1) as
the exact positive root of P(−2) at odd N. `Crack.cs` reads the same count by Descartes'
rule on P(x + 2) and on P(−x − 2), which is exact for a real-rooted polynomial, and meets the law over
N = 3…20 with the odd threshold met at the exact rational and P(−2) = 0 compared as zero.
The typed claim's tests meet the factor forms against the unfactored polynomial at x = ±2
and the count against Descartes' rule on the shifted integer polynomial, which sees a root
on the edge as a vanishing constant term. The experiment's gate E7c-theorem checks the
band-bottom asymptote, k = π − ε: −ε(1+u)[(1−u)N + (1+u)] at odd N and +ε(1−u)[(1+u)N + (1−u)]
at even N, which is Theorem F's P(−2) read as a limit; the bracket vanishing at
u = (N+1)/(N−1) in the odd case and the prefactor (1 − u) in the even case is the same
parity, seen from the curve.

The first draft of this count, on the experiment page and in the sibling's strengthening
convention u = 1 + δ, read *"two levels out of the band at every δ > 0"*, which is true at
even N only; it was caught before commit and is the CAUGHT_ERRORS entry of 2026-08-31.
That is why the parity is in the statement of Theorem F and not in a remark under it.

## (h) Theorem G, the velocity at the chain end

**Theorem G.** At u = 0 the roots of P are the chain comb x_k = 2cos(kπ/(N+1)),
k = 1..N, and

    dx_k/du |_{u=0} = (−1)^{k+1}·(4/(N+1))·sin²(kπ/(N+1)).

Its magnitude is F65's endpoint rate comb α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)).

*Proof, from the polynomial.* Implicit differentiation of P(x_k(u); u) = 0 gives
dx_k/du = −(∂P/∂u)/(∂P/∂x). At u = 0, ∂P/∂u = −2 (the u²-term contributes nothing at u = 0)
and ∂P/∂x = (d/dx)U_N(x/2). With x = 2cos θ, U_N(x/2) = sin((N+1)θ)/sin θ, so
(d/dx)U_N = [(N+1)cos((N+1)θ)·sin θ − sin((N+1)θ)·cos θ]/sin²θ · (−1/(2sin θ)), and at
θ_k = kπ/(N+1), where sin((N+1)θ_k) = 0 and cos((N+1)θ_k) = (−1)^k, this is
−(N+1)(−1)^k/(2sin²θ_k). Hence dx_k/du = 2/P′(x_k) = (−1)^{k+1}·4sin²θ_k/(N+1). ∎

*Proof, from the eigenvectors.* The crack at first order is the perturbation
V = e_0e_{N−1}ᵀ + e_{N−1}e_0ᵀ, the chain spectrum is simple (Lemma J of the blind-seat proof,
or directly from the comb), so dE_k/du = ⟨ψ_k|V|ψ_k⟩ = 2ψ_k(0)ψ_k(N−1) with F2b's
ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)). The chain reflection gives ψ_k(N−1) =
√(2/(N+1))·sin(πk − πk/(N+1)) = (−1)^{k+1}ψ_k(0), so dE_k/du = (−1)^{k+1}·2|ψ_k(N−1)|², and
F65 reads α_k = 2γ₀·|ψ_k(N−1)|². ∎

Two routes to one number: the polynomial's derivative and the endpoint amplitude. The
sign is the chain reflection's (F75's mirror sign η), the magnitude is the rate the same
level pays under one dephased endpoint, twice the Absorption Theorem's light there.
*Velocity* means dE/du along the road, not F2b's group velocity dE/dk on the same chain.
Gate P7 meets both routes against the closed form for N = 3..9 and every k, exactly (the
derivative of U_N as a trigonometric identity in θ, then sin((N+1)θ) = 0 and
cos((N+1)θ) = (−1)^k substituted); the typed claim's test meets the closed form against a
finite difference of the polynomial's roots at u = 1e-3 and 1e-4 by a decade law. What
Theorem G is then good for, F129's collisions moving
along the road by the signed F65 sums, is the subject of
[The Comb on the Road](../../experiments/THE_COMB_ON_THE_ROAD.md).

## (i) The gates, in one table

| Gate | What it meets | Route | Where |
|------|---------------|-------|-------|
| P1 | Theorem A, N = 3..9 | symbolic determinant vs recursion, polynomials in (x, u) | `cracked_ring_exact_curve_proof.py` |
| P2 | Cassini, n = 1..12 | symbolic | same |
| P3 | Corollary B, both forms, N = 3..9 | symbolic, exact by reduction modulo s² + c² − 1 | same |
| P4 | Theorem D, N = 3..9 | symbolic, exact by the same reduction (the Re step as an identity in real A, B) | same |
| P5 | Theorem E | the series in (N, s, c), difference identically 0 | same |
| P6 | Theorem F, N = 3..12 | symbolic; the odd threshold as the exact root | same |
| P7 | Theorem G, N = 3..9, every k | both routes vs the closed form, exact (the derivative identity by reduction, then the substitution at θ_k symbolic in k) | same |
| P8 | Corollary B's simplicity clause: G = 2AB symbolic in N, and the Bézout identity (cos a + u cos b)A + (sin a − u sin b)B = 1 − u² symbolic in N | symbolic | same |
| P9 | the prefactor relations behind the simplicity clause: A = cos(k/2)·χ_e(2cos k), B = sin(k/2)·χ_o(2cos k) at even N and A = ½·χ_e, B = sin k·χ_o at odd N, with χ_e, χ_o the folded blocks' characteristic polynomials, N = 3..9, u symbolic | symbolic (the blocks built in sympy, the identities by reduction in the half angle) | same |
| CrackTests | Theorem A over ℤ, mod two primes to N = 1001; Theorem F by Descartes N = 3..20; the c_m decade law; the 0.971754 pin | exact integers, GF(p), roots as a reading | `compute/MirrorWorld.Tests/CrackTests.cs` |
| CrackedRingExactCurveClaimTests | Theorem A vs Bareiss at eight (N, u); the ring end vs the first-kind recursion and the chain end vs the path pencil's Bareiss determinant; Theorem F's factors vs the unfactored P(±2) and the count vs Descartes on P(y ± 2) over ℤ; c_m by a decade law at five (N, m); Theorem G vs a finite difference by a decade law | exact integers; two float rows with an error model | `compute/RCPsiSquared.Core.Tests/Symmetry/` |
| E1a-c, E2a-d (d in three parts), E3a-c, E4a-b, E5a-e, E6a-b, E7a-e with E7c twice (27 gates) | the curve, the ends, the split as truncation, the join, the split's law and the 0.9705 pin, the road past u = 1 | eigensolver, except E1b, E2a, E2b, E6b and E7c-theorem (the curve alone), E5d (a matrix element of the crack between plane waves), E1c (a determinant), and E3a-b, E4a (identities between displayed formulas) | `cracked_bell_gate.py` stage E |

## (j) Scope and fences

- The XY adjacency book, single-excitation block, uniform bonds off the crack, u ≥ 0.
  The Heisenberg Laplacian is a different curve; this file says nothing about it. Where γ is
  mentioned it is the uniform Z-dephasing of the Absorption Theorem's statement.
- γ appears nowhere in G. The Liouvillian rate of the coherences this spectrum labels is
  the Absorption Theorem's −2γ, and the beat on top of the split, the reversal and the
  visibility wall of the experiment page are first-order statements this file does not
  carry. No time is owned here; the clocks are the page's and `Warble`'s.
- *Departure* is a level leaving the band [−2J, 2J], not the departure from normality of
  the F2b corollary and F89 (a conditioning number of small non-normal matrices; F160's
  own parenthetical filed that word under F86, which never uses it, and is corrected in
  the same change); *band* is that cosine band, not F3's rate band and not the glossary's pre-registration
  constant; *crack* is a weakened wrap bond, and the strengthened one of the sibling is
  the same road past u = 1.
- Not the blind seat: THE_SEAT_THAT_CUTS's open item asks for a detuned bond under a seat
  cut, a different object. What DOES cross between them is Corollary B's fold, and only
  it: [PROOF_BLIND_SEAT_TWO_AXES](PROOF_BLIND_SEAT_TWO_AXES.md) folds the
  end-pair-anisotropic open chain the same way and finds the same two blocks, the
  anisotropy at +Δ giving the crack's even block and at −Δ its odd one. It spends that fold
  twice over, on the blind seat's LOCUS and, in its section (g), on its COUNT. That file
  deliberately does NOT use the simplicity clause for its own degeneracy statement,
  because this file states that clause for u ≥ 0 and half of its u axis is negative; it
  proves the u² = 1 condition from the boundary system instead. Nothing of this file's
  own scope moves.
- Theorem G is a first-order statement at u = 0; what a finite u does to a comb
  coincidence is measured, not derived, on the comb page. The ORDERS BETWEEN the two,
  d₂ to d₅ of the same series, are derived in
  [PROOF_COLLISION_GAP_ODD_ORDERS](PROOF_COLLISION_GAP_ODD_ORDERS.md) (F161), which
  takes Theorem A's polynomial as its whole input and Theorem G as its first term.
- No hardware claim.
