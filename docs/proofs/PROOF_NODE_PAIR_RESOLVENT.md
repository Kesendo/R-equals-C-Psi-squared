# The Node Pair: a resolvent theorem, and which bond can end a blind mode

**Status:** Tier 1 derived for §4 (Theorem 1 and Corollary A), §5 (Corollary B, at all orders in the knob), and the peripheral lower bound in §8. Cited from the parent, not re-derived: the fold in §3. Exact-verified but not derived here: the coefficient in §6 and the arm law in §7; the peripheral exhaustion in §8 is verified by rank and derivable from an owned lemma, as §8 says. Read, not settled: §9's window. No F number is claimed.
**Date:** 2026-09-12
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5)
**Producer:** [`node_pair_resolvent.py`](../../simulations/node_pair_resolvent.py), gates G1-G9 with mutations M1-M4 and controls G1b, G6, G7, G9-control.
**Builds on:** [F157 and the blind-seat node lemma](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), whose Lemma J supplies both ingredients of §4 and the route of §5; [F64](../ANALYTICAL_FORMULAS.md), which is the law both regimes here obey; [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) Theorem 2; [F2b](../ANALYTICAL_FORMULAS.md) for the modes.

---

## 1. What this is

A seat under Z-dephasing sees an excitation only through its amplitude there. A standing wave
with a node at that seat is never charged anything and never decays. It is not protected by a
symmetry or by a gap; it is protected by not being there, which is why moving one bond can end
it. F157 counts these **blind** modes; on the uniform open chain watched at its centre there are
m = (N−1)/2 of them.

Detune one bond. Which blind modes survive, and what do the others pay?

**One case is already owned and is trivial.** A bond incident on the watched seat is struck away
together with the seat, so neither principal block depends on the knob and the blind count cannot
change: [PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) says so
in its setup, calls those cells trivial and excludes them from its gates. Gate G7 confirms it
exactly over ℚ at ε = 1/3, 7/5 and 9, a bond stretched tenfold.

What is here is three things. **A theorem** about the reduced resolvent of a Jacobi matrix
between two nodes of an eigenvector. **A criterion** that follows for the remaining bonds, at
every value of the knob, which is the sufficient direction of a statement that same node lemma
verified over 570 cells and left as "a READING and not a theorem here", with a proof "within
reach through (J3) … not attempted here". And **a fold**: the watched seat cuts the chain in
half, the blind modes are the half-chain's own standing waves, and everything the surviving
modes pay is a half-chain quantity.

## 2. The sweep this document stands on

**F64** is the law underneath all of it: `docs/ANALYTICAL_FORMULAS.md` §F64 gives the effective
rate of a mode watched at one site as 2γ times its squared amplitude there, so rate = 0 exactly
when the mode has a node at the seat. F157's own registry entry routes to it in those words, and
names **F66**, one of whose scope sentences the interior seats correct, and **F152**. That entry
also fences F152 apart: the count is for the single-excitation sector, "NOT F152's (0,1)
coherence block", and the same fence is repeated in
[PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE](PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md): the two
blocks have different dimensions and different rates and must not be merged. This document works
in the (1,1) block throughout. F64 itself is a first-order statement and is named here as the law
that turns an amplitude into a rate, not as a result being extended; the rate step below is the
Absorption Theorem's Theorem 2.

**F2b** owns the sine modes and the band `E_k = 2J cos(πk/(N+1))` used in §3. **F65** owns the
single-excitation rate spectrum for endpoint dephasing; §6's coefficient is numerically
`½·α_{2k}/γ₀`, but that reads F65 outside its declared index range and the identification is not
made here. **F157** owns the blind count `gcd(j+1, N+1) − 1` and, in its
registry entry, the bond criterion §5 proves half of. **F124** owns two scalar contractions of
the bond matrix `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩`, not the entries. **F161**'s proof owns a
resolvent-weighted overlap `R_k = Σ_{l ≠ k, l ≡ k (mod 2)} a_l²/(E_k − E_l)` for the ring's wrap
bond, on the level rather than the light, and says there that it is not an independent claim.
**F163** already sandwiches a reduced resolvent between two perturbations to reach a second
order, and **D6** already carries a second-order rate with a 1/γ; there is no blanket novelty
here and none is claimed.

In the typed layer's Diagnostics half, `KPartnerSelectionRuleClaim` (Tier 1 derived,
machine-exact N = 3…8) owns `⟨ψ_N|V_b|ψ_1⟩ = 0` for every bond defect: the same two-term bond
matrix element §5 kills, in the lit sector, killed by the same staggering, and
[ORTHOGONALITY_SELECTION_FAMILY](../../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) is the
repo's home for that genre. The Core half holds `SeatCutBlindnessClaim` (F157) and
`AbsorptionTheoremClaim`. Also in the Diagnostics half, and also `Claim`s rather than witnesses,
are `JDefectLightMigrationClaim`, which owns §8's rate-is-light identity at every defect value,
and `MissingPhaseRelaxationScaleClaim`, beside the live witness `MissingPhaseRelaxationScaleWitness`
(`inspect --root missingphasescale`) composes the N = 7 second-order coefficients and labels the
γ-rate direction `|d₊⟩⟨d₊| − |d₋⟩⟨d₋|`, which is why §8 is scoped as it is.

Two of the proofs are neighbours and neither is superseded here.
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3 already
derives, for the **end** bond at every r² ≠ 1, that one blind zero ray survives at N ≡ 3 (mod 4)
and none at N ≡ 1 (mod 4): the all-orders criterion for that bond, both directions; §5
generalises the sufficient direction to every non-incident bond by a different route.
[PROOF_DIFFUSION_RAYLEIGH_CLOSURE](PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md) (F123) answers the rate
question for a mode that already carries light, first order in the bond, for the half-filling
survivor under uniform dephasing: a different system, named as the lit counterpart and not
compared term by term.

None of the machinery below is new mathematics outside this repo. The Green's function
factorization in §4 is the standard cofactor formula for the inverse of a tridiagonal matrix; the
node-splits-the-chain fact is Sturm oscillation theory, owned here as Lemma J; a level sitting on
a transmission zero is the textbook nodally protected bound state in the continuum; and §8's
`End(D) ⊕ ⟨I_E⟩` is the usual form of a Lindbladian's peripheral algebra.

The **OpenArcs registry** holds `the_forced_and_the_met`, **open**, whose still-open list names
"a proof of the bond face" and which already carries §5's criterion in §5's own shape, the product
`2·v_b·v_{b+1}`, a node at either end; §5 supplies its sufficient half and the arc is its home.
`relaxation_scale_as_the_defect_vanishes` is retired beside it. **`docs/CAUGHT_ERRORS.md`** carries
this work's own defects, appended the same day, and warned in advance against reading a
characteristic polynomial as a statement about an eigenvector. **`experiments/THE_BLIND_SITE.md`**
§"What is still open" asks for exactly §8's object: the blind **operator** space of the (1,1)
block, and whether it coincides with the stationary count, and §8 answers it. `hypotheses/`,
`reflections/`, `recovered/` and `review/` hold nothing that bears on the theorems; `review/`
holds the counter-evidence §4 now respects. `fw.Confirmations` and the C# `ConfirmationsRegistry`:
nothing, and nothing here has been to hardware. `docs/GLOSSARY.md` fences "dark" as a spent word with four live senses, so the modes
here are **blind**, F157's word. Two further words are avoided on repo usage rather than on a
glossary entry: "level" is the ordinary word for an eigenvalue, is F129's defined object, and is
this arc's own level-versus-light contrast, so §7's steps are **arms**; "profile" names an input
here (the γ assignment, the bond vector), so §7's output is a **bond response**. "Fold
coordinate" is already spent, on `BlindSeatSectorFactorisationClaim`'s `min(j, N−1−j)`, and is
not used for the arm.

## 3. The setting, and the fold

Open chain of N = 2m+1 sites, uniform hopping J, bond b (between sites b and b+1) scaled to
J(1+ε), Z-dephasing at rate γ > 0 on the centre seat c = m alone, single-excitation sector and
within it the (1,1) joint-popcount block, whose generator in row-stack convention is

```text
A(ε) = −i(h_ε ⊗ I − I ⊗ h_εᵀ) + γ(z ⊗ zᵀ − I),      z = I − 2|c⟩⟨c|.
```

At ε = 0, F2b gives `ψ_k(x) = √(2/(N+1))·sin(πk(x+1)/(N+1))` and `E_k = 2J cos(πk/(N+1))`, with
0-based sites and 1-based modes. A mode is **blind** to a seat when it has a node there; at the
centre the blind modes are the even k, m of them, F157's count. Write D for their span, E for its
orthogonal complement. Throughout, the **light** of an operator is the Absorption Theorem's
charged weight for this one watched seat: the fraction of its weight in cells whose bra and ket
differ at site c. On this sector that is a single bit, so for a dyad `u v†` the light is the
exclusive-or, `a + b − 2ab` with a = light(u), b = light(v), MirrorWorld's `Pair` disagreement
restricted to one site. Gate G4 checks the identity exactly on rationals; it is algebra about the
charged cell set and cannot fail for a physical reason.

**The fold (exact-verified).** Since N + 1 = 2(m+1),

```text
sin(π·(2k′)·(x+1)/(N+1)) = sin(π·k′·(x+1)/(m+1)),
```

so the blind modes, restricted to sites 0 … m−1, are exactly the standing waves of the **m-site
uniform chain**: the watched seat is a Dirichlet wall and the blind subspace is the half chain.

This is **derived**, and more generally, in the parent:
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) gives the blind
eigenvectors at any seat in closed form as `v^(c)_l = sin((j − l)·θ_c)` with `θ_c = cπ/h`, and at
the centre seat h = m+1, so substituting l = m−1−x is the display above. Nothing here gates the
fold; it is cited, not measured. Every quantity below is a half-chain quantity, and that is why
the doubled index, the arm law and the peripheral count take the shapes they do.

## 4. The theorem

**Theorem 1 (derived).** Let h be a zero-free open chain, E_k an eigenvalue, ψ_k the eigenvector,
R_k the reduced resolvent at E_k. If ψ_k has a node at x **and** at y, then ⟨x|R_k|y⟩ = 0.

*Proof.* A zero-free Jacobi matrix has simple spectrum, and so does every contiguous principal
submatrix, Lemma J (J2) of
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), so R_k is
defined. For x ≤ y the Green's function factorizes in the standard cofactor form

```text
⟨x|(E − h)⁻¹|y⟩ = u₋(x; E) · u₊(y; E) / W(E),
```

with u₋(x; ·) proportional to the characteristic polynomial of the leading x×x block, u₊(y; ·) to
the trailing block's, and W to the whole characteristic polynomial, so W has a simple zero at
E_k. At E = E_k both solutions are proportional to ψ_k, Lemma J (J3), so each vanishes at its
node, and each does so to order exactly one, its polynomial being squarefree by (J2). The
numerator has a double zero against a simple one, so the quotient vanishes at E_k. That quotient
is the claimed entry: the term subtracted to form the reduced resolvent, ψ_k(x)ψ_k(y)/(E − E_k),
is identically zero at a pair of nodes, so no limit is taken. For x > y exchange the roles. ∎

**Corollary A.** For any seat j and any mode blind to it, ⟨j|R_k|j⟩ = 0. With (J3) this is two
lines from the Cramer identity `χ(H struck at j)/χ(H) = [(xI − H)⁻¹]_{jj}` of
[THE_SEAT_THAT_CUTS](../../experiments/THE_SEAT_THAT_CUTS.md), read at a blind energy.

**Gates.** G1 certifies Theorem 1 on every node-node pair at N = 7, 9, 11, 13, exactly, by
minimal polynomial over ℚ. G1b runs the same sum on node/non-node pairs, where it does not
vanish, which the order count predicts, only one factor then vanishing. G1c certifies a second
vanishing that comes free of the chain's bipartite structure: the zero-energy mode
k₀ = (N+1)/2 has ⟨x|R_{k₀}|y⟩ = 0 for **every** same-parity pair, node or not, because
ψ_{N+1−l}(x) = (−1)^x ψ_l(x) pairs the sum term by term. That is chiral (class BDI) symmetry,
owned here as `ChiralKClaim`. One consequence reaches past this chain and is stated narrowly on
purpose: a single-bond detuning **moves** such a mode: at N = 7 the chain's zero mode goes
(1,0,−1,0,1,0,−1) → (1,0,−r,0,r,0,−r), and `review/EMERGING_QUESTIONS.md` records the chiral
fixed point as the mode a bond defect reorganises most, but it cannot give it weight on the
opposite sublattice, so **a seat on that sublattice stays blind to it**. That is what §5 and §7's
zeros use, and it is not gated beyond the chain.

## 5. Which bond can end a blind mode, at every knob value

**Corollary B (derived, all orders).** Let the bond b not be incident on the watched seat, and let
ψ_k be blind to that seat. If ψ_k also has a node at b or at b+1, then E_k stays in the spectrum
and a blind eigenvector at E_k survives at **every** ε, not merely to first order.

The surviving vector is **not** ψ_k. Detuning re-weights the arm that carries the moved bond: at
N = 7, bond 0, the zero mode has (h_ε ψ₄)₁ = Jε ≠ 0, while the re-weighted
(1, 0, −r, 0, r, 0, −r) is the eigenvector, which is
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3's `v_r`. What
the proof below establishes is the spectral statement, and nothing more.

*Proof.* Blindness at the seat means, by (J3), that E_k is an eigenvalue of both principal blocks
the seat leaves. Expand the left block's characteristic polynomial along the moved bond:

```text
χ_L = P·Q − (J(1+ε))²·P′·Q′,
P = χ(0..b),  Q = χ(b+1..c−1),  P′ = χ(0..b−1),  Q′ = χ(b+2..c−1).
```

A node at b splits the left block there, so by (J3) again E_k is a root of χ(0..b−1) = P′ and of
χ(b+1..c−1) = Q: the first term carries Q, the second carries P′, and both vanish. A node at b+1
gives P(E_k) = Q′(E_k) = 0 and kills the same two terms. Either way χ_L(E_k; ε) ≡ 0 identically
in ε, so E_k stays an eigenvalue of the detuned left block, and of the ε-free right block, at
every ε. By (J3)(⇐) the two eigenvectors glue to an eigenvector of h_ε vanishing at the seat. ∎

The expansion is written for a bond in the left half; a bond in the right half is the same
statement after relabelling the chain end to end, which fixes the centre seat. The one ε the
argument does not reach is ε = −1, where the bond is zero, the chain is cut and Lemma J's
hypotheses lapse: the parent flags the same point for its own family. The conclusion survives
there, E_k sitting in the surviving sub-block, but not by this proof.

This is the sufficient direction of the bond face
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) left open, by the
(J3) route it named. Physically the two vanishing terms are two distinct blocks: shaking a bond is
a source at each of its ends, and a node there silences the **source**, while Theorem 1 severs the
**path**: one zero doing two jobs, which is why the criterion is "a node at either end" rather
than a condition at every site.

**Gates.** G9 certifies the identity in ε symbolically, over every (N, non-incident bond,
node-carrying blind mode) cell at N = 5 … 13, with a control requiring the no-node cells not to be
identically zero. G6 certifies the resulting count by F157's own instrument in exact arithmetic,
the blind count of the **detuned** chain as N minus the rank of the seat's Krylov matrix over ℚ,
at two rational ε, for the centre seat and every non-incident bond of N = 5, 7, 9, 11.

**The converse is not derived.** That a mode with no node at either end does acquire light is
measured here and in the node lemma's 570-cell reading, and proved nowhere, except for the end
bond, where [PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3
settles both directions.

## 6. What a surviving mode pays

Write c_k for the second-order light coefficient: blind mode k acquires light c_k·ε² + O(ε³), and
F64 turns light into rate.

**Statement (exact-verified, not derived).** On every bond of arm 1 (§7): the end bond among
them,

```text
c_k = |ψ_k(1)|² = (2/(N+1))·sin²(2πk/(N+1)).
```

The blind mode pays with its own squared amplitude at one end of the bond that moved: the end
facing the **longer** of the two pieces the bond cuts the half chain into. At the end bond that is
site 1, the inner end; at the other arm-1 bond, b = m−2, it is site b, the outer end, and the two
are exchanged by the half chain's own bond reflection of §7. Gate G2 certifies the value per mode,
exactly, at N = 7, 9, 11; mutation M1 substitutes ψ_k(0)² and fires, and does not test which end
of the moved bond pays.

Read on the fold: a blind mode splits its weight evenly between the two halves, the restriction's
norm² is exactly ½ at every N and every k, so the fold is an isometry up to √2, and

```text
c_k = ½ · |φ_{k′}(1)|²,
```

half the squared amplitude that the corresponding **normalized half-chain mode** φ_{k′} carries at
its own site 1. Measured ratio exactly 0.500000 at N = 7, 9, 11, 13. The zero mode is the case
that matters downstream: sin(π) = 0, so
c_{k₀} = 0, and §5 upgrades that to blindness at every ε. It is blind to the centre exactly when
(N+1)/2 is even, N ≡ 3 (mod 4), which is the one question "is the centre on the zero mode's
sublattice", and the same congruence
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3 reaches from
the two halves' shared root.

## 7. The bond response is the shorter arm

**Statement (exact-verified, not derived).** With the right half mirrored onto the left by
b ↦ N−2−b, write **arm(b) = min(b+1, m−1−b)**. Then

```text
Σ_k c_k(b) = ½ · arm(b),
```

and the arm fixes the whole coefficient vector, not only its sum: two bonds of equal arm give
identical c-vectors. Gate G3 certifies both at N = 5, 7, 9, 11: the sum exactly over ℚ, the
vectors by minimal polynomial, since the individual c_k are irrational ((5±√5)/40 at N = 9),
with mutation M2 shifting the formula by one. The responses:

```text
N=5    ½  0  0  ½
N=7    ½  ½  0  0  ½  ½
N=9    ½  1  ½  0  0  ½  1  ½
N=11   ½  1  1  ½  0  0  ½  1  1  ½
```

The two clauses of the arm add to m, so on the fold the statement is one clause: **the arm is the
shorter of the two pieces the bond cuts the half chain into.** Its two zeros are the seat's own
bonds, §1's trivial cells, where the shorter piece is empty. The symmetry that makes equal arms
give equal vectors is the half chain's own bond reflection b ↦ m−2−b, which is not a symmetry of
h; the cells of G3 that compare a bond with its whole-chain mirror image cannot fail, since that
reflection does commute with h, and at N = 5 every class is such a pair.

## 8. The peripheral algebra, and where the dyad picture is the rate

This section answers an item [THE_BLIND_SITE](../../experiments/THE_BLIND_SITE.md) leaves open in
those words: the blind **operator** space of the (1,1) block, and whether it coincides with the
stationary count, "one killed by the dissipator, the other by the dissipator and the commutator
together".

`End(D) ⊕ ⟨I_E⟩` has no charged cell and is A-invariant. On such a subspace
⟨v, Av⟩ + ⟨Av, v⟩ = −4γ‖Π_c v‖² = 0, so the restriction is anti-Hermitian and the whole subspace
is peripheral: the peripheral dimension is at least **m² + 1**, which is derived, and on the fold
that reads `dim End(half chain) ⊕ ℂ`, the usual shape of a Lindbladian's peripheral algebra. Its
zero-frequency part, the m diagonal dyads plus `I_E`, is the node lemma's `dim ker L_SE(j) =
1 + blind`, which is its Corollary C with Corollary B supplying the zero-free chain: one algebra
read twice, once for stationarity and once for mere non-decay, which is the coincidence
THE_BLIND_SITE asked about.

That the dimension is exactly m² + 1 is exact-verified here and **derivable elsewhere**. Gate G5
computes the largest A-invariant subspace inside the kernel of the charged-cell projector by an
exact rank over GF(p) at two primes with i a square root of −1, no eigensolver, finding 2, 5, 10,
17 at N = 3, 5, 7, 9; mutation M4 moves the charged-cell set to another seat and the dimension
collapses to 1. One prime suffices for the logic: a mod-p computation can only over-count the
subspace and the derived lower bound closes it, so the second guards the implementation. The
derivation this replaces is the **window-edge lemma** of
[PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md): the Hermitian part of A here is the
dephasing dissipator, diagonal with entries in {0, −2γ}, so λ_max = 0 and every peripheral
eigenvalue sits on the window edge; the lemma then gives a joint eigenvector of both parts, and
semisimplicity, in one step. Closing to m² + 1 from there is a few lines in the h-eigenbasis. That
lemma's stated context is a different block and a uniform profile, so porting it needs a scope
sentence this document does not write.

For the rates, `2γ(c_i + c_j)ε²` is the **diagonal** of the second-order effective operator in
the dyad basis. On every block of nonzero frequency it is also the eigenvalue, including the
degenerate ones, because dyad-to-dyad couplings enter at O(ε⁴) while only `I_E` couples at
O(ε²). Gate G8b checks that against the generator at N = 7, 9, 11, at one ε and a float
tolerance: it is the weakest gate in the file and the statement is a **read**, not an
exact-verified one.

**On the zero-frequency block the diagonal is not the rate.** That block holds the m diagonal
dyads together with `I_E`, and the operator mixes them. The reason is exact and holds at every ε:
Z-dephasing is unital, z² = I, so `A[I] = 0` and `I = Σ_i d_i d_i† + I_E` lies in that block,
which therefore always carries an eigenvalue exactly zero. Its stationary space has dimension
1 + blind(ε), which G7 measures. At N = 7 the block's spectrum is {0, 0, −γ, −3γ/2}·ε², the two
nonzero eigenvectors being `|d₊⟩⟨d₊| − |d₋⟩⟨d₋|` and `|d₋⟩⟨d₋| + |d₊⟩⟨d₊| − I_E/2`; that is
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §7's `F_0^(2)`,
whose witness already carries the antisymmetric label. Gates G8a and G8c certify `‖A·vec(I)‖`
exactly 0.0 and that the measured zero-frequency multiset differs from the full diagonal one.

Consequently the slowest nonstationary second-order rate is `2γ ε²·min{c_i + c_j > 0}`, since
c_k = c_{N+1−k} puts every minimal positive pair-sum on an **off-diagonal** dyad where the
diagonal is the eigenvalue. This is a read: it rests on G8b. At N = 7 none of it is needed,
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §6 and §7 derive
that row exactly, the ten-dimensional census, the full second-order split and the (γ/2)ε² slowest
rate, and this section's N = 7 evidence is the weaker of the two.

## 9. What is not settled

- **The converse of §5**, except at the end bond, where the neighbouring proof settles it.
- **The coefficient and the arm law** are exact at the N gated and derived nowhere. The fold is
  not among them: it is the parent's closed form for the blind eigenvectors, cited in §3.
- **Everything about rates is a read**, gated once, at one ε, against a float tolerance that is
  not an error model. That includes §8's last paragraph.
- **The window.** The second-order law has one: at fixed N there is a region in (γ, ε, J) beyond
  which the rate leaves it. A saturation was measured at N = 5, 7 and, over part of the range, at
  N = 11, and no saturation at N = 9 or 13; a closed form for the plateau and for the crossing
  was proposed and fails at two of four N. Neither the window nor what replaces the law beyond it
  is established here, and this document states nothing about either. The resemblance to the
  Absorption Theorem's Zeno asymptote is not pursued: that asymptote is for uniform dephasing and
  its J is the Pauli coupling, whose hopping element is 2J, so a comparison has to cross a
  convention the repo has been caught by before.
- Nothing about the non-peripheral modes; no off-centre seat beyond §4 and §5, which are
  seat-general; no even N; no other topology; no observable lifetime; no hardware; no F number.
