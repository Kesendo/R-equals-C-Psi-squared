# The Node Pair: a resolvent theorem, and which bond can end a blind mode

**Status:** Tier 1 derived for §4 (Theorem 1 and Corollary A), §5 (Corollary B on the zero-free knob domain and Corollary C's pointwise iff for the uniform centre-watched family away from `r=0,+1,-1`), and the peripheral lower bound in §8. Cited from the parent, not re-derived: the fold in §3. Exact-verified but not derived here: the coefficient in §6 and the arm law in §7; the peripheral exhaustion in §8 is verified by rank and derivable from an owned lemma, as §8 says. Read, not settled: the off-centre/nonuniform pointwise converse and §9's window. No F number is claimed.
**Date:** 2026-09-12
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5)
**Producer:** [`node_pair_resolvent.py`](../../simulations/node_pair_resolvent.py), gates G1-G9 (73 current checks) with mutations M1-M9 and controls G1b, G6, G7, G9-control.
**Builds on:** [F157 and the blind-seat node lemma](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), whose Lemma J supplies both ingredients of §4 and the route of §5; [F64](../ANALYTICAL_FORMULAS.md), which is the law both regimes here obey; [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) Theorem 2; [F2b](../ANALYTICAL_FORMULAS.md) for the modes.

---

## 1. What this is

A seat under Z-dephasing sees an excitation only through its amplitude there. A standing wave
with a node at that seat is never charged anything and never decays. In the uniform open chain
watched at its centre, reflection symmetry supplies the whole odd-reflection subspace and forces
its m = (N−1)/2 modes to vanish at the centre. F157 counts these **blind** modes. Theorem 1 below
is more general: once a zero-free Jacobi eigenvector has two nodes, its resolvent conclusion uses
no reflection symmetry.

Detune one bond. Which blind modes survive, and what do the others pay?

**One case is already owned and is trivial.** A bond incident on the watched seat is struck away
together with the seat, so neither principal block depends on the knob and the blind count cannot
change: [PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) says so
in its setup, calls those cells trivial and excludes them from its gates. Gate G7 confirms it
exactly over ℚ at ε = 1/3, 7/5 and 9, a bond stretched tenfold.

What is here is three things. **A theorem** about the reduced resolvent of a Jacobi matrix
between two nodes of an eigenvector. **A criterion** that follows for the remaining bonds, at
every value of the knob. On the uniform centre-watched family the half-chain determinant also
closes the pointwise converse away from three exceptional bond ratios; off-centre and nonuniform
converses remain open. The parent node lemma had verified the same endpoint-node criterion over
570 all-knob cells and left it as "a READING and not a theorem here", with a proof "within reach
through (J3) … not attempted here". And **a fold**: the watched seat cuts the chain in
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

**F2b** owns the sine modes and the band `E_k = 2J cos(πk/(N+1))` used in §3. The coefficient in
§6 is evaluated directly from those sine modes. **F157** owns the blind count
`gcd(j+1, N+1) − 1` and, in its
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
`NodePairResolventWitness` independently recomputes the canonical exact resolvent zero, its
nonzero control, Corollary C's determinant factorization and theorem/control cases, including the
three exceptional ratios; run `inspect --root nodepair`. The typed Core breadcrumb names this
root as a string, preserving the Core-to-Diagnostics dependency direction.

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
node-splits-the-chain fact is Sturm oscillation theory, owned here as Lemma J; and §8's
`End(D) ⊕ ⟨I_E⟩` is the usual form of a Lindbladian's peripheral algebra. The finite chain result
establishes exact nodal decoupling from one watched site; no continuum or scattering embedding is
part of the statement.

The **OpenArcs registry** holds `the_forced_and_the_met`, **open**. It records §5's proved
sufficient bond direction in the product shape `2·v_b·v_{b+1}`, a node at either end, and closes
the pointwise iff only for the uniform centre-watched family at `r != 0,+1,-1`. The off-centre
and nonuniform pointwise converse, equality branch and straddling mechanism remain open.
`relaxation_scale_as_the_defect_vanishes` is retired beside it. **`docs/CAUGHT_ERRORS.md`** carries
this work's own defects, appended the same day, and warned in advance against reading a
characteristic polynomial as a statement about an eigenvector. **`experiments/THE_BLIND_SITE.md`**
§"What is still open" asks whether the dissipator-killed operator space coincides with the
stationary space; §8 answers **no** and also separates the invariant peripheral space from both.
`hypotheses/`, `reflections/` and `recovered/` hold nothing that bears on the theorems; `review/`
holds the chiral counter-evidence §4 respects. `fw.Confirmations` and the C# `ConfirmationsRegistry`:
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

## 5. Which bond can end a blind mode on the zero-free knob domain

**Corollary B (derived, all orders on the zero-free knob domain).** Let `h` be an arbitrary
zero-free real symmetric Jacobi chain, let `c` be any watched seat, and let bond `b` not be incident
on `c`. Replace its nonzero coupling `t_b` by `t_b(1+ε)`, with `ε ≠ −1`. If an eigenvector
`ψ_k` of `h` is blind to `c` and also has a node at `b` or at `b+1`, then `E_k` stays in the
spectrum and a blind eigenvector at `E_k` survives, not merely to first order.

The surviving vector is **not** ψ_k. Detuning re-weights the arm that carries the moved bond: at
N = 7, bond 0, the zero mode has (h_ε ψ₄)₁ = Jε ≠ 0, while the re-weighted
(1, 0, −r, 0, r, 0, −r) is the eigenvector, which is
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3's `v_r`. What
the proof below establishes is the spectral statement, and nothing more.

*Proof.* Blindness at the seat means, by (J3), that E_k is an eigenvalue of both principal blocks
the seat leaves. Expand the left block's characteristic polynomial along the moved bond:

```text
χ_L = P·Q − (t_b(1+ε))²·P′·Q′,
P = χ(0..b),  Q = χ(b+1..c−1),  P′ = χ(0..b−1),  Q′ = χ(b+2..c−1).
```

A node at b splits the left block there, so by (J3) again E_k is a root of χ(0..b−1) = P′ and of
χ(b+1..c−1) = Q: the first term carries Q, the second carries P′, and both vanish. A node at b+1
gives P(E_k) = Q′(E_k) = 0 and kills the same two terms. Either way χ_L(E_k; ε) ≡ 0 identically
in ε, so E_k stays an eigenvalue of the detuned left block, and of the ε-free right block, for
every `ε ≠ −1`. The chain is then zero-free, and (J3)(⇐) glues the two eigenvectors to an
eigenvector of h_ε vanishing at the seat. ∎

The expansion is written for a bond in the principal block left of `c`. For a bond right of an
arbitrary `c`, apply a local order reversal of the right principal block only, putting its free
end first and the watched-seat boundary last. This is a permutation of that block, not a global
reflection of `h`; it assumes neither a centre seat nor reflection symmetry. The one ε the
statement excludes is ε = −1, where the bond is zero, the chain is cut and Lemma J's hypotheses
lapse. The characteristic-polynomial identity remains true there, but this proof makes no
eigenvector claim at the cut point.

**Corollary C (derived, pointwise iff for the uniform centre-watched family).** Retain §3's
uniform odd chain, centre seat and a non-incident bond. Fix the bond ratio `r=1+ε` with
`r != 0,+1,-1`. A baseline blind energy remains a blind energy of the detuned chain **if and
only if** its baseline eigenvector has a node at either endpoint of the moved bond.

*Proof.* It is enough by reflection to put the bond in the left half, `0 <= b <= m-2`. Let
`P_n(E)` be the characteristic polynomial of the uniform n-site path, with
`P_0=1`, `P_1=E`, and `P_n=E P_(n-1)-J²P_(n-2)`. The untouched right half has
`chi_R(E)=P_m(E)`. Expanding the detuned left half across bond b and subtracting its `r=1`
value gives the exact identity

```text
chi_L(E;r) = P_m(E) - J²(r²-1) P_b(E) P_(m-b-2)(E).
```

At a baseline blind energy `E_k`, `P_m(E_k)=0`. Since `r²-1 != 0`, the same energy is also a
root of `chi_L` exactly when `P_b(E_k)P_(m-b-2)(E_k)=0`. Lemma J (J3), applied to the uniform
half path, identifies those two alternatives with a node at b or b+1. The moved chain is
zero-free because `r != 0`, so the common left/right root glues to a unique blind eigenvector.
Corollary B supplies the reverse implication, completing the iff. ∎

All three excluded ratios are structural, not a tolerance fence. At `r=+1` the chain is
unperturbed and every one of its m centre-blind modes remains. At `r=-1` a diagonal sign gauge
returns the same spectrum and blind count. At `r=0` the bond is cut and the zero-free hypothesis
behind Lemma J lapses. None of these exceptional fibres refutes the punctured statement.

This is the sufficient direction of the bond face
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) left open, by the
(J3) route it named. Physically the two vanishing terms are two distinct blocks: shaking a bond is
a source at each of its ends, and a node there silences the **source**, while Theorem 1 severs the
**path**: one zero doing two jobs, which is why the criterion is "a node at either end" rather
than a condition at every site.

**Gates.** G9 certifies the sufficient identity in ε symbolically, over every (N, non-incident bond,
node-carrying blind mode) uniform centre-watched cell at N = 5 … 13, with a control requiring the
no-node cells not to be identically zero. Its coefficient-wise oracle calls `exact_zero` on every
coefficient rather than trusting a residual SymPy form: the seven-root zero identity is its
fail-open control, paired with a genuinely nonzero trigonometric sequence. G9b certifies
Corollary C's determinant factorization exactly over 15 left-half
`(m,b)` cells, and M8 replaces `r²-1` by `r-1`, failing in all 15. G9c compares exact Krylov blind
counts with the endpoint-node count at `r=2/3,3/2,-2`, including both zero-node and positive-node
cells; G9d records the three exceptional fibres at N=7. `NodePairResolventWitness` independently
recomputes canonical cases of the same identities in C# exact arithmetic. G6 certifies the
resulting node-count lower bound by F157's own instrument in
exact arithmetic, the blind count of the **detuned** chain as N minus the rank of the seat's
Krylov matrix over ℚ, at two rational ε, for the centre seat and every non-incident bond of
N = 5, 7, 9, 11. Equality with that lower bound on those sampled knobs is a reading, not a
converse. Its exact `N=5, c=2, b=0, ε=−2` control exposes the sign-gauge return: the node lower
bound is 0 while the Krylov blind count returns to 2.

**Only the uniform centre-watched pointwise converse is derived.** G6's two generic rational
knobs are a sampled count reading on uniform centre-watched chains only; G6 furnishes no evidence
for the off-centre or nonuniform zero-free families. The node lemma's separate 570-cell
resultant reading asks whether blindness persists for **every** knob
value and is not evidence for that pointwise question. Corollary C closes the uniform-centre
family by using its shared half-path polynomial; no such shared polynomial has been established
for the adjacent open families. At the excluded sign-gauge point `r=-1` (`ε=-2`), G6 supplies
the exact counterexample to an unrestricted iff above.

There is also a simple all-knob necessity with a different quantifier. If the simple level E_k
stays fixed throughout a neighbourhood of `ε=0`, Hellmann-Feynman gives
`E'_k(0) = 2J·ψ_k(b)·ψ_k(b+1) = 0`, so the reference eigenvector has a node at one endpoint.
This fixed-energy derivative test neither determines the full blind count at one specified knob
nor settles the OpenArc's equality or straddling branches.

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
c_{k₀} = 0, and §5 upgrades that to blindness throughout the zero-free knob domain. It is blind
to the centre exactly when
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

This section answers an item [THE_BLIND_SITE](../../experiments/THE_BLIND_SITE.md) leaves open:
does the operator space killed by the dissipator alone coincide with the stationary space killed
by the dissipator and commutator together? **No.** Three spaces must be kept separate. In the
N = 2m+1 centre-watched (1,1) block at `ε=0`,

```text
dim ker L_D                                      = (N−1)² + 1,
dim(maximal A-invariant subspace inside ker L_D) = m² + 1,
dim ker A                                        = m + 1.
```

The first count is the uncharged cell space: both indices avoid the watched seat, or both equal
it. It is not invariant under the Hamiltonian; for example `|0⟩⟨0|` is dissipator-killed but does
not commute with the zero-free hopping matrix.

The middle space is `End(D) ⊕ ⟨I_E⟩`. It has no charged cell and is A-invariant. On such a subspace
⟨v, Av⟩ + ⟨Av, v⟩ = −4γ‖Π_c v‖² = 0, so the restriction is anti-Hermitian and the whole subspace
is peripheral: the peripheral dimension is at least **m² + 1**, which is derived, and on the fold
that reads `dim End(half chain) ⊕ ℂ`, the usual shape of a Lindbladian's peripheral algebra. Its
zero-frequency part, the m diagonal dyads plus `I_E`, is the parent node lemma's
`dim ker L_SE(j) = 1 + blind = m+1`, its own Corollary C on the zero-free chain. Thus the stationary space is the
zero-frequency part of the peripheral space, not the whole peripheral space and not the raw
dissipator kernel.

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
O(ε²). Gate G8b checks that against the generator at N = 7, 9, 11 with three acquisitions spanning
two decades, `ε = 10^-2, 10^-3, 10^-4`. Every acquisition must be complete and finite. Its
normalized-error budget is `2|ε| + 100 eps_machine ||A_0||/(γ ε²)`; after subtracting the stated
rounding reserve, successive maximum errors must fall at least linearly with `|ε|`, with a named
1.25 finite-interval slack. This remains a float **read**, not an exact verification or an error
model for the theorem. Mutations M7a-c replace one acquired nonzero-frequency rate with `NaN`,
`+Inf`, or `-Inf`, retain the full length, and must all be rejected through that same multi-decade
path. M9 permanently scales every predicted coefficient by `1.003`; it too goes through the same
three acquisitions and is rejected because a constant coefficient error does not close as `ε`
shrinks.

**On the zero-frequency block the diagonal is not the rate.** That block holds the m diagonal
dyads together with `I_E`, and the operator mixes them. The reason is exact and holds at every ε:
Z-dephasing is unital, z² = I, so `A[I] = 0` and `I = Σ_i d_i d_i† + I_E` lies in that block,
which therefore always carries an eigenvalue exactly zero. On the zero-free chain, F157 gives its
stationary-space dimension as `1 + blind(ε)`; G6 and G7 read that blind count exactly on sampled
non-incident and incident knobs. At N = 7 the block's spectrum is {0, 0, −γ, −3γ/2}·ε², the two
nonzero eigenvectors being `|d₊⟩⟨d₊| − |d₋⟩⟨d₋|` and `|d₋⟩⟨d₋| + |d₊⟩⟨d₊| − I_E/2`; that is
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §7's `F_0^(2)`,
whose witness already carries the antisymmetric label. G8a is a float implementation control:
the NumPy construction gives `‖A·vec(I)‖ = 0.0` bit-for-bit on its three inputs, while the
underlying identity `A[I]=0` is algebraic. G8c first requires all `m+1` zero-frequency branches,
requires every measured and predicted value to be finite, then checks that their measured
multiset differs from the full diagonal one. Mutation M5 retains only the stationary branch and
must be rejected before the spectral comparison. Mutations M6a-c keep the acquisition full-length
but insert `NaN`, `+Inf`, or `-Inf`, once on each side of the comparison; every case must be
rejected even though a separate finite branch still differs.

Consequently the slowest nonstationary second-order rate is `2γ ε²·min{c_i + c_j > 0}`, since
c_k = c_{N+1−k} puts every minimal positive pair-sum on an **off-diagonal** dyad where the
diagonal is the eigenvalue. This is a read: it rests on G8b. At N = 7 none of it is needed,
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §6 and §7 derive
that row exactly, the ten-dimensional census, the full second-order split and the (γ/2)ε² slowest
rate, and this section's N = 7 evidence is the weaker of the two.

## 9. What is not settled

- **The off-centre and nonuniform pointwise fixed-knob converses of §5.** Corollary C settles the
  uniform centre-watched family at `r != 0,+1,-1`; the all-knob fixed-energy derivative test has
  a different quantifier and does not settle either adjacent family.
- **The coefficient and the arm law** are exact at the N gated and derived nowhere. The fold is
  not among them: it is the parent's closed form for the blind eigenvectors, cited in §3.
- **Everything about rates is a read.** G8b now has a complete finite three-acquisition germ gate
  and explicit truncation/rounding budgets, but those finite float checks are not an exact proof
  or a global error model. That includes §8's last paragraph.
- **The window.** The second-order law has one: at fixed N there is a region in (γ, ε, J) beyond
  which the rate leaves it. A saturation was measured at N = 5, 7 and, over part of the range, at
  N = 11, and no saturation at N = 9 or 13; a closed form for the plateau and for the crossing
  was proposed and fails at two of four N. Neither the window nor what replaces the law beyond it
  is established here, and this document states nothing about either. The resemblance to the
  Absorption Theorem's Zeno asymptote is not pursued: that asymptote is for uniform dephasing and
  its J is the Pauli coupling, whose hopping element is 2J, so a comparison has to cross a
  convention the repo has been caught by before.
- Nothing about the non-peripheral modes. Theorem 1 and Corollaries A-B are parity-free and
  seat-general. Corollary C is restricted to the uniform odd chain watched at its centre; there
  is no Corollary-C claim for an even chain, an off-centre seat, or a nonuniform profile. The
  fold, coefficient, arm and rate claims are likewise odd uniform centre-only. No other topology;
  no observable lifetime; no hardware; no F number.
