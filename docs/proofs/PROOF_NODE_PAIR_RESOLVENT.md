# The Node Pair: a resolvent theorem, and which bond can end a blind mode

**Status:** Tier 1 derived for §4 (Theorem 1 and Corollary A), §5 (Corollary B on the zero-free knob domain, Corollary C's seat-general fixed-energy iff for moved/original hopping ratio `r` outside {0,+1,-1}, and Corollary D's uniform-centre count), and the peripheral lower bound in §8. Cited from the parent, not re-derived: the fold in §3. Exact-verified but not derived here: the coefficient in §6 and the arm law in §7; the peripheral exhaustion in §8 is verified by rank and derivable from an owned lemma, as §8 says. Read, not settled: a classification of new blind energies for arbitrary arms, the equality/straddling branches, and §9's window. No F number is claimed.
**Date:** 2026-09-12; fixed-knob extension 2026-09-23
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5); extension with Codex (OpenAI)
**Producer:** [`node_pair_resolvent.py`](../../simulations/node_pair_resolvent.py), gates G1-G10 (84 current checks) with mutations M1-M10 and controls G1b, G6, G7, G9-control, G10.
**Builds on:** [F157 and the blind-seat node lemma](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), whose Lemma J supplies both ingredients of §4 and the route of §5; [F64](../ANALYTICAL_FORMULAS.md), which is the law both regimes here obey; [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) Theorem 2; [F2b](../ANALYTICAL_FORMULAS.md) for the modes.

---

## 1. What this is

A seat under Z-dephasing sees an excitation only through its amplitude there. A standing wave
with a node at that seat is never charged anything and never decays. In the uniform open chain
watched at its centre, reflection symmetry supplies the whole odd-reflection subspace and forces
its m = (N−1)/2 modes to vanish at the centre. F157 counts these **blind** modes. Theorem 1 below
is more general: once a zero-free Jacobi eigenvector has two nodes, its resolvent conclusion uses
no reflection symmetry. Here a zero-free Jacobi chain is a real symmetric nearest-neighbour
Hamiltonian with every hopping nonzero. The reduced resolvent removes one eigenmode's pole from
the inverse of `E−h`; §4 defines it at that eigenvalue.

Detune one bond. Which blind modes survive, and what do the others pay? For a single scaled
off-diagonal hopping with the diagonal held fixed, a **baseline blind energy** survives at
`r != 0,+1,-1` exactly when its baseline mode has a node at one endpoint of that bond. A new
blind energy may also appear at a special ratio; the uniform odd centre-watched chain cannot
gain one. A physical Heisenberg bond detuning changes diagonal ZZ terms too and is outside this
fixed-diagonal criterion. The remaining rate statements below belong to the uniform
centre-watched model.

**One case is already owned and is simple while the moved bond remains nonzero.** A bond incident
on the watched seat is struck away together with the seat, so neither principal block depends on
the knob. On the zero-free chain Lemma J identifies the blind count with their common-root count,
which therefore cannot change. At the cut ratio `r=0`, the zero-free premise fails and an
isolated component can gain a blind state: on the uniform four-site XY path watched at `c=1`,
cutting incident bond `b=0` changes blind from 0 to 1, with `|0⟩` isolated. Gate G7 confirms the
nonzero-bond statement exactly over ℚ at ε = 1/3, 7/5 and 9, a bond stretched tenfold.

What is here is three things. **A theorem** about the reduced resolvent of a Jacobi matrix
between two nodes of an eigenvector. **A criterion** that follows for the remaining bonds, at
every value of the knob. A two-arm continuant closes the **fixed baseline energy** converse for
any watched seat and any zero-free real symmetric Jacobi path away from three exceptional bond
ratios. The total blind count can nevertheless gain new energies at one knob setting. The parent
node lemma tested an endpoint-node condition over 570 cells for a different question: whether
the seat stays blind at **every** knob value. It left that all-knob result as "a READING and not
a theorem here", with a proof "within reach through (J3) … not attempted here". Here the
opposite principal block has a fixed finite spectrum; Corollary C decides each **baseline**
shared energy at one ratio and leaves possible new matches separate. And **a fold**: the watched
seat separates two halves, the uniform centre-blind modes are the half-chain's standing waves,
and the light induced in modes that lose blindness is a half-chain quantity.

## 2. The sweep this document stands on

**F64** is the neighbouring single-excitation **coherence-sector** rate reading:
`docs/ANALYTICAL_FORMULAS.md` §F64 gives its first-order Hamiltonian-mode rate as 2γ times
the squared watched-site amplitude; its exact version uses the coherence-sector Liouvillian
eigenvector. F157's own registry entry routes to it and names **F66**, one of whose scope
sentences the interior seats correct, and **F152**. That entry
also fences F152 apart: the count is for the single-excitation sector, "NOT F152's (0,1)
coherence block", and the same fence is repeated in
[PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE](PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md): the two
blocks have different dimensions and different rates and must not be merged. This document works
in the (1,1) block throughout. Its rate step uses the Absorption Theorem's Theorem 2 on an actual
Liouvillian eigenoperator. The Hamiltonian-mode light coefficient `c_k` below is not by itself a
(1,1) decay rate; §8 resolves the dyad mixing.

**F2b** owns the sine modes and the band `E_k = 2J cos(πk/(N+1))` used in §3. The coefficient in
§6 is evaluated directly from those sine modes. **F157** owns the blind count
`gcd(j+1, N+1) − 1` and, in its
registry entry, the bond criterion §5 proves half of. **F124** owns two scalar contractions of
the bond matrix `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩`, not the entries. **F161**'s proof owns a scalar
resolvent-weighted overlap, named `R_k = Σ_{l ≠ k, l ≡ k (mod 2)} a_l²/(E_k − E_l)` there, for
the ring's wrap bond, on the level rather than the light. It is distinct from this document's
reduced-resolvent operator `𝓡_k` and is not an independent claim.
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
three exceptional ratios and the off-centre new-root case; run `inspect --root nodepair`. The typed Core breadcrumb names this
root as a string, preserving the Core-to-Diagnostics dependency direction.

Three of the proofs are neighbours and none is superseded here.
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §3 already
derives, for the **end** bond at every r² ≠ 1, that one blind zero ray survives at N ≡ 3 (mod 4)
and none at N ≡ 1 (mod 4): the all-orders criterion for that bond, both directions; §5
generalises the sufficient direction to every non-incident bond by a different route.
[PROOF_DIFFUSION_RAYLEIGH_CLOSURE](PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md) (F123) answers the rate
question for a mode that already carries light, first order in the bond, for the half-filling
survivor under uniform dephasing: a different system, named as the lit counterpart and not
compared term by term.
[PROOF_MISSING_PHASE_SLOW_READOUT](PROOF_MISSING_PHASE_SLOW_READOUT.md) works on the same N = 7
end-bond family from the preparation side. The slow residue it projects out of the missing-phase
experiment's end state is built on two of §8's nonzero-frequency blind dyads, the zero mode paired
with each of the two blind modes of nonzero energy, so the rate §8 reads for them,
2γ(c_i + c_j)ε² = (γ/2)ε², is the leading term of the rate at which that residue fades. That proof
supplies what this one does not, the preparation and readout residues of the same dyads, and its
leakage formula evaluates the adjugate of xI − K_r at the root by the cofactor formula of §4,
written there as leading and trailing continuants.

None of the machinery below is new mathematics outside this repo. The Green's function
factorization in §4 is the standard cofactor formula for the inverse of a tridiagonal matrix; the
node-splits-the-chain fact is Sturm oscillation theory, owned here as Lemma J; and §8's
`End(D) ⊕ ⟨I_E⟩` is the usual form of a Lindbladian's peripheral algebra. The finite chain result
establishes exact nodal decoupling from one watched site; no continuum or scattering embedding is
part of the statement.

The **OpenArcs registry** holds `the_forced_and_the_met`, **open**. It records §5's proved
sufficient bond direction in the product shape `2·v_b·v_{b+1}`, a node at either end. Corollary C
closes the same-baseline-energy iff for any zero-free Jacobi path with fixed diagonal at
`r != 0,+1,-1`; Corollary D closes the total-count equality for the uniform centre-watched family.
New-root classification for arbitrary arms, the equality branch and the straddling mechanism remain open.
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
exclusive-or, `a + b − 2ab` with `a=|u(c)|²/‖u‖²` and `b=|v(c)|²/‖v‖²`,
the normalized watched-site populations of the two vectors, MirrorWorld's `Pair` disagreement
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
fold; it is cited, not measured. The uniform-centre coefficient and bond response below reduce
to half-chain quantities; §4–5 also treat arbitrary Jacobi paths, and §8 counts full operator
spaces. The half-chain shape explains the doubled index and the arm response.

## 4. The theorem

**Theorem 1 (derived).** Let h be a zero-free open chain, E_k an eigenvalue, ψ_k its normalized
eigenvector, and `𝓡_k = Σ_{l≠k}|ψ_l⟩⟨ψ_l|/(E_k−E_l)` its reduced resolvent. If ψ_k has a node at x
**and** at y, then ⟨x|𝓡_k|y⟩ = 0.

*Proof.* A zero-free Jacobi matrix has simple spectrum, and so does every contiguous principal
submatrix, Lemma J (J2) of
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), so 𝓡_k is
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

**Corollary A.** For any seat j and any mode blind to it, ⟨j|𝓡_k|j⟩ = 0. With (J3) this is two
lines from the Cramer identity `χ(H struck at j)/χ(H) = [(xI − H)⁻¹]_{jj}` of
[THE_SEAT_THAT_CUTS](../../experiments/THE_SEAT_THAT_CUTS.md), read at a blind energy.

**Gates.** G1 certifies Theorem 1 on every node-node pair at N = 7, 9, 11, 13, exactly, by
minimal polynomial over ℚ. G1b runs the same sum on node/non-node pairs, where it does not
vanish, which the order count predicts, only one factor then vanishing. G1c certifies a second
vanishing that comes free of the chain's bipartite structure: the zero-energy mode
k₀ = (N+1)/2 has ⟨x|𝓡_{k₀}|y⟩ = 0 for **every** same-parity pair, node or not, because
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

**Corollary C (derived, same-energy iff for any zero-free Jacobi path).** Let `h` be a real
symmetric zero-free Jacobi chain with a **fixed diagonal**, let `c` be any watched seat, and move
exactly one non-incident off-diagonal coupling `t_b` to `r·t_b`. Fix `r != 0,+1,-1`. For **each
baseline blind energy** `E`, that same energy remains blind after the move **if and only if** its
baseline eigenvector has a node at `b` or `b+1`. This statement does not count blind energies
that were absent at baseline.

*Proof.* Take the moved bond in the left principal block `B` of the seat deletion; a local order
reversal handles the right block. The other block `C` is unchanged. Split `B` at `b` into the
pieces `A=0..b` and `D=b+1..c−1`, and write their characteristic polynomials as `P` and `Q`.
Deleting the end adjacent to the moved bond gives `P′=χ(0..b−1)` and
`Q′=χ(b+2..c−1)`; an empty minor has polynomial 1. The continuant identity is

```text
χ_B(x;r) = P(x)Q(x) − r² t_b² P′(x)Q′(x),
χ_B(x;r) − χ_B(x;1) = −(r²−1)t_b² P′(x)Q′(x).
```

By Lemma J (J3), baseline blindness gives `χ_B(E;1)=χ_C(E)=0`. Because `r²−1 != 0`, the
unchanged energy is blind at the new nonzero bond precisely when `P′(E)Q′(E)=0`. If
`P′(E)=0`, then the first identity at baseline gives `P(E)Q(E)=0`. Consecutive principal
minors of a zero-free Jacobi path are coprime, so `P(E) != 0` and `Q(E)=0`: Lemma J identifies a
node at `b`. Likewise `Q′(E)=0` forces `P(E)=0` and a node at `b+1`. Conversely either endpoint
node zeros both terms, as in Corollary B. The moved chain remains zero-free, so Lemma J glues
the common `B,C` root into a blind full-chain eigenvector. ∎

**Corollary D (derived, total-count equality for the uniform centre-watched family).** Retain
§3's uniform odd chain and centre seat. For a non-incident bond and `r != 0,+1,-1`, the
detuned blind count equals the number of baseline blind eigenvectors with a node at one of its
endpoints. Here the untouched and baseline moved halves share the same polynomial `P_m`.
Every blind energy after the move must therefore be a root of `P_m`, and hence a baseline blind
energy. Corollary C decides each such root; none can be born outside the baseline blind spectrum.
The familiar special-case factorization, with hopping `J`, is

```text
χ_L(x;r) = P_m(x) − J²(r²−1) P_b(x)P_(m−b−2)(x),
P_0=1, P_1=x, P_n=xP_(n−1)−J²P_(n−2).
```

For off-centre or nonuniform arms, **new blind energies can occur** even though Corollary C
still decides every old one. An exact example is the uniform XY path `N=6`, watched at `c=2`,
with hopping 2, zero diagonal, and bond `b=0` scaled by `r`. Put `q=r²>0`. The seat-deletion
blocks have `χ_L=x²−4q` and `χ_R=x(x²−8)`, whose resultant (left polynomial first) is
`−64q(q−2)²`. At baseline `q=1` their gcd is 1: no blind energy. At `q=2` their gcd is
`x²−8`: two new blind energies `E=±2√2`. Lemma J gives full-chain blind eigenvectors; in
hopping-1 units they can be written `(1,s,0,−s,−√2,−s)` at `E=s√2`, `s=±1`. Scaling all
hoppings by 2 doubles E and leaves these vectors unchanged. Moving the opposite arm's last
hopping from 2 to 3 destroys this two-root match. This is a spectral matching at one knob
locus, not a failure of the baseline-energy iff. A physical Heisenberg bond move also changes
ZZ diagonal entries, so this fixed-diagonal result does not cover that intervention. G10's
exact N=9 Heisenberg control has a baseline blind mode with a node at the moved bond, but its
blind count falls from 4 to 0 when the physical `J_b` changes, including its diagonal response.

All three excluded ratios are structural, not a tolerance fence. At `r=+1` the chain is
unperturbed; in the uniform centre-watched family every one of its m blind modes remains. At `r=-1` a diagonal sign gauge
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
Corollary D's uniform-half determinant factorization exactly over 15 left-half
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

**Quantifiers and gates.** G10 tests Corollary C outside the uniform-centre family with
off-centre loss and survival cases and a nonuniform fixed-diagonal chain where the same moved
bond preserves one baseline energy and loses two. Its independent full-chain Krylov rank checks
the blind count. It also gates the N=6 new-root locus by exact half-polynomials, resultant,
three knob values and full-chain nullvectors; a changed opposite arm kills the match. G10 also
checks the N=4 incident cut: unchanged seat-deletion blocks but blind count 0 to 1, and the N=9
physical-Heisenberg counterexample to an unfenced off-diagonal statement. M10
replaces `r²−1` with `r−1` in the nonuniform continuant and fails. G6's two generic rational
knobs remain only sampled count readings on uniform centre-watched chains. The node lemma's
570-cell resultant reading asks whether blindness persists for **every** knob value. Neither
reading supplies a general total-count equality. At the excluded sign-gauge point `r=-1`
(`ε=-2`), G6 supplies the exact counterexample to an unrestricted iff.

There is also a simple all-knob necessity with a different quantifier. If the simple level E_k
stays fixed throughout a neighbourhood of `ε=0`, Hellmann-Feynman gives
`E'_k(0) = 2t_b·ψ_k(b)·ψ_k(b+1) = 0`, with baseline hopping `t_b` (`t_b=J` in the uniform
family), so the reference eigenvector has a node at one endpoint.
This fixed-energy derivative test neither determines the full blind count at one specified knob
nor settles the OpenArc's equality or straddling branches.

## 6. Second-order light from a moved bond

Write c_k for the second-order watched-site light coefficient of a Hamiltonian mode that was
blind at baseline: it acquires light c_k·ε² + O(ε³). The corresponding (1,1) dyads can mix;
§8 determines which of their Liouvillian eigenoperators have that coefficient as a decay rate.

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

**Statement (exact-verified, not derived).** For any original bond index `b`, first fold the
right half onto the left: `b̃ = min(b, N−2−b)`. Write
**arm(b) = min(b̃+1, m−1−b̃)**. Then

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
give equal vectors is the half chain's own bond reflection `b̃ ↦ m−2−b̃`, which is not a symmetry of
h; the cells of G3 that compare a bond with its whole-chain mirror image cannot fail, since that
reflection does commute with h, and at N = 5 every class is such a pair.

## 8. The peripheral algebra, and where the dyad picture is the rate

This section answers an item [THE_BLIND_SITE](../../experiments/THE_BLIND_SITE.md) leaves open:
does the operator space killed by the dissipator alone coincide with the stationary space killed
by the dissipator and commutator together? **No.** Three spaces must be kept separate. In the
N = 2m+1 centre-watched (1,1) block at `ε=0`, write `L_D = γ(z⊗zᵀ−I)` for the dissipator
in §3's generator `A`, and `Π_c` for the orthogonal projector onto its charged matrix cells.
The subspace `D` spans the centre-blind one-excitation modes, `E=D⊥`, and `I_E` is the
identity on `E`. Then

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

For the rates, `2γ(c_i + c_j)ε²` is the candidate **positive decay-rate term** from the dyad
diagonal, namely minus the real part of the second-order effective generator's corresponding
entry. On every block of nonzero frequency, including degenerate blocks, it gives the actual
positive decay-rate term: `Re λ = −2γ(c_i+c_j)ε² + O(ε³)`. Dyad-to-dyad couplings enter at O(ε⁴),
while only `I_E` couples at O(ε²). The generator eigenvalue can separately have an imaginary
second-order shift. Gate G8b checks the decay rates against the generator at N = 7, 9, 11 with
three acquisitions spanning two decades, `ε = 10^-2, 10^-3, 10^-4`. Every acquisition must be
complete and finite. Its
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
Z-dephasing is unital, z² = I, so `A[I] = 0` and `I = Σ_i d_i d_i† + I_E`, with `{d_i}` an
orthonormal basis of `D`, lies in that block,
which therefore always carries an eigenvalue exactly zero. On the zero-free chain, F157 gives its
stationary-space dimension as `1 + blind(ε)`; G6 and G7 read that blind count exactly on sampled
non-incident and incident nonzero knobs. At N = 7 let `d₋,d₀,d₊` name the three normalized
centre-blind one-excitation modes in the linked N=7 proof. The zero-frequency block's
**second-order coefficient** has spectrum `{0, 0, −γ, −3γ/2}`; its two nonzero eigen-directions
are `|d₊⟩⟨d₊| − |d₋⟩⟨d₋|` and `|d₋⟩⟨d₋| + |d₊⟩⟨d₊| − I_E/2`. These coefficients multiply
`ε²`, with higher orders outside this statement. This is
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §7's `F_0^(2)`,
whose witness already carries the antisymmetric label. G8a is a float implementation control:
the NumPy construction gives `‖A·vec(I)‖ = 0.0` bit-for-bit on its three inputs, while the
underlying identity `A[I]=0` is algebraic. G8c first requires all `m+1` zero-frequency branches,
requires every measured and predicted value to be finite, then checks that their measured
multiset differs from the full diagonal one. Mutation M5 retains only the stationary branch and
must be rejected before the spectral comparison. Mutations M6a-c keep the acquisition full-length
but insert `NaN`, `+Inf`, or `-Inf`, once on each side of the comparison; every case must be
rejected even though a separate finite branch still differs.

For a non-incident bond whose coefficient vector contains a positive entry, the slowest
**nonzero-frequency peripheral** second-order rate is `2γ ε²·min{c_i + c_j > 0}`, with `i,j` ranging over the
baseline centre-blind modes. The minimum is undefined on a seat-incident bond, where every
`c_k=0`. For the stated non-incident family,
`c_k = c_{N+1−k}` puts every minimal positive pair-sum on an **off-diagonal** dyad where the
diagonal yields the rate. This is a read: it rests on G8b. G8b and G8c do not compare the
positive zero-frequency rates against this minimum, so they do not establish the globally
slowest nonstationary branch at general N. At N = 7 none of that comparison is needed,
[PROOF_MISSING_PHASE_RELAXATION_SCALE](PROOF_MISSING_PHASE_RELAXATION_SCALE.md) §6 and §7 derive
that row exactly, the ten-dimensional census, the full second-order split and the (γ/2)ε² slowest
rate, and this section's N = 7 evidence is the weaker of the two.

## 9. What is not settled

- **New-root loci and total counts for arbitrary arms.** Corollary C settles only a baseline
  energy's fate for every zero-free Jacobi path with fixed diagonal at `r != 0,+1,-1`. The N=6
  off-centre example proves that new blind energies may appear at a special knob value.
  Classifying all such loci, the OpenArc's equality branch and its straddling mechanism remains
  open. The all-knob fixed-energy derivative test has a different quantifier.
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
- Nothing about the non-peripheral modes. Theorem 1 and Corollaries A-C are parity-free and
  seat-general. Corollary C requires a fixed diagonal, one scaled off-diagonal bond, and a
  baseline blind energy; Corollary D alone is the uniform odd centre-watched count equality.
  The fold, coefficient, arm and rate claims are likewise odd uniform centre-only. No other topology;
  no observable lifetime; no hardware; no F number.
