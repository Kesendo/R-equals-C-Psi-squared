# PROOF: A physical readout of the N=7 missing-phase slow cluster

**Status:** Tier 1 derived, local in the one-end defect at fixed γ > 0.
**Date:** 2026-09-22
**Authors:** Thomas Wicht, Codex (OpenAI), Claude (Anthropic)
**Typed claim:** [`MissingPhaseSlowReadoutClaim`](../../compute/RCPsiSquared.Diagnostics/Foundation/MissingPhaseSlowReadoutClaim.cs)
**Live exact witness:** [`MissingPhaseSlowReadoutWitness`](../../compute/RCPsiSquared.Diagnostics/Foundation/MissingPhaseSlowReadoutWitness.cs), §2 to §7.1
**Symbolic three-spin check:** [`missing_phase_three_spin_readout.py`](../../simulations/missing_phase_three_spin_readout.py)
**Owning experiment:** [The Motion and the Missing Phase](../../experiments/THE_MOTION_AND_THE_MISSING_PHASE.md)

## What this is about

Seven quantum spins form a chain. Neighbouring spins can exchange an
excitation. In [the project's light picture](../GLOSSARY.md#parameters),
light falls only on the middle spin. This is our way of picturing local
Z-dephasing, with strength γ: phase relationships fade between alternatives
that differ at that spin.

With equal couplings, some patterns of internal motion have a node at the
middle site: their amplitude is zero there, even while the excitation moves
elsewhere. This internal motion stays untouched by the light. Changing one
end coupling slightly reshapes the motion and exposes part of it to the lit
centre. The illumination γ stays fixed, but that part now slowly fades. The
smaller the change, the longer this fading takes.

The equations already tell us that such a slow motion exists. Here we ask
how to see it. We use a starting state built from an excitation shared
between the two ends with opposite signs, coherently combined with its
all-spins-flipped copy. Then we read the two end spins together, comparing
how often their up-or-down outcomes agree or disagree. Their correlation
contains a slowly fading oscillation whose amplitude stays finite as the
coupling change becomes small. Making this contribution longer lived does
not make its starting amplitude disappear. The flipped copy is not what
this end reading needs: the shared excitation alone gives the same
correlation. What the copy changes is which small groups of spins can see
the direction of the internal motion, the subject of the last paragraph
below.

Other measurements can miss the same contribution. This gives us three
things to keep track of: which motion the dynamics allows, whether our
starting state excites it, and whether the chosen measurement can read it.
That connection is what we want to learn from this small system. Sections
1–3 work through it; the later sections explore the missing signals and give
an exact calculation we can run ourselves.

At equal couplings we can also compare two moments that look identical to
every two-spin measurement. The two copies carry the internal current in
opposite directions, so on any pair of spins they cancel. Any third spin
tells the copies apart, because it points up in one and down in the other;
read together with a neighbouring pair that carries the current, the end
pair (0,1) for instance, it reveals opposite internal currents at the two
moments. Following the end correlation over time gives
the same information through its slope. Section 7.1 connects these two ways
of reading motion: a snapshot and a time series.

## Abstract

We connect a slow decay rate to a physical readout in an open seven-spin XY
chain with Z-dephasing of strength γ > 0 only at its centre. The left end
coupling is J₀ = 1+ε and the others are 1. For the coherent end preparation
defined in §1, we project the initial state onto the two conjugate slow
eigenspaces of the internal motion block A, each of dimension two. The
endpoint correlation Z₀Z₆ has a complex residue tending to −1/2, so the corresponding real
oscillation has an amplitude tending to 1 as ε → 0. Its envelope decays on
the already derived scale τ = 1/Δ_A ∼ 2/(γε²).

The same projection explains different outcomes for other readouts. Three
odd-site spin pairs miss this slow contribution exactly throughout the local
branch. At ε = 0 four even-site pairs miss it as well, and away from that
point they see it only at order ε; the leakage residue starts at order ε².
At equal couplings, every two-spin readout loses the sine component of this
slow contribution while the ideal decoder retains it. At that uniform point
Y₀X₁Z_k recovers the sine component for any spin k outside the pair, the
centre k = 3 and the neighbour k = 2 alike; three sites are minimal for an
instantaneous readout before decoding, and the flipped copy is what makes
them so. On the full uniform trajectory the expectation of Y₀X₁Z₃ is one
eighth of the time derivative of ⟨Z₀Z₆⟩.
The derivation and exact checks connect these outputs to the prepared state. The result
is local in ε at fixed γ; how it governs the nonlinear distances d_out and
d₂ from a continuing unitary reference at γ = 0 remains a further question.

The named-store search of `docs/ANALYTICAL_FORMULAS.md` returned F70's
partial-trace selection rule, F157's blind seat and F158's two-end count.
`docs/proofs/` returned the [N=7 relaxation proof](PROOF_MISSING_PHASE_RELAXATION_SCALE.md),
the [blind-seat proof](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md), the
[chiral trajectory proof](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md), and the
[node-pair resolvent proof](PROOF_NODE_PAIR_RESOLVENT.md), whose §8 reads the
second-order rates of the same blind dyads this residue is built on.
`experiments/`, including null results and prior hardware flights, returned
the owning preparation and physical maps, the separate N=4 readout problem
in [Route-B virtual readout](../../experiments/ROUTE_B_N4_VIRTUAL_READOUT.md),
and the [price-pair flight](../../experiments/PRICE_PAIR_HARDWARE_PREDICTION.md).
Neither that sweep nor searches of both `simulations/framework/confirmations.py`
and `compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs` returned
a hardware confirmation of this N=7 readout. Their price-pair and F120 flights
are other uses of multi-spin correlations, not tests of this residue.
`docs/GLOSSARY.md` supplies the distinctions between blindness, a site and
light content; `docs/CAUGHT_ERRORS.md` supplies the fixed-point preparation
and Riesz-versus-eigenvector warnings. The OpenArcs registry's
`relaxation_scale_as_the_defect_vanishes` is retired for its spectral question
and identifies physical coupling as a separate continuation.

Both typed stores were searched. Core owns
[`F70DeltaNSelectionRulePi2Inheritance`](../../compute/RCPsiSquared.Core/Symmetry/F70DeltaNSelectionRulePi2Inheritance.cs)
and [`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs);
Diagnostics owns
[`MissingPhaseRelaxationScaleClaim`](../../compute/RCPsiSquared.Diagnostics/Foundation/MissingPhaseRelaxationScaleClaim.cs),
its spectral/light witness, and
[`MissingPhaseOnsetWitness`](../../compute/RCPsiSquared.Diagnostics/Foundation/MissingPhaseOnsetWitness.cs).
The adjacency check, including complete reads of the relaxation proof and
owning experiment, joins the former's rank-two Riesz space to the latter's
preparation and physical maps. It also joins this proof to the node-pair
proof, which works on the same N = 7 end-bond family: there the two dyads
carry their decay rates, here their preparation and readout residues (§2),
and both use the same tridiagonal cofactor formula, there for a Green's
function, here for an adjugate at the root (§5).
The other readout and hardware entries do not own that composition. The
new claim has the relaxation claim, F70 and `ChiralKClaim` as its three
direct typed parents. The generic cofactor
derivation belongs to this proof; the live witness reconstructs the uniform
residue and its completeness, the physical maps with their null tiers, the
three-spin reading and the leakage derivative, as §8 describes.

## 1. System and physical preparation

Sites are zero-based. Take N = 7, H = Σⱼ Jⱼ(XⱼXⱼ₊₁ + YⱼYⱼ₊₁),
J₀ = r = 1+ε and J₁ = ⋯ = J₅ = 1, with Z-dephasing only at c = 3.
In the one-excitation basis |0⟩,…,|6⟩ the real Hamiltonian hᵣ has zero
diagonal, h₀₁ = h₁₀ = 2r and every other neighbouring hopping equal to 2.
Set

    P_c = |3⟩⟨3|,                 z = I − 2P_c,
    L_A(X) = −i[hᵣ,X] + γ(zXz−X),
    d = (|0⟩−|6⟩)/√2,           A_init = B_init = dd†.

With F = X⊗⁷ and W mapping the two seven-dimensional copies to |j⟩ and
F|j⟩, the prepared physical state is

    |Ψ_init⟩ = (|d⟩+F|d⟩)/√2,
    ρ(t) = ½ W [[A(t),B(t)],[B(t),A(t)]] W†.

The two copies are orthogonal. B evolves under
L_B(X) = −i[hᵣ,X] − γ(zXz+X). F70 removes B from every two-site marginal
because its excitation-number difference is 5. B must still be retained for
the full decoder and for comparisons with the unitary reference.

The relaxation proof owns

    vᵣ = (1,0,−r,0,r,0,−r)ᵀ/√(1+3r²),     P_v = vᵣvᵣ†,
    Kᵣ = −ihᵣ − 2γP_c,
    hᵣvᵣ = 0,                              zvᵣ = vᵣ.

For each fixed γ > 0 there is a neighbourhood of ε = 0 in which the simple
K eigenvalue λ₋ continuing −2√2i and its conjugate λ₊ each embed twice in A.
In a sufficiently small punctured neighbourhood they give its slowest
nonstationary rate,

    Δ_A = −Re λ₋ = (γ/2)ε² + (γ/2)ε³ + O_γ(ε⁴).

The radius need not be uniform in γ. Neither this result nor the readout
below is a claim about the complete 4⁷ Liouvillian or all odd N.

## 2. Project the preparation onto the complete cluster

The decomposition span(vᵣ) ⊕ span(vᵣ)⊥ is preserved separately by hᵣ and z,
so all four operator blocks are L_A-invariant. Let E₋ and E₊ denote the
simple K Riesz projectors at λ₋ and λ₊. The two cross blocks give the full
rank-two A Riesz projection at λ₋:

    ℘₋ᴬ(X) = E₋ X P_v + P_v X E₊†,
    M₋ = ℘₋ᴬ(A_init),                       M₊ = M₋†.

The relaxation proof's complete peripheral census ensures there is no
additional direction in this isolated cluster. Its two exact embeddings
also exclude an internal Jordan term. Thus its trajectory contribution is
e^(λ₋t)M₋ + e^(λ₊t)M₋†, with no polynomial factor in t.
For pᵣ(x) = det(xI−Kᵣ),

    E₋ = adj(λ₋I−Kᵣ)/pᵣ′(λ₋).

This is generally an oblique projector. The orthogonal right/left projectors
used by the spectral witness to read centre light do not project a prepared
state onto its spectral residue.

At r = 1 use the reflection-odd basis bⱼ = (|j⟩−|6−j⟩)/√2, j = 0,1,2.
Write v_* = vᵣ|ᵣ₌₁, avoiding a collision with the prepared d. Then

    v_* = (b₀−b₂)/√2,
    e₊ = (b₀+√2b₁+b₂)/2,      e₋ = (b₀−√2b₁+b₂)/2,
    h₁e₊ = 2√2 e₊,            h₁e₋ = −2√2 e₋,
    d = b₀ = e₊/2 + v_*/√2 + e₋/2.

The preparation therefore gives

    M₋(0) = (e₊v_*† + v_*e₋†)/(2√2),
    ‖M₋(0)‖²_F = 1/4.

Each of the two orthogonal dyads carries coefficient 1/(2√2). Keeping one
dyad is a rank-one mutation of this projection, not an alternative
normalization of the same prepared state.

The vectors e₊, v_* and e₋ are the three centre-blind modes of the uniform
chain, at energies 2√2, 0 and −2√2, so both dyads belong to the block of
blind dyads at frequency 2√2 that
[the node-pair proof](PROOF_NODE_PAIR_RESOLVENT.md) §8 treats. Its end-bond
light coefficients (§6 there) are 1/4 for e± and 0 for
v_*, and its second-order rate 2γ(c_i + c_j)ε² is therefore (γ/2)ε² for
both dyads, the leading term of the Δ_A of §1, which the relaxation proof
derives exactly at N = 7. Those two proofs give these dyads their rate; this
proof gives the same dyads their preparation and readout residues.

## 3. Endpoint correlations see an order-one amplitude

For the physical state of §1,

    ⟨Z₀Z₆⟩ = 1 − 2(A₀₀+A₆₆).

The product is unchanged by flipping both spins. Its linear action on a
traceless residue is f_ZZ(M) = −2(M₀₀+M₆₆). Since both diagonal entries of
M₋(0) are 1/8,

    f_ZZ(M₋(0)) = −1/2.

The isolated slow contribution to this expectation value is exactly

    2 Re[f_ZZ(M₋(ε)) e^(λ₋(ε)t)],

with envelope amplitude 2|f_ZZ(M₋(ε))| e^(−Δ_A t). Analyticity of the
isolated Riesz projectors and vᵣ makes the residue analytic in ε. Its nonzero
limit gives a smaller neighbourhood in which it never vanishes, and
2|f_ZZ(M₋(ε))| → 1. This proves the envelope scale for this prepared linear
readout. The complete signal also has a stationary part and other A modes;
the oscillation has zero crossings even when its envelope is nonzero.

The physical X₀X₆ map gives

    f_XX(M) = M₀₆+M₆₀,         f_XX(M₋(0)) = −1/4,

so its paired amplitude tends to 1/2. Single-site Z expectations vanish by
the equal flip copies. The global F expectation is the separate exact
e^(−2γt) phase readout, not this slow motion.

## 4. Complex-linear pair maps and the null pairs

For a complex traceless operator M the physical two-site map is

    s = M_aa + M_bb,            v = (M_ab+M_ba)/2,
    Φ_ab(M) = [[−s/2, 0,   0,   0],
               [0,    s/2, v,   0],
               [0,    v,   s/2, 0],
               [0,    0,   0,  −s/2]].

Tracing the two physical copies gives this complex-linear map. The familiar
v = Re M_ab is valid only after forming a Hermitian contribution or state
difference. Applying that shortcut to one complex residue would destroy
linearity and can fabricate or erase a quadrature.

Let C = diag((−1)ʲ). The zero-diagonal XY hopping gives ChᵣC = −hᵣ and
CP_cC = P_c, hence CKᵣC = conjugate(Kᵣ). For real ε this pairs the simple
branches: E₊ = C conjugate(E₋) C. Since Cd = d and Cvᵣ = vᵣ, with
w = E₋d and a = vᵣᵀd the complete residue is

    M₋ = a(wvᵣᵀ + vᵣwᵀC).

Its mixed even/odd cells are antisymmetric, so on a pair of one even and one
odd site the symmetric coherence (M_ab+M_ba)/2 vanishes and only the
diagonal sum reads the residue.

The odd/odd block needs no pairing. Every entry (a,b) of E₋XP_v + P_vXE₊†
is a sum of a term carrying (vᵣ)_b and a term carrying (vᵣ)_a, and vᵣ has
only even-site entries: it is the zero mode of the bipartite path, which
lives on the larger sublattice, the four even sites (the node-pair proof §4
records the same fact for a single-bond detuning). Hence

    Φ₁₃(M₋) = Φ₁₅(M₋) = Φ₃₅(M₋) = 0

exactly throughout the local isolated branch. These three physical pages
miss this whole slow contribution; this says nothing about their other
modes. F70 removes B first, while the blind vector's sublattice support
removes this A residue. These are distinct cancellations.

At the uniform point four more pages are silent. There e₊ and e₋ agree on
the even sites, both equal to u/√2 with u = (b₀+b₂)/√2, so the even/even
block of M₋(0) is (uv_*ᵀ + v_*uᵀ)/4. On sites 0, 2, 4, 6 the vector u is
proportional to (1, 1, −1, −1) and v_* to (1, −1, 1, −1). A pair map of
this block vanishes exactly when one of the two vectors takes equal values
on the pair and the other opposite values, which selects (0,2), (0,4),
(2,6) and (4,6); the other fourteen, the twelve mixed pairs and (0,6) and
(2,4), are not null (the witness computes all 21). At ε = 0 the null set is
therefore the seven pairs

    (0,2), (0,4), (1,3), (1,5), (2,6), (3,5), (4,6),

and the same seven annihilate the cosine quadrature of §7. Away from the
uniform point the four even pairs see the residue at order ε, with the
coefficients of §6, none of which vanishes at any γ. For every sufficiently
small ε ≠ 0 the three odd pairs are therefore the only null pairs. The
fourteen other pairs, the endpoints among them, read the residue at order
one.

## 5. Leakage residue from the endpoint cofactors

Let R|j⟩ = |6−j⟩ and Q_R = (I+R)/2. The experiment's leakage is
Tr(Q_R A), a population outside the original reflection-odd space. Put

    T(x) = x⁴ + 2γx³ + 12x² + 8γx + 16.

Its slow residue is

    Tr(Q_R M₋)
      = 2(1−r²)² T(λ₋)/[(1+3r²)pᵣ′(λ₋)].

Here is a finite polynomial route to the numerator, without an inverse at
a singular root. In D(x) = xI−Kᵣ the diagonal entries are dⱼ = x + 2γδⱼ₃,
and the adjacent entries are t₀ = 2ir, t₁ = ⋯ = t₅ = 2i. Define leading and
trailing continuants

    L₀ = 1,  L₁ = d₀,
    Lₖ₊₁ = dₖLₖ − tₖ₋₁²Lₖ₋₁                 (k = 1,…,6),
    U₇ = 1,  U₆ = d₆,
    Uⱼ = dⱼUⱼ₊₁ − tⱼ²Uⱼ₊₂                   (j = 5,…,0).

Then pᵣ = L₇ = U₀ and the symmetric adjugate entries for i ≤ j are

    adj(D)_ij = (−1)^(i+j) Lᵢ Uⱼ₊₁ ∏ₖ₌ᵢ^(j−1) tₖ.

The empty product is 1. This is the tridiagonal cofactor formula that
[the node-pair proof](PROOF_NODE_PAIR_RESOLVENT.md) §4 uses for the Green's
function of a Jacobi chain, applied here to the complex symmetric Kᵣ.
Substitution in the endpoint diagonal entries gives

    adj(D)₀₀ − adj(D)₆₆ = U₁ − L₆ = 4(1−r²)T(x).

For the raw blind vector b = (1,0,−r,0,r,0,−r)ᵀ set y = adj(D)d. Applying
the same entry formula, now retaining both dyads, gives the polynomial
contraction

    Tr[Q_R (dᵀb)/(bᵀb) · (ybᵀ + byᵀC)]
      = 2(1−r²)²T(x)/(1+3r²).

Division by pᵣ′(λ₋) at the simple root proves the residue formula. This
last contraction fixes the normalization of the full rank-two residue;
the endpoint-minor difference alone would not fix its preparation overlap.

At r = 1 and λ_* = −2√2i the relaxation polynomial gives

    p₁′(λ_*) = 512,             T(λ_*) = −16 + 16√2iγ.

Since (1−r²)² = 4ε² + O(ε³), the result is

    Tr(Q_R M₋(ε)) = [(−1+i√2γ)/16] ε² + O_γ(ε³).

The leakage's slow complex amplitude is quadratic even though the endpoint
correlation's amplitude tends to one. Q_R measures a different component of
the same moving residue. This is an expansion in the defect, not the
experiment's short-time expansion of a leakage difference between two runs.

## 6. A derivative route to the leakage coefficient and the first-order pairs

An independent exact route uses only an eigenvector derivative at the uniform
point. Let e₊(ε) be the eigenvector of K at λ₋(ε) that continues e₊ (so
e₊(ε) spans the range of E₋), so e₊(0) = e₊, normalized by e₊ᵀe₊(ε) = 1. The bordered system

    [[K₁−λ_*I, −e₊], [e₊ᵀ, 0]] [e₊′, λ′]ᵀ = [−K′e₊, 0]ᵀ,
    K′ = −2i(|0⟩⟨1|+|1⟩⟨0|),

is invertible because λ_* is simple. For rational γ its entries and solution
lie in ℚ(√2,i). At the uniform point K and the preparation's blind energy
vectors give w(0) = E₋(0)d = e₊/2. The reflection-even parts satisfy

    Q_R vᵣ′ = −(|0⟩+|6⟩)/4,       Q_R w′ = Q_R e₊′/2.

Both unperturbed even parts vanish. Since C fixes Q_Rvᵣ, the chiral residue
in §4 gives Tr(Q_RM₋) = 2a vᵣᵀQ_Rw, and its quadratic coefficient is

    (Q_Rvᵣ′)ᵀ(Q_Re₊′)/√2 = −((e₊′)₀+(e₊′)₆)/(4√2).

The bordered solve yields (e₊′)₀+(e₊′)₆ = √2/4 − iγ/2, recovering
(−1+i√2γ)/16. Terms differentiating a or the scalar coefficient multiplying
e₊(ε) do not contribute, because both even parts start at order ε. This route
checks the leakage germ without recomputing the generic adjugate family.

The live witness first evaluates λ′ = e₊ᵀK′e₊ and solves the equivalent
bordered system with right-hand side (λ′I−K′)e₊ and border column +e₊.
Its extra compatibility multiplier must be zero.

The same derivatives give the first-order pair maps that §4 announces. With
s and v the diagonal sum and the symmetric coherence of §4, differentiating
M₋ at ε = 0 gives, exactly and for every γ,

    (0,2):  s′ = 1/8,              v′ = −3/16,
    (4,6):  s′ = −1/8,             v′ = 1/16,
    (0,4):  s′ = −(1+i√2γ)/8,      v′ = (1−i√2γ)/16,
    (2,6):  s′ = (1+i√2γ)/8,       v′ = (1+i√2γ)/16,

while s′ and v′ vanish on the three odd pairs. Since ‖Φ_ab‖²_F = |s|² +
2|v|², these four pages read the residue with Frobenius norms √22/16, √6/16
and √(12γ²+6)/16 per |ε| at first order, the last for each of (0,4) and
(2,6); none vanishes at any γ.

For these coefficients the witness does not assume the pairing of §4. It
also solves the bordered system at λ₊ with border e₋ for e₋′, the derivative
of the eigenvector e₋(ε) that continues e₋ and spans the range of E₊, and
assembles the derivative of E₋A_initP_v + P_vA_initE₊† from e₊′, e₋′ and vᵣ′
directly. The pairing then appears as a result, e₋′ = C·conjugate(e₊′)
exactly, the first-order face of E₊ = C·conjugate(E₋)·C. A real on-site
energy on the centre breaks ChC = −h while leaving e±, v_* and their energies
untouched. In the witness's tests it makes that check and the mixed-cell
antisymmetry fail, while the three odd pairs stay null at first order: the
odd-pair null rests on vᵣ's support, not on the pairing.

The table holds at every γ, not only at the rational γ the witness
evaluates. γ enters each of the two bordered matrices only through its
centre diagonal entry, and the cofactor of that entry vanishes: with the
centre struck out, both three-site halves carry the eigenvalue, so the
reduced bordered matrix is singular. Each bordered determinant is therefore
γ-free: for the system displayed above it is e₊ᵀ adj(K₁ − λ_*I) e₊ = p₁′(λ_*)
= 512 of §5 (−512 in the witness's form with border column +e₊, and likewise
at λ₊). The right-hand sides carry no γ either, so by Cramer's rule every
entry of e₊′ and e₋′, hence every s′ and v′ and the leakage coefficient, is affine in γ. An
affine function that matches an affine closed form at two values of γ
matches it at all of them, and the witness's tests compare at five.

## 7. The decoder retains the missing quadrature

For the uniform residue form its Hermitian quadrature

    C_θ = e^(−iθ)M₋(0) + e^(iθ)M₋(0)†
        = (cos θ/2)(b₀b₀†−b₂b₂†)
          − (i sin θ/2)(b₁v_*†−v_*b₁†).

At cos θ = 0 all diagonals and all symmetric off-diagonal sums vanish.
Every one of the 21 physical two-site maps therefore annihilates C_θ.
The ideal decoder U = ∏ⱼ≠₃ CNOT(3→j), followed by tracing the centre, instead
has the complex-linear map from experiment §8.3:

    D(A,B) = Σᵢ,ⱼ≠₃ A_ij |i⟩_out⟨j| + A₃₃|f⟩_out⟨f|
             + Σⱼ≠₃ (B₃ⱼ|f⟩_out⟨j| + Bⱼ₃|j⟩_out⟨f|).

Here |f⟩_out is the all-ones outer state. Every vector in C_θ has zero
centre component, so D(C_θ,0) preserves its whole nonzero block. In the
three-dimensional blind space its characteristic polynomial is

    x[x²−(cos²θ+sin²θ)/4] = x(x²−1/4).

Thus its eigenvalues are +1/2, −1/2 and zeros, and

    ½‖D(C_θ,0)‖₁ = 1/2

for every θ. C_θ is a traceless operator contribution, not a density
matrix. This isolated uniform-limit norm is not the full d_out of
experiment §7, the distance between a γ > 0 run and its matched γ = 0 run.
That distance subtracts a moving unitary reference before taking a norm.
Likewise d₂ takes norms and then a maximum over pages. Their nonlinear
approach-to-plateau scales remain separate.

### 7.1 Three spins read the direction of motion

The cancellation in §7 is between the two flipped copies, and the flipped
copy is the whole reason for it. On the single copy |d⟩ the pair current
Y₀X₁ reads the sine quadrature directly, on two sites: removing the flipped
copy exposes the sine quadrature on pair (0,1), a control the live witness
runs. The endpoint readings do not need the copy at all. ZZ, XX and the
leakage read only the A block, so they take the same values on one copy or
two. What the copy adds is the cancellation, and with it the need for a
third site.

Any spin outside the pair supplies the missing distinction. On the (0,1)
coherence copy one has every other spin up and copy two has every other
spin down, so for each k ∈ {2, …, 6} the value of Z_k is +1 on one copy and
−1 on the other. Take the centre,

    O = Y₀X₁Z₃.

For n₀ = (I−Z₀)/2 and the Hamiltonian of §1 at J = 1, the outward bond
current is j₀→₁ = −i[H,n₀] = −(X₀Y₁−Y₀X₁). This is twice the normalized
adjacent-bond [current already used in the repository](../carbon/BENZENE_THREE_DEPHASE_LETTERS.md#y-axis-selected-single-site-model-axis-tier-4-candidate).
Both the current and Z₃ change sign under the global flip. Their product
does not, so its reading adds the two copies instead of cancelling them.
The unconditioned current remains zero. Y₀X₁Z_k gives exactly the reading
of O for every k = 2, …, 6, so the neighbouring triple Y₀X₁Z₂ serves as well
as the centre.

The residue from §2 gives

    f_O(M₋(0)) = Tr[O · ½W diag(M₋(0),M₋(0))W†] = i√2/8,
    Tr[O · ½W diag(C_θ,C_θ)W†] = sin θ/(2√2).

Thus O reads the uniform sine component and gives zero on the cosine
component. Since every one- or two-spin reduction of C_{π/2} is zero by §7,
three sites are necessary and sufficient for an instantaneous readout of
this component before decoding. Of the 35 three-site reductions, 24 are
nonzero: exactly those holding site 1 or 5 together with an even site.
C_{π/2} couples {1, 5}, where b₁ lives, to {0, 2, 4, 6}, where v_* lives,
and on three sites the two copies land on different local basis states
instead of cancelling.

The same distinction occurs between physical states of the original
preparation, not only between traceless operator contributions. At ε = 0 set

    u = (b₀+b₂)/√2,
    d(θ) = (v_* + cos θ · u − i sin θ · b₁)/√2,
    θ = 2√2t,               A(t) = d(θ)d(θ)†,
    B(t) = e^(−2γt) A(t).

Here d(0) = b₀ is the original end preparation. The hopping matrix obeys
i∂ₜd = h₁d, and z d = d, so the centre dephasing leaves A untouched.
The state ρ = ½W[[A,B],[B,A]]W† is a convex mixture of the symmetric and
antisymmetric coherent copies of d, with weights (1±e^(−2γt))/2.
F70 removes B from every two- and three-spin readout because its excitation
number difference is five. For every γ ≥ 0 the full trajectory therefore has

    ⟨O⟩ = sin θ(1+cos θ)/(2√2)
         = sin θ/(2√2) + sin(2θ)/(4√2),
    ⟨Z₃j₀→₁⟩ = 2⟨O⟩,       ⟨Y₀X₁⟩ = 0,
    ⟨Z₀Z₆⟩ = 1 − (1+cos θ)²/2,
    d⟨Z₀Z₆⟩/dt = 8⟨O⟩.

The full O signal contains a second harmonic in addition to the isolated
slow contribution. The slope has an operator reason. The centre dephasing
commutes with Z₀Z₆, and

    i[H, Z₀Z₆] = 2(Y₀X₁ − X₀Y₁)Z₆ + 2(X₅Y₆ − Y₅X₆)Z₀,

two end-bond currents, each tagged by the far end's Z. On this trajectory
Z₆ tags the copies as Z₃ does, X₀Y₁Z₆ reads −⟨O⟩, and the reflection
symmetry of d(θ) makes the right-end term equal the left one, which gives
8⟨O⟩. The slope equality therefore holds on this trajectory with J = 1; it
is not an operator identity for arbitrary states.

At the two positive times t₊ = π/(4√2) and t₋ = 3π/(4√2),
d(t₊) = (v_*−ib₁)/√2 and d(t₋) = (v_*+ib₁)/√2. Their A blocks differ
by 2C_{π/2}. Their B blocks have different decay factors, but are invisible
to all pairs by F70. Hence all 21 pair density matrices agree between these
times, while the 24 three-site reductions named above differ and the
three-spin reading changes sign:

| Reading | t₊ | t₋ |
|---|---:|---:|
| ⟨Z₀Z₆⟩ | 1/2 | 1/2 |
| ⟨Y₀X₁⟩ | 0 | 0 |
| ⟨Y₀X₁Z₃⟩ | √2/4 | −√2/4 |
| ⟨Y₀X₁Z₂⟩ | √2/4 | −√2/4 |

![Equal two-spin snapshots and opposite three-spin motion readings](../../simulations/results/missing_phase_three_spin_readout.png)

Three spins give an instantaneous reading of this internal-current
correlation; a two-spin time series supplies it through the slope. To read O,
prepare afresh at the chosen time, measure Y on site 0, X on site 1 and Z on
site 3, and average the product of the three outcomes; Z on the neighbouring
site 2 gives the same average. This does not read
the separately fading phase between the copies: that B phase retains its
five-spin threshold. The minimal-support statement here is at ε = 0; no
finite-defect or all-N minimum is asserted.

## 8. Live reconstruction and controls

From the repository root:

```powershell
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root missingphasereadout
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release --filter FullyQualifiedName~MissingPhase --nologo
```

The live witness uses exact ℚ(√2,i) arithmetic. At N = 7 it reconstructs
both uniform dyads and certifies that they are the whole cluster: the
shifted 49 × 49 generator L_A − λ_*I and its square both have rank 47, and
both dyads are also eigenoperators of L_A† at the conjugate value, so the
orthogonal projection onto them is the Riesz projection. It projects
`MissingPhaseOnsetWitness.InitialBlock(7)` and reads the endpoint residues
−1/2 and −1/4 and squared norm 1/4. Its physical-copy partial traces
compare the pair map on all 21 pairs, including complex inputs, and it
computes the seven null pairs of §4. For finite defect it takes the kernel
of hᵣ at r = 4/3 by exact elimination and finds one ray with no odd-site
entry. Its quadrature checks compare exact matrices; its bordered 8 × 8
solves reconstruct the leakage germ and the first-order pair coefficients
of §6. It evaluates §7.1 bit by bit in the 2⁷ space: the residue under
every tag Z_k, k = 2, …, 6, the one-copy controls, the trajectory at seven
exact points of the phase circle (a polynomial of degree two in cos θ and sin θ
that vanishes at five points of the circle vanishes on all of it), the
commutator identity on all 128 basis words, the two snapshots with their
24 separating three-site reductions, and the F70 threshold of the
intercopy block. The CLI's canonical instance uses γ = 3/10 and lists every
exact check with its result; the tests repeat the reconstruction at γ from
1/100 to 100. Equality checks have residual exactly zero, with no floating
tolerance.

The separating controls change the preparation, the projection or the
readout, and the witness runs each of them:

- Preparing P_v at the uniform point makes the slow residue exactly zero.
- Omitting one dyad changes the endpoint ZZ residue from −1/2 to −1/4.
- The sine quadrature kills every pair map while the decoder retains norm
  1/2. Using Re M_ab instead of (M_ab+M_ba)/2 fails the linear-map
  comparison in 168 entries, two for each of the 84 off-diagonal inputs of
  the complex basis.
- Removing the Z tag, or replacing it by X, reads zero on the residue;
  removing the flipped copy lets the untagged pair read what O reads.

One more changes the claim rather than the object: putting the endpoints
(0,6) into the claimed null set in place of the odd pair (1,3) is rejected,
which the computed null set already implies.

Three mutations run in the witness's tests, through two entry points that
are not physical knobs. The witness's constructor accepts a real on-site
energy on the centre, which no physical caller sets: it breaks ChC = −h
while leaving the uniform dyads intact, and the first-order pairing check
and the mixed-cell antisymmetry fail while the odd pairs stay null at first
order. The blind-ray gates, which the witness runs on the physical chain,
also run on a mutated one. A bond between two even sites leaves hᵣ at
r = 4/3 with no kernel. Hopping entries h₀₂ = 1 and h₄₆ = −3/4 leave one
ray, (−3/4, −3/8, 1, 3/4, −1, −3/8, 1), which the dimension gate passes and
the odd-site gate rejects.

A separate [SymPy check](../../simulations/missing_phase_three_spin_readout.py)
keeps cos θ and sin θ symbolic for the centre-tagged reading, and a
[formula plot](../../simulations/plot_missing_phase_three_spin_readout.py)
draws its formulas:

```powershell
python simulations/missing_phase_three_spin_readout.py
python simulations/plot_missing_phase_three_spin_readout.py
```

The check writes its [exact results](../../simulations/results/missing_phase_three_spin_readout.json).
It covers the trajectory's Schrödinger identity and its two decay rules, the
Y₀X₁Z₃ residue and quadratures, the trajectory reading ⟨O⟩, the reversed
current X₀Y₁Z₃, the untagged, X-tag and one-copy
controls, the 21 pair snapshots and one separating triple, (0,1,3),
positivity by the explicit convex decomposition, the current normalization,
and the slope through the full commutator i[H, Z₀Z₆]. The plot draws the
derived formulas. The tags Z_k with k = 2, 4, 5, 6, the 24 separating triples, the
two-term form of i[H, Z₀Z₆] and the one-copy ZZ and XX readings are the
live witness's alone.

The proof owns the local analytic quantifiers, the chiral branch identity
at every ε and the generic cofactor formula. The live checks reconstruct
exact finite objects and the first ε-derivatives; they do not compute a
generic symbolic adjugate family or certify a finite ε interval. The original
spectral arc remains retired. All-odd-N
extensions, γ-uniform neighbourhood control and the nonlinear d_out/d₂
lifetimes remain outside this result.
