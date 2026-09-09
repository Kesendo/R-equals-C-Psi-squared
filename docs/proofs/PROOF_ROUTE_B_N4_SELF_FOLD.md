# Real EP branches from the N=4 self-fold

The full 24-dimensional (1,2) coherence block of the four-site open XY
chain has two individually self-mirrored EP branches under each of the
three end-bond profiles below. For sufficiently small real nonzero ε,
their couplings are real and their merged eigenvalues obey Re λ = −4.
This is a local theorem; no explicit ε radius is asserted.

The repository survey found the operator identity in
[F89d](../../experiments/F89_BRANCH_LOCUS_PALINDROME.md), the seed and its
independent semisimple character in
[F89Path3OcticEpClaim](../../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs),
and the complete response transport in
[F163 §6](PROOF_ROUTE_B_N6_UNFOLDING.md#6-the-palindrome-transports-the-entire-response).
[F131](PROOF_MIRROR_ORDER_SORTING.md) already distinguishes a multiset
symmetry from a tracked-branch statement. The F-registry, typed Core and
Diagnostics, OpenArcs, GLOSSARY, CAUGHT_ERRORS, experiments, hypotheses,
reflections, recovered, framework and MirrorWorld were searched by these
primitives. The previous N=4 end-response coefficients were numerical;
MirrorWorld owns the fold, not these EP paths. The Confirmations registry
contains a dissipator-half pairing measurement, not this EP result.

## 1. Model and exact local data

Use γ = 1, Δ = 0 and q = J/γ in H = J Σ wₛ(XX+YY), with hopping 2q.
The bond weights are (1+aε, 1, 1+bε):

| Profile | a | b | η |
|---|---:|---:|---|
| One end | 1 | 0 | ε |
| Equal ends | 1/2 | 1/2 | ε |
| Opposite ends | 1/2 | −1/2 | ε² |

The seed is

    q₀ = √((√13−1)/6),     λ₀ = −4+2iq₀.

The [exact producer](../../simulations/route_b_n4_self_fold.py) works in
Q(t), 3t⁴−t²−1 = 0, in the embedding t = iq₀. It assembles the full block
from the existing framework-backed builder, checks its integer entries,
and compares an asymmetric bond profile against an independent full
Hilbert-space hopping construction. All subsequent matrix operations are
exact number-field operations; no eigensolver or rank tolerance is used.

Write L₀ = D+tT. The producer finds nullity(L₀−λ₀I) =
nullity((L₀−λ₀I)²) = 2. Thus the eigenvalue has algebraic and geometric
multiplicity two. Its kernel basis U has invertible UᵀU. Define

    W = (UᵀU)⁻¹Uᵀ,  P = UW,  Q = I−P,
    S = Q(λ₀I−L₀+P)⁻¹Q.

Spatial reflection acts as +I on U. With E = tV for the chosen bond
direction, set

    A = WtTU,
    B = WEU                     (one or equal ends),
    B = W·E·S·E·U               (opposite ends).

Here A differentiates a *relative* q displacement, q = q₀(1+rη).
For opposite ends WEU = 0 exactly, by reflection.

Let δ(X) = 2tr(X²)−tr(X)². The normalized leading discriminant is

    δ(rA+B)/δ(A) = r²+b r+g.

The exact result is:

| Profile | b | g | b²−4g |
|---|---|---|---|
| One or equal | (663−131√13)/299 | (87−25√13)/26 | (−41072+13088√13)/6877 |
| Opposite | −(11323+557√13)/1196 | −(12373651+3475859√13)/9568 | (36210033+10053738√13)/6877 |

Also δ(A) = (7319728−9469200t²)/196249 ≠ 0. Both discriminants are
strictly positive: for the first, √13 > 7/2 already suffices; all terms
in the second numerator are positive. Hence the two roots r₋ and r₊ are
distinct and real. Closed forms for the leading q coefficients are

    c± = q₀(−b ± √(b²−4g))/2.

They are approximately −0.520880800087, +0.100646300073 for one/equal
ends and −30.148714776875, +37.494111387670 for opposite ends. The signs
are exact too: g < 0 in each case. Each profile therefore sends one
branch to lower q and the other to higher q near the seed.

## 2. Existence and fixation of individual branches

The isolated semisimple plane admits an analytic reduced 2×2 matrix.
For the odd profile it can be chosen even in ε, because reflection fixes
the plane and sends ε to −ε. Its discriminant after q = q₀(1+rη) is

    η² [δ(A)(r²+b r+g)+O(η)].

The two simple roots at η = 0 give unique holomorphic roots r±(η) by
the implicit-function theorem. Thus q± = q₀+c±η+O(η²). At nonzero small
η the discriminant zero is simple in q. The reduced matrix cannot be
scalar there, since the derivative of δ vanishes at any scalar matrix.
It is consequently a defective EP2, isolated from the other 22 modes.
This is the local reduction used in F163, with its hypotheses checked
here at N=4 rather than inherited from the N=6 family certificate.

Bra complement B_f acts inside this block and gives, bond by bond,

    L(q̄,ε̄) = −B_f conj(L(q,ε)) B_f⁻¹−8I.

It fixes the seed (q₀,λ₀) and sends the germ with leading coefficient c
to a germ with leading coefficient c̄. Since c₋ and c₊ are real and
distinct, the transformed germ must be the same germ, by uniqueness.
Therefore, individually for both branches,

    q±(ε̄) = conj(q±(ε)),
    λ±(ε̄) = −conj(λ±(ε))−8.

For real sufficiently small ε this proves q± ∈ R and Re λ± = −4.
Every q Taylor coefficient is real, and every nonconstant λ Taylor
coefficient is purely imaginary. Spatial reflection additionally makes
the opposite-end branches even in ε. A staggered site-sign gauge gives
L(q̄,ε̄) = G conj(L(q,ε)) G, carrying the statement to the conjugate
seed −4−2iq₀ with the same q branches; the producer checks this identity
for D and every bond direction as well.

## 3. What is held, and what is not

For an EP eigenvector v at these real parameters, the absorption identity
Re λ = v†Dv/(v†v) applies because the Hamiltonian commutator is
anti-Hermitian. This block has only disagreement rungs k = 1 and 3, with
rates −2 and −6. Re λ = −4 therefore forces equal squared-norm weights
on the two rungs. This concerns the eigenvector, not a positive density
matrix, and does not assert that D stays scalar on the whole two-plane
or on a Jordan chain. The eigenvalue's decay rate is pinned while the
frequency and q location may move.

Self-fold alone is insufficient. The complex-symmetric two-dimensional
family

    L(q,ε) = −4I+i [[q−1, ε], [ε, −(q−1)]]

obeys the same antiunitary identity with B_f = I. Its discriminant is
−4((q−1)²+ε²), so its two EP branches are q = 1±iε, λ = −4. They are
exchanged by conjugation and leave real q immediately. The producer
checks this control through the same discriminant formula. A separate
one-cell dissipator mutation breaks the full-block fold equation and is
rejected. The positive discriminants above are essential instance data,
not a consequence of the word “self-fold”.

The theorem does not provide a maximal real interval, a certified
numerical ε radius, a guarantee for other directions or ZZ perturbations,
or an experimental visibility claim. N=5 has real-q crossings without an
internal (1,2) bra fold; the present sufficient mechanism is not necessary.

The [range companion](../../experiments/ROUTE_B_N4_RANGE.md) tracks these
germs farther and certifies seven boundary points independently. It finds
both EP3 collisions and EP2 turning points, including the exact equal-end
turn at ε = √2−2, q = 2, λ = −4+2i. The selected negative-slope
equal-end germ has a separate
[validated forward continuation](../../experiments/ROUTE_B_N4_SECOND_ARM.md)
through this turn to infinity. Its R = +1 mode stays EP2 for all q > q₀;
two cross-parity coincidences have full-block character J₂ ⊕ J₁.
The other boundary connections remain numerical; point certificates
alone are not maximal-interval proofs.

## Reproduction

```powershell
python simulations/route_b_n4_self_fold.py
```

The [JSON](../../simulations/results/route_b_n4_self_fold.json) stores the
exact coefficients, decimal readings, control outcomes and source hashes.
