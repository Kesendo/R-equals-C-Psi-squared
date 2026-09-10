# The crossing opens in two orders

**Result:** F163, exact local unfolding of the N=6 Route-B A₂ layer.

The system is an open six-site XY chain under uniform local Z dephasing. We
read its generator on the 90 coherences |a⟩⟨b| with one ket excitation and two
bra excitations. At a diabolic crossing, two independent eigenvectors share
one eigenvalue. The question is what happens when the end bonds change.

Every one of the 266 parity-labelled A₂ loci is semisimple by the algebraic
argument below. At each locus, each of the three end profiles opens two
distinct exceptional-point branches for sufficiently small nonzero ε.
Changing the ends equally already opens the crossing. Changing them
oppositely enters only at second order. Spatial reflection sorts the order
of the response; it does not preserve this coincidence.

The starting point is the [Route-B atlas](../../experiments/ROUTE_B_N6_A2_LOCUS_ATLAS.md).
The prior-work sweep searched the F-registry, proof corpus, experiments
including null results and hardware flights, glossary, caught-errors ledger,
Core claims, OpenArcs, Diagnostics and Python confirmations. F131 and
[PROOF_MIRROR_ORDER_SORTING](PROOF_MIRROR_ORDER_SORTING.md) own the reflection
selection rule. [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md)
owns the conditional residual twin-scalar discussion; it does not certify
single-Hamiltonian-multiplet descent. The reduced-resolvent and Schur-complement
methods also appear in [PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md) and
[PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md). The caught-errors ledger
warns against reading Jordan character from discriminant order alone.
Core/OpenArcs and Diagnostics supply the A₂ source and its numerical character
readings. No matching end-profile flight or confirmation was found. This is
an application of algebraic conjugacy and standard analytic perturbation
theory to those existing objects, not a new perturbation method.

## 1. The exact objects

Use γ=1, Δ=0 and the atlas parameter t=i q, q=qCSharp. Physical eigenvalues
are λ; the exported integer polynomials use Λ=2λ. The analytically continued
block is

    L(t,ε) = D + t(T + ε V),
    D[a,b;a,b] = −2 popcount(a xor b).

T has real entries −2 for a ket hop and +2 for a bra hop. Thus L(t,0) is real
symmetric at real t. For complex t this is a continued pencil, not a physical
Lindblad evolution. Site 0 is the least significant bit. Let T_L,T_R be the
contributions of bonds (0,1) and (4,5), and define

    V₊ = (T_L + T_R)/2,     V₋ = (T_L − T_R)/2.

The three end profiles, with all three interior bonds fixed at one, are

| V | Left end | Right end |
|---|---|---|
| T_L | 1+ε | 1 |
| V₊ | 1+ε/2 | 1+ε/2 |
| V₋ | 1+ε/2 | 1−ε/2 |

R reverses the six sites on both ket and bra. It commutes with D,T, preserves
V₊ and negates V₋. All reflection orbits of this coherence basis have size two,
so the R-even and R-odd blocks each have dimension 45.

Let a(t) be the exact degree-133 R-even A₂ polynomial produced by
[`route_b_a2_n6.py`](../../simulations/route_b_a2_n6.py), and let F(Λ,t) be the
degree-32 residual factor of the even characteristic polynomial. The remaining
AT factor has degree 13. Both come from the strict
[integer fixture](../../simulations/tests/fixtures/route_b_a2_n6_residual.json).
The coefficients of a are stored in the
[exact certificate](../../simulations/results/route_b_n6_exact_unfolding.json).

## 2. One Hermitian conjugate certifies all 133 even loci

The [gate](../../simulations/route_b_n6_exact_unfolding.py) establishes the
following exact premises.

1. The characteristic polynomials of the two integer 45D sectors equal their
   exported residual times AT factors. Their coefficients have t-degree at
   most 45 on both sides. Equality at the 46 integers 0,…,45 is consequently
   a polynomial identity. The matrices are formed directly from spin hops.
2. The existing bounded-CRT producer independently reconstructs and proves
   `disc_Λ F = constant · t^536 · A₁(t) · a(t)^2` over the integers.
3. a is irreducible modulo 367 and retains degree 133. The Rabin certificate
   checks `x^(367^133)=x mod a` and the two proper Frobenius gcds, for divisors
   7 and 19 of 133. Thus a is irreducible over ℚ. Its odd degree supplies a
   real root.
4. At the good reduction p=101, t=31, the gcd of F and F_Λ is linear, with
   Λ=48. The full 90D characteristic polynomial has multiplicity exactly two
   at physical λ=24. The A₂ derivative is 65, and the complementary
   characteristic factor evaluated at λ is 39, both nonzero modulo 101.

Here is why these premises prove semisimplicity, rather than merely suggest
it. Work over the field K=ℚ[t]/(a). The exact discriminant identity implies
`gcd_K(F,F_Λ)` has degree at least one. Its good modular specialization has
degree one, so its degree over K is exactly one. Equivalently, its first
subresultant is u(t)Λ+v(t), with u nonzero at this prime. The repeated
physical eigenvalue is therefore the well-defined element

    h = −v/(2u) ∈ K.

The full characteristic polynomial contains (λ−h)². It cannot contain its
cube: after localization at the good prime, such divisibility would survive
specialization and contradict the modular multiplicity two. This also
excludes another-sector or AT eigenvalue from joining the pair.

Embed K into ℝ through a real root of a. The corresponding real symmetric
matrix has a double eigenvalue h, hence a two-dimensional eigenspace. Matrix
rank over a number field is preserved by every embedding: a minor is zero
or nonzero before and after the embedding. Thus `dim ker(L−hI)=2` at **every**
root of a, real or complex. Algebraic multiplicity is exactly two everywhere
as well. All 133 even loci are therefore semisimple.

This reasoning does not infer semisimplicity from a squared discriminant.
The independent input is the Hermitian embedding, carried by irreducibility.

## 3. Exact algebraic coefficients without expanding the number field

Fix an embedding, write t₀ for its root and λ₀=h(t₀), and let P be the
spectral projector onto the double eigenspace, Q=I−P. P is not generally an
orthogonal projector. Put

    S = Q(λ₀I−L₀|ran Q)⁻¹Q,       L₀ = L(t₀,0),
    A = (PTP)|ran P,
    B₊ = (P t₀V₊ P)|ran P,
    B₂ = (P t₀V₋ S t₀V₋ P)|ran P.

These are exact matrices over K. A basis-free rational construction is
available: if `χ(x)=(x−λ₀)² g(x)`, then

    P = g(L₀)/g(λ₀),
    S = Q(λ₀I−L₀+P)⁻¹Q.

The complement need not be diagonalizable for this construction. Its full
algebraic multiplicities remain in g. The two-dimensional matrices can be
held as 90D operators zero on ran Q; their traces are then unchanged.

For B=B₊ or B=B₂ define

    δ(M) = 2 tr(M²) − tr(M)²,
    α = δ(A),
    β = 4 tr(AB) − 2 tr(A) tr(B),
    γ₀ = δ(B),
    Ω = β² − 4αγ₀.

Here γ₀ is a polynomial coefficient, not the dephasing rate. The leading
location coefficients have the closed algebraic form

    c± = (−β ± √Ω)/(2α).                           (1)

This is a quadratic extension of the already specified degree-133 field.
The trace expressions specify its coefficients exactly; a 133-coefficient
expanded representation of each field element is not required.

The modular gate evaluates this rational construction at its good prime:

| Profile | α | β | γ₀ | Ω |
|---|---:|---:|---:|---:|
| One end | 19 | 4 | 0 | 16 |
| Equal ends | 19 | 4 | 0 | 16 |
| Opposite ends | 19 | 52 | 18 | 23 |

These are residues modulo 101, not approximate real coefficients. All
inverses used by P and S exist at this prime. Hence α and Ω are nonzero
elements of K and remain nonzero at every embedding. The zero residue for
γ₀ in the first two rows makes no characteristic-zero zero claim.

## 4. Why the closed coefficients give two EPs

On the entire original plane R acts as +I, so `PV₋P=0`. One-end and equal-end
changes have the same first-order matrix because `T_L=V₊+V₋`.

The isolated two-dimensional spectral cluster admits a holomorphic graph
over the fixed plane ran P. In these fixed coordinates, with x=t−t₀, its
effective matrix K_eff has the expansions

    K_eff−λ₀I = xA + ε B₊ + O(x²,xε,ε²),
    K_eff−λ₀I = xA + ε²B₂ + O(x²,xε²,ε⁴),          (2)

for equal ends and opposite ends respectively. The one-end profile shares
the first line's leading terms. The second line is even in ε: conjugation
by R changes ε to −ε and acts identically on the reference plane. The
second-order term is the excursion into the opposite-parity complement
and its return, with the positive reduced-resolvent sign used above.

Set η=ε for the first line and η=ε² for the second, and x=ηc. The reduced
discriminant divided by η² tends to

    δ(cA+B) = αc² + βc + γ₀.

Because αΩ≠0, its two finite roots are simple. The holomorphic implicit
function theorem gives two distinct local branches

    t±(ε) = t₀ + c± ε + O(ε²)       (one end or equal ends),
    t±(ε) = t₀ + c± ε² + O(ε⁴)      (opposite ends).       (3)

Along each branch, for sufficiently small nonzero ε,

    ∂t δ(K_eff) = η(2αc±+β) + O(η²) ≠ 0.

The differential of δ vanishes at every scalar 2×2 matrix. Thus the matrix
at these simple discriminant zeros is not scalar: its repeated eigenvalue
has geometric multiplicity one. Each branch consists of EP2s, with the
usual square-root eigenvalue splitting transverse to the branch in t.
The EP **locations** themselves are analytic in ε, or in ε² for opposite
ends. In the atlas q convention, replace `t₀,c±` by `−it₀,−ic±`.

The common eigenvalue along either branch is

    λ±(ε) = λ₀ + ½ tr(c±A+B) ε^m + O(ε^(2m)),

where m=1 or 2. One-end and equal-end locations agree to first order;
their higher coefficients need not agree.

## 4a. The next coefficient links the three profiles

Write the one/equal-end branches as
`t_j=t0+c_j epsilon+d_j epsilon^2+O(epsilon^3)`, matched by their
common first coefficient c_j. Write the opposite-end branches as
`t_j=t0+a_j epsilon^2+O(epsilon^4)`. Then, as an unordered pair,

    {d_one,+ - d_equal,+, d_one,- - d_equal,-} = {a_+, a_-}.       (4)

This is an exact corollary on the same A2 loci. Multiplying all parameter
coefficients by -i carries it unchanged to the q convention. The
[next-coefficient study](../../experiments/ROUTE_B_N6_LOCAL_VALIDITY.md#the-next-coefficient-from-the-full-block)
gives a graph recurrence evaluating the coefficients without fitting paths.

Define the polar form `b(X,Y)=4 tr(XY)-2 tr(X)tr(Y)` of the discriminant.
Along `t=t0+c_j epsilon`, let the fixed-plane effective matrix be
`lambda0 I+epsilon H_j+epsilon^2 J_j+...`, with `H_j=c_j A+B_+`.
Adding the unknown next coefficient adds `d_j A` at second order. Thus

    d_j = -b(H_j,J_j)/b(H_j,A).

The denominator is nonzero because c_j is a simple discriminant root.
Since the one-end direction is V_++V_-, reflection removes every mixed
second-order term and the direct compression of c_j V_-. The odd graph
feedback has zero compression as well. Consequently
`J_one,j-J_equal,j=B_2` and

    d_one,j-d_equal,j = -b(H_j,B_2)/b(H_j,A).                     (5)

A, B_+ and B_2 are self-adjoint for the same nondegenerate **complex
bilinear** form on the spectral plane: the full pencil, spectral projector
and reduced resolvent are transpose symmetric. In coordinates U this form
is G=U^T U. Its nondegeneracy follows from spectral separation and
semisimplicity of the selected eigenvalue of the complex-symmetric L0.
A change of basis makes G=I. Traceless symmetric 2x2 matrices then have form

    [[x,y],[y,-x]],     delta=4(x+i y)(x-i y).

There are two null lines, with linear forms ell_+ and ell_-. Because alpha
and Omega are nonzero, the two nonzero traceless H_j occupy different null
lines and neither line kills A. On the line containing H_j,
`b(H_j,Y)/b(H_j,A)=ell_j(Y)/ell_j(A)`. Thus (5) is a root of
`delta(zA+B_2)=0`; the two j give its two roots, proving (4).

The common bilinear symmetry is essential. For `A=diag(1,-1)`,
`B_+=[[0,1],[1,0]]`, and `B_2=[[0,1],[-1,0]]`, both quotients in (5)
vanish whereas `delta(zA+B_2)=4(z^2-1)`. The antisymmetric last matrix
violates the shared bilinear symmetry. Reflection alone does not imply (4).

**Dimension-independent conditional form.** The argument for (4) uses no
particular matrix dimension or value of N. It applies to a complex-symmetric
affine pencil with an orthogonal reflection commuting with the unperturbed
pencil, reflected end directions V_+ and V_-, and an isolated semisimple
two-plane on which reflection is scalar, either +I or -I. Nonzero alpha and
simple even-profile leading roots are required. Nonzero odd Omega is the
additional condition for two distinct quadratic odd EP branches. These
hypotheses must be established for each new instance; the N=6 family
certificate is not an all-N proof of them. The
[N=4/5 comparison](../../experiments/ROUTE_B_OTHER_N_UNFOLDING.md) checks two
N=4 crossings and all 58 stored N=5 q loci numerically.

## 5. The odd sector and the boundary of the result

The diagonal sign operator G on |a⟩⟨b| is
`(−1)^(sum of occupations on sites 1,3,5 on both sides)`.
Direct integer identities give

    GDG=D,   GT_sG=−T_s for every bond s,   GR=−RG.

Thus `G L(t,ε) G = L(−t,ε)` for each of the three profiles, and G exchanges
the two reflection sectors. The 133 odd loci and their local unfoldings
are the images of the even ones under t↦−t. This proves the stated result
for all 266 parity-labelled loci of this particular A₂ layer.

The family theorem is local in ε and supplies no radius uniform over its loci.
For `N6-E-A2-T-007`, the [ball-arithmetic companion](PROOF_ROUTE_B_N6_REMAINDER_BOUND.md)
certifies explicit small disks and a Taylor remainder for all three profiles.
These bounds do not reach ε=0.05. The [numerical continuation](../../simulations/results/route_b_n6_validity_probe.json)
reaches ε=0.10 only for that selected crossing. The theorem does not cover
another N, a ZZ term, an arbitrary bond direction, or physical-real-q
accessibility. All three profiles retain a quadratic XY Hamiltonian.
The result therefore does not diagnose a loss of free-fermion additivity.

F131 alone supplies evenness, not a nonzero quadratic term. The nonzero Ω
for B₂ is the extra information that makes the second order actually open.
The gate distinguishes this: scalar or commuting effective responses give
Ω=0, and adding an even component destroys the odd first-order suppression.

## 6. The palindrome transports the entire response

Let B be the bra-complement permutation from the (1,2) block to (1,N-2).
For these XY end profiles and uniform gamma=1, F89d gives

    L_partner(conj(q),conj(epsilon))
       = -B conj(L_source(q,epsilon)) B^(-1) - 2N I.          (6)

Both parameters must conjugate off the real axis. On a basis coherence,
complementing the bra changes disagreement k to N-k, hence the two dephasing
entries sum to -2N. Global spin flip on the bra commutes with its XY
Hamiltonian for every bond profile. This proves the Hamiltonian part of (6)
as well; the full operator identity, not just the diagonal pairing, is used.

For a selected semisimple pair, with the reduced-resolvent convention of
this proof, the corresponding objects obey

    lambda_partner = -conj(lambda_source)-2N,
    P_partner = B conj(P_source) B^(-1),
    S_partner = -B conj(S_source) B^(-1).

The last minus sign follows by inverting the complementary restriction of
`lambda_partner I-L_partner = -B conj(lambda_source I-L_source) B^(-1)`.
The q derivative and the end perturbation each also minus-conjugate. Write
A_q=P(partial_q L)P=i A, where A in sections 1-4 uses t=i q. In matched
plane coordinates, A_q, B_+ and B_2 therefore all minus-conjugate;
the two perturbations and one resolvent in B_2 contribute three minus signs.
Their q-coordinate discriminant coefficients conjugate, since delta is quadratic.

More generally, uniqueness of the local branches gives

    q_partner(epsilon)=conj(q_source(conj(epsilon))),
    lambda_partner(epsilon)=-conj(lambda_source(conj(epsilon)))-2N.

Thus every q Taylor coefficient conjugates, whereas every nonconstant
lambda coefficient minus-conjugates. Labels must be matched by their
leading coefficients; an unordered numerical pair is not an imposed order.
In t coordinates the coefficient rule is instead minus conjugation,
because t_partner=-conj(t_source).
Each individual complementary-mode residue and its contribution to B_2
obeys the same transport. Permutation/conjugation preserve singular values:
spectral gaps, projector/resolvent norms and response magnitudes agree in
matched orthonormal coordinates. The fold does not make a large response small.

The bra complement commutes with spatial reflection. It preserves the
spatial parity of both the selected plane and its intermediate modes.
Spatial reflection instead changes the sign of the odd end parameter
inside a block. Equation (6) pairs whole responses between blocks; these
are different consequences of different operations.

Because conjugation preserves disks and absolute errors, the selected N=6
[remainder certificate](PROOF_ROUTE_B_N6_REMAINDER_BOUND.md) transfers to its
(1,4) partner with the same epsilon radii 2^-20, 2^-18, 2^-10 for one,
equal and opposite ends, and the same relative remainder bounds. This is
a consequence of the certified source and exact identity (6), not a second
independent ball computation. At N=4 the partner block is the same block;
at N=5 it is (1,3), and at N=6 it is (1,4).

The [response probe](../../simulations/route_b_mirror_response.py) checks
dyadic matrix identities and independently rebuilt partner planes at the
selected N=4/5/6 seeds. It also rejects a complex-epsilon control that
conjugates q but incorrectly leaves epsilon unchanged. Its coefficient
readings are numerical checks of this algebraic transport, not premises
of the proof or a new enumeration of partner loci.

At the N=4 self-fold, the [exact companion](PROOF_ROUTE_B_N4_SELF_FOLD.md)
checks that the two leading q coefficients are distinct and real for
each profile. Local uniqueness then fixes each EP branch individually:
real sufficiently small ε gives real q and Re λ = −4 to all orders.
This extra sign information is not implied by self-fold alone, and the
N=4 statement supplies no explicit ε radius.

## Reproduction

```powershell
python simulations/route_b_n6_exact_unfolding.py
```

The run rebuilds integer matrices, proves the source-polynomial identities,
reruns the bounded-CRT discriminant proof, executes the Rabin checks and
the modular nonzero witnesses, and writes the exact JSON certificate.
The certificate binds the wrapper, the exact source producer, every transitive
local Python dependency discovered from that producer, and the integer fixture
by LF-normalized SHA-256; the C# consumer requires the complete dependency set.
It invokes no eigenvalue solver and no floating-point rank tolerance.
The companion numerical probe evaluates (1) and checks the local paths;
its decimal values are not premises of this proof.

The typed layer consumes that certificate. `RouteBN6A2UnfoldingClaim` (F163)
carries the statement, `RouteBN6A2UnfoldingWitness` recomputes the reading at
inspect time, and both are gated under one category:

```powershell
dotnet run --project compute/RCPsiSquared.Cli -- inspect --root n6unfolding
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release --filter "Category=ROUTE_B_A2_N6_UNFOLDING"
dotnet test compute/RCPsiSquared.Runtime.Tests -c Release --filter "Category=ROUTE_B_A2_N6_UNFOLDING"
```

The witness reads the certificate rather than the proof text, so a producer
rerun that changes any bound hash fails the witness before it reaches a number.
