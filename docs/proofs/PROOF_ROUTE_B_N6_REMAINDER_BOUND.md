# A certified local remainder at the N6 crossing

The selected F163 crossing `N6-E-A2-T-007` now has an explicit, conservative
complex neighborhood on which both EP-location branches are holomorphic and
their corrected Taylor remainder is bounded. This is a local computation
for one crossing, not a radius uniform over the 266-locus family.

## What the repo already held

The stores were searched by this proof's primitives: a certified parameter disk, a
Banach contraction with a fixed preconditioner, and a Taylor remainder bound at one
Route-B crossing.

The contraction design is not new here, and its source is in the same arc. The N=4
certificates already choose a fixed nonsingular matrix Y, bound the row sums of
I − Y·∂E/∂x below one, keep the Newton image strictly inside the chosen radius, and
close with Banach's theorem in outward-rounded high-precision balls:
[`ROUTE_B_N4_SECOND_ARM`](../../experiments/ROUTE_B_N4_SECOND_ARM.md) over 197 tiles
on z, and [`ROUTE_B_N4_RANGE`](../../experiments/ROUTE_B_N4_RANGE.md), whose
inclusion inequality is the one restated below. What is new is the object the design
is applied to, and it is why nothing could be imported: the N=4 certificates run real
balls on a scalar polynomial system, while the crossing here is a complex-ball problem
on the 90-dimensional Schur-complement pencil with an eta = epsilon^2 triangular
desingularization for the opposite-end profile.

[F163](../ANALYTICAL_FORMULAS.md) states the gap this proof fills, that "the algebraic
corollary alone supplies no remainder bound or explicit ε radius". Two stores already
carry this result's radii by value, the registry and the OpenArcs arc
`diabolic_over_higher_n`; the parent proof [F163](PROOF_ROUTE_B_N6_UNFOLDING.md) and the
locus atlas name the disks without their sizes. All four say the family theorem supplies
no uniform radius, and nothing in the repo claims one.

The remaining stores return nothing. `docs/proofs/` holds no other ball or interval
arithmetic anywhere. There is no NULL result on uniform radii or all-N extension; what
exists instead is a chain of stated non-extensions, and this proof's fences agree with
every one of them. `docs/GLOSSARY.md` defines no term for a certified parameter radius;
its one "radius" entry, the F135 record radius, is a different object. `docs/CAUGHT_ERRORS.md` holds nothing on radius arithmetic. `fw.Confirmations`
holds no hardware measurement touching this.

The typed layer holds this result on its claim side only.
`RouteBN6RemainderBoundClaim` carries the statement with F163 as its single typed parent,
and `RouteBN6RemainderBoundClaimTests` reads the certified radii, the η orders, ρ, the
worst-branch remainders, every Newton image and every producer and premise digest back out
of the committed certificate, so a producer rerun that moves a radius, or a premise that
moves underneath it, turns the gate red instead of leaving this page ahead of its evidence:

```powershell
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release --filter "Category=ROUTE_B_A2_N6_REMAINDER"
dotnet test compute/RCPsiSquared.Runtime.Tests -c Release --filter "Category=ROUTE_B_A2_N6_REMAINDER"
```

No witness
recomputes the balls at inspect time, which is the honest state: the ball arithmetic lives
in the Python producer, and under the repository's C# witness rule that producer is a
witness waiting to be ported. `RouteBN6ResidualProducerTests` gates the even sector's
dimension 45 and the residual's Λ-degree 32, the algebraic base this proof starts from.

## Statement and conventions

Use the full 90D (ket weight 1, bra weight 2) block, N=6, XY (Delta=0),
uniform gamma=1, q=qCSharp, t=i q and physical lambda. The three internal
bonds are 1. End strengths are one=(1+epsilon,1),
even=(1+epsilon/2,1+epsilon/2), odd=(1+epsilon/2,1-epsilon/2).
Put eta=epsilon for one/even and eta=epsilon^2 for odd.

The [ball producer](../../simulations/route_b_n6_remainder_ball.py), using
512-bit complex ball arithmetic, certifies the following eta radii R_j.
Labels match the leading slopes in the existing numerical study.

| Profile | Branch | Certified R_j in eta | Contraction factor, upper rounded | Newton image radius, upper rounded |
|---|---:|---:|---:|---:|
| One | 0 | 2^-20 | 0.285019 | 0.001874215 |
| One | 1 | 2^-18 | 0.286145 | 0.001908126 |
| Even | 0 | 2^-18 | 0.289475 | 0.003865486 |
| Even | 1 | 2^-16 | 0.290469 | 0.002947518 |
| Odd | 0 | 2^-20 | 0.197320 | 0.001756499 |
| Odd | 1 | 2^-20 | 0.197354 | 0.002383162 |

Every image radius is strictly less than rho=2^-8=0.00390625. The certificate
uses outward enclosures, not the rounded table entries. In each branch there
are holomorphic functions z(eta), h(eta) satisfying

    q(eta)=q0+eta z(eta),     lambda(eta)=lambda0+eta h(eta),
    |z(eta)-z_center| <= rho, |h(eta)-h_center| <= rho.

For nonzero eta these are EP2s of the full block. At zero the original
crossing remains semisimple, as proved in F163. The first two q coefficients
are the exact c=z(0), d=z'(0) of the preceding graph derivation.

For 0<|eta|<R_j, set s=|eta|/R_j. The explicit relative remainder bound is

    |q(eta)-q0-c eta-d eta^2| / |q(eta)-q0|
       <= rho/(|z_center|-rho) * s^2/(1-s).                 (1)

The denominator measures the actual displacement, not the leading term or
the total q. All six lower bounds |z_center|-rho are positive. The formula
uses exact Taylor coefficients; rounding c and d adds their rounding errors.

A common radius for both branches is 2^-20 (one), 2^-18 (even), and
2^-20 in eta, hence 2^-10=1/1024 in epsilon (odd). At half the common eta
radius, conservative upper bounds for both branches' relative errors are
0.03630%, 0.03630%, and 0.04582%, respectively. In the odd case this means
|epsilon| <= 1/(1024 sqrt(2)), not half of 1/1024.

These are sufficient neighborhoods. A failed larger candidate means this
enclosure failed; it does not locate a singularity or a Taylor radius.

## Algebraic base, not a rounded replacement model

The [base helper](../../simulations/route_b_n6_ball_base.py) imports the exact
rank, semisimplicity, characteristic-factor and unique repeated-root results
of [F163](PROOF_ROUTE_B_N6_UNFOLDING.md). Its fixture and producer hashes and
the stored CRT bound are checked. This helper does not rerun the CRT proof;
the exact theorem is an explicit premise of this computation.

It isolates all 133 A2 roots with complex balls and verifies that the chosen
t0 lies strictly inside the inventory's exact rational box. At this t0 it
isolates the 31 roots of the residual polynomial's Lambda derivative. Thirty
are excluded from being roots of the residual itself. F163 supplies a unique
repeated root, so the remaining critical root is that root; lambda0=Lambda0/2.
Merely finding a ball containing zero would not establish this equality.

An invertible 43x43 minor in the even 45D sector, together with the exact
nullity-two premise, gives an enclosed basis U for the exact kernel. Floating
QR chooses the minor only; a successful ball solve establishes invertibility.
After lifting U to 90D, set W=(U^T U)^(-1)U^T, P=UW, Q=I-P and

    S=Q (lambda0 I-L0+P)^(-1) Q.

Complex symmetry and semisimplicity justify the transpose dual. The reported
zero-containing residuals are consistency checks; exact identities come
from this construction and the imported theorem.

## Desingularized EP equations

In q coordinates C=iT and C_left/right=iT_left/right from F163.
For one/even write L=D+q(C+eta V), q=q0+eta z and lambda=lambda0+eta h.
Then L-L0=eta Ehat, where

    Ehat=z C+q0 V+eta z V,
    R=[lambda0 I-L0+P+eta(h I-Q Ehat Q)]^(-1),
    M=h I_2-W Ehat U-eta W Ehat Q R Q Ehat U.

R is enclosed and invertible over the entire parameter/variable box, not
just at a collection of points. Its Q block is the complementary resolvent;
the extra invertible P block is immaterial under the Q sandwiches. The
normalized Schur matrix M is regular at eta=0. Solve the two equations

    F(z,h;eta) = (det M, partial_h det M) = (0,0).

The producer differentiates this expression analytically, including the
dependence of R on h and z. Constants such as W C U and Q V Q are multiplied
before interval parameter substitution to reduce repeated-variable inflation.
No off-block floating residual is silently rounded to zero.

## The odd profile is analytic in eta without a square-root cut

Let E and O be the reflection-even and reflection-odd projections and
V=(C_left-C_right)/2. Set V_low=O V E and V_high=E V O.
For epsilon nonzero the similarity with T_epsilon=E+epsilon O gives

    T_epsilon^(-1) L(q,epsilon) T_epsilon
      = D+q C+q V_low+eta q V_high,       eta=epsilon^2.

This transformed pencil is analytic at eta=0 even though the similarity
itself is singular there. Its base is Lbar0=L0+q0 V_low, with

    Ubar=U+q0 S V U,  Wbar=W,  Pbar=Ubar Wbar,  Qbar=I-Pbar.

The odd diagonal sector excludes lambda0 by F163, so this triangular base
has the same semisimple multiplicity two. The displayed Ubar and Wbar are
its right and left spectral bases. Use the preceding equations with
Cbar=C+V_low, V=V_high and the corresponding base_matrix/Pbar/Qbar.
For nonzero epsilon similarity carries its EP2s back to the original block.
No transpose-Gram reconstruction is used for this nonsymmetric transformed base.

## Uniform contraction and holomorphy

Let x=(z,h), choose exact binary centers x_center near the known leading
roots, and let Y be the fixed midpoint of an enclosed inverse Jacobian at
eta=0. The code checks Y is invertible. It evaluates F(x_center;eta) and
J=partial_x F on |eta|<=R_j and the polydisk |x_i-x_center,i|<=rho.
Complex rectangles enclose these disks during arithmetic; the norms below
are complex moduli and hence apply to the disks themselves.

For each row i it verifies, with outward bounds,

    k_i=sum_j |(I-Y J)_ij| < 1,
    |(Y F(x_center;eta))_i|+k_i rho < rho.

Thus x -> x-Y F(x;eta) is a uniform contraction mapping the polydisk strictly
inside itself. It supplies a unique root there for every eta in the closed
disk. Analyticity of F and invertibility of J give a holomorphic root branch;
the strict margins allow continuation to a neighborhood of the closed disk.

At the root F_2=partial_h det M=0. Invertibility of J then implies
F_1,z and F_2,h are nonzero. The first excludes M=0; the second makes the
lambda multiplicity exactly two. Hence M has rank one, and invertibility of
the complementary block gives geometric multiplicity one in the full block
for eta nonzero. This establishes EP2, not merely a repeated eigenvalue.

## The remaining Taylor tail

Expand z(eta)=c+d eta+sum_(n>=2) a_n eta^n. Since z-z_center is holomorphic
and bounded by rho on the circle, Cauchy's inequality gives
|a_n|<=rho/R_j^n for n>=1. Summing n>=2 and multiplying by |eta| yields

    |q-q0-c eta-d eta^2| <= rho |eta| s^2/(1-s).

Meanwhile |q-q0|=|eta| |z(eta)| >= |eta|(|z_center|-rho). This proves (1).
It controls the entire remaining series, not just the next observed order.

## Reproduction, controls and scope

PowerShell from the repository root:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python simulations/route_b_n6_remainder_ball.py
```

The [JSON](../../simulations/results/route_b_n6_remainder_ball.json) records
successful and failed candidate radii, outward bounds, centers, source hashes,
precision, and the Cauchy reading. The program raises if any branch has no
successful radius. Moving the first z center by +1 fails the same inclusion
test; it is retained as a negative control. Independent differentiation checks
the analytic Jacobian, and separate reviews check the interval and Schur logic.

The broader [complex-circle scout](../../simulations/route_b_n6_complex_remainder.py)
samples much larger circles and observes closing branches and finite scaled
remainders. Those observations are not premises of this proof. Nor do its
small certified radii contradict the accurate real-positive continuation at
epsilon=0.1: interval overestimation and the fixed contraction coordinates
are conservative. No maximal convergence radius is claimed. This result
concerns the selected complex-q crossing and these profiles; it provides no
physical-real-q operating point, different-N result, or all-locus radius.

The exact [F89d response transport](PROOF_ROUTE_B_N6_UNFOLDING.md#6-the-palindrome-transports-the-entire-response)
does carry this certificate to the selected (1,4) partner at conjugate q0
and lambda_partner=-conj(lambda0)-12. Its epsilon radii and relative
remainder bounds are identical. This is inherited through the exact mirror,
not a separately executed ball solve and not a transfer to a different N.
