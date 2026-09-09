# How far the local crossing formula carries

**Date:** 2026-09-09

## What this is about

We have a small-change formula for how a crossing moves when the end bonds
change. How long can we trust its first prediction? This note compares it
with the full six-spin calculation. Strengthening one end while weakening
the other makes the leading formula stay accurate farther in the tested
example than changing both ends together. The question is accuracy along
this particular path, not whether the formula works everywhere. Part of
the path also uses complex coupling, so it is a mathematical exploration
rather than a sequence of laboratory settings.
We then add the next correction and ask a stricter question: within what
small region can we bound everything the truncated formula leaves out?

## Abstract

At `N6-E-A2-T-007`, the leading F163 formula stays accurate appreciably farther
along the antisymmetric end direction than along the one-sided or symmetric
ones. This is a numerical accuracy reading for one crossing, not the convergence
radius of the Taylor series and not a bound uniform over all 266 crossings.
The document also derives the next coefficient from the full-block
reduction, certifies both branches on conservative complex parameter
disks, and bounds the remaining Taylor tail there. Those certified disks
are distinct from the wider numerical accuracy scan. The relation between
the three profiles' coefficients is exact under F163 and applies to all
266 loci covered by its certificate.

## Object and measurement

N=6, full 90-dimensional (ket weight 1, bra weight 2) coherence block;
XY (Delta=0), uniform gamma=1, q=J/gamma=qCSharp and t=i q.
The complex-q continuation is not a physical Lindblad parameter scan.
The three profiles are (left,right)=(1+epsilon,1),
(1+epsilon/2,1+epsilon/2), and (1+epsilon/2,1-epsilon/2);
all three internal bonds remain 1. Only positive epsilon is scanned here.

The reference is the full-block ordered-Schur discriminant continuation in
[the producer](../simulations/route_b_n6_validity_probe.py), using the existing
[end-profile construction](../simulations/route_b_n6_end_profile_probe.py).
The leading q slopes are evaluated from the F163 reduced matrices, not fitted
to these finite-epsilon rows. For each branch the error is

    abs(q_EP(epsilon) - q0 - c_q epsilon^m) / abs(q_EP(epsilon) - q0),

with m=1 for one/even and m=2 for odd. Thus percentages measure the error in
the movement itself. The table takes the larger error of the two branches.

| epsilon | One-sided | Symmetric | Antisymmetric |
|---:|---:|---:|---:|
| 0.01 | 3.552% | 2.823% | 0.0518% |
| 0.02 | 7.179% | 5.682% | 0.207% |
| 0.03 | 10.874% | 8.573% | 0.466% |
| 0.05 | 18.439% | 14.440% | 1.298% |
| 0.075 | 28.161% | 21.914% | 2.934% |
| 0.10 | 38.097% | 29.517% | 5.249% |

For a 5% tolerance, the sampled crossing brackets are (0.01,0.02) for
one/symmetric and (0.075,0.10) for antisymmetric. These are grid brackets,
not a certified first threshold or a guarantee between samples. Individual
branches need not have monotone error: the symmetric second branch improves
again after epsilon=0.05. No global monotonicity is assumed.

The odd direction removes the linear effective perturbation. F163 gives an
absolute remainder O(epsilon^4) after its quadratic term, hence relative error
O(epsilon^2) when the leading branch slope is nonzero. The one/even remainder
is O(epsilon^2), hence relative error O(epsilon). The readings agree with
these local orders. The different profiles also move the ends by different
amounts at equal epsilon; epsilon=0.10 in the odd profile means end bonds
1.05 and 0.95.

## Isolation is a different question

At epsilon=0.10 the minimum distance from the repeated eigenvalue to the
remaining 88 eigenvalues is 0.4760 (one), 0.4461 (even), and 0.6597 (odd),
taking the smaller of each profile's two branch values. Loss of leading-term
accuracy therefore precedes any observed loss of eigenvalue separation on
this grid. These distances alone do not control a nonnormal resolvent.

For an independent sufficient isolation criterion, fix a circle about lambda0
and write E=L(q,epsilon)-L(q0,0). If

    ||E||_2 < inf_on_circle sigma_min(z I - L(q0,0)),

the Neumann argument preserves the rank-two Riesz cluster throughout the
straight operator homotopy L0+sE, 0<=s<=1. It does not prove an error bound on
the leading EP formula, uniqueness of EP locations, or a Taylor radius.
For n equally spaced nodes on a circle of radius r, singular-value Lipschitz
continuity gives the exact-arithmetic lower bound

    min_nodes sigma_min(z I-L0) - 2r sin(pi/(2n)).

Among five tested radii, r=0.4115136800 with n=512 gives 0.1513117509 after
the angular correction. The SVDs and seed are floating point, so this is a
continuum lower *estimate*, not an interval-certified bound. For both branches, the resulting
norm ratio is below one at the sampled one/even epsilon=0.001 and odd
epsilon=0.02. The test no longer covers both branches at one/even 0.01
and odd 0.03, although the numerical cluster stays isolated. Failure of this sufficient test means
only that this fixed-contour norm estimate is too conservative.

## Reproduction and limits

PowerShell, from the repository root:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python simulations/route_b_n6_validity_probe.py
```

[Machine-readable rows](../simulations/results/route_b_n6_validity_probe.json)
carry source/script hashes, both branches, residuals, isolation distances,
perturbation norms and contour estimates. The inherited contour/signature
checks at epsilon=0.01 and 0.05 are retained. The maximum scaled discriminant
residual in this run is 1.50e-9; EP numerical gaps are roundoff-sensitive.
This continuation uses a fixed eigenvalue selection disk and is not a global
branch tracker; no conclusion is drawn beyond epsilon=0.10.

The certified companion below encloses the algebraic seed and uses a uniform
contraction of the normalized EP equations. It supplies a remainder bound on
smaller disks without relying on this floating-point contour estimate.

Related: [F163 proof](../docs/proofs/PROOF_ROUTE_B_N6_UNFOLDING.md),
[initial atlas and end-profile study](ROUTE_B_N6_A2_LOCUS_ATLAS.md).


## The next coefficient from the full block

The next term can be computed at the crossing without fitting any finite-epsilon
EP position. Write

    q(epsilon)=q0+c epsilon^m+d epsilon^(2m)+O(epsilon^(3m)),

with m=1 for one/even and m=2 for odd. Here c and d are in the q convention;
no extra t=i q conversion is applied to the numerical coefficients below.
The existing F163 holomorphic branches justify this Taylor expansion locally.

Choose a basis U of the semisimple spectral plane and its spectral dual W,
WU=I, P=UW, Q=I-P. For this complex-symmetric pencil the implementation uses
W=(U^T U)^(-1)U^T, not the Hermitian adjoint. Let

    S=Q (lambda0 I-L0+P)^(-1) Q,
    L(q,epsilon)=D+q(C+epsilon V).

Start with the known path q=q0+c epsilon^m, omitting d. Its perturbation
E=sum E_n epsilon^n has E_1=q0 V, E_m+=c C, E_(m+1)+=c V;
coincident indices are added. Express the invariant subspace as U+X, WX=0,
and its effective matrix as lambda0 I+K. The graph equation gives

    F_n = E_n U + sum_(j=1)^(n-1) E_j X_(n-j),
    K_n = W F_n,
    X_n = S [Q F_n - sum_(j=1)^(n-1) X_j K_(n-j)].

Compute these terms through n=2m. Put H=K_m, J=K_(2m), A=W C U and

    B(X,Y)=4 tr(XY)-2 tr(X) tr(Y).

For m=2, reflection makes K_1=K_3=0 exactly in this fixed-plane gauge.
Inserting the omitted d adds dA to K_(2m). The coefficient of epsilon^(3m)
in the discriminant is consequently B(H,J)+d B(H,A), so

    d = -B(H,J)/B(H,A).

The denominator is the derivative of the leading discriminant polynomial at
its simple root c; F163 proves it nonzero. These are algebraic operations on
the original pencil and its spectral projector. The displayed decimals are
floating-point evaluations of that formula, not new exact decimal constants.

| Profile | Branch | d in q convention |
|---|---:|---:|
| One-sided | 0 | 18.586518485 - 3.486687067 i |
| One-sided | 1 | -2.588941305 - 4.488950741 i |
| Symmetric | 0 | 15.064990263 - 1.076414495 i |
| Symmetric | 1 | -0.354298811 - 0.612167558 i |
| Antisymmetric | 0 | 16.076044049 + 8.033093998 i |
| Antisymmetric | 1 | -20.827887182 + 7.346441782 i |

Branch labels refer to the stored leading slopes, not to a real-part sorting
at each epsilon. Using q0+c epsilon^m+d epsilon^(2m), the larger of the two
relative displacement errors becomes:

| epsilon | One-sided | Symmetric | Antisymmetric |
|---:|---:|---:|---:|
| 0.01 | 0.2013% | 0.0672% | 0.0000282% |
| 0.05 | 4.1582% | 1.7661% | 0.01755% |
| 0.10 | 13.7830% | 7.3878% | 0.27855% |

The expected relative remainders are now O(epsilon^2) and O(epsilon^4).
The improvement is not uniform over every individual row: the symmetric
branch 1 had a small leading-only error at epsilon=0.10 (0.464%), whereas
its corrected error is 4.076%. Asymptotic truncations need not improve at
every finite argument, and the larger-branch table must not conceal that.

[Coefficient producer](../simulations/route_b_n6_next_coefficient.py) and
[its output](../simulations/results/route_b_n6_next_coefficient.json) retain
both branches, observed next-coefficient quotients, scaled remainders and
source hashes. Reproduce with the same single-thread environment as above,
then `python simulations/route_b_n6_next_coefficient.py`. The graph recurrence
uses only matrices at epsilon=0; the comparison rows come from the separate
full-block continuation. No fitting is performed.

### A certified bound on the remaining tail

The [ball-arithmetic proof](../docs/proofs/PROOF_ROUTE_B_N6_REMAINDER_BOUND.md)
now certifies holomorphic EP branches for this crossing on explicit complex
disks. A normalized Schur-complement system, evaluated with outward complex
balls at 512 bits, contracts uniformly on the entire parameter disk. The odd
profile is made analytic in eta=epsilon^2 by an exact triangular similarity.
The full 90D complement remains in the equations.

Both branches are certified on common epsilon radii 2^-20 (one), 2^-18
(even), and 2^-10=1/1024 (odd). These are sufficient radii, much smaller than
the numerical scan, not measured points where the branches fail.

Write q=q0+eta z(eta), with |z(eta)-z_center|<=rho=2^-8 on |eta|<=R.
Cauchy's inequality applied to z-z_center bounds the entire remaining tail:

    |q-q0-c eta-d eta^2| / |q-q0|
       <= rho/(|z_center|-rho) * s^2/(1-s),  s=|eta|/R < 1.

Here c=z(0) and d=z'(0) are exact Taylor coefficients. R is the certified
branch-specific eta radius recorded by the producer. At half the common eta
radius the worst-branch relative bounds are 0.03630% (one), 0.03630% (even),
and 0.04582% (odd). For odd, this corresponds to epsilon=1/(1024 sqrt(2)).
Rounding the coefficients introduces additional rounding error.

Reproduce with `python simulations/route_b_n6_remainder_ball.py` in the same
single-thread environment; the [certificate](../simulations/results/route_b_n6_remainder_ball.json)
records failed and successful attempts and rejects a displaced-center control.
The program fails if any branch has no certified radius. The F163 exact
algebraic certificate is an imported premise, with source hashes checked.

The [complex-circle scout](../simulations/route_b_n6_complex_remainder.py)
also explores larger circles (up to |epsilon|=0.05 one/even and 0.10 odd).
Its sampled scaled-remainder maxima, branch closures and eigenvalue distances
are exploratory readings, not premises of the interval proof. Sharpening the
conservative certified region toward those circles remains open.


### A relation among the three profiles

The differences between the one-sided and symmetric second coefficients are
3.521528222 - 2.410272571 i and -2.234642494 - 3.876783183 i. They equal the
two leading antisymmetric coefficients as an unordered pair. The numerical
set residual is below 2.75e-13. This is an exact relation, proved in
[F163 section 4a](../docs/proofs/PROOF_ROUTE_B_N6_UNFOLDING.md#4a-the-next-coefficient-links-the-three-profiles),
for all its 266 loci. One/even branches must be matched by their common
leading coefficient before subtracting; the odd pair has no imposed ordering.

Reflection makes the difference of the second-order effective matrices the
same opposite-parity excursion that opens the odd profile. The common complex
bilinear symmetry then turns that response into the same two scalar roots.
The producer retains a nonsymmetric 2x2 counterexample: reflection-style
cancellation alone would not imply the identity.
