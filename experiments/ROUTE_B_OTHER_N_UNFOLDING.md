# The end-bond crossing at N=4 and N=5

The local coefficient relation of F163 also appears at the smaller chains:
two N=4 crossings at one positive-real q, and all 58 stored N=5 A2 q loci.
The new instance checks are numerical. The relation itself has a
dimension-independent conditional proof; this survey does not supply an
exact nonzero/rank certificate for every N=5 locus or a theorem at every N.

## Same object, different block sizes

All readings use the full (ket weight 1, bra weight 2) coherence block:
dimension N*binomial(N,2), namely 24, 50, and 90 for N=4,5,6.
The model is the open XY chain, Delta=0, uniform gamma=1,
q=J/gamma=qCSharp and hopping amplitude 2q. The end profiles are
one=(1+epsilon,1), even=(1+epsilon/2,1+epsilon/2), and
odd=(1+epsilon/2,1-epsilon/2), with all interior bonds fixed at 1.

The [producer](../simulations/route_b_other_n_unfolding.py) reuses the
framework's weight-block primitive and the existing graph recurrence.
An independent full Hilbert-space weighted-XY construction, followed by
leg restriction, agrees entrywise with each tested end direction.

N=4 uses the known algebraic coordinates

    q0=sqrt((-1+sqrt(13))/6)=0.6589829632931832,
    lambda0=-4 +/- 2 i q0.

Both lambda signs are included in the coefficient survey; only the positive
imaginary sign is followed in epsilon. This is not a census of all N=4 loci.
The seed is inherited from the F89 path-3 octic study, including
[the existing full-block Jordan probe](../simulations/f89_jordan_definitive.py).

N=5 consumes every entry of
[the stored inventory](../simulations/results/route_b_a2_n5.json): 29 A2(w)
roots, each with its two q square roots. Its convention is
qPhysicalCSharp=qUnitHop/2, while lambdaSeed is already physical lambda.
The representative for finite-epsilon continuation is
`N5-O-A2-W-006-Q-plus`, q0=1.1292509708747671,
lambda0=-4.791960365179641. It lies on the positive-real q axis.

| Selected crossing | Reflection on the plane | Nearest third eigenvalue distance | Reduced resolvent norm, 2-norm |
|---|---:|---:|---:|
| N=4, positive Im(lambda) | +I | 0.0235478589 | 113.113797 |
| N=5, O-W-006-Q-plus | -I | 1.1614749462 | 4.058028 |

These distances are in physical-lambda units at gamma=1. Comparing different
selected crossings is not a monotone scaling law in N.

## What survives across the tested loci

At all 60 surveyed crossings the full-block two-plane is numerically scalar
under L0, reflection acts as +I or -I, and the odd first compression vanishes
to numerical precision. The largest scalar-plane residual is below 2e-14;
the largest odd compression norm is below 6.5e-13. All computed alpha and
Omega values are nonzero; among the N=5 rows the smallest |alpha| is 0.15876
and the smallest |Omega| across the three profiles is 9.9150e-5.
These are floating-point readings, not algebraic nonzero certificates.

The one/even branches share their first coefficients c_j. For the next
coefficients d_j, their differences reproduce the leading odd coefficients
a_j as an unordered pair:

    {d_one,j-d_even,j : j=0,1} = {a_0,a_1}.

One/even branches are matched by their common c_j before subtraction.
The maximum numerical mismatch, divided by max(1,max_j |a_j|), is below
1.63e-13 across these 60 rows. The matching never presumes an ordering of
the odd roots.

The proof in [F163 section 4a](../docs/proofs/PROOF_ROUTE_B_N6_UNFOLDING.md#4a-the-next-coefficient-links-the-three-profiles)
uses an isolated semisimple rank-two spectral plane in one reflection parity,
a common complex bilinear symmetry, and nonzero alpha/Omega for the even
leading polynomial. None of these arguments depends on N=6 or dimension90.
Odd Omega must additionally be nonzero to conclude two distinct quadratic
odd EP branches. Reflection alone guarantees evenness, not a nonzero second
order. Nor does this prove the hypotheses at every N or every crossing.

## The coefficients need not be small

For the selected representatives, write q=q0+c epsilon+... for one/even and
q=q0+a epsilon^2+... for odd. The coefficients below are real to numerical
precision; their small imaginary parts are retained in the JSON.

| Selected crossing | c, one/even | a, odd |
|---|---|---|
| N=4 | -0.520880800, +0.100646300 | +37.494111388, -30.148714777 |
| N=5 | -0.766036465, -0.393734713 | +3.729208560, +2.123889178 |

Thus the odd response is suppressed in order, not necessarily in coefficient.
At the selected N=5 point, the one/even leading movements both go toward
smaller q, whereas the odd movements both go toward larger q.

For N=4 the next odd epsilon^4 coefficients are approximately -207461 and
+323853. Their size is a warning against assigning the N=6 approximation's
useful epsilon range to a different crossing without checking it.

## Which complementary modes supply the response?

For a simple complementary mode use its complex-symmetric spectral residue

    P_j=v_j v_j^T/(v_j^T v_j),   S_j=P_j/(lambda0-lambda_j).

The quadratic effective response is B2=W(q0 V_odd)S(q0 V_odd)U.
Recompute it with S_j alone and with S-S_j. This is a decomposition of the
response formula, not a physically modified Hamiltonian.

At N=4 the closest complementary mode is
lambda_j=-4+1.2944180677094614 i. It is reflection odd, opposite the even
seed plane. Its projector norm is about 2.6818. Its B2 contribution has
Frobenius norm 403.549, compared with 400.883 for the full B2 in the same
orthonormal seed basis. The coefficient readings are:

| N=4 resolvent contribution | Odd coefficients |
|---|---|
| Full complement | +37.49411139, -30.14871478 |
| Nearest mode alone | +38.08641535, -30.13578375 |
| Nearest mode removed | -0.59230397, -0.01293103 |

The nearby allowed mode therefore dominates the large response in this
example. The conclusion comes from recomputing its contribution, not merely
from inspecting the eigenvalue gap.

At N=5 the closest modes are the conjugate reflection-even pair
-3.77789159523160 +/- 0.56629372288187 i. Both are retained together in this
control. Their contribution alone gives +6.24903949 and -0.76275776;
removing the pair gives +2.88664694 and -2.51983093, compared with the full
+3.72920856 and +2.12388918. Several contributions combine here. A smallest
gap by itself does not explain the signs or magnitudes.

## Finite-epsilon checks and their reach

Both branches of every profile were continued with the full-block Schur
discriminant at each representative. N=4 uses epsilon=0.0001,0.0003,0.001
for one/even and 0.001,0.002,0.003 for odd; N=5 uses 0.001,0.003,0.01
for each profile. Both branches remain separated from each other and from
the complementary spectrum on these grids. An independent full-block SVD
check at all 36 points finds a single numerical null direction, with the
second-smallest singular value well above roundoff: numerical EP2 evidence.

At epsilon=0.01 the selected N=5 q locations are:

| Profile | Two continued q locations |
|---|---|
| One | 1.122166502879, 1.125590939225 |
| Even | 1.121779246490, 1.125365138337 |
| Odd | 1.129623370701, 1.129463268225 |

Their imaginary parts are below 8e-15 on the grid. These are numerical
real-q EP locations in the stated model, not a certified interval of real-q
branches. No observation strength or experimental visibility is inferred.

For odd perturbations the discriminant is divided by epsilon^4 in the
solver, so its reported scaled residual is not a direct error in q. At the
N=4 odd points it reaches about 6.84e-6 while the leading scaled derivative
is about sqrt(|Omega|)=9073.6. Solving farther from zero reduces roundoff
amplification; a raw scaled residual alone should not label a locus as failed.

Reproduce with:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python simulations/route_b_other_n_unfolding.py
```

The [N=4 self-fold companion](../docs/proofs/PROOF_ROUTE_B_N4_SELF_FOLD.md)
proves the N=4 leading coefficients exactly and establishes individual
real-q EP branches with Re λ = −4 for sufficiently small real ε, using
their distinct real slopes and the internal bra fold. It supplies no
explicit radius and does not certify the finite grids in this note.
The N=5 statements here remain numerical.

The [JSON](../simulations/results/route_b_other_n_unfolding.json) retains the
60 coefficient rows, 36 continuation rows, projector-residue controls,
imaginary components and source/dependency hashes. Independent reviews
recomputed the coefficients and representative Jordan signatures. This work
does not transfer the N=6 ball-radius certificate to N=4/5 and does not survey
N=7. At odd N the staggered gauge preserves reflection parity, so N=6's
even-to-odd-sector transport is not an unchanged shortcut for N=5.
