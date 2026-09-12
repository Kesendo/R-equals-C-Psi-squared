# PROOF: The N=7 Missing-Phase Relaxation Scale

**Status:** Tier 1 derived for the local N = 7 A-sector theorem stated below.
**Date:** 2026-09-11
**Authors:** Thomas Wicht, Codex (OpenAI)
**Exact producer:** [`missing_phase_relaxation_scale.py`](../../simulations/missing_phase_relaxation_scale.py)
**Certificate:** [`missing_phase_relaxation_scale.json`](../../simulations/results/missing_phase_relaxation_scale.json)
**Tests:** [`test_missing_phase_relaxation_scale.py`](../../simulations/tests/test_missing_phase_relaxation_scale.py)
**Owning experiment:** [The Motion and the Missing Phase](../../experiments/THE_MOTION_AND_THE_MISSING_PHASE.md)

## 1. Statement and fences

Fix `gamma > 0`. Put `r = 1 + epsilon` in the seven-site path defined in
Section 2. There is a radius `delta(gamma) > 0` such that, for
`0 < |epsilon| < delta(gamma)`, the slowest nonzero decay rate of the
49-dimensional A generator in the joint-popcount sector `(1,1)` is

```text
Delta_A(epsilon)
  = (gamma/2) epsilon^2 + (gamma/2) epsilon^3 + O_gamma(epsilon^4).
```

Equivalently, its reciprocal spectral scale satisfies

```text
tau_gap := 1/Delta_A  ~  2/(gamma epsilon^2).
```

The quantifier is **for fixed gamma>0, in a punctured neighbourhood of
epsilon=0 whose radius is not claimed uniform in gamma**. The
`O(epsilon^4)` remainder belongs to the analytic eigenvalue branch and hence
to its gap. `tau_gap` is only the reciprocal spectral scale. It is not a
proved relaxation time for `d_out`, `d_2`, leakage, or any other observable.
This is not the gap of the full `4^7` Liouvillian and not an all-N theorem.
No F-number is assigned.

Before this proof was written, the named stores were read as one object.
`docs/ANALYTICAL_FORMULAS.md` contained the absorption and blind-seat
ingredients but no epsilon-to-zero law; the absorption and blind-seat proofs
contained the two structural inputs but not their local composition; the
owning experiment contained the reduction, the historical `+/-0.1` readings,
and this question; neither Confirmations registry nor the hardware pages
contained a relevant flight; the glossary fixed the meanings of seat,
blindness, stationary span, and light content; the OpenArcs registry held the
question under `relaxation_scale_as_the_defect_vanishes`; and
`docs/CAUGHT_ERRORS.md` warned against substituting a short-time exact witness,
a sorted eigenvalue index, or a characteristic polynomial alone for the
object required here. The new content is therefore the local composition and
its exact perturbative calculation, not a re-derivation of an already owned
claim.

## 2. The N=7 A generator

Use the one-excitation site basis `|0>,...,|6>`, the watched centre
`c = 3`, and the real path Hamiltonian

```text
h_r[j,j+1] = h_r[j+1,j] = 2r  for j = 0,
                              2   for j = 1,...,5,
h_r[j,j] = 0.
```

Let

```text
P_c = |3><3|,                 z = I - 2 P_c,
L_A(X) = -i[h_r,X] + gamma(zXz-X).
```

This is the A block of the historical missing-phase reduction in the `(1,1)`
joint-popcount sector. In row-stack convention,

```text
vec_C(L_A)
  = -i(h_r tensor I - I tensor h_r^T)
    + gamma(z tensor z^T - I_49).
```

The factor `2r` is the physical hopping convention of the experiment: one
end bond moves, all other bonds remain at `2`. A symmetric two-end detuning is
a control, not this path.

## 3. Exact blind vector and 7 by 7 invariant reduction

For every `r`, define

```text
v_r = (1,0,-r,0,r,0,-r)^T / sqrt(1+3r^2).
```

Direct multiplication gives

```text
h_r v_r = 0,       P_c v_r = 0,       z v_r = v_r.
```

For every column vector `u`, therefore,

```text
L_A(u v_r^dagger)
  = (-i h_r - 2 gamma P_c)u v_r^dagger
  = K_r u v_r^dagger,

K_r := -i h_r - 2 gamma P_c.
```

The producer proves this identity exactly over `Q(i)` using the unnormalised
vector `(1,0,-r,0,r,0,-r)^T`; normalization adds no content and would only
introduce a square root. Thus a seven-dimensional invariant family produces
the candidate eigenvalue branch. The 49-dimensional generator is constructed
separately and is not defined by this reduction.

## 4. Two-sided row-stack embedding and doubled clusters

The adjoint-side family is invariant as well. If

```text
K_r u_- = lambda_- u_-,
K_r u_+ = lambda_+ u_+,
lambda_+ = conjugate(lambda_-),
```

then the following are exact A eigenoperators:

```text
u_- v_r^dagger  and  v_r u_+^dagger   at lambda_-,
u_+ v_r^dagger  and  v_r u_-^dagger   at lambda_+.
```

For a nonzero root the two operators on each line are independent. Hence each
of the conjugate A clusters has rank two and cannot split internally. In the
repository's row-major convention,

```text
vec_C(u_- v_r^dagger) = u_- tensor conjugate(v_r),
vec_C(v_r u_+^dagger) = v_r tensor conjugate(u_+).
```

The independently assembled 49 by 49 action verifies all four eigenoperators.
This two-sided construction is load-bearing: a column-stack substitution such
as `v_r tensor u_-` fails against that independently assembled matrix.

## 5. Characteristic polynomial and branch expansion

Exact elimination in the seven-dimensional reduction gives

```text
det(lambda I-K_r) = lambda Q_r(lambda),

Q_r(lambda) = lambda^6 + 2 gamma lambda^5
  + (4r^2+20) lambda^4
  + (8 gamma r^2+24 gamma) lambda^3
  + (64r^2+96) lambda^2
  + (64 gamma r^2+64 gamma) lambda
  + 192r^2+64.
```

At the uniform point `r=1`,

```text
Q_1(lambda)
  = (lambda^2+8)
    (lambda^4+2 gamma lambda^3+16 lambda^2+16 gamma lambda+32).
```

The roots `lambda_0 = +/-2 sqrt(2)i` are simple for every `gamma>0`.
The implicit-function theorem therefore supplies analytic branches rather
than a Puiseux or absolute-value cusp. Substituting `r=1+epsilon` and solving
the branch from `-2 sqrt(2)i` order by order gives

```text
lambda_-(epsilon) = -2 sqrt(2)i
  - i epsilon/sqrt(2)
  - (gamma/2 + 3 sqrt(2)i/16) epsilon^2
  - (gamma/2 - 19 sqrt(2)i/64) epsilon^3
  + O_gamma(epsilon^4).
```

The conjugate branch has the same real part. Consequently

```text
-Re lambda_-(epsilon)
  = (gamma/2) epsilon^2 + (gamma/2) epsilon^3
    + O_gamma(epsilon^4).
```

The cubic term matters: the positive and negative defects need not have equal
finite-epsilon rates even though both share the same quadratic leading law.

## 6. The complete epsilon=0 peripheral space

At `epsilon=0`, reflection splits the one-excitation Hilbert space into the
three-dimensional blind space `D` and its four-dimensional complement `E`.
The restriction of `h_1` to `D` has the three simple energies
`-2 sqrt(2), 0, +2 sqrt(2)`. Every operator in `End(D)` pays zero dissipative
cost, and the scalar `I_E` supplies one further zero-cost invariant direction.
These ten independent operators exhaust the A peripheral space:

| A eigenvalue | exact multiplicity |
|---|---:|
| `-4 sqrt(2)i` | 1 |
| `-2 sqrt(2)i` | 2 |
| `0` | 4 |
| `+2 sqrt(2)i` | 2 |
| `+4 sqrt(2)i` | 1 |

Thus the peripheral dimension is 10 while the kernel dimension is 4. The
producer obtains `D` as the nullspace of the transpose of the centre-seat
Krylov matrix, constructs the ten operators, verifies their frequencies by
the exact 49-dimensional action, and separately computes the exact kernel
nullity. The peripheral census is not inferred from the kernel.

## 7. Kato/Feshbach splitting and why the selected pair is locally minimal

Write the exact generator pencil as `L(epsilon)=L_0+epsilon L_1`. For each
unperturbed peripheral eigenvalue `lambda_0`, let the columns of `R` span its
right Riesz space and let `W^dagger R=I`. With

```text
P = R W^dagger,
Q = I-P,
G = (L_0-lambda_0 I+P)^(-1)-P,
F_1 = W^dagger L_1 R,
F_2 = -W^dagger L_1 G L_1 R,
```

the producer verifies exactly `P^2=P`, `[P,L_0]=0`, the reduced-resolvent
identities, and the nonvanishing of every denominator for `gamma>0`. In the
ordered Riesz bases the first-order operators are

```text
F_(-2sqrt2 i)^(1) = -i/sqrt(2) I_2,
F_(+2sqrt2 i)^(1) = +i/sqrt(2) I_2,
F_(-4sqrt2 i)^(1) = -i sqrt(2) I_1,
F_(+4sqrt2 i)^(1) = +i sqrt(2) I_1,
F_0^(1)            = 0_4.
```

They are purely reactive: no peripheral direction acquires a decay rate at
first order. At second order,

```text
F_(-2sqrt2 i)^(2) = (-gamma/2 - 3 sqrt(2)i/16) I_2,
F_(+2sqrt2 i)^(2) = (-gamma/2 + 3 sqrt(2)i/16) I_2,
F_(-4sqrt2 i)^(2) = (-gamma   - 3 sqrt(2)i/8)  I_1,
F_(+4sqrt2 i)^(2) = (-gamma   + 3 sqrt(2)i/8)  I_1,

F_0^(2) = [ -gamma      0       0       gamma/2 ]
          [    0        0       0          0    ]
          [    0        0    -gamma      gamma/2 ]
          [ gamma/2     0     gamma/2   -gamma/2 ].
```

The zero-frequency basis is
`(|d_-><d_-|, |d_0><d_0|, |d_+><d_+|, I_E/2)`, and the last matrix has
spectrum `{0,0,-gamma,-3gamma/2}`. Therefore the four directions at
`+/-2 sqrt(2)i` acquire rate `(gamma/2)epsilon^2`, the two directions at
`+/-4 sqrt(2)i` acquire `gamma epsilon^2`, and the two departing zero-cluster
directions acquire `gamma epsilon^2` and `(3gamma/2)epsilon^2`. Two zero
directions remain stationary.

Every nonperipheral A eigenvalue is already strictly left of the imaginary
axis at `epsilon=0`: Section 6 exhausts the equality case of the dissipative
Rayleigh quotient. A finite spectrum then gives a strictly positive
complementary separation, and continuity retains it locally. The exact
second-order ordering is strict, so the `+/-2 sqrt(2)i` pair is the slowest
nonstationary A rate for all sufficiently small nonzero `|epsilon|` at fixed
`gamma>0`. This proves local minimality; following one polynomial root alone
would not.

## 8. Absorption/light identity

For centre-only Z-dephasing the Hermitian part of `L_A` charges exactly those
matrix cells whose bra and ket differ at site 3. Let `Delta_c` be their
orthogonal coordinate projector. For any right A eigenoperator `M`, the
site-resolved Absorption Theorem gives

```text
-Re(lambda) = 2 gamma <M,Delta_c M>/<M,M>.
```

The continued eigenvalue has two exact embeddings, so the degeneration-safe
cluster reading uses the orthogonal rank-two projectors separately:

```text
w_c^R = Tr(Pi_R Delta_c)/2,
w_c^L = Tr(Pi_L Delta_c)/2.
```

Here `Pi_R` projects onto the right invariant cluster of `L_A`, while `Pi_L`
projects onto the independently constructed right invariant cluster of
`L_A^dagger`. Applying the Rayleigh identity to an orthonormal basis of each
space yields

```text
Delta_A = 2 gamma w_c^R = 2 gamma w_c^L,

w_c^R = w_c^L
      = epsilon^2/4 + epsilon^3/4 + O_gamma(epsilon^4).
```

This is the physical mechanism: breaking one end bond gives the formerly
centre-blind oscillatory mode centre-charged operator weight only at
quadratic order in `epsilon`. The light is operator-cell weight, not a
probability that an excitation occupies the centre. The raw biorthogonal
Riesz projector is not either orthogonal light, and omitting the division by
rank doubles the reading.

## 9. Punctured kernel and finite complementary separation

For `r^2 != 1` near `r=1`, the left and right principal characteristic
polynomials share only the zero root. Hence `v_r` spans the complete blind
space. Because every bond remains nonzero, the blind-seat commutant theorem
applies on the zero-free path and the complement restriction is simple. The
A kernel is therefore exactly two-dimensional, spanned by

```text
v_r v_r^dagger,
I-v_r v_r^dagger.
```

Kernel completeness is not by itself peripheral completeness. The exact
certificate also constructs the Hamiltonian commutator `A=-i ad(h_r)` and the
charged-cell constraint `D`, then iterates

```text
V_0 = ker D,
V_(j+1) = {x in V_j : A x lies in V_j}
```

until the largest A-invariant subspace inside `ker D` stabilizes. At the
punctured rational anchors it has dimension two, equals the stationary span,
and its restricted characteristic polynomial is `lambda^2`. Thus no nonzero
imaginary-axis A direction survives. This exact invariant-subspace route is
separate from the direct 49 by 49 kernel-nullity gate.

At the limiting point, Section 6 already places every remaining eigenvalue a
positive distance from the imaginary axis. This finite complementary
separation, the exact Kato ordering, and continuity together supply the local
radius used in Section 1. No uniform lower bound on that radius as gamma
varies is asserted.

## 10. B carrier boundary and why its gap is not 2 gamma minus the A gap

The historical second carrier is built directly as

```text
L_B(X) = -i[h_r,X] - gamma(zXz+X).
```

Taking the Hilbert-Schmidt adjoint of A gives the operator identity

```text
L_B = -2 gamma I - L_A^dagger.
```

It follows mode by mode that

```text
lambda_B,k = -2 gamma - conjugate(lambda_A,k),
Delta_B,k  =  2 gamma - Delta_A,k.
```

The minimum B rate is consequently

```text
Delta_B = 2 gamma - max_k Delta_A,k,
```

not `2 gamma-Delta_A`, because `Delta_A` denotes the minimum positive A rate.

It remains to exclude a B mode on the imaginary axis at the limiting point.
Write the uniform Hamiltonian centre-first as
`h=[[0,u^dagger],[u,h_o]]`. Zero B dissipative cost restricts an eigenoperator
to the centre/outer off-diagonal cells. The eigenoperator equations make its
two outer vectors proportional to `u` and then require `h_o u=0`. For this
seven-site path the exact product is

```text
h_o u = (0,4,0,0,4,0)^T != 0.
```

Thus B has no peripheral direction at `epsilon=0`; its finite spectrum has a
positive gap that remains positive locally. Since `Delta_A(epsilon)` tends to
zero, A supplies the gap of the historical invariant A-plus-B carrier in the
punctured neighbourhood. The two global-flip copies are identical. Nothing
here classifies the other many-body sectors.

## 11. Numerical certificate and precision contract

The exact argument above owns the theorem. The numerical certificate is an
independent reconstruction and finite-anchor reading, not the source of the
series. It sweeps exact tokens

```text
gamma in {3/100, 3/10, 3},
epsilon in {+/-1/10, +/-1/2^k for k=3,...,13},
```

for 72 rows. Positive and negative branches are continued independently from
the exact rank-two Riesz seed, through `+/-1/10`, `+/-1/8`, and then down the
dyadic tail. Selection is by principal-angle overlap of two-dimensional right
and left invariant spaces, never by sorted eigenvalue index.

Each attempt independently rebuilds `K`, `L_A`, `L_A^dagger`, `L_B`, and both
global-flip copies at precision pairs `(53,106)`, `(106,212)`, or `(212,424)`
bits. A row is `TRUSTED` only if all of the following survive at both
precisions:

- rank two and unambiguous continuation with principal angles below `pi/2`;
- two-sided residual/backward-error bounds and A-versus-K eigenvalue agreement;
- ordinary complex separation and the true Sylvester singular-value separation;
- the half-radius subspace bound and p/2p orthogonal-projector agreement;
- separately positive, stable right and left light intervals and both
  absorption identities;
- exact stationary projection followed by strict separation of the next
  distinct A rate;
- fresh B/global-flip operator identities and full-spectrum numerical reads;
- authoritative construction and measurement provenance.

Any unresolved row makes the top-level verdict `UNRESOLVED`; an empty sweep
cannot pass. The JSON stores compact diagnostics and source/dependency
provenance, not eigenvectors or full spectra. Its scaling exponents and the
finite reads of the quadratic and cubic coefficients are labelled regression
readings: the exact polynomial and Kato calculation, not finite convergence,
are the theorem gates.

Reproduce the certificate and tests from the repository root with

```powershell
python -m pytest simulations/tests/test_missing_phase_relaxation_scale.py -q
python simulations/missing_phase_relaxation_scale.py
```

## 12. Controls and mutations

The gate set is designed to fail in the directions that would change the
claim:

- commutator-sign, operator-side, watched-seat, missing-identity, and
  row/column-stack mutations break the exact reduction or 49-dimensional
  embedding;
- dropping or changing any displayed branch coefficient breaks the exact
  characteristic-polynomial substitution through cubic order;
- replacing the ten-dimensional peripheral census by the four-dimensional
  kernel, changing a Kato operator, or deleting the `I_E/2` direction changes
  exact multiplicities or characteristic polynomials;
- a boundary with `u != 0` but `h_o u=0` creates a B peripheral direction,
  while missing-adjoint, missing-right-action, and price-sign mutations break
  the A/B operator or full-spectrum map;
- equal detuning of both end bonds preserves the ten peripheral directions,
  moving the watched seat to site 2 opens a finite gap, and changing the path
  to `r=1+epsilon^2` changes the leading gap to fourth order;
- the exact gauge relation is
  `Delta_A(epsilon)=Delta_A(-2-epsilon)`, not an invented evenness under
  `epsilon -> -epsilon`;
- promoting one low-precision matrix is rejected as an independent rebuild;
  a fixed absolute kernel cutoff fails on the dyadic tail; a nonnormal toy
  separates ordinary eigenvalue distance from Sylvester separation; and two
  rotated rank-two projectors with equal light expose why light agreement
  cannot replace subspace agreement;
- the new A construction agrees with the preserved historical `+/-0.1`
  generator without modifying that producer.

These controls separate the physical one-end symmetry breaking from a
reflection-preserving path, the local theorem from a branch-selection
artifact, and an error-bounded certificate from a fixed-tolerance fit.

## 13. What is not proved

This proof does not establish:

- the gap of the full `4^7` Liouvillian;
- a relaxation time of `d_out`, `d_2`, leakage, or any other named observable;
- that the relevant preparation or readout has nonzero overlap with the slow
  cluster;
- the same expansion for every odd N, or for any N other than 7;
- a punctured-neighbourhood radius uniform in gamma;
- an F-number or a new hardware-confirmed result;
- the false formula `Delta_B=2 gamma-Delta_A` for the B gap;
- a statement for symmetric two-end detuning or a different watched seat.

The open continuations are therefore sharply smaller than the retired local
question: derive or refute the law for all odd N, control the neighbourhood as
gamma varies, and determine which physical observables actually couple to the
rank-two slow cluster.
