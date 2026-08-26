# The V-Effect Mechanism: What the Exact-Palindrome Construction Shows

**Status:** No mechanism established. The palindrome is a spectral involution
and forbids no oscillation; the coupled construction is not palindromic. What
stands is a set of frequency counts.
**Result date:** March 27, 2026 (the change history lives in git).
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:** [Palindrome Proof](PROOF_PALINDROME_NEURAL.md),
[Hierarchy of Incompleteness](../../HIERARCHY_OF_INCOMPLETENESS.md)

**What the repo already held, and what it returned. It held the answer, in one
sentence, since April.** `docs/KMS_DETAILED_BALANCE.md:37-39`: *"QDB relates L to its
adjoint L† (real eigenvalues); Π relates L to its negative (palindromic complex
pairs)"*, with the table at `:200-207` and `:262-266` (*"It allows complex eigenvalues.
It pairs decay rates rather than forcing them to be real."*). That is this file's error
and its refutation, written before the error was made, together with the symmetry that
WOULD force reality. `docs/proofs/PROOF_K_PARTNERSHIP.md:164` goes further and names the
neighbouring version of it a standing hazard: *"The recurrent terminology error of
calling KHK = −H 'time reversal'"*, a linear chiral symmetry read as the antiunitary one
that does force reality; `compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs:20-21`
repeats the warning in code. `docs/proofs/MIRROR_SYMMETRY_PROOF.md:453-454` writes the
partner as −d + iω → 2Σγ − d, −ω, carrying the frequency explicitly. For the residual:
`docs/ANALYTICAL_FORMULAS.md:5922` (F137) states the general law this page needed,
*"The centre is an identity and carries no evidence"*, the point being that exactly ONE
candidate centre exists, so that a large distance at it makes "broken" a MEASUREMENT
rather than a failed search. Fitting the centre per seat reverses that. `experiments/`:
`NEURAL_GAMMA_CAVITY.md` convicts a sibling result for reading a matcher's threshold as
a spectral fact. `fw.Confirmations`: nothing neural. `docs/CAUGHT_ERRORS.md`: carries the
normalisation and matcher-threshold shapes and the general "a gate built from its own
measurement cannot fail" (`:131`), but no entry for reading a linear spectral reflection
as reality. The OpenArcs registry: no arc on the neural V-Effect.

None of this was consulted when the file was written. The quantum side of this
repository had the distinction, stated, tested and flagged as recurrent, and the neural
side re-made the error next door.

---

## What this document is about

Coupling two neural networks through a shared mediator changes the number of
oscillation frequencies in the linearised dynamics, non-monotonically in the
coupling strength. That much is measured. No mechanism accounts for it. The
natural one, that each network is individually silent because an exact
palindrome forbids oscillation and the coupling releases what the symmetry held
down, fails twice over: the palindrome forbids no oscillation, and the coupled
construction never satisfied it.

---

## What the palindrome actually gives

If Q·J·Q + 2sI = −J holds with **one scalar** s and Q a permutation
involution, then for every eigenvalue μ of J the value −μ − 2s is also an
eigenvalue, with the same multiplicity: the map v ↦ Qv is a
multiplicity-preserving involution of the spectrum.

That is the whole content. The map μ ↦ −μ − 2s sends the complex plane to
itself and fixes the vertical line Re μ = −s. **It forbids nothing off the
real axis, and it bounds nothing.** A conjugate pair −s ± iω is its own
palindromic partner.

### The reality claim is false

Take the smallest admissible network: one excitatory seat (τ_E = 5) and one
inhibitory seat (τ_I = 10), Q the swap, so s = 0.15. The magnitude condition
W[Q(i),Q(j)] = −(τ_{Q(i)}/τ_i)·W[i,j] reduces here to J[0,1] = −J[1,0], so

```
J = [[-0.2,  -b],
     [   b, -0.1]]      b > 0
```

satisfies the palindrome **exactly** when built from the reciprocals 1/τ:
‖QJQ + J + 2sI‖ = 0.0, not small, at every b. The signs are Dale-legal as
written, the excitatory column carrying +b and the inhibitory column −b. Its
eigenvalues solve λ² + 0.3λ + (0.02 + b²) = 0, discriminant 0.01 − 4b²:

| b | eigenvalues | residual |
|---|---|---|
| 0.02 | −0.1958, −0.1042 | 0.0 |
| 0.05 | −0.15, −0.15 (defective) | 0.0 |
| **0.10** | **−0.15 ± 0.0866 i** | **0.0** |
| 0.30 | −0.15 ± 0.2958 i | 0.0 |

The pair sums to −0.3 = −2s in every row. The palindrome holds *while* the
spectrum is complex.

Nor does the palindrome give stability. It forces Re μ + Re μ′ = −2s, so an
exactly palindromic network is stable only if its whole spectrum sits inside
the fixed-width strip −2s ≤ Re μ ≤ 0. The spectrum outgrows a fixed strip as
the coupling α rises. Over 200 seeds of this file's own generator
([veffect_exact.py](../../../simulations/neural/veffect_exact.py),
`build_exact_palindromic_network`, N = 10, τ_E = 5, τ_I = 10, density 0.3),
counting only draws whose scalar-s residual is exactly 0.0:

| α | exactly palindromic | complex spectrum | unstable |
|---|---|---|---|
| 0.5 | 200/200 | 24 | 0 |
| 1.5 | 200/200 | 110 | 1 |
| 3.0 | 200/200 | 149 | 15 |
| 5.0 | 200/200 | 159 | 24 |
| 10.0 | 200/200 | 167 | 45 |

Every draw is exactly palindromic and most of them oscillate. A sweep over
N = 10, 20 and 30 at seed 42 does give real spectra throughout, which is a fact
about seed 42: at the same setting, 24 of its 200 siblings do not. One draw is
not a witness for a theorem.

---

## The residual the scripts measure is not the residual the statement needs

The neural residual routines in the repo
([veffect_exact.py](../../../simulations/neural/veffect_exact.py) lines 90-106,
[validation_checks.py](../../../simulations/neural/validation_checks.py) lines 76-84,
and the siblings) compute

```
S_diag = -(diag(QJQ) + diag(J)) / 2          # S fitted PER SEAT
R_off  = offdiag(QJQ + J + 2*diag(S_diag))   # diagonal then DISCARDED
```

The statement above needs 2S = 2sI, one scalar. A per-seat S absorbs any
mismatch in the self-decay rates into the fit, and dropping the diagonal
discards what is left. So residual = 0 as measured does **not** imply the
hypothesis, and where the seats do not all pair it demonstrably does not.

`celegans_trichotomy.py` keeps the diagonal in its norm, but fits S the same
way, so the diagonal it keeps is zero by construction and the value is the same;
its `palindrome_residual_norm` docstring says so ("S is chosen to zero the
diagonal"). The error is not in the code; it is in reading that number as the
algebraic condition.

---

## The coupled system was never palindromic

The construction couples two N-neuron networks through one mediator, so
N_c = 2N + 1 is **odd**. The Q the script builds pairs E-seats with I-seats
and leaves the mediator as a **fixed point**, at every coupling strength
including zero. By Step 3 of the [Palindrome Proof](PROOF_PALINDROME_NEURAL.md),
a fixed seat requires 2/τ_M = 1/τ_E + 1/τ_I; here that is 0.4 against 0.3,
so condition (a) fails at the mediator unconditionally.

Measured on the N = 10 coupled system at coupling c = 0:

| quantity | value |
|---|---|
| Q fixed points | seat 20, the mediator |
| fitted per-seat S | takes two values, 0.15 and 0.20 |
| residual as the script measures it | **0.0** |
| residual with the scalar s = 0.15 the statement needs | **0.132** |
| spectrum against μ ↦ −μ − 2s, matched to 1e-8 | 10 matched pairs, **1 unpaired** |

The unpaired eigenvalue is the mediator's own −0.2, whose partner would have
to be −0.1 and is not there. So the row the coupling table prints as
0.00 / 0 / 0 / 0.000, the row that carries the whole "two silent palindromes
before we couple them" reading, is not a palindromic system. Nothing was
holding the oscillation down, because there was no exact symmetry to break.

This also disposes of the "mediator cannot serve two masters" argument. That
argument requires Q_A·J·Q_A + J + 2S = 0 and Q_B·J·Q_B + J + 2S = 0
separately, with Q_A the A-swap extended by the identity on B. The B-block of
the first expression is 2(J_B + S), which is not zero at any coupling, with
or without a mediator. The code forms neither equation; it uses
Q = Q_A ⊕ Q_B ⊕ id_M. The frustration story describes no computed object.

---

## What the coupling sweep does measure

The frequency counts are real measurements of the linearised spectrum, and
the non-monotonicity in the coupling is real. What is not established is that
either is a symmetry effect.

**N=10 (5E+5I per network), correlation space:**

| Coupling | Frequencies (activity) | Frequencies (correlation) | Residual |
|----------|----------------------|--------------------------|----------|
| 0.00 | 0 | 0 | 0.000 |
| 0.01 | 2 | 6 | 0.005 |
| 0.05 | 3 | 12 | 0.025 |
| 0.10 | 3 | 12 | 0.050 |
| 0.30 | 1 | 2 | 0.150 |
| 1.00 | 0 | 0 | 0.474 |

**N=20 (10E+10I per network), correlation space:**

| Coupling | Frequencies (activity) | Frequencies (correlation) | Residual |
|----------|----------------------|--------------------------|----------|
| 0.00 | 0 | 0 | 0.000 |
| 0.01 | 6 | 48 | 0.003 |
| 0.05 | 7 | 62 | 0.017 |
| 0.10 | 6 | 47 | 0.035 |
| 1.00 | 5 | 31 | 0.335 |

**The residual column is not independent evidence: its numerator is exactly
linear in the coupling.** The printed column bends (at N=10 the ratio to the
coupling runs 0.503123, 0.503047, 0.502808, 0.500285, 0.495350, 0.474022; at
N=20 it runs 0.345163, 0.345139, 0.345062, 0.344243, 0.342623, 0.335322), and
the bend is entirely in the denominator. The reported quantity is ‖R_off‖/‖J‖,
and instrumenting
[veffect_exact.py](../../../simulations/neural/veffect_exact.py) prints the two
parts separately: the numerator is 0.38078865529·c to every digit across the
whole range at N=10 and 0.36055512755·c at N=20, while ‖J(c)‖ grows from 0.75685
to 0.80331 and from 1.04459 to 1.07525 respectively. So the residual
measures the coupling that was dialled in, divided by a norm that the same
coupling inflates; it is not a measure of how much structure the coupling broke.
Two readings this rules out: there is no saturation of the symmetry breaking
(the numerator never bends), and the two chain lengths do not differ in kind
(N=20 shows the same bend, 0.345 to 0.335). This is the same signature the
sibling result carries, where the palindrome residual collapses to a multiple of
‖W_eff‖/‖J‖ whenever no Q-partner pair of edges is present
([Algebraic Palindrome Neural](../ALGEBRAIC_PALINDROME_NEURAL.md)).

The script prints three further rows: at N=10, c = 0.50 gives residual 0.248
with 1 activity and 2 correlation frequencies; at N=20, c = 0.30 gives 0.103
with 5 and 34, and c = 0.50 gives 0.171 with 5 and 33.

### The coupling window is an N=10 statement

At N=10 the count rises to 12 at c = 0.05-0.10 and falls to 0 by c = 1.00. At
N=20 it rises to 62 and then **settles**, not collapses: 34 frequencies at
c = 0.30, 33 at c = 0.50 and 31 at the strongest coupling tested. "Strong
coupling destroys the structure" is a statement about N = 10.

The mediator coupling also breaks Dale's Law, the stated source of the sign
antisymmetry: the script writes W[2N, offset] = +c regardless of the source
seat's type, and for seed 42 the seat at offset 0 is **inhibitory** and sends
a positive weight. Whatever the sweep varies, it is not a Dale-legal
perturbation of a palindromic network.

---

## Correlation space: what the counts are counting

The correlation Liouvillian L_C = J⊗I + I⊗Jᵀ has eigenvalues λ_i + λ_j,
which is correct. The ceiling built on it was not.

With K distinct activity frequencies, each from a conjugate pair ±iω, the
attainable correlation frequencies are the sums (K(K+1)/2), the differences
(K(K−1)/2) and the complex-plus-real singles (K): at most **K(K+1)**, not
K(K+1)/2.

| N | K_activity | K_correlation | claimed ceiling K(K+1)/2 | true ceiling K(K+1) |
|---|-----------|--------------|--------------------------|---------------------|
| 10 | 3 | 12 | 6, violated | 12, attained exactly |
| 20 | 7 | 62 | 28, violated | 56, apparently exceeded |

At N=10 the ceiling is attained exactly: three genuine frequencies, twelve
symbolic combinations, twelve bins, nothing merged.

**The N=20 row does not exceed anything; the two counts are taken at
incompatible thresholds.** `J_c` there has **ten** conjugate pairs, not seven.
Three carry ω = 2.49e-9, 4.37e-9 and 3.67e-7, all below `count_freqs`'s
`tol = 1e-6`, so K_activity drops them and reports 7; `count_corr_freqs`
applies the same cut to the SUM, where ω_a ± 3.67e-7 survives. Exactly six of
the 62 bins are one of the seven large ω shifted by that 3.67e-7. Counting the
pairs as they are, K = 10, the ceiling is 110 and 62 sits well inside it.

The three sub-threshold ω are numerical, not physical: J_c is real, so its
transpose must have the same spectrum, and it returns a different sub-threshold
set (1.44e-9, 5.09e-7). They are a near-defective real cluster resolving
differently. **So the 62 should be recomputed with one threshold before it is
used**, and any amplification ratio built on it inherits the mismatch. The
binning itself is not the culprit: it only merges, seven of the 62 bins holding
more than one combination, and splits nothing.

---

## Heat: the section's own run contradicts it

No drive P creates oscillation here: zero modes at every P in the table below.
The reason once given, that the drive shifts the operating point without
touching the operator structure, is not the reason.

The drive enters as J_nl[i,j] = α·W[i,j]·dS_i/τ_i, a **per-row** gain. Row
scaling of W is exactly the magnitude condition's failure mode, so the drive
does change the operator structure, and the script's own residual says so at
every P, P = 0 included:

| P | oscillatory modes | residual | row gains, E seats | row gains, I seats |
|---|---|---|---|---|
| 0.0 | 0 | 1.38e-03 | 0.0071 – 0.0071 | 0.0012 – 0.0012 |
| 2.0 | 0 | 4.39e-03 | 0.0830 – 0.0851 | 0.0626 – 0.0669 |
| 3.5 | 0 | 4.43e-02 | 0.2809 – 0.2929 | 0.4616 – 0.4990 |
| 4.0 | 0 | 3.31e-02 | 0.3218 – 0.3250 | 0.3835 – 0.4919 |
| 8.0 | 0 | 1.74e-03 | 0.0070 – 0.0095 | 0.0001 – 0.0004 |

The two populations carry different sigmoid parameters, so the gain splits by
type. At P = 0 the E seats sit at 0.0071 and the I seats at 0.0012, a factor of
six, which is where the residual at zero drive comes from.

The residual at P = 3.5 is 4.4e-2 and no modes appear, while the coupling
sweep credits 3.45e-2 with releasing six activity frequencies. The two numbers
are not directly comparable, being ‖R_off‖/‖J‖ on different objects with
different denominators, but they are the same order, and one releases modes
while the larger one releases none. So breaking the symmetry is not what
releases oscillation here. The heat result is a null, and
it is a null against the mechanism, not for it.

---

## What is left

- The palindrome is a spectral involution μ ↦ −μ − 2s. That is proved
  ([Palindrome Proof](PROOF_PALINDROME_NEURAL.md), Step 6) and it is all.
- Exactly palindromic networks oscillate, and at large enough coupling they
  go unstable. "Neural noble gases, perfectly stable, perfectly dead" was a
  reading of one seed.
- Coupling through a mediator changes the frequency count non-monotonically.
  The counts are real; at N=20 the correlation count is at the resolution of
  its own binning and should be recomputed before it is used.
- No mechanism links the two. The construction never held an exact symmetry
  for the coupling to break.

The open question the sweep leaves is worth keeping: **what does the
non-monotonic frequency count in the coupling actually track**, if not
symmetry breaking? A first thing to rule out is the trivial one, that the
mediator's own two edges move eigenvalues across the real axis independently
of any palindrome.

---

## Scripts

| Script | What it computes |
|--------|-----------------|
| [veffect_exact.py](../../../simulations/neural/veffect_exact.py) | The construction above: single networks, coupling sweep, drive sweep |
| [veffect_and_heat.py](../../../simulations/neural/veffect_and_heat.py) | Approximate networks, thermal window, 2× law |

---

*See also:*
[Palindrome Proof](PROOF_PALINDROME_NEURAL.md) (the algebraic condition),
[Hierarchy of Incompleteness](../../HIERARCHY_OF_INCOMPLETENESS.md) (C=1 as dead end),
[V-Effect Palindrome](../../../experiments/V_EFFECT_PALINDROME.md) (quantum V-Effect)
