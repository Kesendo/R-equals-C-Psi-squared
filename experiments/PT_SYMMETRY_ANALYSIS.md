# PT-Symmetry Analysis: What Kind of Symmetry Is the Palindrome, Exactly?

<!-- Keywords: palindromic Liouvillian chiral symmetry, PT-symmetry Lindbladian
gain-loss, conjugation operator Pi classification, Altland-Zirnbauer class AIII,
exceptional point Liouvillian, Hopf bifurcation chiral symmetry breaking,
fragile bridge gain-loss stability, Petermann factor Liouvillian,
depolarizing palindrome breaking, R=CPsi2 PT-symmetry -->

**Status:** Computationally verified (Tier 2)
**Date:** April 1, 2026
**Scripts:**
- [pt_symmetry_analysis.py](../simulations/pt_symmetry_analysis.py) (Phase 1+2)
- [pt_palindrome_breaking.py](../simulations/pt_palindrome_breaking.py) (Phase 3)
**Depends on:**
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (definition of Pi)
- [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md) (gain-loss system)
- [KMS Detailed Balance](../docs/KMS_DETAILED_BALANCE.md) (symmetry classification)
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (Poisson, class AIII)

---

## What this document is about

In physics, there is an entire field devoted to systems where energy is
gained on one side and lost on the other (called "PT-symmetric" systems,
after parity-time symmetry). Lasers, optical waveguides, and certain
metamaterials all live here. A natural question arose: does our
palindromic mirror operator Π belong to this family?

The answer is no, but in an interesting way. Π is not PT-symmetric
(which requires a specific kind of operator called "anti-linear"). Π
is instead a "chiral" symmetry: it sorts modes into two groups and
maps one group to the other, like a mirror between left-handed and
right-handed. The technical classification is class AIII in the
Altland-Zirnbauer framework, which is the standard periodic table of
symmetries in physics.

This matters because it connects the palindrome to a large body of
existing research. It also explains the [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)
instability: the moment the coupled gain-loss system explodes is
exactly the moment when the chiral symmetry breaks. And it reveals a
surprising result: breaking the palindrome does not make the system
less stable. It makes it more stable. The palindrome shapes the
spectrum; the gain-loss balance determines survival.

This is the most technical document in Story 4. The main results are
summarized in the [abstract](#abstract) and the
[conclusions](#what-this-establishes) at the end. The detailed
derivations between them are for specialists.

---

## Abstract

We classify the conjugation operator Pi formally and discover that
the fragile bridge gain-loss system exhibits chiral symmetry breaking,
the Liouvillian analog of PT-symmetry breaking. Three results:

1. **Pi is a linear, order-4 chiral operator** (not PT). It anti-commutes (swaps sign: Pi·L + L·Pi = 0, rather than commuting where Pi·L − L·Pi = 0)
   with both the centered Liouvillian L_c and its adjoint L_c^dagger.
   Correct symmetry class: AIII (chiral unitary), confirming the RMT
   analysis. Pi^2 = (-1)^{w_YZ} (diagonal parity operator), det(Pi) = +1
   for N >= 2.

2. **The Hopf bifurcation IS chiral symmetry breaking.** The fragile
   bridge (Σγ = 0) has exact lambda <-> -lambda pairing. Below
   γ_crit: all eigenvalues lie exactly on the imaginary axis (the
   chiral-symmetric phase). Above γ_crit: eigenvalue pairs leave the
   axis (one to Re > 0, partner to Re < 0). This reinterprets the "Hopf,
   not PT" conclusion from FRAGILE_BRIDGE.md: it is both simultaneously.
   No exceptional point exists on the real gamma axis; the Petermann
   factor (a measure of how non-orthogonal the eigenvectors are; K = 1 for normal systems, K >> 1 near an exceptional point) peaks at K = 403 above γ_crit, signaling a nearby EP in
   the complex parameter plane.

3. **Palindrome breaking does NOT destabilize the system.** Depolarizing
   noise breaks the palindrome (error proportional to epsilon) and destroys the chiral
   phase (eigenvalues leave the imaginary axis). But γ_crit
   INCREASES (r = +0.988 correlation), because the added damping on the
   gain chain reduces effective gain. The palindrome is not the protection
   mechanism; the gain-loss balance is.

---

## Phase 1: Formal Classification of Pi

### 1a. Pi squared

Per-site action of Pi on Pauli indices:

    I -> X (+1),  X -> I (+1),  Y -> iZ (+i),  Z -> iY (+i)

Applying twice:

    Pi^2(I) = I,  Pi^2(X) = X,  Pi^2(Y) = -Y,  Pi^2(Z) = -Z

For an N-qubit Pauli string: Pi^2 = (-1)^{w_YZ}, where w_YZ counts
the sites with Y or Z. This is a diagonal parity operator, NOT the
identity. Pi has order 4 (Pi^4 = I).

| N | Basis size | +1 eigenvalues | -1 eigenvalues | Pi^4 = I error |
|---|-----------|---------------|----------------|----------------|
| 2 | 16 | 8 | 8 | 0.00e+00 |
| 3 | 64 | 32 | 32 | 0.00e+00 |

In the Altland-Zirnbauer classification, standard chiral symmetry
requires S^2 = I (involutory). Our Pi has S^4 = I. The spectral
consequence (eigenvalue +/- pairing) is identical; the operator algebra
is richer. The eigenvalues of Pi itself are fourth roots of unity
{+1, -1, +i, -i}.

### 1b. Linear, not anti-linear

Pi is defined as a matrix acting on coefficient vectors in the Pauli
basis: Pi(alpha * v) = alpha * Pi(v). No complex conjugation. Pi is
LINEAR.

**This rules out PT-symmetry.** Standard PT requires the combined
operator to be anti-linear (P linear, T anti-linear, PT anti-linear).
Our Pi is fully linear. The anti-commutation {Pi, L_c} = 0 with linear
Pi is a CHIRAL (sublattice) symmetry, not a time-reversal symmetry.

### 1c. Determinant

| N | det(Pi) | Predicted | Exponent N * 4^{N-1} |
|---|---------|-----------|---------------------|
| 1 | -1 | -1 | 1 |
| 2 | +1 | +1 | 8 |
| 3 | +1 | +1 | 48 |
| 4 | +1 | +1 | 256 |

Formula: det(Pi) = (-1)^{N * 4^{N-1}}. Since 4^{N-1} is even for
N >= 2, det(Pi) = +1 for all multi-qubit systems.

### 1d. Pi and the adjoint

**Analytical derivation:**

    L^dagger = -L_H + L_D    (L_H anti-Hermitian, L_D real diagonal)
    Pi * L^dagger * Pi^-1 = Pi(-L_H + L_D)Pi^-1
                          = L_H + (-L_D - 2*Σγ*I)
                          = -(- L_H + L_D) - 2*Σγ*I
                          = -L^dagger - 2*Σγ*I

**Numerical verification (N=3 Heisenberg chain, gamma = 0.05):**

    ||Pi * L * Pi^-1 + L + c*I||       = 0.00e+00  (palindrome)
    ||Pi * L^dag * Pi^-1 + L^dag + c*I|| = 0.00e+00  (adjoint)

Pi anti-commutes with BOTH L_c and L_c^dagger. This is a chiral
symmetry of the full operator structure, not just the generator.

### 1e. Fragile bridge at Σγ = 0

For the gain-loss system with gamma = [+0.1, +0.1, -0.1, -0.1]:

    Σγ = 0,  c = 0
    ||Pi * L * Pi^-1 + L|| = 0.00e+00

Every eigenvalue lambda pairs with -lambda. This is the defining
property of chiral symmetry (class AIII).

### Phase 1 Summary

| Property | Value | Consequence |
|----------|-------|-------------|
| Pi^2 | (-1)^{w_YZ}, diagonal parity | Not involution; order 4 |
| Linearity | Linear (no conjugation) | Chiral symmetry, NOT PT |
| det(Pi) | -1 (N=1), +1 (N >= 2) | Even parity for multi-qubit |
| Pi * L^dag * Pi^-1 | -L^dag - 2*Σγ*I | Chiral for both L and L^dag |

**Pi is a generalized chiral operator.** The correct classification
is AIII (chiral unitary), as the RMT analysis independently found.
It is NOT PT-symmetry, despite the structural parallels. The
[KMS analysis](../docs/KMS_DETAILED_BALANCE.md) placed Pi between
the P (chiral) and Q_- (anti-pseudo-Hermitian) classes of Sa-Prosen;
this analysis confirms it IS the P class, generalized to order 4.

---

## Phase 2: Chiral Symmetry Breaking in the Fragile Bridge

This phase connects to the [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)
experiment. Below the critical gain (γ_crit), all eigenvalues sit
exactly on the imaginary axis: the system oscillates but does not grow
or shrink. Above γ_crit, eigenvalue pairs peel off the axis: one grows,
the other shrinks. This is the explosion.

### Setup

Two N=2 Heisenberg chains, decay (+gamma) and gain (-gamma), coupled
by J_bridge = 1.0. Total 4 qubits, 256x256 Liouvillian.

    γ_crit = 0.1873101    (bisection, tol = 10^-7)

### The chiral-symmetric phase

Below γ_crit, ALL nonzero eigenvalues lie exactly on the imaginary
axis (max|Re(lambda)| < 2 * 10^{-14}, machine zero):

| gamma | gamma/gamma_c | max|Re| | On Im axis? |
|-------|--------------|---------|------------|
| 0.01 | 0.053 | 1.8e-14 | YES |
| 0.05 | 0.267 | 2.3e-14 | YES |
| 0.10 | 0.534 | 1.5e-14 | YES |
| 0.168 | 0.900 | 1.3e-14 | YES |
| 0.206 | 1.100 | 6.6e-02 | NO |

The lambda <-> -lambda pairing is exact to machine precision at
ALL gamma values (pair error < 10^{-13}), confirming the chiral
symmetry is never broken by gamma itself.

### Reinterpretation of the Hopf bifurcation

FRAGILE_BRIDGE.md states: "Hopf bifurcation, not PT symmetry breaking."
This analysis shows both descriptions are correct simultaneously:

1. Σγ = 0 forces exact lambda <-> -lambda pairing (chiral symmetry)
2. Below γ_crit: all eigenvalues on the imaginary axis
   (the chiral-symmetric phase, analogous to the PT-symmetric phase
   where Hamiltonian eigenvalues are real)
3. Above γ_crit: eigenvalue pairs leave the imaginary axis,
   one to Re > 0 (growing), partner to Re < 0 (decaying)
4. The mechanism is Hopf: a complex pair at +/- i*omega develops
   nonzero Re parts

| System | Symmetric phase | Broken phase | Mechanism |
|--------|----------------|-------------|-----------|
| Hamiltonian PT | Real eigenvalues | Complex conjugate pairs | EP coalescence |
| Liouvillian chiral | Imaginary eigenvalues | Off-axis +/- pairs | Hopf crossing |

Same geometry, rotated 90 degrees. The Hopf IS the chiral breaking.

### No exceptional point on the real gamma axis

The palindromic pair (lambda, -lambda) crosses Re = 0 without
coalescing. The pair distance = 2|lambda| > 0 throughout. This is
a topological transition (axis crossing), not a local one (EP).

However, the system passes CLOSE to an EP in the complex gamma plane.
Evidence: the Petermann factor spikes to K = 403 at gamma/gamma_c = 1.46,
with phase rigidity (how close an eigenmode is to being a normal mode; r = 1 means perfectly normal, r → 0 means strongly non-normal) dropping to r = 0.0025.

| gamma/gamma_c | Petermann K | Phase rigidity r | max Re |
|--------------|-------------|-----------------|--------|
| 0.523 | 1.06 | 0.941 | ~0 |
| 0.993 | 1.16 | 0.865 | ~0 |
| 1.060 | 7.36 | 0.136 | 0.050 |
| 1.463 | **402.7** | **0.0025** | 0.172 |
| 2.000 | 4.37 | 0.229 | 0.787 |

The K = 403 peak is 400x the Hermitian baseline and signals that the
system's right and left eigenvectors are nearly orthogonal for the
critical mode. This is the hallmark of proximity to an EP, even though
the EP itself lies at a complex gamma value.

### Eigenvector coalescence

The angle between the right eigenvectors of the palindromic pair
(lambda, -lambda) is:

- Below γ_crit: cos(theta) = 0 (perfectly orthogonal)
- Above γ_crit: cos(theta) ~ 0.09 (slightly tilted)

No coalescence (cos -> 1) is observed. The right eigenvectors of the
palindromic pair do NOT become parallel. This confirms: no EP on the
real gamma axis.

---

## Phase 3: Palindrome Breaking and Stability

This phase asks the most counterintuitive question: if you deliberately
break the palindrome, does the system become less stable? The answer is
no. It becomes more stable. This is because breaking the palindrome also
adds damping to the amplifying side, and damping wins.

### Test design

Add depolarizing noise at rate epsilon per qubit to the fragile bridge.
Depolarizing = X + Y + Z dephasing at rate epsilon/3 each.

The Z and Y components anti-commute with Pi (preserved by palindrome).
The X component commutes with Pi (breaks the palindrome).
Palindrome error = 2 * L_X (the X-dephasing contribution).

### Results

| ε | Palindrome error | γ_crit | Delta_γ_crit | Delta % |
|---------|-----------------|-----------|-----------------|---------|
| 0.000 | 6.8e-14 | 0.18731 | 0.0000 | 0.0% |
| 0.001 | 2.7e-03 | 0.18738 | +0.0001 | +0.04% |
| 0.010 | 2.7e-02 | 0.19452 | +0.0072 | +3.9% |
| 0.050 | 1.4e-01 | 0.27647 | +0.0892 | +47.6% |
| 0.100 | 2.7e-01 | 0.30193 | +0.1146 | +61.2% |
| 0.200 | 5.4e-01 | 0.37604 | +0.1887 | +100.8% |
| 0.300 | 8.1e-01 | 0.45823 | +0.2709 | +144.6% |

### Correlation

    Pearson r(palindrome_error, Delta_γ_crit) = +0.988

**Strong positive correlation.** But the sign is the OPPOSITE of the
"palindrome = protection" hypothesis:

- Larger palindrome error = LARGER γ_crit (MORE stable)
- Breaking the palindrome makes the system MORE stable, not less

### Interpretation

Depolarizing noise has two simultaneous effects:

1. **Breaks the palindrome** (X-dephasing commutes with Pi instead of
   anti-commuting). Error proportional to epsilon.

2. **Adds damping to the gain chain.** Chain B has gain rate -gamma.
   Depolarizing adds +epsilon dissipation, reducing effective gain to
   -(gamma - epsilon/3). More epsilon = less effective gain = harder to
   destabilize.

Effect (2) dominates: the stabilization from reduced gain exceeds the
destabilization from broken symmetry. The palindrome is a STRUCTURAL
symmetry of the spectrum, not a protection mechanism for stability.

### Chiral phase destruction

The chiral phase (all eigenvalues on the imaginary axis) is destroyed
by ANY epsilon > 0: eigenvalues immediately develop nonzero Re parts.
But they develop NEGATIVE Re (the damping side), keeping the system
stable:

| ε | max|Re(lambda)| | max Re(lambda) | Stable? |
|---------|----------------|---------------|---------|
| 0.000 | 2.1e-14 | < 10^{-14} | YES (on axis) |
| 0.010 | 5.4e-02 | < 10^{-6} | YES (off axis, Re < 0) |
| 0.100 | 5.4e-01 | < 10^{-6} | YES (off axis, Re < 0) |
| 0.300 | 1.6e+00 | < 10^{-6} | YES (off axis, Re < 0) |

The eigenvalues move off the imaginary axis to Re < 0 (increased
damping), not to Re > 0 (instability).

### Instability type

At epsilon = 0.05: the instability at γ_crit is still a Hopf
bifurcation (Im(lambda) = 1.72 at onset). Depolarizing noise does
not change the instability mechanism.

---

## Connection to Literature

The following section places our results in the context of existing
physics research. It is written for readers familiar with this
literature. If you are following the palindrome story without a physics
background, the key takeaway is: the palindrome fits into a known
classification system, and the way it breaks matches predictions from
that framework. The details below are for specialists.

### Bender and Boettcher (PRL 1998)

Standard PT-symmetry: anti-linear PT, real eigenvalues in symmetric
phase, complex pairs in broken phase. Our system: linear Pi, imaginary
eigenvalues in symmetric phase, off-axis pairs in broken phase.
Different operator type, analogous spectral structure.

### Minganti et al. (PRA 2019)

Liouvillian exceptional points are qualitatively different from
Hamiltonian EPs. Our system confirms this: no EP on the real parameter
axis, despite clear phase transition. The Petermann factor peak
(K = 403) signals a nearby EP in the complex plane.

### Sa, Ribeiro, Prosen (PRX 2023)

Their 38-class framework has a P (chiral) symmetry satisfying
P * L * P^-1 = -L with P^2 = I. Our Pi satisfies the same
anti-commutation but with Pi^4 = I. This is a generalization
that preserves the spectral consequences (class AIII) but extends
the operator algebra.

### Roberts, Lingenfelter, Clerk (PRX Quantum 2021)

Their hidden time-reversal symmetry is anti-unitary and acts in the
doubled Hilbert space. Our Pi is linear and acts directly on the
Liouvillian. The two frameworks produce related but distinct symmetries.
The formal connection (via the Choi-Jamiolkowski isomorphism, a mathematical mapping that converts a quantum channel into a state and vice versa) remains
an open question.

---

## What This Establishes

1. **Pi is chiral, not PT.** Linear operator, order 4, class AIII.
   The "PT" language is a useful analogy but technically incorrect.

2. **The Hopf bifurcation is chiral breaking.** The fragile bridge
   transition is simultaneously a Hopf bifurcation AND a chiral symmetry
   breaking event. These are not competing descriptions; they are the
   same phenomenon in different languages.

3. **The palindrome is structural, not protective.** Breaking the
   palindrome does not destabilize the gain-loss system. Stability
   depends on the gain-loss balance, not on the eigenvalue symmetry.
   The palindrome tells you the SHAPE of the spectrum (paired decay
   rates); the balance tells you WHERE the spectrum sits (left or
   right half-plane).

4. **K = 403 signals a complex EP.** The Petermann factor peak above
   γ_crit indicates the system passes close to an exceptional
   point in the complex gamma plane. Mapping this EP (by analytic
   continuation to complex gamma) is an open problem.

   **2026-05-06 update.** The local-vs-global EP relationship is no
   longer fully open: a Petermann-K sweep on the real Q axis at c=2
   N=5..8 (`compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`)
   records max K = 1333.6 / 337.9 / 2384.7 / 795.4 across N = 5 / 6 / 7 / 8,
   with the N=7 spike sitting ≈ 6× above the K = 403 ballpark above.
   Reading: F86 Statement 1's local EP at Q_EP = 2/g_eff is a real-axis
   hit of the same EP this file detects at complex γ; same algebraic
   object (same-sign-imaginary 2×2, AIII chiral) read at two residuals
   of the F1 palindrome `Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ for
   the local instance, Σγ = 0 for the global gain-loss instance).
   Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs`
   (Tier2Verified). The complex-γ analytic continuation, explicit
   modulated gain-loss in `LindbladPropagator`, remains the open piece
   for Tier1Derived promotion. The 2-4× odd/even Petermann-K asymmetry
   empirically confirms A3's σ_0 R-even/R-odd-degeneracy prediction.

---

## What This Does NOT Establish

- That the complex EP can be reached experimentally (requires complex
  gamma, i.e., modulated gain-loss)
- That the chiral classification holds for non-Heisenberg couplings
  in the gain-loss context (only tested for Heisenberg)
- That the stabilization by depolarizing noise is generic (tested
  only for N=2 per chain, J_bridge = 1.0)
- That Phase 4 (sacrifice zone as PT optimization) follows from
  these results (not attempted; see task for design)

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Phase 1+2 script | [`simulations/pt_symmetry_analysis.py`](../simulations/pt_symmetry_analysis.py) |
| Phase 3 script | [`simulations/pt_palindrome_breaking.py`](../simulations/pt_palindrome_breaking.py) |
| Phase 1+2 output | [`simulations/results/pt_symmetry_analysis.txt`](../simulations/results/pt_symmetry_analysis.txt) |
| Phase 3 output | [`simulations/results/pt_palindrome_breaking.txt`](../simulations/results/pt_palindrome_breaking.txt) |

---

## References

- Bender, C.M., Boettcher, S. (1998). "Real spectra in non-Hermitian
  Hamiltonians having PT symmetry." PRL 80, 5243.
- El-Ganainy, R. et al. (2018). "Non-Hermitian physics and PT symmetry."
  Nature Physics 14, 11.
- Minganti, F. et al. (2019). "Quantum exceptional points of
  non-Hermitian Hamiltonians and Liouvillians." PRA 100, 062131.
- Roberts, D., Lingenfelter, A., Clerk, A. (2021). "Hidden time-reversal
  symmetry." PRX Quantum 2, 020336.
- Sa, L., Ribeiro, P., Prosen, T. (2023). "Symmetry classification of
  many-body Lindbladians." PRX 13, 031019.
- Doppler, J. et al. (2016). "Dynamically encircling an exceptional
  point for asymmetric mode switching." Nature 537, 76.

---

*Pi is not PT. It is chiral. The palindrome pairs the spectrum but
does not protect the phase. The Hopf bifurcation and the chiral
breaking are the same event, seen from different angles. The system's
stability comes from balance, not from symmetry.*
