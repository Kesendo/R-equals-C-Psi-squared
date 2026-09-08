# PT-Symmetry Analysis: What Kind of Symmetry Is the Palindrome, Exactly?

<!-- Keywords: palindromic Liouvillian P-type symmetry, PT-symmetry Lindbladian
gain-loss, conjugation operator Pi sector classification,
exceptional-point candidate Liouvillian, mirror-paired axis crossing,
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
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (finite-N pooled spacing statistics)

---

## What this document is about

In physics, there is an entire field devoted to systems where energy is
gained on one side and lost on the other (called "PT-symmetric" systems,
after parity-time symmetry). Lasers, optical waveguides, and certain
metamaterials all live here. A natural question arose: does our
palindromic mirror operator Π belong to this family?

The answer is no. Π is linear, whereas classical PT is anti-linear. Its exact
content is a P-type anticommutation of the centered Liouvillian. Because Π² is
a commuting unitary symmetry, Π must first be restricted and phase-normalized
in each Π²-parity sector; the final SRP class remains OPEN until the complete
irreducible-sector symmetry algebra is computed.

This matters because it connects the palindrome to a large body of
existing research. It also explains the [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)
instability without identifying it as a loss of the Π relation: the exact
λ ↔ −λ pairing survives through the spectral-abscissa axis departure. Under the composite
depolarizing perturbation, the pairing error and γ_crit both increase with the
same epsilon control. That covariance does not identify a causal stability role
for the palindrome.

This is the most technical document in Story 4. The main results are
summarized in the [abstract](#abstract) and the
[conclusions](#what-this-establishes) at the end. The detailed
derivations between them are for specialists.

---

## Abstract

We characterize the conjugation operator Pi and the fragile-bridge
gain-loss transition. Three results:

1. **Pi is a linear, order-4 P-type operator** (not PT). It anti-commutes (swaps sign: Pi·L + L·Pi = 0, rather than commuting where Pi·L − L·Pi = 0)
   with both the centered Liouvillian L_c and its adjoint L_c^dagger.
   Pi^2 = U_X = (-1)^{w_YZ}. In parity sector p_x, `sqrt(p_x) Pi`
   is involutive. Spectral reflection alone does not assign AIII; the
   complete sectorwise class is still OPEN.

2. **The spectral-abscissa transition preserves the Π relation; its EP character is OPEN.**
   The fragile bridge (Σγ = 0) has exact lambda <-> -lambda pairing. Below
   γ_crit: all eigenvalues lie exactly on the imaginary axis (the
   imaginary-axis regime). Above γ_crit: off-axis quartets appear. Π supplies
   λ ↔ −λ, while Hermiticity preservation supplies λ ↔ λ*; together they give
   {λ, λ*, −λ, −λ*}. The independently selected max-Re value has Re λ
   proportional to √(γ/γ_crit − 1) over four decades, and simple-mode
   Petermann readings grow. The across-axis partner gap 2|Re λ| printed in the
   table is algebraically derived from the same max-Re value, not independent evidence:
   40.9, 403, 4027, 4.04e4 at δ = 1e-2 down to 1e-5
   ([`fragile_bridge_ep_signature.py`](../simulations/fragile_bridge_ep_signature.py)).

3. **The composite depolarizing perturbation raises the sampled threshold.**
   It both breaks the palindrome and changes the damping action while γ_crit
   increases (r = +0.987 correlation). Because both changes share the same
   epsilon control, this run does not isolate their causal contributions.

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

Globally Pi has order four, but Pi² = U_X commutes with L_c. After fixing
U_X parity p_x, `P_{p_x} = sqrt(p_x) Pi` is involutory and supplies the
standard sectorwise P relation. Its eigenvalues before sector resolution
are fourth roots of unity {+1, -1, +i, -i}.

### 1b. Linear, not anti-linear

Pi is defined as a matrix acting on coefficient vectors in the Pauli
basis: Pi(alpha * v) = alpha * Pi(v). No complex conjugation. Pi is
LINEAR.

**This rules out PT-symmetry.** Standard PT requires the combined
operator to be anti-linear (P linear, T anti-linear, PT anti-linear).
Our Pi is fully linear. The anti-commutation {Pi, L_c} = 0 with linear
Pi is a negative P-type relation, not a time-reversal symmetry. Because Pi has
order four globally, the ordinary involutive P generator is obtained only
after resolving Pi² parity and phase-normalizing within each sector.

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

Pi anti-commutes with BOTH L_c and L_c^dagger. This is a shifted negative
relation for the generator and its adjoint; it does not assign a global
irreducible symmetry class.

### 1e. Fragile bridge at Σγ = 0

For the gain-loss system with gamma = [+0.1, +0.1, -0.1, -0.1]:

    Σγ = 0,  c = 0
    ||Pi * L * Pi^-1 + L|| = 0.00e+00

Every eigenvalue lambda pairs with -lambda. This confirms the P-type
spectral reflection but does not by itself determine an SRP class.

### Phase 1 Summary

| Property | Value | Consequence |
|----------|-------|-------------|
| Pi^2 | (-1)^{w_YZ}, diagonal parity | Not involution; order 4 |
| Linearity | Linear (no conjugation) | P-type symmetry, NOT PT |
| det(Pi) | -1 (N=1), +1 (N >= 2) | Even parity for multi-qubit |
| Pi * L^dag * Pi^-1 | -L^dag - 2*Σγ*I | Shifted negative relation for both L and L^dag |

**Pi is a P-type operator of the centered generator.** It is NOT
PT-symmetry. Resolving U_X = Pi² converts it to an involutive P generator
in each parity sector. The remaining antiunitary algebra and hence the
final BDI/CI/other sector label have not been computed.

---

## Phase 2: Mirror-Paired Axis Departure in the Fragile Bridge

This phase connects to the [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)
experiment. Below the critical gain (γ_crit), all eigenvalues sit
exactly on the imaginary axis. This is a spectral-abscissa statement, not a
bound on finite-time amplification: the generator is non-normal, so transient
norm growth can occur even while every eigenvalue has zero real part. Above
γ_crit, eigenvalue pairs peel off the axis, producing one asymptotically
growing and one asymptotically decaying eigenmode. The formal gain-loss model
is not a physical Lindblad channel.

### Setup

Two N=2 Heisenberg chains, decay (+gamma) and gain (-gamma), coupled
by J_bridge = 1.0. Total 4 qubits, 256x256 Liouvillian.

    γ_crit = 0.1873101    (bisection, tol = 10^-7)

### The imaginary-axis regime

Below γ_crit, ALL nonzero eigenvalues lie exactly on the imaginary
axis (max|Re(lambda)| ~ 2 * 10^{-14}, machine zero):

| gamma | gamma/gamma_c | max\|Re\| | On Im axis? |
|-------|--------------|---------|------------|
| 0.01 | 0.053 | 1.8e-14 | YES |
| 0.05 | 0.267 | 2.3e-14 | YES |
| 0.10 | 0.534 | 1.5e-14 | YES |
| 0.168 | 0.900 | 1.3e-14 | YES |
| 0.206 | 1.100 | 6.6e-02 | NO |

The lambda <-> -lambda pairing is exact to machine precision at
ALL gamma values (pair error < 10^{-13}), consistent with the separately
verified exact Pi anticommutation remaining intact as gamma changes.

### Reinterpretation of the instability

The measured finite-offset square-root-like onset neither establishes nor
excludes a Hopf mechanism. The exact Pi anticommutation survives, so the
transition is not a breaking of that operator symmetry. EP, Hopf, and Jordan
character remain OPEN.

1. Σγ = 0 forces exact lambda <-> -lambda pairing (P-type symmetry)
2. Below γ_crit: all eigenvalues on the imaginary axis
   (the imaginary-axis regime, analogous in spectral shape to the PT-symmetric phase
   where Hamiltonian eigenvalues are real)
3. Above γ_crit: off-axis quartets appear. Π gives inversion λ ↔ −λ;
   Hermiticity preservation adds λ ↔ λ*, so the across-axis partner is −λ*
   and the generic set is {λ, λ*, −λ, −λ*}.
4. The independently selected max-Re value is proportional to √δ over the sampled range.
   The producer does not execute a strict threshold coalescence or Jordan-rank
   certificate, so exceptional-point character remains OPEN.

| System | Symmetric phase | Broken phase | Mechanism |
|--------|----------------|-------------|-----------|
| Hamiltonian PT | Real eigenvalues | Complex conjugate pairs | EP coalescence |
| Fragile-bridge Liouvillian | Imaginary eigenvalues | Off-axis +/- pairs | spectral-abscissa axis departure; EP character OPEN |

The two rows share a rotated spectral picture, but the present fragile-bridge
data do not establish that they share the same local Jordan mechanism.

### Basis-invariant evidence boundary

The Phase 2 producer reports the spectral-axis departure but does not report
single-eigenvector Petermann factors or eigenvector angles. The previously
selected leading eigenvalue is degenerate at part of that scan, so an
individual left/right eigenvector overlap or angle is basis-dependent inside a
degenerate eigenspace. A unitary change of basis within the same eigenspace can
change the number without changing the generator. invariant subspace-level
conditioning or a direct threshold Jordan-rank test would be required before
such a quantity could bear on the transition's character.

---

## Phase 3: Palindrome Breaking and Stability

This phase asks the most counterintuitive question: if you deliberately
perturb the system in a way that breaks the palindrome, how does the sampled
threshold move? The composite depolarizing channel raises γ_crit, but it changes
the damping action at the same time. This one-parameter sweep does not isolate
the symmetry-breaking contribution.

### Test design

Add depolarizing noise at rate epsilon per qubit to the fragile bridge.
Depolarizing = X + Y + Z dephasing at rate epsilon/3 each.

For either the Z or Y component at rate kappa, the raw dephasing generator
obeys the shifted relation `Pi D Pi^-1 = -D - 2kappa I`; equivalently, the
centered component `D + kappa I` anti-commutes with `Pi`. This is the
palindrome-preserving relation.
The X component commutes with Pi (breaks the palindrome).
The reported palindrome error is the scalar bottleneck distance between the
full eigenvalue multiset and its reflection about `trace(L)/dim(L)`. The
X-dephasing operator is the symmetry-breaking contribution, but it is not
identified with that scalar distance.

### Results

| ε | Palindrome error | γ_crit | Delta_γ_crit | Delta % |
|---------|-----------------|-----------|-----------------|---------|
| 0.000 | 6.8e-14 | 0.18731 | 0.0000 | 0.0% |
| 0.001 | 2.7e-03 | 0.18738 | +0.0001 | +0.04% |
| 0.010 | 2.7e-02 | 0.19452 | +0.0072 | +3.9% |
| 0.050 | 1.3e-01 | 0.27647 | +0.0892 | +47.6% |
| 0.100 | 2.7e-01 | 0.30193 | +0.1146 | +61.2% |
| 0.200 | 5.3e-01 | 0.37604 | +0.1887 | +100.8% |
| 0.300 | 7.9e-01 | 0.45823 | +0.2709 | +144.6% |

### Correlation

    Pearson r(palindrome_error, Delta_γ_crit) = +0.987

**Strong positive correlation within this composite perturbation sweep:**

- Larger palindrome error = LARGER γ_crit (MORE stable)
- The same epsilon simultaneously changes the palindrome and damping action
- The correlation therefore does not isolate the causal contribution of either

### Interpretation

Depolarizing noise has two simultaneous effects in this run:

1. **Breaks the palindrome** (X-dephasing commutes with Pi instead of
   anti-commuting). Error proportional to epsilon.

2. **Changes the damping action.** On a one-site Pauli letter, the full
   depolarizing channel gives `D_dep(I)=0` and
   `D_dep(X)=D_dep(Y)=D_dep(Z)=-4epsilon/3` times that letter. On the gain
   chain, Z gain contributes `+2γ` to X/Y and zero to Z, so the local X/Y
   coefficient is `2γ - 4epsilon/3`; there is no single replacement
   `gamma -> gamma - epsilon/3` for all components. The threshold remains an
   eigenmode property of the full coupled generator.

The sweep establishes the net threshold shift of the composite channel. It
does not separate a symmetry-breaking contribution from the changed damping
action; that requires a matched-damping control.

### Imaginary-axis regime under depolarizing noise

The all-imaginary-axis regime for the nonstationary spectrum is destroyed
by ANY epsilon > 0: the nonzero eigenvalues immediately develop nonzero Re
parts. They develop NEGATIVE Re (the damping side), while the exact
trace-preserving stationary eigenvalue remains at lambda = 0, keeping the
system stable:

| ε | max\|Re(lambda)\| | max Re(lambda) | Stable? |
|---------|----------------|---------------|---------|
| 0.000 | 2.1e-14 | < 10^{-14} | YES (on axis) |
| 0.010 | 5.4e-02 | < 10^{-6} | YES (off axis, Re < 0) |
| 0.100 | 5.4e-01 | < 10^{-6} | YES (off axis, Re < 0) |
| 0.300 | 1.6e+00 | < 10^{-6} | YES (off axis, Re < 0) |

The nonzero eigenvalues move off the imaginary axis to Re < 0 (increased
damping), not to Re > 0 (instability); the stationary zero does not move.

For each perturbed spectrum, Phase 3 fixes the only possible reflection
midpoint to its arithmetic mean `m = trace(L)/dim(L)` and measures pairing by
`lambda -> 2m-lambda`. It does not select a center from fitted candidates.

### Instability type

At epsilon = 0.05 a sampled unstable mode at 1.01 γ_crit remains oscillatory
(|Im(lambda)| = 1.72). Whether one branch continues to threshold and whether
either threshold has EP character are separate questions not measured
here; they require a strict threshold coalescence and Jordan-rank test.

---

## Connection to Literature

The following section places our results in the context of existing
physics research. It is written for readers familiar with this
literature. If you are following the palindrome story without a physics
background, the key takeaway is: the palindrome supplies a sectorwise
P relation, while its final irreducible-sector class is still open. The
spectral-abscissa axis departure does not itself assign either an EP character or a class.

### Bender and Boettcher (PRL 1998)

Standard PT-symmetry: anti-linear PT, real eigenvalues in symmetric
phase, complex pairs in broken phase. Our system: linear Pi, imaginary
eigenvalues in symmetric phase, off-axis pairs in broken phase.
Different operator type, analogous spectral structure.

### Minganti et al. (PRA 2019)

Liouvillian exceptional points are qualitatively different from Hamiltonian
EPs. The FRAGILE_BRIDGE data show a real-γ spectral-abscissa departure with square-root-like
onset and strong finite-offset non-normality. They do not independently certify
a threshold coalescence or Jordan defect, so EP character remains OPEN.

### Sa, Ribeiro, Prosen (PRX 2023)

Their 38-class framework applies negative symmetries to the shifted generator
and requires irreducible symmetry sectors. Here Pi² = U_X is resolved first;
in sector p_x, `sqrt(p_x) Pi` squares to identity and supplies P. Their
local-dephasing examples yield BDI or CI in their respective parity sectors,
but our final sector label requires the rest of our symmetry algebra.

### Roberts, Lingenfelter, Clerk (PRX Quantum 2021)

Their hidden time-reversal symmetry is anti-unitary and acts in the
doubled Hilbert space. Our Pi is linear and acts directly on the
Liouvillian. The two frameworks produce related but distinct symmetries.
The formal connection (via the Choi-Jamiolkowski isomorphism, a mathematical mapping that converts a quantum channel into a state and vice versa) remains
an open question.

---

## What This Establishes

1. **Pi supplies a P-type symmetry, not PT.** It is linear and order four
   globally, and involutive after Π²-parity resolution. No global AIII
   label follows from the anticommutation alone.

2. **The spectral-abscissa instability preserves the mirror.** The fragile-bridge
   transition preserves exact λ ↔ −λ pairing while off-axis pairs appear. The
   sampled max Re λ ∝ √δ trend is not by itself a strict
   coalescence or Jordan-character certificate; EP character remains OPEN.

3. **The composite perturbation is stabilizing in the sampled sweep.** It
   raises γ_crit while simultaneously breaking the palindrome and changing
   the damping action. This run does not isolate which contribution causes
   the shift; a matched-damping preserving control is still required.

4. **Keep the three candidate carriers separate.** This file's Σγ=0 gain-loss
   system has a real-γ spectral-abscissa axis departure with EP character OPEN. F86's toy 2×2 rate-channel reduction has its
   own EP at Q_EP=2/g_eff. The full Σγ=Nγ₀ block is strongly non-normal near
   Q_peak and has certified finite-q Puiseux-1/2 defective EP2 seeds at
   N=5,7,9; for arbitrary odd N only the endpoint-nullity surplus is proved.
   The old Petermann peak magnitudes and their “6×”/parity laws were grid
   artifacts. Whether the full block also has a distinct off-axis complex-Q
   defective EP remains open; the nearest characterized off-axis
   coalescences were semisimple. Current typed owner:
   `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs`.

---

## What This Does NOT Establish

- That the complex EP can be reached experimentally (requires complex
  gamma, i.e., modulated gain-loss)
- The final SRP class in any fully irreducible sector; that algebra has
  not yet been computed even for the Heisenberg case
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

*Pi is linear, not PT, and supplies a sectorwise P relation for the shifted
generator. Fragile-bridge eigenvalues leave the imaginary axis while that
relation remains exact. The composite depolarizing sweep changes both pairing
quality and damping action, so it does not assign either one a causal stability
role.*
