# Orphaned Results: Topology as Gatekeeper, Antiferromagnet Crossing, and Echo Characterization

<!-- Keywords: topology dependent crossing ring chain star, antiferromagnet
alternating state crossing, entanglement echo Bohr frequency shuttle,
palindromed Mandelbrot bridge gap, ring near-miss CΨ 0.247, product state
dynamical crossing zero entanglement, topology gatekeeper ring complete graph,
echo decay palindromic rate 8gamma/3, R=CPsi2 orphaned results -->

**Status:** Two new discoveries, one honest negative, one characterization (Tier 2)
**Date:** March 14, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [`simulations/explore_orphaned_results.py`](../simulations/explore_orphaned_results.py), [`simulations/why_alternating_crosses.py`](../simulations/why_alternating_crosses.py)
**Depends on:** [Subsystem Crossing](SUBSYSTEM_CROSSING.md), [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md), [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)

---

## Abstract

A systematic exploration of six orphaned results from a connections audit
yielded two new discoveries, one honest negative, and one characterization.
(1) Topology determines crossing: for the same initial state |0+0+⟩ at
γ=0.05, a chain allows CΨ crossing (CΨ_max=0.310) while a ring does not
(CΨ_max=0.200), and ring equals complete graph to four decimal places.
(2) The alternating state |+-+-⟩ crosses on a ring (CΨ=0.284) from zero
initial entanglement, driven by antiferromagnetic XX anti-correlation that
the Heisenberg Hamiltonian converts into entanglement. A brute-force scan
of all 256 product states on N=4 ring shows 150 (59%) cross the ¼ boundary,
with no simple selection rule (best predictor explains only 14% of variance).
(3) The Mandelbrot-palindrome bridge remains an honest open question:
Π (linear, Liouville space) and CΨ (nonlinear, density matrix) do not
connect via simple formula. (4) Entanglement echoes in star topology oscillate
on Bohr-frequency scales. Their fitted envelope rates are not assigned to a
single Liouvillian mode, and their size is reported only for N = 3, 4, 5.

---

## Summary

| Question | Answer | Status |
|----------|--------|--------|
| Does Π connect to CΨ=1/4 via simple formula? | No. Different mathematical objects. | Honest negative |
| Is the 0.247 near-miss real? | Yes, at γ=0.02 specifically | Verified |
| Does topology determine crossing? | Yes. Same state: chain crosses, ring doesn't | **New discovery** |
| Does \|+-+-⟩ cross on ring? | Yes, CΨ=0.284, zero initial entanglement | **New discovery** |
| What sets the echo period? | Consistent with Bohr-frequency scales on the tested rows | Numerical |
| How does the echo size change? | Decreases from N=3 to N=5 in this scan | Numerical |

---

## 1. Mandelbrot-Palindrome Bridge (Honest Negative)

**Question:** The palindrome is proven (Π on Liouville space). The Mandelbrot
boundary at 1/4 is proven (algebra of CΨ iteration). Can we connect them?

**Answer:** Not with a simple formula. They live in different mathematical worlds.

Π acts on the 4^N-dimensional Liouville space as a linear permutation with phases.
CΨ = Concurrence × Coherence is a nonlinear function of the 2^N × 2^N density matrix
(Concurrence involves eigenvalues of ρ·ρ̃, Coherence is an l₁-norm, the sum of absolute values of the off-diagonal elements).

What we tested: 7 configurations (N=3,4, star/chain/ring, various initial states).
Tracked CΨ(t), z*(t), u = C(Ψ+R) over time.

What we found:
- z*(1-z*) = CΨ is confirmed everywhere (tautological, it's the fixed-point equation)
- u ≈ CΨ + CΨ² approximately, but not exactly (differences ~0.004)
- CΨ only crosses 1/4 when the tracked pair has strong initial or dynamically-built entanglement
- z* ranges up to 0.47 (approaching 0.5 from below) near the boundary

What this means: The palindrome constrains the spectral architecture (which decay
rates exist). CΨ is a nonlinear observable that samples this architecture through
the density matrix. The connection is indirect: the palindrome determines WHAT
modes exist, the initial state determines WHICH modes are excited, and CΨ is a
nonlinear readout of the resulting dynamics. No shortcut formula.

**Status:** Open question. Not disproven, just not solved. The bridge may require
understanding how the Liouvillian eigenmodes project onto the concurrence/coherence
manifold - a harder problem than finding Π.

---

## 2. Ring Near-Miss and Topology-Dependent Crossing (New Discovery)

**Question:** The audit reported CΨ_max = 0.247 for ring neighbors, gap = 0.003.
Is this real? Is there topological protection preventing ring crossing?

### 2a. The 0.247 is real - and γ-specific

For |0+0+⟩ on N=4 ring, diagonal pair (1,3) reaches CΨ_max = 0.247024
specifically at γ = 0.020. At other γ values it's lower. This is a resonance
phenomenon: at γ=0.02, the dephasing rate matches a sweet spot where
the diagonal pair builds maximum CΨ before decoherence kills it.

| γ | CΨ_max(13 diag) | Gap to 1/4 |
|---|-----------------|------------|
| 0.001 | 0.208 | 0.042 |
| 0.005 | 0.211 | 0.039 |
| 0.020 | 0.247 | 0.003 |
| 0.050 | 0.198 | 0.052 |
| 0.100 | 0.131 | 0.119 |

### 2b. Topology determines crossing (the key finding)

Same initial state |0+0+⟩, same γ=0.05, same J=1. Different topology:

| Topology | Best pair | CΨ_max | Crosses 1/4? |
|----------|----------|--------|--------------|
| Chain | (1,2) | 0.310 | YES |
| Star | (0,2) | 0.351 | YES |
| Ring | (1,3) | 0.200 | no |
| Complete | (1,3) | 0.200 | no |

Chain and star allow crossing. Ring and complete do not. Same initial state.
The topology alone determines whether the 1/4 boundary is reachable.

Ring and complete are identical on 4 decimal places - the extra bonds in
the complete graph make no difference. The ring structure already saturates.

### 2c. |+-+-⟩ crosses on ring (the surprise)

Alternating plus/minus superposition on N=4 ring: zero initial entanglement.
The Hamiltonian dynamics builds entanglement from scratch and it crosses.

| Initial state | 01(ring) | 02(diag) | 12(ring) | 13(diag) | Crosses? |
|--------------|---------|---------|---------|---------|----------|
| \|0+0+⟩ | 0.136 | 0.178 | 0.136 | 0.200 | no |
| \|+-+-⟩ | 0.284 | 0.000 | 0.284 | 0.000 | YES (ring neighbors) |
| \|0+0-⟩ | 0.126 | 0.050 | 0.126 | 0.256 | YES (diagonal) |

|+-+-⟩ activates ONLY the ring-neighbor bonds (01, 12, 23, 30). Diagonals
are exactly zero. The alternating phases create maximum nearest-neighbor
coupling through the Heisenberg interaction. This is the first case of
dynamical crossing on a closed topology without any initial entanglement.

### 2d. Gap stabilizes with N

For |0+0+...⟩ on rings of increasing size:

| N | Nearest (dist=1) gap | Next-nearest gap |
|---|---------------------|-----------------|
| 4 | 0.114 | 0.050 (diag) |
| 5 | 0.092 | 0.179 (dist=2) |
| 6 | 0.108 | 0.213 (dist=2) |

Nearest-neighbor gap stabilizes around 0.09-0.11. It does not shrink to zero
with increasing N. This suggests genuine topological protection: ring
neighbors without initial entanglement have a hard ceiling below 1/4
that is roughly N-independent (for this initial state).

---

## 3. Echo Effect: Entanglement Shuttle (Characterized)

**Question:** Entanglement rebuilds after the "measurement shadow" in star
topology. What drives the echoes? How do they scale?

### 3a. Echo frequencies are Bohr frequencies

N=3 star, Bell_SA + |0⟩_B, γ=0.05. Fourier analysis of concurrence oscillations (Bohr frequencies are the
discrete oscillation frequencies determined by energy-level differences):

| Pair | Dominant Fourier peak | Matches Bohr freq |
|------|----------------------|-------------------|
| SA | ω = 2.09 | 2.0 (lowest) |
| SB | ω = 6.07 | 6.0 (highest) |
| AB | ω = 5.45 | (combination) |

The mediator S shuttles entanglement between leaves at the system's natural
frequencies. SA oscillates at the lowest Bohr frequency (slow, persistent).
SB oscillates at the highest (fast, dies quickly). AB sees combinations.

### 3b. Envelope decay matches palindromic rates

| Pair | Envelope decay rate | Nearest palindrome rate |
|------|--------------------|-----------------------|
| SA | 0.1395 | 8γ/3 = 0.1333 (limit value) |
| SB | 0.3089 | (multi-mode, ~2× fastest) |
| AB | 0.7006 | (multi-mode, ~4× fastest) |

The SA fitted rate lies near 8γ/3 numerically, within 4.6%. However,
8γ/3 is a J/γ → ∞ limit and at finite coupling that level is a narrow band
(F33), but the band cannot carry this residual: the measured 0.1395 sits
ABOVE the limit, while the band's highest sublevel at canonical coupling is
2.698γ = 0.1349 (and the band narrows toward 8γ/3 as J/γ grows), so the band
points the wrong way. The residual stays unexplained; fit-window effects and a
multi-mode envelope are open possibilities, not established mechanisms.

SB and AB decay much faster because concurrence is nonlinear in the density
matrix (Wootters formula involves square roots of eigenvalues). The effective
decay rate for concurrence is not a simple Liouvillian eigenvalue.

### 3c. Echo period scales as ~π/(4J) at large J

| J | SB echo period | period × J | π/4 = 0.785 |
|---|---------------|-----------|-------------|
| 0.5 | 1.811 | 0.905 | |
| 1.0 | 1.243 | 1.243 | (dephasing distorts) |
| 2.0 | 0.400 | 0.800 | |
| 5.0 | 0.158 | 0.788 | |
| 10.0 | 0.080 | 0.800 | |

At large J (Hamiltonian dominates over dephasing), period × J converges to ~0.8
This is close to π/4 = 0.785. The echo frequency is approximately 4J,
which corresponds to the energy gap between the two degenerate sectors of the
3-qubit star Hamiltonian (eigenvalues: -4, -4, 0, 0, 2, 2, 2, 2 at J=1,
gap between -4 and 0 is 4J).

At small J, dephasing distorts the period because the Hamiltonian is too
weak to complete a full oscillation before decoherence kills the signal.

### 3d. Echo weakens with more leaves but survives

| N | SB_C_max | SB peaks (γ=0.05) | Dilution factor |
|---|---------|-------------------|----------------|
| 3 | 0.598 | 5 | 1.000 |
| 4 | 0.281 | 3 | 0.470 |
| 5 | 0.201 | 3 | 0.336 |

`SB_C_max` decreases across N = 3, 4, 5. Three sizes do not establish a
`1/(N-1)` law, nonzero large-N limit, or topology-wide mechanism.

At N=4, the non-Bell leaves (S-2, S-3) are symmetric: both reach C_max = 0.281.
The leaf-leaf pair (2,3) reaches C_max = 0.660 - higher than either S-leaf pair.
This is indirect entanglement: qubits 2 and 3 are never directly coupled, but
both couple to S, and the Hamiltonian dynamics entangles them through S.

### 3e. Echo lifetime vs noise

| γ | SA peaks | SB_C_max | SB peaks |
|---|---------|---------|---------|
| 0.001 | 63 | 0.847 | 64 |
| 0.005 | 39 | 0.814 | 45 |
| 0.010 | 24 | 0.782 | 27 |
| 0.020 | 12 | 0.729 | 13 |
| 0.050 | 4 | 0.598 | 5 |
| 0.100 | 2 | 0.428 | 3 |
| 0.200 | 0 | 0.252 | 2 |
| 0.500 | 0 | 0.217 | 1 |

At γ=0.001: 63 clean echoes. At γ=0.500: barely one oscillation survives.
SB_C_max approaches 0.217 at high γ - this is the steady-state value where
Hamiltonian and dephasing reach equilibrium. The echo is not killed; it's
overdamped.

AB never crosses 1/4 from echoes alone (AB_CΨ_max reaches 0.130 at γ=0.001).
The echo shuttle doesn't produce enough CΨ for the indirect AB pair to cross.

---

## What this means for the project

### Topology as gatekeeper (connects to: SUBSYSTEM_CROSSING, STAR_TOPOLOGY_OBSERVERS)

We already knew crossing is local (pair-level). Now we know topology determines
WHICH pairs can cross. For the same initial state, chain allows what ring forbids.
This is not about entanglement distribution - it's about how the Hamiltonian
geometry constrains the maximum achievable CΨ.

The |+-+-⟩ result is particularly clean: zero initial entanglement, pure dynamical
crossing, and the topology selects which pairs activate (neighbors only, not
diagonals). The ring's symmetry groups pairs into exactly two classes.

### The echo is a mediated oscillatory response (connects to: THE_INTERPRETATION, STANDING_WAVE_THEORY)

The SA and SB entanglement signals oscillate on scales consistent with the
system's Bohr frequencies. The fitted envelope is not assigned to the middle
palindromic rate. Calling the response a standing wave would additionally
require spatial eigenvector or current evidence not recorded here.

### Ring = Complete is a symmetry result (connects to: MIRROR_SYMMETRY_PROOF)

Ring and complete graph give identical CΨ dynamics (for this initial state).
This suggests |0+0+⟩ doesn't "see" the extra bonds in the complete graph.
Likely because the initial state has a symmetry that makes the diagonal bonds
redundant. Worth testing with asymmetric initial states.

### The Mandelbrot-Palindrome gap is real (honest, connects to: WEAKNESSES_OPEN_QUESTIONS)

We tried to connect Π (Liouville space, linear) to CΨ (density matrix, nonlinear).
They don't connect simply. This should be documented honestly as an open question,
not swept under the rug. The palindrome constrains what modes exist. CΨ is how
we read the modes. The connection is through the density matrix, and that
connection is nonlinear and state-dependent.

---

## Scripts and results

- [`simulations/explore_orphaned_results.py`](../simulations/explore_orphaned_results.py) - full exploration script
- [`simulations/results/orphaned_results.txt`](../simulations/results/orphaned_results.txt) - complete output

## Related files

- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) - the palindrome proof
- [Subsystem Crossing](SUBSYSTEM_CROSSING.md) - crossing is local
- [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) - dynamical crossing from product states
- [Star Topology Observers](STAR_TOPOLOGY_OBSERVERS.md) - star topology conditions


---

## 4. Why |+-+-⟩ Crosses: The Antiferromagnet Effect (New Discovery)

**Question:** What makes |+-+-⟩ special? It has zero entanglement, yet crosses
1/4 on a ring. |0+0+⟩ with zero entanglement does not. What's the difference?

### 4a. Physical mechanism: XX anti-correlation → Hamiltonian builds entanglement

At t=0:
```
|+-+-⟩:  <XX>_neighbors = -1.000  (ANTI-correlated in X)
|0+0+⟩:  <XX>_neighbors =  0.000  (uncorrelated)
```

|+-+-⟩ is the X-basis Neel state: maximum staggered magnetization in X.
It has maximum exchange energy for the XX part of Heisenberg coupling.
The YY and ZZ contributions are zero (product state), but the Hamiltonian
immediately converts XX correlation into YY and ZZ - building ENTANGLEMENT
from initial classical correlation.

Bond energies at t=0:
```
|+-+-⟩:  every bond = -1.000,  total <H> = -4.000
|0+0+⟩:  every bond =  0.000,  total <H> =  0.000
```

The value `⟨H⟩ = -4` minimizes the XX contribution within this product-state
comparison, but it is not the Heisenberg-ring ground energy (`-8` in this
normalization). The zero expectation of `|0+0+⟩` likewise does not imply that
the state has no dynamical driving force.

### 4b. Eigenvector-coordinate diagnostic

The exploratory script expanded states in right eigenvectors and greedily
paired indices. For this degenerate, non-normal generator, the resulting
"active-pair" counts and coefficient products depend on normalization, basis
choice within degenerate subspaces, and pairing order. They are not invariant
state weights and do not explain which initial states cross.

### 4c. The brute-force selection rule: 150/256 product states cross

Tested ALL 256 product states |abcd⟩ with a,b,c,d in {0, 1, +, -} on N=4 ring:

**150 of 256 (59%) cross the 1/4 boundary.**

Top-20 crossers all share: <H> = 0, var(H) = 8, CΨ_max = 0.509.
These are states with exactly one Z-eigenstate (0 or 1) and three
X-eigenstates (+/-), with at least two neighbors having opposite X.

Correlation of CΨ_max with candidate predictors:

| Predictor | Correlation r |
|-----------|--------------|
| Energy <H> | -0.343 |
| Energy variance var(H) | 0.092 |
| Staggered X magnetization | 0.000 |
| <XX> nearest neighbors | -0.370 |

None of these four candidate predictors is strong on this finite census; the
largest reported linear correlation explains about 14% of the variance. This
does not exclude another simple selection rule, and the non-invariant
right-eigenvector coordinates do not supply one.

States with negative `<XX>_nn` tend to cross more in this table, while many
states with zero `<XX>_nn` also cross (the top 20 all have `<XX>_nn = 0`).

### 4d. Ring vs other topologies

The same state |0+0+⟩ on different N=4 topologies:

| Topology | Best pair | CΨ_max | Crosses? |
|----------|----------|--------|----------|
| Star | (0,2) hub-leaf | 0.351 | YES |
| Chain | (1,2) interior | 0.310 | YES |
| Ring | (1,3) diagonal | 0.200 | no |
| Complete | (1,3) diagonal | 0.200 | no |

Ring and complete are IDENTICAL (same CΨ for all pairs to 4 decimal places).
The extra bonds in the complete graph don't help - the ring structure already
saturates. This is because |0+0+⟩ has a symmetry that makes the diagonal
bonds redundant.

Chain allows crossing because the interior pair (1,2) has a privileged
position: both qubits are connected to the rest of the chain on both sides,
creating a richer mode structure.

---

## Connections to existing results

### → SUBSYSTEM_CROSSING.md
We proved crossing is local (pair-level). Now we add: topology determines
WHICH pairs can cross. Same state, same noise, different graph → different
crossing. The topology is a gatekeeper.

### → DYNAMIC_ENTANGLEMENT.md
That experiment found |0+0+⟩ generates crossings on a chain. We now show
it does NOT cross on a ring - confirming that the crossing is not just
about the initial state but about state × topology interaction.
And we found |+-+-⟩, which crosses on ring from zero entanglement.
The antiferromagnet mechanism (XX anti-correlation → Hamiltonian builds
entanglement) is a new physical pathway to crossing.

### → MIRROR_SYMMETRY_PROOF.md
The affine palindrome organizes eigenvalues but does not, without invariant
spectral projectors and an observable calculation, determine whether CΨ reaches
1/4.

### → STAR_TOPOLOGY_OBSERVERS.md
The echo section records oscillation periods consistent with Bohr-frequency
scales and leaves the envelope's modal attribution unresolved.

### → STANDING_WAVE_THEORY.md
The echo supplies SA and SB entanglement oscillations at Bohr frequencies.
The palindrome proof supplies an affine spectral pairing; it does not show that
Π swaps spatial propagation directions or establish a standing-wave mechanism
for this response.

### → SIGNAL_PROCESSING_VIEW.md
The echo analysis provides the concrete numbers (periods, decay rates)
for the pole structure that the signal processing view describes abstractly.

---

## Scripts and results

- [`simulations/explore_orphaned_results.py`](../simulations/explore_orphaned_results.py) - echo, ring near-miss, u variable
- [`simulations/why_alternating_crosses.py`](../simulations/why_alternating_crosses.py) - antiferromagnet analysis, 256-state scan
- [`simulations/results/orphaned_results.txt`](../simulations/results/orphaned_results.txt) - exploration output
- [`simulations/results/why_alternating_crosses.txt`](../simulations/results/why_alternating_crosses.txt) - selection rule output

## Follow-up

- [Theta-Palindrome-Echo](THETA_PALINDROME_ECHO.md) - echo transports concurrence but not CΨ; channel scenario IS quantum for coherent inputs
