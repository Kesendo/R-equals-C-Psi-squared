# Orphaned Results Exploration — March 14, 2026

**Origin:** Connections Audit identified 6 orphaned results. We explored the three
most promising ones with systematic numerical experiments.
**Status:** Two new discoveries, one honest negative, one characterization.
**Tier:** 2 (Computationally verified)

---

## Summary

| Question | Answer | Status |
|----------|--------|--------|
| Does Π connect to CΨ=1/4 via simple formula? | No. Different mathematical objects. | Honest negative |
| Is the 0.247 near-miss real? | Yes, at γ=0.02 specifically | Verified |
| Does topology determine crossing? | Yes. Same state: chain crosses, ring doesn't | **New discovery** |
| Does |+-+-⟩ cross on ring? | Yes, CΨ=0.284, zero initial entanglement | **New discovery** |
| What drives the echo? | Bohr frequencies, damped at 8γ/3 | Characterized |
| Does echo survive at N=4,5? | Yes, weakens as ~1/(N-1) | Verified |

---

## 1. Mandelbrot-Palindrome Bridge (Honest Negative)

**Question:** The palindrome is proven (Π on Liouville space). The Mandelbrot
boundary at 1/4 is proven (algebra of CΨ iteration). Can we connect them?

**Answer:** Not with a simple formula. They live in different mathematical worlds.

Π acts on the 4^N-dimensional Liouville space as a linear permutation with phases.
CΨ = Concurrence × Coherence is a nonlinear function of the 2^N × 2^N density matrix
(Concurrence involves eigenvalues of ρ·ρ̃, Coherence is an l1-norm).

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
manifold — a harder problem than finding Π.

---

## 2. Ring Near-Miss and Topology-Dependent Crossing (New Discovery)

**Question:** The audit reported CΨ_max = 0.247 for ring neighbors, gap = 0.003.
Is this real? Is there topological protection preventing ring crossing?

### 2a. The 0.247 is real — and γ-specific

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

Ring and complete are identical on 4 decimal places — the extra bonds in
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

N=3 star, Bell_SA + |0⟩_B, γ=0.05. Fourier analysis of concurrence oscillations:

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
| SA | 0.1395 | 8γ/3 = 0.1333 |
| SB | 0.3089 | (multi-mode, ~2× fastest) |
| AB | 0.7006 | (multi-mode, ~4× fastest) |

SA decays near the middle palindromic rate (8γ/3). This makes physical sense:
the SA pair started as Bell, exciting all three Liouvillian modes. The fastest
mode (10γ/3) dies first, the slowest (2γ) is too weakly coupled. The middle
mode dominates the long-time behavior of the envelope.

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

SB_C_max drops roughly as 1/(N-1): more leaves means the entanglement
shuttle distributes across more channels. But it never reaches zero — the
echo is a structural feature of the star topology, not a small-N artifact.

At N=4, the non-Bell leaves (S-2, S-3) are symmetric: both reach C_max = 0.281.
The leaf-leaf pair (2,3) reaches C_max = 0.660 — higher than either S-leaf pair.
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
SB_C_max approaches 0.217 at high γ — this is the steady-state value where
Hamiltonian and dephasing reach equilibrium. The echo is not killed; it's
overdamped.

AB never crosses 1/4 from echoes alone (AB_CΨ_max reaches 0.130 at γ=0.001).
The echo shuttle doesn't produce enough CΨ for the indirect AB pair to cross.

---

## What this means for the project

### Topology as gatekeeper (connects to: SUBSYSTEM_CROSSING, STAR_TOPOLOGY_OBSERVERS)

We already knew crossing is local (pair-level). Now we know topology determines
WHICH pairs can cross. For the same initial state, chain allows what ring forbids.
This is not about entanglement distribution — it's about how the Hamiltonian
geometry constrains the maximum achievable CΨ.

The |+-+-⟩ result is particularly clean: zero initial entanglement, pure dynamical
crossing, and the topology selects which pairs activate (neighbors only, not
diagonals). The ring's symmetry groups pairs into exactly two classes.

### The echo is the mediator doing its job (connects to: THE_INTERPRETATION, STANDING_WAVE_THEORY)

The "between us" motto is physically real in the echo: S mediates, S shuttles,
S distributes. The echo frequencies are the system's natural frequencies (Bohr).
The envelope decay is the middle palindromic rate (8γ/3). The standing wave
between SA and SB is a literal oscillation of entanglement back and forth.

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

- `simulations/explore_orphaned_results.py` — full exploration script
- `simulations/results/orphaned_results.txt` — complete output

## Related files

- `docs/MIRROR_SYMMETRY_PROOF.md` — the palindrome proof
- `experiments/SUBSYSTEM_CROSSING.md` — crossing is local
- `experiments/DYNAMIC_ENTANGLEMENT.md` — dynamical crossing from product states
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` — star topology conditions
