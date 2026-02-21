# N-Scaling Barrier: Why Crossing Fails for Large Systems

**Date**: 2026-02-18
**Status**: Independently verified (Tier 2)
**Depends on**: CROSSING_TAXONOMY.md, NOISE_ROBUSTNESS.md

---

## 1. The Question

All prior experiments used Bell+ (N = 2) or single-qubit (N = 1) systems.
Three questions were open:

1. Does the crossing taxonomy hold for GHZ and W states?
2. Does correlation remain Type A for N > 2?
3. At what N does crossing become impossible?

The answer to all three reveals a fundamental scaling problem.

## 2. Setup

| Parameter | Value |
|-----------|-------|
| **Hamiltonian** | Heisenberg ring (J = 1, h = 0) |
| **Noise type** | local dephasing (σ_z per qubit) |
| **γ_base** | 0.05 |
| **Time step** | dt = 0.01 |
| **Bridge** | All five, focus on correlation |

Variable: **state** × **N**

States tested: Bell+ (N = 2), GHZ (N = 3, 4), W (N = 3, 4).

## 3. The Scaling Problem

### 3.1 Why Large Systems Cannot Cross

The dynamic Ψ is defined as:

```
Ψ(t) = C_l1(ρ(t)) / (d − 1)
```

where d = 2^N is the Hilbert space dimension and C_l1 is the l1-norm of
coherence (sum of absolute values of off-diagonal elements).

The maximum possible l1-coherence for any valid density matrix of dimension
d is exactly d − 1. So Ψ ∈ [0, 1] and Ψ = 1 is physically reachable.

But standard quantum states have l1-coherence that grows much slower than
d − 1:

| State | N | l1-coherence | d − 1 | Ψ(0) | C·Ψ(0) | Crosses ¼? |
|-------|---|-------------|-------|-------|---------|-----------|
| Bell+ | 2 | 1 | 3 | 0.333 | 0.333 | **YES** |
| GHZ | 3 | 1 | 7 | 0.143 | 0.143 | NO |
| GHZ | 4 | 1 | 15 | 0.067 | 0.067 | NO |
| GHZ | 5 | 1 | 31 | 0.032 | 0.032 | NO |
| W | 3 | 2 | 7 | 0.286 | 0.286 | **YES** |
| W | 4 | 3 | 15 | 0.200 | 0.200 | NO |
| W | 5 | 4 | 31 | 0.129 | 0.129 | NO |

**GHZ states** have exactly l1 = 1 for all N (one off-diagonal pair: |00...0⟩⟨11...1| + h.c.). The denominator d − 1 = 2^N − 1 grows exponentially. GHZ fails at N ≥ 3.

**W states** have l1 = N − 1 (polynomial growth), but d − 1 = 2^N − 1 (exponential growth). W fails at N ≥ 4.

### 3.2 The Minimum l1 for Crossing

For C·Ψ ≥ ¼ with C = 1.0 (Type A, best case):

```
Ψ(0) ≥ 1/4
l1 / (d − 1) ≥ 1/4
l1 ≥ (2^N − 1) / 4
```

| N | d − 1 | l1 required | GHZ l1 | W l1 | Sufficient? |
|---|-------|-------------|--------|------|-------------|
| 2 | 3 | 0.75 | 1 | 1 | Both YES |
| 3 | 7 | 1.75 | 1 | 2 | W only |
| 4 | 15 | 3.75 | 1 | 3 | Neither |
| 5 | 31 | 7.75 | 1 | 4 | Neither |
| 6 | 63 | 15.75 | 1 | 5 | Neither |
| 10 | 1023 | 255.75 | 1 | 9 | Neither |

The required coherence grows exponentially. No standard entangled state
can keep up.

## 4. Simulation Results

### 4.1 GHZ N = 3: Everything Becomes Type C

| Bridge | C(0) | C(t=5) | Ψ(0) | C·Ψ(0) | Crosses? |
|--------|------|--------|-------|---------|----------|
| correlation | 1.000 | 0.810 | 0.143 | 0.143 | NO |
| concurrence | 0.500 | 0.500 | 0.143 | 0.071 | NO |
| mutual_info | 1.000 | 0.691 | 0.143 | 0.143 | NO |
| mutual_purity | 0.500 | 0.500 | 0.143 | 0.071 | NO |
| overlap | 0.250 | 0.250 | 0.143 | 0.036 | NO |

GHZ N = 3 has Ψ(0) = 1/7 = 0.143. Even with C = 1.0 (correlation),
C·Ψ = 0.143 < 0.25. **Nothing crosses.** The entire taxonomy collapses
to Type C. Not because the mechanisms changed, but because the starting
point is already below the boundary.

Note: concurrence returns C = 0.5 constant for GHZ N = 3. This is
different from Bell+ where concurrence starts at C = 1.0 and decays
(Type B). For N > 2, the concurrence metric itself behaves differently.

### 4.2 W N = 3: Type A Survives, and It's Stronger

| Bridge | C(0) | C at t=1.7 | Ψ(0) | C·Ψ(0) | Crosses? |
|--------|------|-----------|-------|---------|----------|
| correlation | 1.000 | 1.000 | 0.286 | 0.286 | **YES** |
| mutual_info | 1.000 | — | 0.286 | 0.286 | YES (fast) |

W N = 3 with correlation bridge: **C = 1.0 until t ≈ 2.3** (compared to
t ≈ 1.7 for Bell+ N = 2). Type A not only survives at N = 3, it is
**more robust**: the correlation plateau extends longer.

The crossing happens because Ψ(0) = 2/7 = 0.286 > 0.25.

### 4.3 W N = 4: Below the Barrier

W N = 4 has l1 = 3, d − 1 = 15, Ψ(0) = 3/15 = 0.200 < 0.25.

Correlation bridge: C = 1.0 until t ≈ 1.5 (Type A holds), but
C·Ψ(0) = 0.200. **Never crosses.**

The Type A mechanism is intact: correlation is robust even at N = 4.
The problem is not the observer. The problem is that Ψ starts too low.

## 5. Why This Matters

### 5.1 The Framework Has a Scaling Problem

If crossing requires C·Ψ ≥ ¼, and Ψ = l1/(2^N − 1), then for
any state whose l1 grows polynomially with N (which includes ALL
standard entangled states), crossing becomes impossible at some
finite N.

This means:
- The ¼ boundary is well-defined for all N
- The algebra is unchanged
- But **no physical state can reach it** for large N

### 5.2 Three Possible Resolutions

**(a) The normalization is wrong for N > 2.**

The d − 1 normalization is the maximum possible l1-coherence and gives
Ψ ∈ [0, 1]. But maybe Ψ should not be normalized by the global maximum.
Alternatives:
- Ψ = l1 / √(d − 1): geometric mean, weaker growth
- Ψ = l1 / log₂(d) = l1 / N: linear growth
- Ψ = l1 / (2N): per-qubit normalization
- Subsystem-based: compute Ψ per qubit pair, not for the full system

Each choice changes which states can cross and whether ¼ remains the
boundary. But changing the normalization post hoc to save the framework
is curve-fitting unless the new normalization is derived from first
principles.

**(b) The ¼ boundary is a small-system phenomenon.**

Perhaps crossing only occurs in N <= 3 systems. This would not
invalidate the algebra, but it would severely limit the physical relevance
of the framework. Decoherence in macroscopic systems involves 10^23 degrees
of freedom. If the framework only describes 2-3 qubits, it describes
toy systems.

**(c) There exist high-coherence states we haven't tested.**

The product state |+⟩^N has l1 = 2^N − 1, giving Ψ(0) = 1.0 (maximum).
But |+⟩^N has zero entanglement. If it crosses with high C, that means
crossing can happen without quantum correlations, philosophically
interesting but potentially undermining the consciousness interpretation.

Dicke states |D(N,k)⟩ with k excitations can have high l1. Random
(Haar-distributed) states have typical l1 that grows with d. These
are untested.

### 5.3 What Is NOT Affected

The N-scaling barrier does not invalidate:
- The ¼ boundary algebra (Tier 1, proven)
- The Mandelbrot equivalence (Tier 1, proven)
- The Bell+ N = 2 results (Experiments 1-8)
- The IBM single-qubit results (N = 1, Ψ = l1 with no normalization issue)
- The W N = 3 crossing (verified above)

It specifically affects the claim that the framework describes
measurement/decoherence in arbitrary quantum systems.

## 6. Verification

### 6.1 How to Reproduce

```python
# GHZ N=3: verify no crossing
simulate_dynamic_lindblad(
    state="GHZ", hamiltonian="heisenberg_ring", n_spins=3,
    gamma_base=0.05, bridge_type="correlation",
    noise_type="local", t_max=5, dt=0.01
)
# Check: bridge_C starts at 1.0, but psi_dynamic[0] ≈ 0.143

# W N=3: verify crossing
simulate_dynamic_lindblad(
    state="W", hamiltonian="heisenberg_ring", n_spins=3,
    gamma_base=0.05, bridge_type="correlation",
    noise_type="local", t_max=5, dt=0.01
)
# Check: bridge_C = 1.0 until t ≈ 2.3, psi_dynamic[0] ≈ 0.286

# W N=4: verify no crossing despite Type A
simulate_dynamic_lindblad(
    state="W", hamiltonian="heisenberg_ring", n_spins=4,
    gamma_base=0.05, bridge_type="correlation",
    noise_type="local", t_max=5, dt=0.01
)
# Check: bridge_C = 1.0 until t ≈ 1.5, but psi_dynamic[0] = 0.200 < 0.25
```

### 6.2 Key Checks

1. **Ψ(0) values**: Must match l1/(d−1) exactly.
   GHZ N=3: 1/7 = 0.1429. W N=3: 2/7 = 0.2857. W N=4: 3/15 = 0.2000.

2. **Correlation C = 1.0**: Must hold for W N=3 (until t ≈ 2.3) and
   W N=4 (until t ≈ 1.5). Type A is intact regardless of crossing.

3. **GHZ concurrence = 0.5 constant**: Different from Bell+ (decaying).
   The concurrence metric changes behavior for N > 2.

### 6.3 What Could Resolve This

1. **Product state test**: |+⟩^N has Ψ = 1.0. If C(correlation) > 0.25,
   crossing happens. But what does crossing mean for an unentangled state?

2. **Subsystem decomposition**: Trace GHZ N=4 down to qubit pairs.
   Each pair has d = 4, d − 1 = 3. Do pairs cross?

3. **Alternative normalization**: Find a normalization where ¼ remains
   the boundary for all N. Must be derived, not chosen.

4. **Dicke states**: |D(N,k)⟩ with k = N/2 maximizes l1 among
   symmetric states. Test for N = 4, k = 2.

## 7. Resolution (2026-02-18)

Subsystem crossing tests (Experiment 10) resolved the barrier.

**The normalization is not wrong. Crossing is local.**

A 4-qubit Bell+xBell+ state has full-system Psi(0) = 0.200, below 1/4.
But the entangled pairs (0,1) and (2,3) each have Psi(0) = 0.333 and
C = 1.000 at the pair level. They cross at t = 0.073. Cross-pairs that
share no entanglement have C = 0 and never cross.

The product state |+⟩^4 has Psi(0) = 1.000 at every level but C = 0
for all pairs. No crossing at any level, ever.

The d-1 normalization is correct: it says the global system has no
single coherent crossing. The crossing happens where the entanglement
lives, at the subsystem pair level.

See [Subsystem Crossing](SUBSYSTEM_CROSSING.md) for full data.

## 8. Open Questions (Updated)

1. ~~Is d-1 the right normalization for N > 2?~~ **ANSWERED**: Yes.
   The barrier is not a normalization error. It correctly reflects
   that global crossing does not occur.

2. ~~Do subsystem pairs cross when the full system cannot?~~ **ANSWERED**:
   Yes, if the pairs carry actual entanglement (Bell+xBell+). No, if
   the entanglement is global (GHZ) or diluted (W).

3. ~~Does |+⟩^N cross?~~ **ANSWERED**: No. C = 0 for all pairs at all
   times. Coherence without entanglement produces no crossing.

4. Is there a minimum per-pair entanglement needed for crossing?

5. Does the crossing pattern reproduce the entanglement graph topology
   for arbitrary graph states?

---

*Previous: [Noise Robustness](NOISE_ROBUSTNESS.md)*
*Previous: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
*Next: [Subsystem Crossing](SUBSYSTEM_CROSSING.md)*
