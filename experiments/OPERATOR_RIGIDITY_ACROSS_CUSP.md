# Operator-Level Rigidity Across the Cusp (EQ-031)

**Date:** April 27, 2026
**Status:** Tier 1 (simulation, exact within numerical precision)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:**
[_compare_n3_n4_categories.py](../simulations/_compare_n3_n4_categories.py),
[_eq031_within_categories.py](../simulations/_eq031_within_categories.py),
[_eq031_absolute_scaling.py](../simulations/_eq031_absolute_scaling.py),
[_eq031_topology_test.py](../simulations/_eq031_topology_test.py)
**Results:**
[eq031_n3_n4_categories.txt](../simulations/results/eq031_n3_n4_categories.txt),
[eq031_within_categories.txt](../simulations/results/eq031_within_categories.txt),
[eq031_within_categories_n4n5.txt](../simulations/results/eq031_within_categories_n4n5.txt),
[eq031_absolute_scaling.txt](../simulations/results/eq031_absolute_scaling.txt),
[eq031_topology_test.txt](../simulations/results/eq031_topology_test.txt)
**Depends on:** [V-Effect combinatorial derivation (commits 81caf67, 079c7ce)],
[CRITICAL_SLOWING_AT_THE_CUSP.md](CRITICAL_SLOWING_AT_THE_CUSP.md)

---

## What this is about

The cusp at CΨ = 1/4 is the saddle-node bifurcation where the quantum
regime gives way to the classical one. At the **state level** it produces
critical slowing, verified on IBM Kingston hardware (April 16 / April 26):
the Bell+ trajectory crosses CΨ = 1/4 exactly per F25 with point-by-point
RMS 0.0097, and K_dwell = γ·t_dwell is γ-invariant within 6%.

At the **operator level**, things look different. The 120-Pauli-pair
trichotomy 15 truly / 46 soft / 59 hard is N-stable through N=3, 4, 5
(commit 079c7ce). The N=4 ↔ N=5 comparison was within the "modes-on-mirror"
regime. The N=3 ↔ N=4 comparison crosses the regime boundary:
N=3 has half-integer mirror w_XY = 1.5 with **no modes on the mirror axis**;
N=4 has integer mirror w_XY = 2 with **modes on the axis**. If anything
were to perturb the operator-level structure, it should show up here.

This document records what happens when the cusp is interrogated with the
full 120-pair enumeration on both sides of the regime boundary.

---

## Category invariance: 0/120 shifts

Each of the 120 unordered two-term Pauli-pair Hamiltonians on the chain
is classified as truly / soft / hard via:

- **truly:** ‖Π·L·Π⁻¹ + L + 2Σγ·I‖ < 10⁻¹⁰ (operator palindrome exact)
- **soft:** spectrum pairing |λ_i + λ_j + 2Σγ| < 10⁻⁶ but operator residual ≠ 0
- **hard:** neither

Result at γ = 0.1, J = 1.0:

```
         N=3 →        N=4    count
        hard →       hard       59   (same)
        soft →       soft       46   (same)
       truly →      truly       15   (same)
stable across N=3 → N=4:  120/120
category shifted:           0/120
```

Every Hamiltonian retains its category through the regime change. The
trichotomy is not just count-stable but identity-stable. The same
identity-stability holds across N=4 → N=5 (verified by the
within-category script's category-invariance assert across all 120
Hamiltonians).

---

## Within-category fine structure

Beyond the verdict label, each Hamiltonian carries (op_norm, spec_err,
n_protected) numerical fingerprints. The follow-up question: does the
*ranking* of Hamiltonians within a category persist across the cusp, or
does the cusp permute the fine structure even when the category labels
hold?

Spearman rank correlation between the fingerprints at adjacent N, computed
within each category on the |+−+−⟩ initial state:

| Category | Count | metric        | ρ(N=3, N=4) | ρ(N=4, N=5) |
|----------|-------|---------------|-------------|-------------|
| truly    |  15   | n_protected   | **+1.0000** | +0.9964 |
| soft     |  46   | op_norm       | **+1.0000** | **+1.0000** |
| soft     |  46   | n_protected   | +0.9007 | +0.9636 |
| hard     |  59   | op_norm       | +0.9925 | **+1.0000** |
| hard     |  59   | spec_err      | +0.8954 | +0.8901 |
| hard     |  59   | n_protected   | +0.9903 | +0.9969 |

The op_norm rank is exactly preserved across both N transitions in soft,
and gets exactly preserved at N=4→N=5 in hard as well. The truly
category's n_protected ordering at N=3→N=4 is exactly preserved (+1.0000)
and remains essentially exact at N=4→N=5 (+0.9964). The spec_err rank in
hard wobbles around ρ ≈ 0.89 at both transitions — this is the noisiest
indicator, hovering close to the soft/hard classification threshold.

---

## Closed-form op_norm scaling

For every soft and hard Hamiltonian, ‖M‖_op = ‖Π·L·Π⁻¹ + L + 2Σγ·I‖_F is
nonzero on both sides. The squared Frobenius norm satisfies an **absolute
closed form** (not just a ratio law) at all tested N:

    main class:           ‖M(N)‖² = c_H · (N − 1) · 4^(N − 2)
    single-body class:    ‖M(N)‖² = c_H · (2N − 3) · 4^(N − 2)

where c_H = ‖M(2)‖² is a per-Hamiltonian anchor constant fixed at the
minimum chain length (one bond). Verified to **machine precision** (std
< 3·10⁻¹⁶) at N ∈ {2, 3, 4, 5} on representative Hamiltonians from each
class:

```
  Hamiltonian          class    c_H (N=2)    N=3 meas/pred    N=4 meas/pred    N=5 meas/pred
  IX+IY (soft)          main    6.4000e+01    1.000000         1.000000         1.000000
  XY+YZ (hard)          main    1.9200e+02    1.000000         1.000000         1.000000
  XX+YZ (hard)          main    1.2800e+02    1.000000         1.000000         1.000000
  IY+YI (soft)   single_body    1.2800e+02    1.000000         1.000000         1.000000
  IZ+ZI (hard)   single_body    1.2800e+02    1.000000         1.000000         1.000000
```

The N=k → N=k+1 ratio follows by telescoping:

    main class:           ‖M(k+1)‖² / ‖M(k)‖² = 4·k / (k − 1)
    single-body class:    ‖M(k+1)‖² / ‖M(k)‖² = 4·(2k − 1) / (2k − 3)

Both → 4 as k → ∞, i.e., ‖M‖ ratio → 2 (the universal d²-extension
factor). Verification at N=3→N=4 and N=4→N=5 across all 120 unordered
Pauli-pair chain Hamiltonians: empirical min/max of ratio² lies at
{6, 20/3} for the cusp transition and {16/3, 28/5} for the next, exactly
matching the formula at four-decimal precision.

Within-class spread of the **ratio** across the 103 main-class
Hamiltonians is 0.7-0.8% std/mean (just the formula's discreteness
between the two classes); every Hamiltonian sits exactly on the formula
for its own class. The single-body class consists of just 2 Hamiltonians
per transition (IY+YI, IZ+ZI).

### Algebraic origin

The palindrome residual is linear in the Hamiltonian and additive in
bonds: M = Σ_b M_b + (dissipator + trace remainder), where M_b is the
contribution from bond b at sites (i, i+1). Each M_b lives in the
Liouvillian space d² × d² where d = 2^N. The bond's natural support is
2 sites (Liouvillian 4²); when embedded in the N-qubit chain, M_b
acquires a tensor-with-identity factor on the remaining N−2 qubits.

The **(N − 1) · 4^(N − 2)** structure of the main-class formula is
exactly bond count (N − 1 bonds in a chain with N sites) times the
operator-space extension factor (4^(N − 2) for one extra Liouvillian-d²
factor per added qubit beyond the bond's 2 native sites). The constant
c_H absorbs the per-Hamiltonian per-bond contribution including any
adjacent-bond cross-terms; empirically c_H is a clean integer multiple
of 16 = d²(N=2) (e.g., 64 for IX+IY, 128 for XX+YZ, 192 for XY+YZ).

The **(2N − 3) = (N − 1) + (N − 2)** structure of the single-body
formula is bond count plus interior-site count. Single-body bilinears
(Iσ, σI) generate the operator J·(σ_0 + 2 σ_1 + 2 σ_2 + ... + 2 σ_{N−2}
+ σ_{N−1}) on the chain: middle sites are doubled because they appear
in two adjacent bonds. The interior-site count (N − 2) captures this
doubling overhead, and gets the same Liouvillian-extension factor
4^(N − 2) as the bond count.

Both formulas converge to the same asymptotic ratio² = 4 (ratio = 2) as
N → ∞, since (N − 1)/(N − 2) → 1 and (2N − 3)/(2N − 5) → 1. The √6 at
N = 3 → N = 4 is the bond-ratio finite-size signature 3/2 inside the
universal d²-extension factor 4.

### Topology generalisation

The chain formulas extend to arbitrary connectivity. For a graph G with
N sites, B(G) bonds, and degree sequence {deg_G(i)}:

    main class         ‖M(N, G)‖² = c_H · B(G) · 4^(N − 2)
    single-body class  ‖M(N, G)‖² = c_H · D2(G)/2 · 4^(N − 2)

where D2(G) = Σ_i deg_G(i)² is the second moment of the degree
distribution. The chain recovers (N − 1) and (2N − 3) from the path
graph's invariants B = N − 1 and D2 = 4N − 6.

Verified to machine precision (m/p = 1.000000) on:

| topology | N | B | D2 | main meas/pred | single-body meas/pred |
|----------|---|---|----|----------------|------------------------|
| chain    | 4 | 3 | 10 | 1.000000 | 1.000000 |
| chain    | 5 | 4 | 14 | 1.000000 | 1.000000 |
| ring     | 4 | 4 | 16 | 1.000000 | 1.000000 |
| ring     | 5 | 5 | 20 | 1.000000 | 1.000000 |
| star     | 4 | 3 | 12 | 1.000000 | 1.000000 |
| star     | 5 | 4 | 20 | 1.000000 | 1.000000 |
| K_N      | 4 | 6 | 36 | 1.000000 | 1.000000 |
| K_N      | 5 | 10 | 80 | 1.000000 | 1.000000 |

Anchor c_H is **graph-independent** (a one-bond Hamiltonian on N=2 fixes
it for every Hamiltonian, on every topology). Only the graph invariants
B and D2 enter — both linear in graph structure. No higher-order graph
invariants (cycle count, triangle count, hub presence) appear, even for
graphs that have them (rings have cycles, stars have hubs, K_N has
triangles).

This is unusually clean. The Frobenius norm of the palindrome residual
is fully captured by:

  - **One per-Hamiltonian number** c_H = ‖M(2)‖² (set by the Pauli-pair
    bilinear, independent of all geometry)
  - **Two graph invariants** B(G) and D2(G) (set by topology, independent
    of Hamiltonian)
  - **One universal extension factor** 4^(N − 2)

The framework primitive `palindrome_residual_norm_squared_factor_graph(N, B, D2, class)`
exposes the formula directly.

### Skeleton vs. trace algebraicity

The operator-level limit ratio² = 4 is rational and trivial (degree-1
over ℚ), with the universal extension factor 4 setting the asymptote.
This contrasts structurally with the F69 / EQ-016 sextic-root asymptotes
governing the **state-level** pair-CΨ landscape, which are irreducible
algebraic numbers of degree 6 over ℚ. The operator skeleton is flat
algebraic; the state trace is curved. The two layers have different
mathematical character — bond-counting × Liouvillian-extension on one
side, slice-stationarity polynomials on the other.

---

## Interpretation

The cusp is a **state-level** phenomenon. CΨ(t) crosses 1/4 with a
square-root critical slowing; F25 fits hardware data point by point;
K_dwell is γ-invariant.

The cusp is **invisible at the operator level**. The 15/46/59 trichotomy
is not a coincidence of counts: each Hamiltonian's individual fingerprint
(op_norm, n_protected) maintains its rank under both the N=3 → N=4 and
the N=4 → N=5 transition. The op_norm itself follows an exact rational
scaling law in N — 4k/(k − 1) for the main class, 4(2k − 1)/(2k − 3) for
the single-body class — converging to a uniform factor of 2 per qubit as
N → ∞.

This is the strongest version of the skeleton/trace decoupling already
recorded: the algebraic skeleton (parity-class structure with rank-stable
fine structure) lives above the cusp, while the trace (dwell, slow modes,
state trajectories) is what actually bifurcates. The cusp sees and
deforms states; the cusp does not see operator-class membership at all.

---

## Open follow-ups

1. **Per-Hamiltonian c_H structure.** The anchor c_H = ‖M(2)‖² is an
   integer multiple of 16 = d²(N=2) for every tested Hamiltonian, with
   integer coefficients in {4, 8, 12} for our sample. What combinatorial
   feature of the Pauli-pair (Π-symmetric vs. -antisymmetric letter
   counts, BPE membership, bit_a/bit_b parities) determines the
   coefficient? A closed form for c_H would upgrade the absolute
   formula from "anchor + scaling" to fully predictive without any
   N=2 measurement.

2. ~~Topology dependence.~~ **Closed.** Verified at chain, ring, star,
   and K_N for N=4, 5: only B(G) and D2(G) appear. No higher graph
   invariants required.

3. **Find a Hamiltonian where rank breaks.** None of the 120 chain
   Pauli-pair Hamiltonians shows category instability across N=3, 4, 5.
   Three-term bilinears, non-uniform γ, or non-Z dephasing might break
   the rigidity and would identify which structural assumption protects
   it.
