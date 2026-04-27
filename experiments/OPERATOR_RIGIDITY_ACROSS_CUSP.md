# Operator-Level Rigidity Across the Cusp (EQ-031)

**Date:** April 27, 2026
**Status:** Tier 1 (simulation, exact within numerical precision)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:**
[_compare_n3_n4_categories.py](../simulations/_compare_n3_n4_categories.py),
[_eq031_within_categories.py](../simulations/_eq031_within_categories.py)
**Results:**
[eq031_n3_n4_categories.txt](../simulations/results/eq031_n3_n4_categories.txt),
[eq031_within_categories.txt](../simulations/results/eq031_within_categories.txt),
[eq031_within_categories_n4n5.txt](../simulations/results/eq031_within_categories_n4n5.txt)
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
nonzero on both sides. The N=k → N=k+1 amplification factor follows two
clean rational laws:

    main class:           ‖M(k+1)‖² / ‖M(k)‖² = 4·k / (k − 1)
    single-body class:    ‖M(k+1)‖² / ‖M(k)‖² = 4·(2k − 1) / (2k − 3)

Both → 4 as k → ∞, i.e., ‖M‖ ratio → 2.

Verification:

| transition | main ratio² (formula) | main ratio² (measured) | single-body ratio² (formula) | single-body ratio² (measured) |
|------------|-----------------------|------------------------|------------------------------|-------------------------------|
| N=3 → N=4  | 4·3/2 = **6** | 6.000 | 4·5/3 = **20/3** ≈ 6.667 | 6.667 |
| N=4 → N=5  | 4·4/3 = **16/3** ≈ 5.333 | 5.333 | 4·7/5 = **28/5** = 5.600 | 5.598 |

Within-class spread across the 103 main-class Hamiltonians is 0.7-0.8%
std/mean, i.e., effectively zero — every Hamiltonian in the main class
sits exactly on the formula. The single-body class consists of only 2
Hamiltonians per transition (IY+YI, IZ+ZI), which sit exactly on the
single-body formula.

### Algebraic origin

In Pauli basis, Π acts diagonally with sign (−1)^{nXY(α)} on the Pauli
string σ_α (counting X and Y letters). The palindrome residual is
linear in the Hamiltonian: M = Π·L·Π⁻¹ + L + 2Σγ·I splits as a sum of
per-bond contributions plus a dissipator-correction piece.

For a chain with N sites and N−1 bonds, the **main class** consists of
Hamiltonians whose bilinear bond term has at least one non-trivial 2-body
component (e.g., XY+YZ, IX+ZZ). For these, the residual is dominated by
the bond commutator, and:

    ‖M(N)‖² ∝ (N − 1) · 4^(N − 2)

The (N − 1) is the bond count; the 4^(N − 2) is the Liouvillian
extension factor — adding one qubit beyond the bond's native 2 amplifies
the Frobenius norm by sqrt(4) = 2 per qubit. The ratio between adjacent
N is therefore 4·(N − 1)/(N − 2), evaluated at source N = k:
4·k/(k − 1).

The **single-body class** consists of Pauli pairs of the form (Iσ, σI)
with σ ∈ {X, Y, Z}. The bond bilinear σ·I + I·σ on each bond reduces to a
sum of single-site Pauli terms — there is no genuine 2-body interaction.
For X and Y, ΠσΠ⁻¹ = −σ on each site, so the H-commutator part of M
cancels (i[H, ·] − i[H, ·] = 0) and only the dissipator + trace remainder
contributes. This produces a different N-dependence yielding
4·(2k − 1)/(2k − 3). The ratio is slightly larger than the main-class
formula and converges to it from above as k → ∞.

The IZ+ZI case lands in **hard**, not soft, because Π Z Π⁻¹ = +Z (Z
commutes with the chiral conjugation): the H-commutator does not cancel,
but the resulting residual still respects the same single-body algebraic
structure that gives ratio² = 4(2k − 1)/(2k − 3).

### Asymptotic behaviour

Both classes converge to ratio² = 4 (ratio = 2) as N → ∞. This is the
pure d²-extension factor: an additional qubit doubles the operator-space
dimension. The "structural" content — bond count vs. dissipator-only
piece — washes out at large N. The √6 at N=3→N=4 is a finite-size
signature of the bond ratio 3/2 inside the universal extension factor 4.

This is structurally different from the F69 / EQ-016 sextic-root
asymptotes that govern the **state-level** pair-CΨ landscape. The
state-level limit is irrational (degree-6 over ℚ); the operator-level
limit is the trivial integer 4 (degree-1 over ℚ). The skeleton is flat
algebraic; the trace is curved.

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

1. **Rigorous proof of the closed-form scaling.** The rational laws
   4k/(k − 1) and 4(2k − 1)/(2k − 3) are heuristically derived from
   bond counting × Liouvillian extension factor, and verified at
   N=3→N=4 and N=4→N=5 to six significant figures. A first-principles
   derivation that produces both formulas from the Pauli-basis
   structure of M would close this from "verified" to "theorem".

2. **Topology dependence.** The chain has a special bond structure
   (open boundary, sequential connectivity). On a ring, complete graph,
   or star, the bond count differs and the formula should change.
   Mapping how the rational law generalises to other topologies would
   identify which factors are topological vs. universal.

3. **Find a Hamiltonian where rank breaks.** None of the 120 chain
   Pauli-pair Hamiltonians shows category instability across N=3, 4, 5.
   Three-term bilinears, non-uniform γ, or non-Z dephasing might break
   the rigidity and would identify which structural assumption protects
   it.
