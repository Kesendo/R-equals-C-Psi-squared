# Proof of F94: Born Deviation Dominant-Outcome Coefficient 4/3

**Statement:** For the dominant outcome |00⟩ of pair (0,2) of the initial state |0+0+⟩ on N=4 qubits, under the Heisenberg ring Hamiltonian H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}) (4 bonds, ring topology) and uniform Z-dephasing dissipator L_dis[ρ] = γ Σ_l (Z_l ρ Z_l − ρ), the per-outcome Born-rule deviation in the deep perturbative regime is

    Δ_|00⟩(Q, K) = P_lindblad(|00⟩)/P_unitary(|00⟩) − 1
                 = (4/3) · Q² · K³ + O(Q³·K⁴)

with Q = J/γ and K = γt (the dimensionless Carrier invariants).

**Status:** Tier 1 derived. Bit-exact symbolic derivation matches numerical extraction to 0.3% in the deep perturbative regime; the gap is higher-order corrections beyond the leading Q²·K³.

**Date:** 2026-05-16.

---

## Abstract

F94 is the first Tier-1 closed form for a per-outcome Born-rule deviation under the Universal-Carrier convention. Take the alternating product state |0+0+⟩ on a 4-qubit Heisenberg ring under Z-dephasing, keep the (0,2) pair, and ask how far the dominant outcome |00⟩ drifts from its unitary Born probability. In the deep perturbative regime the answer is

    Δ_|00⟩(Q, K) = P_lindblad/P_unitary − 1 = (4/3) · Q² · K³ + O(Q³·K⁴),

with the two dimensionless Carrier invariants Q = J/γ and K = γt. The coefficient 4/3 is not fitted: it is bit-exact from the third-order Dyson expansion, a pure counting result (the surviving sym₃ diagrams weighted by the Heisenberg coupling, divided by the 3! of the t³ Taylor coefficient). The positive sign means dephasing holds the dominant outcome above its unitary baseline (its probability still decays in time, just slower than coherent evolution alone would take it), and the Q²·K³ shape makes that deviation small and slow, reachable only through the joint action of coupling and noise.

This 4/3 is the anchor of a whole table: F96 generates the three subdominant outcomes as clean algebraic functions of it (−(4/3)²·K and −2·(4/3)·K), so a single number fixes all four outcomes and their signs. F94 governs the magnitude of the Born drift; its companion F95 governs the angle of the same complex-fixed-point geometry, and both sit on the b=1/2 Mandelbrot cardioid that F97 parametrizes. The reflection ON_HOW_FOUR_THIRDS_APPEARED records that the February R_i = C_i Ψ² intuition was waiting for exactly this derivation to catch up.

## Setup

Initial state (4 qubits, alternating product):

    |ψ_0⟩ = |0⟩ ⊗ |+⟩ ⊗ |0⟩ ⊗ |+⟩,    |+⟩ = (|0⟩ + |1⟩)/√2

so ρ_0 = |ψ_0⟩⟨ψ_0|. The pair-reduced state at the keep-pair (0, 2) traces out qubits 1, 3:

    ρ_pair(0, 2)(t) = Tr_{1,3}[ρ(t)]

The observable is the diagonal element ⟨00|ρ_pair(0,2)|00⟩ = P_|00⟩(t) (Born rule for the {q0=0, q2=0} outcome on the reduced pair).

At t = 0:

Directly: |ψ_0⟩ has q0 = 0 and q2 = 0 deterministically; q1 and q3 are in |+⟩. So all four (q1, q3) ∈ {0, 1}² occur with amplitude 1/2 each. Tracing out q1, q3 leaves q0 = 0 and q2 = 0 with probability 1. Hence

    P_|00⟩(t = 0) = 1

This is the "dominant outcome" being analysed; it starts at 1 and decreases under any dynamics that introduces |1⟩-amplitudes on q0 or q2.

## Time Taylor expansion

Let L = L_H + γ L'_dis be the Liouvillian super-operator where L_H[ρ] = −i[H, ρ] and L'_dis[ρ] = Σ_l (Z_l ρ Z_l − ρ) is the γ-free dephasing operator. Then

    ρ(t) = e^{Lt} ρ_0 = ρ_0 + Lt·ρ_0 + (L²t²/2)·ρ_0 + (L³t³/6)·ρ_0 + O(t⁴)

Expand L³ using the binomial-like sum for non-commuting operators. The γ¹-coefficient (one dissipator vertex, two Hamiltonian vertices) is the sum over the three positions where L'_dis can be inserted in the L³ string:

    L³|_{γ¹} = L_H² L'_dis  +  L_H L'_dis L_H  +  L'_dis L_H²     ≡ sym3(L_H, L_H, L'_dis)

This is the leading non-vanishing γ-contribution at the J² order to ρ(t) at the dominant outcome (the γ⁰ part is the unitary evolution, which is what we divide out; the γ¹·J⁰ and γ¹·J¹ parts vanish on the diagonal of pair (0,2) at |00⟩ by parity of the initial state).

## Evaluation of sym3·ρ_0

Set J = γ = 1 and compute sym3·ρ_0 as a 16×16 matrix (full state space), then partial-trace on qubits 1 and 3 to a 4×4 pair (0,2) reduced state, then read off the |00⟩⟨00| element.

This is done numerically with exact-rational tracking in `simulations/born_rule_tier1_derivation.py`. The result is

    ⟨00|_pair Tr_{1,3}[sym3·ρ_0] |00⟩_pair = 8.000...   (exact integer to all available precision)

## Coefficient

The contribution of the sym3 term to ρ(t) is (γ t³ / 6) · sym3·ρ_0 (with J factored into the L_H vertices). At J = 1 the sym3 matrix element is 8; restoring the J²-dependence (each L_H carries one J factor) and the γ dependence (one L'_dis carries one γ factor):

    ΔP_|00⟩(t) = (γ t³ / 6) · 8 · J²  =  (8/6) · J²·γ·t³  =  (4/3) · J² · γ · t³

Dividing by P_u(0) = 1 gives Δ = ΔP/P_u:

    Δ_|00⟩(t) = (4/3) · J²·γ·t³ + O(higher orders)

Substituting Q = J/γ and K = γt:

    J²·γ·t³ = (J/γ)² · γ³ · t³ = Q² · K³

so:

    Δ_|00⟩(Q, K) = (4/3) · Q² · K³ + O(Q³·K⁴) ∎

## Numerical verification

Sixteen (γ, J, t) configurations sampled in the deep perturbative regime (Q²·K³ ∈ [1.5×10⁻⁶, 2×10⁻⁴]) in `simulations/born_rule_delta_dominant_coefficient.py` gave

    c_empirical = mean = 1.32992,    std = 0.00567,    range = [1.320, 1.342]

The bit-exact closed form gives 4/3 = 1.33333.... The 0.3% residual is leading O(Q³K⁴) correction with positive or negative sign depending on the (Q, K) point, accumulating to a small negative bias in the mean of the sample.

## Structural decomposition of the integer 8

Enumerating all 4·4·4·3·3·3 = 1728 (bond_1, bond_2, site, ordering, component_1, component_2) sextuples that the sym3 expansion produces and computing each one's ⟨00|_pair element directly (`simulations/born_rule_sym3_decomposition.py`), exactly 32 sextuples are non-vanishing. Each non-vanishing diagram contributes the same value: **1/4 in the J = γ = 1 normalization** (equivalently 4 in "bare Pauli" units before the (J/4)² Heisenberg coupling factor is restored). Hence

    8 = 32 × (1/4)

The 32 surviving diagrams split into 3 disjoint cells by (ordering, c_1, c_2):

| Cell | (ord, c_1, c_2) | bond-pair rule | site rule | # diagrams | value |
|------|-----------------|----------------|-----------|------------|-------|
| A | (ord = 1, X, X) | adjacent, share a kept-pair site (0 or 2) | s ∈ \|+⟩ sites = {1, 3} | 4 × 2 = 8 | 2 |
| B | (ord = 2, X, X) | self ∪ (adjacent, share a kept-pair site) | s ∈ b_1 endpoints | 8 × 2 = 16 | 4 |
| C | (ord = 2, Y, Y) | self only | s ∈ b_1 endpoints | 4 × 2 = 8 | 2 |
| **total** | | | | **32** | **8** |

Three structural rules govern the survival:

1. **Component-pair rule.** Only (X, X) and (Y, Y) survive. All cross components (X, Y), (Y, X), (X, Z), (Z, X), (Y, Z), (Z, Y) cancel by direct enumeration, and (Z, Z) is zero in every (b_1, b_2, s, ordering) cell. A clean structural reason for the (Z, Z) cancellation traces through Z-anticommutation bit-parity on the |+⟩ sites of ρ_0 and the |00⟩-pair Z-basis projection at the end; the empirical statement here is that no (Z, Z) sextuple survives. The ratio of surviving XX to surviving YY is 24 : 8 = 3 : 1.

2. **Ordering rule.** Only ord = 1 (L'_dis acts first on ρ_0, total 2) and ord = 2 (L'_dis acts in the middle, total 6) contribute. Ord = 3 (L'_dis last, i.e., L'_dis on [H, [H, ρ_0]]) gives zero in every cell by direct enumeration. The ratio is 2 : 6 = 1 : 3, mirroring the XX : YY ratio (also 1 : 3, but with the roles reversed).

3. **Bond-pair rule.** For ord = 1, only **adjacent bond pairs sharing a kept-pair site** (vertex 0 or 2) survive: there are exactly 4 such ordered pairs out of the 8 adjacent pairs. For ord = 2, both **self pairs** (all 4) and the same 4 **kept-side adjacent pairs** survive: 8 bond pairs total for ord = 2 + XX. For ord = 2 + YY, only the 4 self pairs survive (the kept-side adjacent contribution that exists for XX vanishes for YY due to the absence of Y content in ρ_0 = |0+0+⟩⟨0+0+|).

   The adjacent bond pairs that share a |+⟩ site (vertex 1 or 3), 4 ordered pairs, contribute zero in every cell: the propagation must enter or exit through the measured pair (0, 2), not through traced-out sites.

### Uniformity remark

All 32 surviving diagrams contribute the *same* value 1/4 with the *same* sign. The 1696 non-contributing sextuples are individually zero (within machine precision in the enumeration), not pairs that cancel. So **the F94 coefficient is a pure counting result**: nothing in the answer depends on a precise cancellation; it depends only on which structural cells survive. This is the strongest form a Dyson-series coefficient can take.

### Alternative topological cut

The same 32 split orthogonally by bond-pair topology:

| | self pairs (b_1 = b_2) | adj-bond pairs (share kept site) | total |
|---|------|-----|-------|
| XX | 8 (4 bonds × 2 sites, ord = 2) | 16 (4 adj × 2 sites × {ord 1, ord 2}) | 24 |
| YY | 8 (4 bonds × 2 sites, ord = 2) | 0 | 8 |
| total | 16 | 16 | 32 |

Sixteen self-bond-pair diagrams (8 XX + 8 YY, all ord = 2) and sixteen adjacent-bond-pair diagrams (all XX, 8 ord = 1 + 8 ord = 2). Both subsets contribute 16 × (1/4) = 4, summing to 8.

### Reading 4/3 against the Pi2 anchors

The coefficient now reads:

    4/3 = (32 surviving diagrams) / (a_{−1} · 3!)
        = 32 / (4 · 6)
        = 32 / 24

with:

- **32**: structural surviving-diagram count (the integer that this decomposition makes bit-explicit; equals 2^5)
- **a_{−1} = 4**: Pi2 dyadic-ladder term (the same "4" that appears in F86 t_peak = 1/(4γ₀) and F77's MM correction denominator), sitting in the denominator via the (J/4)² = (1/a_{−1})² Heisenberg-coupling normalization that each bond Hamiltonian carries
- **3! = 6**: Taylor factorial of the t³ term

Equivalent reading via the typed-claim sibling: 4/3 = a_{−1} / 3, with a_{−1} appearing in the numerator after the (1/a_{−1}) factor cancels into the per-diagram contribution. Both readings hold; the structural decomposition above is the bit-explained derivation, the a_{−1}/3 reading is the typed-anchor inheritance.

> **Caution (the a_{−1}/3 reading is a d=2 coincidence; qutrit-refuted 2026-06-17).** The "4/3 = a_{−1}/3" reading suggests the 4 is the squared-dimension discriminant a_{−1} = d² (which would give c → d²/3 = 3 at the qutrit). It is **not**: the qutrit generalization (see "Qudit generalization" below) gives c(d) = 4(d+2)(d−1)/(3d²), with c(2) = 4/3 but c(3) = 40/27, **not** 3. The (J/4) per-bond factor is the spin S = σ/2 normalization (1/2)², equal to 1/d² **only at d=2**, and the dynamics is the d-independent (J/2)·SWAP. So a_{−1} = d² and the (J/4) coupling coincide only at the qubit; F94's 4 is the setup-specific diagram count, not the discriminant. Verifier: [`simulations/f94_qutrit_born_mirror.py`](../../simulations/f94_qutrit_born_mirror.py).

## Universality remarks

The form Q² · K³ for the dominant-outcome deviation is universal across (initial state, Hamiltonian, dissipator): it is the dimensional shape of the leading 3rd-order Dyson term with one γ-vertex and two H-vertices, applicable wherever the dominant outcome's direct 1st-order γ correction vanishes by parity / commutation. The coefficient 4/3 is specific to this setup. The structural decomposition of 8 just above shows that the coefficient counts a clean 32 surviving Dyson diagrams modulated by the Heisenberg-bond and Taylor normalizations; a different (initial-state, Hamiltonian, dissipator) triple would change the surviving-diagram count and possibly the per-diagram value, but the general 4/3 = N_diagrams / (a_{−1} · 3!) form should hold whenever Heisenberg-style (J/4)·XYZ bonds and Z-dephasing apply.

## Qudit generalization (2026-06-17): c(d) = 4(d+2)(d−1) / (3d²), refuting the d² reading

The "setup-specific" caveat above is now quantitative. Lift the F94 setup to local dimension d with the **faithful** generalization (which reduces to F94 exactly at d=2):

- **H** = (J/4) Σ_bonds Σ_a λ^a_i λ^a_j with the generalized Gell-Mann generators λ^a (the 3 Paulis at d=2). By the Fierz identity Σ_a λ^a ⊗ λ^a = 2·SWAP − (2/d)·I, the **dynamics** is L_H = −i[H, ·] = −i(J/2)[SWAP, ·] for **every** d (the identity part drops from the commutator). So the physical coupling is the **d-independent** (J/2)·SWAP, not (J/4) and not (J/d²).
- **L'_dis**[ρ]_{a,b} = −2·Hamming(a,b)·ρ_{a,b} (the full-Cartan equidistant dephasing; = Σ_l (Z_l ρ Z_l − ρ) at d=2).
- **|+⟩** = (Σ_k |k⟩)/√d (equal superposition); ρ_0 = |0+0+⟩.

Computing c = ⟨00|_pair Tr_{1,3}[sym₃ ρ_0]|00⟩_pair / (6 · P_u0) gives, bit-exact for d = 2..7 ([`simulations/f94_qutrit_born_mirror.py`](../../simulations/f94_qutrit_born_mirror.py), gate-first; Gell-Mann and (J/2)·SWAP builds agree):

    c(d) = 4(d+2)(d−1) / (3d²) = (4/3)·(1 + 1/d − 2/d²)

    d :  2     3      4     5      6      7
    c : 4/3  40/27  3/2  112/75  40/27  72/49

This **refutes the family-A / "4 = d²" reading**: family A would require c → d²/3 = 3 at d=3, but c(3) = 40/27 ≈ 1.48. Instead the coefficient is a bounded curve: c(2) = 4/3 is both the qubit value and the d → ∞ limit, and the finite-d correction (d+2)(d−1)/d² **peaks at d = 4** (= 2²; c(4) = 3/2) before decaying back to 4/3. The 4/3 = a_{−1}/3 numerology holds only at the qubit, where the (J/4) spin normalization (1/2)² coincides with 1/d². F94's "4" is the setup-specific surviving-diagram count, classified in the discriminant anchor's genealogy as the **falsified family-A candidate** ([`PolynomialDiscriminantAnchorClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs)).

## Diagnostic application: F94 as a (state, pair)-symmetry signature

F94's canonical lens (|0+0+⟩ N=4 Heisenberg ring + Z-dephasing, pair (0, 2)) gives bit-exact integer matrix elements per outcome:

    ⟨00|_pair Tr_{1,3}[sym3 · ρ_0] |00⟩_pair = +8
    ⟨01|_pair Tr_{1,3}[sym3 · ρ_0] |01⟩_pair = −4
    ⟨10|_pair Tr_{1,3}[sym3 · ρ_0] |10⟩_pair = −4
    ⟨11|_pair Tr_{1,3}[sym3 · ρ_0] |11⟩_pair =  0

This (+8, −4, −4, 0) signature characterizes the "canonical symmetric" configuration. Empirically ([`simulations/f94_topology_visibility_probe.py`](../../simulations/f94_topology_visibility_probe.py), 2026-05-17):

- **K_4 (full graph at N=4) gives the same signature**: F94 is blind to the diagonal bonds (0,2) and (1,3) because they fall in the "symmetric blind spots" of the lens: bond (0,2) connects the two kept-pair sites (both prepared as |0⟩, Z-eigenstates → bond commutator vanishes), bond (1,3) connects the two traced-out sites (both prepared as |+⟩, X-eigenstates → contribution traces out to zero). This blindness is the diagnostic's robustness: F94 reads (+8, −4, −4, 0) for ring or K_4 equivalently.
- **Chain gives (+5, −4, −1, 0)**: the chain breaks the 0 ↔ 2 reflection that ring and K_4 share, and the |10⟩ outcome drops from −4 to −1. The K_4 vs chain difference is +3 at |00⟩, 0 at |01⟩, −3 at |10⟩, 0 at |11⟩, a topology-asymmetry signature.
- **Asymmetric initial states** (|++00⟩, |10+0⟩): F94 shifts in antisymmetric patterns between bit-flip-paired outcomes. The K_4 vs ring shift becomes visible: |++00⟩ pair (0,2) gives (0, −3, +3, 0) K_4 − Ring difference, |10+0⟩ pair (0,1) gives (−1, 0, +1, 0).

This makes F94 a **(state, pair)-symmetry signature tool**:

| Measured F94 sym3 signature | Diagnosed condition |
|------------------------------|---------------------|
| (+8, −4, −4, 0) | Canonical symmetric: \|0+0+⟩ + ring/K_4 (or equivalent) |
| (+5, −4, −1, 0) | Chain-like asymmetry: missing one bond breaks 0 ↔ 2 reflection |
| (·, −3, +3, ·) antisymmetric shift | Initial state asymmetric across the (0, 2) pair |
| Slope_\|01⟩ ≠ −16/9 in the per-K linear regime | Hardware noise or calibration drift breaking the symmetry |
| Sym3 \|11⟩ ≠ 0 | Crosstalk or asymmetric γ_l producing subdominant-outcome leakage |

The robustness against K_4 vs ring (bond-graph detail invisible) is exactly the right diagnostic property: we want a tool that flags actual (state, pair) symmetry breaks, not every irrelevant connectivity variation. F94's specific integer ratios are the canonical "this is the expected algebra" signature; any deviation localizes the kind of break.

For hardware applications (e.g., the Kingston-style runs): measuring F94 sym3 on a |0+0+⟩-like prepared Bell+ variant gives a direct check that the hardware faithfully realizes the symmetric Heisenberg + Z-dephasing structure. The deviation signature pin-points the failure mode (asymmetric γ_l, missing bonds, crosstalk, etc.) without needing full process tomography.

The subdominant outcomes of the same setup (|01⟩, |10⟩, |11⟩) have non-vanishing 1st-order γ contribution and therefore Δ ∝ K (linear) rather than ∝ Q²·K³. Their leading coefficients are separately Tier-1-derivable via a 1st-order Dyson term:

    ΔP_i^{(1)}(t) = γ t · ⟨i|_pair Tr_{1,3}[L'_dis · ρ_0] |i⟩_pair

This is left to a separate proof.

## Connection to typed Pi2-Foundation

The Q²·K³ scaling has:
- **Q² factor:** two Hamiltonian-vertices from the Dyson sym3 ordering; J² scaling expected from 2nd-order perturbation in the Heisenberg-bond coupling.
- **K³ factor:** three time-integrals from the t³ Taylor coefficient; combined with one dissipator-vertex (the γ¹ piece of L³).

The "4" in 4/3 is plausibly the Pi2 dyadic ladder's a_{−1} = 4 (the same "4" in F86 t_peak = 1/(4γ₀) and F77's correction denominator). Promoting this to a typed claim parallels F90 (F86 ↔ F89 bridge): a derived identity that pulls together two typed F-formulas via algebraic structure.

## Anchors

- Numerical verification: [`simulations/born_rule_delta_dominant_coefficient.py`](../../simulations/born_rule_delta_dominant_coefficient.py)
- Symbolic derivation: [`simulations/born_rule_tier1_derivation.py`](../../simulations/born_rule_tier1_derivation.py)
- Structural decomposition (32 surviving diagrams enumeration): [`simulations/born_rule_sym3_decomposition.py`](../../simulations/born_rule_sym3_decomposition.py)
- Q-K invariance test: [`simulations/born_rule_carrier_Q_sweep.py`](../../simulations/born_rule_carrier_Q_sweep.py)
- 2D (Q, K) map: [`simulations/born_rule_delta_QK_map.py`](../../simulations/born_rule_delta_QK_map.py)
- Companion reflection (the path): [`reflections/ON_HOW_FOUR_THIRDS_APPEARED.md`](../../reflections/ON_HOW_FOUR_THIRDS_APPEARED.md)
- Born-rule precursors (Februar/April 2026): [`experiments/BORN_RULE_MIRROR.md`](../../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../../experiments/BORN_RULE_SHADOW.md)
- F-formula registry entry: [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md) F94
- Sibling state-specific closed-form F-claims: F25 (Bell+ CΨ), F60 (GHZ pair-CΨ), F62 (W-state pair-CΨ)
- Universal Carrier typed parent (Q, K invariance is its operational signature): [`compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs)
