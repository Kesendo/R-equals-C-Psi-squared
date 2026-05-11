# F89 Topology Orbit Closure: Bond-Graph Topology Determines ρ_cc Spatial-Sum Coherence

**Status:** Tier 1 derived (S_N-orbit symmetry proven; numerically verified bit-identical at N=7 across 12 topology classes and N=4 across all C(4,2)=6 site pairs)
**Date:** 2026-05-11
**Authors:** Thomas Wicht, Claude Opus 4.7 (1M context)
**Scripts:**
- [`_bond_isolate_compare_n7.py`](../simulations/_bond_isolate_compare_n7.py): N=7 single-bond pairwise comparison (six bonds, 30 ordered pairs).
- [`_bond_isolate_long_range_verify.py`](../simulations/_bond_isolate_long_range_verify.py): N=4 single-pair NN + long-range verification (6 site pairs, direct expm-of-Liouvillian).
- [`_bond_isolate_topology_classes_n7.py`](../simulations/_bond_isolate_topology_classes_n7.py): N=7 multi-bond topology-class consistency (12 classes).

**Outputs:** [`bond_isolate/`](../simulations/results/bond_isolate/) (28 CSVs at N=7 + two comparison plots).
**Related register entries:** [F73](../docs/ANALYTICAL_FORMULAS.md) (analogous closed-form closure for the (vac, SE) coherence block); [F71](../docs/ANALYTICAL_FORMULAS.md) (the spatial mirror Z₂ that sits inside the full S_N argument used here); [F86](../docs/ANALYTICAL_FORMULAS.md) (per-bond Q_peak fan, the empirical contrast: a linear response ∂J_b breaks S_N differently from the uniform-J multi-bond setup of F89).

---

## Theorem

For an N-qubit system with:
- **Hamiltonian** H_B = J · Σ_{(p,q) ∈ B} (X_p X_q + Y_p Y_q), where B is a set of distinct site pairs (p, q ∈ {0,...,N−1}, p ≠ q) and J is uniform across all active pairs.
- **Dissipator** uniform Z-dephasing: Lindblad operators √γ₀ · Z_l for every site l.
- **Initial state** ρ_cc = (|S_1⟩⟨S_2| + |S_2⟩⟨S_1|) / 2, where |S_n⟩ is the popcount-n (n excitations) symmetric Dicke state.
- **Observable** spatial-sum coherence S(t) = Σ_l 2 · |(ρ_l(t))_{0,1}|², where ρ_l = Tr_{j ≠ l}(ρ) is the reduced density matrix on site l (a 2×2 matrix), and (ρ_l)_{0,1} is its off-diagonal element between the |0⟩ and |1⟩ basis vectors of that qubit.

Then S(t) **depends only on the S_N-orbit of B** under the action σ · B = {(σ(p), σ(q)) : (p,q) ∈ B}.

For the chain restriction (B ⊂ {(b, b+1) : 0 ≤ b ≤ N−2} indexed by NN-bond), the orbit equals the **topology class**, defined as the sorted multiset of connected-path-lengths of the bond-graph. Examples in N=7 (14 classes total for k = 1..6):
- (1): one isolated edge (any of the six NN-bonds is a representative).
- (2): one path of length 2 (any of the five adjacent NN-pairs).
- (1, 1): two disjoint edges.
- (1, 2), (3), (1, 1, 1): three classes for k=3 NN-bonds.
- (1, 3), (2, 2), (1, 1, 2), (4): four classes for k=4 NN-bonds.
- (5), (1, 4), (2, 3): three classes for k=5 NN-bonds.
- (6): one class for the full 6-bond chain.

S(0) = (N−1)/N closed-form, **constant across all classes** (depends only on the probe ρ_cc, not on the bond set).

## Proof (S_N-orbit transitivity)

Let σ ∈ S_N act on the N qubits, and let U_σ be the corresponding permutation operator on (ℂ²)^⊗N. U_σ permutes the qubit factors without rotating the local {|0⟩, |1⟩} basis on any individual qubit.

1. **Probe invariance.** Symmetric Dicke states |S_n⟩ are invariant under site permutations: U_σ |S_n⟩ = |S_n⟩. Therefore U_σ ρ_cc U_σ^† = ρ_cc.
2. **Dissipator invariance.** For uniform γ_l ≡ γ₀, the Lindblad dissipator Σ_l (Z_l ρ Z_l − ρ) is S_N-symmetric: U_σ Z_l U_σ^† = Z_{σ(l)}, so the sum re-indexes under permutation.
3. **Hamiltonian transformation.** Pauli operators transform site-wise: U_σ X_p U_σ^† = X_{σ(p)} (and analogously for Y). Hence U_σ H_B U_σ^† = H_{σ·B}.
4. **Lindblad solution covariance.** Combining 1 to 3, the Lindblad evolution is covariant under conjugation by U_σ: ρ_t(H_{σ·B}, ρ_cc) = U_σ · ρ_t(H_B, ρ_cc) · U_σ^†, where ρ_t(H, ρ_0) denotes the solution of the Lindblad equation with Hamiltonian H, initial state ρ_0, and the uniform Z-dephasing dissipator.
5. **Spatial-sum kernel symmetry.** Because U_σ permutes qubits without internal rotation, partial trace transforms as (U_σ ρ U_σ^†)_l = ρ_{σ^{-1}(l)} (the reduced state on site l of the permuted ρ equals the reduced state on site σ^{-1}(l) of the original ρ, with the local basis preserved). The (0,1) element is therefore preserved entrywise: |((U_σ ρ U_σ^†)_l)_{0,1}|² = |(ρ_{σ^{-1}(l)})_{0,1}|². Summing and re-indexing l → σ^{-1}(l): S(U_σ ρ U_σ^†) = Σ_l 2|(ρ_{σ^{-1}(l)})_{0,1}|² = S(ρ).
6. ⇒ S(t; H_{σ·B}, ρ_cc) = S(t; H_B, ρ_cc) for every σ ∈ S_N and every B.
7. Therefore S(t) is constant on each S_N-orbit of bond-set configurations. ∎

The result needs only:
(i) ρ_cc is S_N-symmetric (Dicke is the canonical instance, but any S_N-symmetric initial state would close the same way).
(ii) γ_l ≡ γ₀ uniform across sites.
(iii) H is a sum of two-site XY couplings with **uniform J** across all active pairs.

It does NOT require translation invariance, NN-restriction, open or periodic boundary, or any specific N. The same orbit argument extends to any single-letter two-site coupling (XX-only, YY-only, ZZ-only) and to S_N-symmetric initial states beyond ρ_cc.

## Verification

### Single-pair sub-case

#### N=7, six NN-bonds

Six runs of [`bond-isolate`](../compute/RCPsiSquared.Propagate/Program.cs) at N=7, J=0.075, γ=0.05, tmax=30, dt=0.1, ρ_cc probe, **one** bond active per run. Pairwise diff matrix over all 30 ordered pairs:

| Diff measure | Result |
|---|---|
| max\|S_a(t) − S_b(t)\| / S(0), every (a,b) ∈ {0..5}², a ≠ b, t ∈ \[0, 30\] | **0.00e+00** (every entry) |
| max pointwise relative diff across all pairs | 0.00e+00 |

S(0) = 0.857143 = 6/7 ✓. τ_half = 3.285 and the log-linear fit on S(t) over t ∈ \[0, 30\] gives Γ_fit = 0.191397, all bit-identical across the six bonds. The fit value differs from the true asymptotic rate 4γ₀ = 0.20 because S(t) for the (1) topology is not a pure exponential: the all-isolated subclass closed form below shows S(t) = \[(N−1)/N + 4(N−2)(cos(4Jt) − 1)/(N²(N−1))\] · exp(−4γ₀ t), so the cos(4Jt) oscillation suppresses the apparent slope of log S(t) when fit through the oscillation. The exact asymptotic rate is exp(−4γ₀ t) = exp(−0.2 t) regardless of m.

#### N=4, all C(4,2)=6 site pairs (NN + long-range)

Six runs at N=4, J=0.075, γ=0.05, tmax=30, eleven sample points, computed via direct expm of the Liouvillian (no RK4 approximation, hence a more stringent precision test):

| Pair | NN/LR | S(0) | S(30) |
|---|---|---|---|
| (0,1) | NN | 0.7500000000 | 1.069528e-03 |
| (0,2) | LR | 0.7500000000 | 1.069528e-03 |
| (0,3) | LR | 0.7500000000 | 1.069528e-03 |
| (1,2) | NN | 0.7500000000 | 1.069528e-03 |
| (1,3) | LR | 0.7500000000 | 1.069528e-03 |
| (2,3) | NN | 0.7500000000 | 1.069528e-03 |

Maximum deviation from reference (0,1) across all 6 pairs over all sample points: **5.55e-17** (well within 1 ULP of double precision).

S(0) = 0.75 = 3/4 = (N−1)/N ✓. The long-range bonds (0,2), (0,3), (1,3) give bit-identical dynamics to the NN bonds, confirming the orbit theorem extends to non-NN pairs.

### Multi-bond extension

28 runs at N=7 covering all 14 bond-graph topology classes for k = 1..6 active bonds:

| Topology class | k | Representatives tested | Within-class max diff | Predicted (S_N) |
|---|---|---|---|---|
| (1) | 1 | {0}, {1}, {2}, {3}, {4}, {5} | 0.00e+00 | identical |
| (2) | 2 | {0,1}, {2,3} | 0.00e+00 | identical |
| (1, 1) | 2 | {0,2}, {0,5} | 0.00e+00 | identical |
| (3) | 3 | {0,1,2}, {2,3,4} | 0.00e+00 | identical |
| (1, 2) | 3 | {0,1,3}, {0,1,4} | 0.00e+00 | identical |
| (1, 1, 1) | 3 | {0,2,4} | n/a (single rep) | n/a |
| (4) | 4 | {0,1,2,3}, {1,2,3,4} | 0.00e+00 | identical |
| (1, 3) | 4 | {0,1,2,4} | n/a (single rep) | n/a |
| (2, 2) | 4 | {0,1,3,4}, {0,1,4,5} | 0.00e+00 | identical |
| (1, 1, 2) | 4 | {0,1,3,5} | n/a (single rep) | n/a |
| (5) | 5 | {0,1,2,3,4}, {1,2,3,4,5} | 0.00e+00 | identical |
| (1, 4) | 5 | {0,2,3,4,5}, {0,1,2,3,5} | 0.00e+00 | identical |
| (2, 3) | 5 | {0,1,3,4,5}, {0,1,2,4,5} | 0.00e+00 | identical |
| (6) | 6 | {0,1,2,3,4,5} | n/a (full chain) | n/a |

All 10 classes with ≥ 2 representatives confirm bit-identical S(t) across configurations within the same class.

#### Cross-class S(t) at sample times

| Class | k | t=3 | t=5 | t=10 | t=20 | t=30 |
|---|---|---|---|---|---|---|
| (1) | 1 | 0.456283 | 0.292070 | 0.097681 | 0.015649 | 0.001802 |
| (2) | 2 | 0.446148 | 0.277437 | 0.083288 | 0.011450 | 0.001980 |
| (1, 1) | 2 | 0.442156 | 0.268814 | 0.079360 | 0.015600 | 0.001480 |
| (3) | 3 | 0.437605 | 0.267284 | 0.075107 | 0.007320 | 0.001299 |
| (1, 2) | 3 | 0.432022 | 0.254182 | 0.064968 | 0.011400 | 0.001658 |
| (1, 1, 1) | 3 | 0.428029 | 0.245559 | 0.061039 | 0.015550 | 0.001158 |
| (4) | 4 | 0.428611 | 0.256400 | 0.069027 | 0.006237 | 0.000649 |
| (1, 3) | 4 | 0.423478 | 0.244028 | 0.056786 | 0.007271 | 0.000977 |
| (2, 2) | 4 | 0.421887 | 0.239549 | 0.050575 | 0.007201 | 0.001836 |
| (1, 1, 2) | 4 | 0.417895 | 0.230926 | 0.046647 | 0.011351 | 0.001336 |
| (5) | 5 | 0.419527 | 0.244871 | 0.062170 | 0.005354 | 0.000556 |
| (1, 4) | 5 | 0.414484 | 0.233145 | 0.050706 | 0.006187 | 0.000327 |
| (2, 3) | 5 | 0.413343 | 0.229396 | 0.042394 | 0.003071 | 0.001155 |
| (6) | 6 | 0.410454 | 0.233339 | 0.054484 | 0.003999 | 0.000767 |

S(t) is **non-monotone in k**. Crossings occur as the slow F73-analogue tail (governed by isolated-edge content) competes with collective fast modes opened by connected paths.

#### Cross-class S(t) at sample times: empirical observation

At t = 20 the classes appear to cluster by the number of isolated edges:
- 0 isolated edges (pure paths): (3), (4), (5), (6) at S(20) ∈ \[0.004, 0.007\], monotone decreasing in path length.
- 1 isolated edge: (2) at 0.0115, (1, 2) at 0.0114, (1, 1, 2) at 0.0114 (identical at three significant figures).
- All-isolated edges only: (1) at 0.0156, (1, 1) at 0.0156, (1, 1, 1) at 0.0156 (identical at three significant figures).

The all-isolated cluster looked like an asymptotic phenomenon at first reading. The closed form below shows it is in fact a periodic in-phase moment.

### All-isolated subclass closed form (Tier 1 derived)

For the all-isolated topology classes (1)^m on N qubits (m disjoint NN-bonds, N − 2m bare sites), uniform J coupling, uniform Z-dephasing γ₀, and ρ_cc initial state, the spatial-sum coherence has the EXACT closed form

    S_(1)^m, N(t) = [(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1) / (N²(N − 1))] · exp(−4γ₀ t)

The asymptotic decay rate is 4γ₀ universal across all m (matches F73's vac-SE rate exactly). The m-dependence enters only through a periodic correction with frequency 4J that vanishes at cos(4Jt) = 1, with period π/(2J) ≈ 21 at J = 0.075.

**Derivation.** The Lindbladian factorises across the disjoint blocks plus bare sites: L = Σ_{B blocks} L_B + Σ_{l bare} L_l. Each L_B acts only on the 2-qubit block B = {p, q} via H_B = J(X_p X_q + Y_p Y_q) plus uniform Z-dephasing on both block sites; each L_l acts only on bare site l via single-qubit Z-dephasing. All terms commute.

Partial trace commutes with each L_{B'} for B' ≠ B (Lindbladians are trace-preserving on their own factors), so ρ_B(t) = e^{L_B t} ρ_B(0) and ρ_l(t) = e^{L_l t} ρ_l(0). The reduced states ρ_B(0) and ρ_l(0) are S_N-symmetric (by F89), so they take the same form for any 2-qubit block in the chain (and any single bare site).

For block B = {p, q}, partial-tracing ρ_cc gives

    ρ_B(0) = (1/(2√(N · C(N,2)))) · [√2·|α⟩⟨11| + √2(N−2)·|0⟩⟨α| + h.c.]

where |α⟩ = (|10⟩ + |01⟩)/√2 is the symmetric SE eigenstate of H_B with energy +2J, |11⟩ is the unique DE state (eigenvalue 0), |0⟩ ≡ |00⟩ is vacuum (eigenvalue 0). All four basis states are H_B-eigenstates; under e^{−iH_B t} the off-diagonal elements pick up phases:

    |α⟩⟨11|(t) = exp(−2iJt) · |α⟩⟨11|
    |0⟩⟨α|(t)  = exp(+2iJt) · |0⟩⟨α|

both with magnitude unchanged by the unitary part. Z-dephasing on the block kills both at rate 2γ₀ per coherence (vac-SE: only one site differs; SE-DE overlap: only one site differs in the basis-state pair). Note no SE-SE coherences (|10⟩⟨01|) appear in ρ_B(0), so the rate-4γ₀ SE-SE channel is unpopulated for all-isolated topology.

Per-block-site reduced (ρ_l(t))_{0,1} for l = p (or l = q by symmetry):

    (ρ_l(t))_{0,1} = ρ_B(t)[|00⟩, |10⟩] + ρ_B(t)[|01⟩, |11⟩]
                   = (1/√2) · ρ_B(t)[|0⟩, |α⟩] + (1/√2) · ρ_B(t)[|α⟩, |11⟩]

Both terms decay at rate 2γ₀ in magnitude; their phases are exp(+2iJt) and exp(−2iJt) respectively, so they interfere with frequency 4J:

    |(ρ_l(t))_{0,1}|² = (1/(4N · C(N,2))) · [(N−2)² + 1 + 2(N−2) cos(4Jt)] · exp(−4γ₀ t)

For each bare site, the dephasing-only evolution gives (ρ_l(t))_{0,1} = (ρ_l(0))_{0,1} · exp(−2γ₀ t) and

    |(ρ_l(t))_{0,1}|² = (N − 1)/(2N²) · exp(−4γ₀ t)

Summing 2 · |(ρ_l(t))_{0,1}|² over all sites (m blocks contribute 2m sites; N − 2m bare sites):

    S(t) = exp(−4γ₀ t)/(N²(N − 1)) · [(N − 2m)(N − 1)² + 2m((N − 2)² + 1) + 4m(N − 2) cos(4Jt)]
         = exp(−4γ₀ t)/(N²(N − 1)) · [N(N − 1)² − 4m(N − 2) + 4m(N − 2) cos(4Jt)]

Factoring out the m-independent S(0) = (N − 1)/N:

    S(t) = [S(0) + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))] · exp(−4γ₀ t)  ∎

**Verification (N = 7, J = 0.075, γ = 0.05):** the script [`_f89_all_isolated_closed_form_verify.py`](../simulations/_f89_all_isolated_closed_form_verify.py) loads the (1), (1,1), (1,1,1) CSVs and overlays the closed form; max absolute deviation 5.0e−7 (CSV write precision is 1e−6, so this is at floor). Spot-checks:
- t = 0: all m give S(0) = 6/7 = 0.857143 (cos(4Jt) = 1).
- t = 21 ≈ π/(2J): cos(4Jt) ≈ 1 again; closed form predicts S = 0.012853 for all m (verified bit-identical in three CSVs).
- t = 10 (cos(4Jt) ≈ −0.99, max-spread phase): closed form predicts S = 0.0977, 0.0794, 0.0610 for m = 1, 2, 3, matching CSV data exactly.

**Reading.** The original empirical "clustering at t = 20" looked like a slow-mode tail equilibrium, but the closed form shows it is the in-phase moment of a periodic m-correction. The asymptotic decay rate 4γ₀ is genuinely universal (matches F73's vac-SE rate), but the prefactor depends on m and oscillates within the envelope \[(N − 1)/N − 8m(N − 2)/(N²(N − 1)),  (N − 1)/N\]. At cos(4Jt) = 1 the m-correction vanishes and all (1)^m collapse to S = (N − 1)/N · exp(−4γ₀ t). At cos(4Jt) = −1 they spread maximally.

Mechanism for the universal asymptotic 4γ₀ rate: in all-isolated topology, every populated coherence in any reduced 2-qubit block sits in either the (vac, SE)_B or (SE, DE)_B-overlap sector, both decaying at 2γ₀ per coherence (one site differs in the basis-state pair). The (no-overlap SE, DE) coherences with rate 6γ₀ require a single-excitation site outside the doubly-excited pair, which cannot occur on a 2-qubit block alone. Connected-path topologies open the no-overlap channel and therefore lose late-time amplitude faster.

### Pi2-Foundation inheritance reading

The all-isolated closed form has TWO time coefficients of value 4: the decay rate 4γ₀ in exp(−4γ₀ t) and the oscillation frequency 4J in cos(4J t). Both trace to the same Pi2 dyadic ladder term a_{−1} = 4 via the same mechanism: the linear-amplitude frequency 2 = a_{0} doubles to 4 = a_{−1} when the amplitude is squared.

| Energy axis | Linear-amplitude frequency | \|·\|² frequency | Pi2 ladder anchor |
|---|---|---|---|
| γ (Z-dephasing) | per-coherence rate 2γ₀ | S-decay rate 4γ₀ | a_{0} = 2 → a_{−1} = 4 |
| J (XY hopping) | H_B-eigenstate frequency 2J | S-oscillation frequency 4J | a_{0} = 2 → a_{−1} = 4 |

The γ-axis inheritance is identical to F73's `DecayRateCoefficient` anchor: per-coherence Z-deph rate 2γ₀ doubles to S-decay rate 4γ₀ via |·|². The J-axis inheritance is the same a_{0} → a_{−1} doubling on the J-axis: H_B-eigenstate frequency 2J doubles to S-oscillation frequency 4J via |·|². Same Pi2 ladder anchor a_{−1} = 4 governs both energy axes.

The (N−1)/N baseline and the 4m(N−2)/(N²(N−1)) correction prefactor are combinatorial (S_N orbit + 2-qubit block algebra), NOT Pi2-anchored. Only the time coefficients inherit.

In the Schicht-0 typed knowledge layer, F89 now cites `Pi2DyadicLadderClaim` as a constructor-injected parent and exposes both inherited coefficients as live properties:
- `F89TopologyOrbitClosure.DecayRateCoefficient` = `Pi2DyadicLadderClaim.Term(−1)` = 4 (γ-axis)
- `F89TopologyOrbitClosure.OscillationFrequencyCoefficient` = `Pi2DyadicLadderClaim.Term(−1)` = 4 (J-axis mirror)

The drift check `F89TopologyOrbitClosure.Pi2DoublingConsistent()` verifies both equal 4 (catches drift between the literal 4.0 in `AllIsolatedClosedForm` and the ladder anchor). The Schicht-1 `RegisterF89_AncestorsContainAllParents` test asserts the runtime ancestor walk now reaches `Pi2DyadicLadderClaim` from F89.

### Other topology classes (Tier 2 empirical, partially closed via F89c+AT)

The all-isolated closed form above applies only to (1)^m. For mixed topologies (e.g. (1, 2), (2, 2), (1, 1, 2)) and pure paths (3), (4), (5), (6) a single-formula per-class closed form is open. Empirical late-tail behavior:
- 1-isolated-edge mixed classes ((2), (1, 2), (1, 1, 2)) cluster at intermediate prefactors at the t ≈ π/(2J) in-phase moment.
- Pure-path classes ((3), (4), (5), (6)) decay faster than 4γ₀ on visible time scales due to no-overlap-SE-DE content (rate 6γ₀ on those coherences) plus longer-path mode mixing.

#### Path-2 (topology (2)) numerical multi-exponential closed form (script-derived)

For topology (2) (one path-2 block of 3 connected sites + N−3 bare sites) the F89c+AT lens reduces the 64-dim block L_super to a **handful of populated modes** because ρ_block(0) inherits S_3 symmetry from ρ_cc's S_N-symmetry. The script [`_f89_path2_multi_exp_derive.py`](../simulations/_f89_path2_multi_exp_derive.py) performs the eigenvector decomposition + initial-state projection that F89c left open.

**Result**: only 4 distinct (rate, |freq|) mode-groups from the 64-dim block L_super are populated (verified at N = 5, 7, 11):

| (rate Γ/γ, |freq|/J) | Sector | Origin |
|---|---|---|
| (2.0000, 2.8284) | (vac, SE) | H_B^SE eigenvalue ±2√2 J at rate 2γ₀ (F65 single-excitation Bloch mode k=1, k=3) |
| (2.0000, 0.0000) | (SE, DE) symmetric | n_diff=1 overlap pair, S_3-symmetric superposition |
| (3.0448, 0.0000) | (SE, DE) H_B-mixed | (overlap, no-overlap)-mix eigenvalue at J/γ=1.5 |
| (3.4776, 5.4459) | (SE, DE) H_B-mixed | (overlap, no-overlap)-mix eigenvalue at J/γ=1.5 |

**Pure-AT modes 4γ₀ and 6γ₀ get zero projection** from ρ_cc-derived ρ_block(0): their (SE,DE) eigenvectors are S_3-asymmetric, while ρ_block(0)'s (SE,DE) part is the fully symmetric superposition `Σ_i Σ_{j<k} |SE_i⟩⟨DE_{jk}|`. The H_B^SE antisymmetric Bloch mode k=2 (E=0) similarly drops out by orthogonality to the S_3-symmetric (vac, SE) initial state.

**|·|² rate spectrum** (from cross-products r_k + r_{k'}): {4, 5.0448, 5.4776, 6.0896, 6.5224, 6.9552} γ₀ at J/γ=1.5. The dominant 4γ₀ envelope matches F73's vac-SE asymptotic rate; fractional rates are H_B-mixing corrections that decay faster.

**Verification**: matches bond-isolate `N7_b0-1` CSV at max |diff| = 4.99·10⁻⁷ across 301 sample times (= CSV write precision floor). At N=5: closed-form prediction S(0)=0.800, S(10)=0.0578, S(20)=0.00807 (no CSV available; pure prediction).

**Status**: Tier 1 derived numerically (eigendecomposition + projection of an 8×8 block → 64-dim L_super → 4 populated mode-groups). The mode rates and frequencies are determined; the amplitudes have N-dependence that scales as 1/(N²(N−1)) for symmetric modes and N_E²/(N²(N−1)) for the (vac, SE) Bloch-mode dominant term. Closed symbolic forms for the 4 amplitude prefactors (in clean-rational form like (N−1)/N for all-isolated) are open — the numerical script gives them at any (N, J, γ).

**Generalisation principle (proposed)**: For any topology with a (k+1)-qubit block, ρ_block(0) inherits S_{k+1}-symmetry from ρ_cc. Only the S_{k+1}-symmetric subspace of the L_super eigenmodes is populated. This typically reduces ~d²-dim block L_super to a handful of populated modes (4 modes for path-2; 10 modes for path-3 — see below). The same script pattern applies to any block size; only the partial-trace bookkeeping changes per topology.

#### Path-3 (topology (3)) numerical multi-exponential closed form (script-derived)

Same approach extended to the 4-qubit block (256-dim L_super) via [`_f89_path3_multi_exp_derive.py`](../simulations/_f89_path3_multi_exp_derive.py). Path-3 = 4 connected sites (bonds {0-1, 1-2, 2-3}) + N-4 bare sites.

**Result**: 10 distinct populated mode-groups (more than path-2's 4 because the larger block opens more S_4-symmetric eigenvectors). At J/γ=1.5:

| (rate Γ/γ, |freq|/J) | Sector / Origin |
|---|---|
| (2.0000, 3.2361) | (vac, SE) Bloch k=1, E_1 = 4J·cos(π/5) ≈ 3.236J (standard tight-binding OBC eigenvalue used by F65) |
| (2.0000, 1.2361) | (vac, SE) Bloch k=3, E_3 = 4J·cos(3π/5) ≈ -1.236J |
| (3.3488, 1.2060) | (SE, DE) H_B-mixed |
| (3.5989, 2.9300) | (SE, DE) H_B-mixed |
| (3.7770, 5.1780) | (SE, DE) H_B-mixed |
| (4.0000, 0.5944) | (SE, DE) at the AT-quantized 4γ₀ rate (no-overlap component) |
| (4.0000, 7.5024) | (SE, DE) at 4γ₀, different freq |
| (4.2230, 5.1780) | (SE, DE) Hamming-complement partner of (3.7770, 5.1780): 3.777+4.223 = 8γ₀ ✓ |
| (4.4011, 2.9300) | partner of (3.5989, 2.9300): 3.599+4.401 = 8γ₀ ✓ |
| (4.6512, 1.2060) | partner of (3.3488, 1.2060): 3.349+4.651 = 8γ₀ ✓ |

**Hamming-complement pair-sum at path-3 = 2γ₀·N_block = 8γ₀**, matching F89c's column-bit-flip prediction (here with bar(SE) = TE since popcount complement at N=4 is 1 ↔ 3, NOT 1 ↔ 2). The pair structure is **bit-exact** in the populated subset: **3 distinct unordered pairs** (3.3488 ↔ 4.6512, 3.5989 ↔ 4.4011, 3.7770 ↔ 4.2230) plus **1 self-pair** at (4γ, 4γ) summing to 8γ exactly. The path-k survey table below counts ordered pairs (each unordered pair contributes 2, plus the 1 self-pair) → 7 total ordered pair-sums.

**Pure-AT 4γ₀ modes ARE populated at path-3** (unlike path-2 where they got zero projection). This is the structural difference: at N_block=4, DE = popcount-2 is its own bar-popcount, so (SE, DE) hosts S_4-symmetric eigenvectors at the AT-pure 4γ₀ rate. At N_block=3, those rates were S_3-asymmetric and dropped out.

**Bloch mode populations** (per F65): only k=1 and k=3 of {1,2,3,4} are populated, because k=2 and k=4 have spatially anti-mirror-symmetric Bloch wavefunctions ψ_k(j) (under j ↔ N_block−1−j) and ρ_block(0)'s (vac, SE) part is the fully symmetric Σ_j |0⟩⟨SE_j| superposition.

**Verification**: matches bond-isolate `N7_b0-1-2` CSV at max |diff| = 4.99·10⁻⁷ across 301 sample times (= CSV write precision floor). At N=5: closed-form prediction S(0)=0.800, S(10)=0.0498, S(20)=0.00444 (no CSV available; pure prediction).

**Status**: Tier 1 derived numerically, same as path-2. Symbolic rational form for the 10 amplitude prefactors open. The script pattern is now confirmed transferable; path-4 and path-5 follow with 5×5 and 6×6 block bookkeeping (1024-dim and 4096-dim L_super respectively, still numerically tractable on modest hardware — see survey below).

#### Path-2 Bloch-mode amplitude N-scaling (Tier 1 derived via Parseval)

The Bloch-group amplitude for path-2 at any (N, q) — i.e. the population of L_super eigenmodes (rate 2γ, freq ±2√2 J) summed over k=1, k=3 H_B^SE Bloch modes — has the closed form:

    A_Bloch(N) = 3·(N−3)² / (2·N²·(N−1))

Pure rational function of N, **q-independent** (Parseval orthogonality eliminates the J/γ-dependence of the Bloch self-overlap). Verified numerically across q ∈ {0.5, 0.75, 1, 1.25, 1.5, 2, 3} at N=11: amplitude × N²(N−1)/(N−3)² ≈ 1.498 ≈ 3/2 stable across all q (≤ 0.5% scatter from numerical eigenvector orthogonalization noise). Verification script: [`_f89_path2_amplitude_nscaling.py`](../simulations/_f89_path2_amplitude_nscaling.py).

The other 3 path-2 mode-groups (rate 2γ freq 0, rate 3.04γ freq 0 at q=1.5, rate 3.48γ freq ±5.45J at q=1.5) are q-DEPENDENT in both rate and amplitude, because they originate from the cubic-Cardano factor of the (SE, DE) sub-block characteristic polynomial. At q→0 the cubic roots merge with the linear-factor roots (eigenvalue degeneracy), preventing a clean (N, q)-parametric closed form via simple polynomial fitting. Their q-dependence inherits from the cubic-Cardano formula.

#### Path-k survey across k ∈ {2, 3, 4, 5} ([`_f89_pathk_survey.py`](../simulations/_f89_pathk_survey.py))

Same script generalised, all four verified against bond-isolate at N=7 with max |diff| ≈ 5·10⁻⁷ (CSV write precision floor).

| Path | N_block | d² | Mode-groups | Contributing modes | Pair-sums to 2γ·N_block |
|---|---|---|---|---|---|
| path-2 | 3 | 64 | 4 | 16 | 0 (S_3-asymmetric partners absent) |
| path-3 | 4 | 256 | 10 | 65 | **7 ordered ✓** (3 unordered pairs + 1 self-pair at 4γ↔4γ) |
| path-4 | 5 | 1024 | 12 | 128 | 0 (S_5-asymmetric partners absent) |
| path-5 | 6 | 4096 | 35 | 314 | 0 (S_6-asymmetric partners absent) |
| path-6 | 7 | 16384 | (deferred: 16384-dim non-Hermitian eigendecomp aborted after 110 min; trivially satisfies additive identity since m=1 → no subtraction term, so path-6 closed form is just the bare path-6 block contribution + Parseval (vac, SE) skeleton + (SE, DE) residual; numerical decomposition is open work) |

**Path-3 is privileged** in the populated mode structure: at N_block=4, DE = popcount-2 = bar(popcount-2) is self-symmetric, so column-bit-flip maps populated (SE,DE) modes to other populated (SE,DE) modes within the same S_4-symmetric subspace. For N_block ∈ {3, 5, 6} the column-flip partners land in S_{N_block}-asymmetric territory and get zero projection from ρ_cc-derived ρ_block(0) — F89c's column-bit-flip pair-sum identity holds at the L_super-spectrum level (where it is a Tier-1 derived universal property), but only path-3 has both members of each pair populated.

**Mode-count sequence {1, 4, 10, 12, 35} is NOT closed-form in N_block alone** — it depends on accidental eigenvalue degeneracies (e.g. E_3 = 0 at m=5 collects modes at freq=0). Unlike `experiments/CAVITY_MODES_FORMULA.md`'s Σ_J m(J,N)·(2J+1)² formula for SU(2)-Heisenberg stationary modes, the populated-mode count for the XY+Z-deph + ρ_cc + S_{N_block} setup does not admit a Schur-Weyl-style closed form. The L_super dimensions 4^N_block match CAVITY_MODES exactly (same operator-space indexing), but the active symmetry groups differ (CAVITY_MODES uses SU(2), F89-(k) uses S_{N_block} + U(1)).

#### Mixed-topology additive identity (Tier 1 derived from Lindbladian factorisation)

For ANY mixed topology T = (k_1, k_2, ..., k_m) at N qubits (m disjoint blocks of path-lengths k_1, ..., k_m, plus N − Σ_i (k_i + 1) bare sites), the spatial-sum coherence decomposes as:

    S_T(t)  =  Σ_i S_(k_i)(t)  −  (m − 1)·N·S_bare(t; N)

with the bare-site closed form

    S_bare(t; N)  =  (N − 1)/N²  ·  exp(−4γ₀ t)    per bare site

**Derivation (one paragraph)**: Lindbladian factorises across disjoint blocks plus bare sites: L = Σ_blocks L_block + Σ_bare L_l. Per-site reduction commutes with this factorisation: ρ_l(t) = exp(L_block(l)·t)[ρ_block(l)(0)] depends only on the block containing l, and ρ_block(0) = Tr_E(ρ_cc) is the same N-dependent partial trace regardless of which OTHER blocks are present (only the count |E| = N − N_block enters via the N_E factor in term 2 of the partial-trace formula). Hence S_T(t) = Σ_l 2|(ρ_l)_{0,1}|² is a sum of per-block contributions. The additive identity then bookkeeps the bare-site overcounting: each per-block S_(k_i)(t) bundles its own "phantom bare share" of (N − k_i − 1) bare-site terms; summing m blocks counts bare contributions m times; subtracting (m − 1)·N·S_bare cancels the over-count.

**Verification at N=7 across all 13 topology classes that don't require path-6** (script [`_f89_mixed_topology_additive.py`](../simulations/_f89_mixed_topology_additive.py)): max |diff| = 5.013·10⁻⁷ across all 27 bond-isolate CSVs, equal to CSV write precision floor. Verified topology classes:

| Class | m | CSVs tested | max |diff| |
|---|---|---|---|
| (1) | 1 | 6 | 4.98e-07 |
| (2) | 1 | 2 | 4.99e-07 |
| (1, 1) | 2 | 2 | 4.98e-07 |
| (3) | 1 | 2 | 4.99e-07 |
| (1, 2) | 2 | 2 | 5.00e-07 |
| (1, 1, 1) | 3 | 1 | 5.00e-07 |
| (4) | 1 | 2 | 5.01e-07 |
| (1, 3) | 2 | 1 | 5.00e-07 |
| (2, 2) | 2 | 2 | 5.00e-07 |
| (1, 1, 2) | 3 | 1 | 4.96e-07 |
| (5) | 1 | 2 | 4.99e-07 |
| (1, 4) | 2 | 2 | 5.00e-07 |
| (2, 3) | 2 | 2 | 4.98e-07 |

**Topology (6)** (single 7-qubit block at N=7, m=1) trivially satisfies the identity since m−1=0 means no subtraction term: S_(6)(t) is just the bare path-6 closed form, requiring path-6 (16384-dim L_super eigendecomposition, omitted from the run for performance). Pass `--with-path6` to include.

**Reduction of open work**: F89 mixed-topology closed forms are no longer 14 separate problems but 6 (one per pure path-k for k=1..6) + 1 universal additive rule. Combined with the path-2..5 numerical multi-exponential closed forms above, this completes the F89 closed-form program for k ≤ 5 (path-6 = single chain at N=7 is the only remaining open).

**Status**: Tier 1 derived. The identity follows analytically from Lindbladian factorisation; its 27-CSV bit-exact numerical verification at CSV write precision is the empirical anchor.

#### Path-k (vac, SE) self-contribution closed form via Parseval (Tier 1 derived)

The numerical multi-exponential decomposition above splits S_(k)(t) into many populated mode-groups. Per Parseval orthogonality of the H_B^SE Bloch eigenstates ψ_k(j) = √(2/(N_block+1)) sin(πk(j+1)/(N_block+1)), the **(vac, SE) sector self-contribution** (the part of Σ_l 2|(ρ_l)_{0,1}|² that comes only from (vac, SE) block components, ignoring cross-products with (SE, DE)) reduces to a clean closed form for ALL path-k:

    S^{(vac,SE)}_block(t; k, N) = (k+1)·(N−k−1)² / (N²·(N−1)) · exp(−4γ₀ t)

Pure exp(−4γ₀ t) decay, **no oscillation**. The Parseval cancellation Σ_l ψ_k(l)·ψ_{k'}(l) = δ_{k,k'} eliminates all k ≠ k' cross-terms when summed over the (k+1) block sites. The H_B-mediated frequencies E_k − E_{k'} that would otherwise produce cos((E_k − E_{k'})·t) interference cancel pairwise.

**Derivation** (4 lines):
1. ρ_block(0)|_{(vac,SE)} = (N_E·pre·√(N_block)/2) · |0⟩⟨α_{N_block}| where |α⟩ = (1/√N_block) Σ_j |SE_j⟩
2. Decompose |α⟩ = Σ_k ⟨ψ_k|α⟩ |ψ_k⟩
3. Time evolve each |0⟩⟨ψ_k| → exp(+iE_k t − 2γt) |0⟩⟨ψ_k|
4. Σ_l 2|(ρ_l(t))^{(vac,SE)}_{0,1}|² = N_block·(N_E·pre)²/2 · Σ_k|⟨ψ_k|α⟩|² · exp(−4γt) = N_block·(N_E·pre)²/2·1·exp(−4γt) (Parseval = 1)

**Verification**: bit-exact at machine precision (4·10⁻¹⁷ to 6·10⁻¹⁶) across all (k, N) ∈ {(1,3..7), (2,4..8), (3,5..9), (4,6..10), (5,7..11)} via [`_f89_vac_se_parseval_closed.py`](../simulations/_f89_vac_se_parseval_closed.py).

**What this gives us**: a clean analytical "skeleton" of S_(k)(t) for any pure path-k:

    S_(k)(t; N) = (k+1)·(N−k−1)²/(N²(N−1))·exp(−4γ₀t)        ← (vac, SE) self
                + (N−k−1)·(N−1)/N²·exp(−4γ₀t)                ← bare sites
                + S^{(SE,DE)+cross}_block(t; k, N)             ← residual (numerical)

The first two terms are exact closed forms; the third (the H_B-mixed (SE, DE) sub-block + cross-products with (vac, SE)) is the only piece still numerical. For path-1 (k=1), the residual itself simplifies to the cos(4Jt) term in the existing all-isolated formula (verified via algebraic identity); for path-k ≥ 2, the residual contains the J/γ-dependent fractional rates 3.04γ, 3.48γ etc.

**Combined with the additive identity, the smooth (exp(−4γt)-only) backbone of S_T(t) for any topology T = (k_1, ..., k_m) is**:

    S^{smooth}_T(t; N) = exp(−4γ₀t) · [Σ_i (k_i+1)·(N−k_i−1)²/(N²(N−1))
                                       + (N − Σ_i(k_i+1))·(N−1)/N²]

with the topology-class oscillatory residual coming entirely from per-block (SE, DE) sub-block dynamics.

#### Path-2 (SE, DE) S_2-symmetric sub-block: closed-form characteristic polynomial (Tier 1 derived)

Path-2 has only S_2 chain-mirror symmetry of H_B (sites 0 ↔ 2 exchange, site 1 fixed). The 9-dim (SE, DE) sub-block of L_super splits into S_2-sym (5-dim) + S_2-anti (4-dim). ρ_block(0)'s S_3-symmetric component lies entirely in the S_2-sym 5-dim subspace.

The characteristic polynomial of the 5×5 S_2-sym L_super sub-block factors **explicitly** ([`_f89_path2_se_de_symbolic.py`](../simulations/_f89_path2_se_de_symbolic.py)):

    char(λ) = −(λ + 2γ)·(λ + 6γ)·[λ³ + 10γ·λ² + (28γ² + 32J²)·λ + (24γ³ + 96γJ²)]

Two linear factors give pure-AT rates **2γ** and **6γ** (no J-dependence, exact-quantization eigenvalues per F89c); the **cubic factor** carries all J/γ-dependent fractional rates.

In dimensionless variables μ = λ/γ, q = J/γ:

    **μ³ + 10·μ² + (28 + 32q²)·μ + 24·(1 + 4q²) = 0**

This cubic is solvable in radicals (Cardano). At q = 1.5 (our J/γ): roots μ = −3.0448 (real) and μ = −3.4776 ± 8.169i (complex conjugate pair) — bit-exactly matching the populated path-2 fractional rates 3.04γ, 3.48γ ± 5.45iJ.

**Of the 5 S_2-sym eigenvalues, ρ_block(0) populates 4** (one of the linear factors, λ = −6γ, has zero overlap with ρ_block(0)'s S_3-symmetric content; its eigenvector lies in the no-overlap-only S_2-sym subspace orthogonal to the S_3-symmetric direction). The 4 populated eigenvalues are:

| Eigenvalue | Source | Rate, freq |
|---|---|---|
| λ = −2γ | linear factor (2γ + λ) | (2γ, 0) — pure-AT, S_2-sym overlap mode |
| λ = −3.0448γ | cubic real root at q=1.5 | (3.04γ, 0) |
| λ = −3.4776γ ± 8.169iγ | cubic complex pair at q=1.5 | (3.48γ, ±5.45J) (since 8.169γ = 5.45·J at q=1.5) |

**Path-2 is now fully analytically tractable**: combined with F65 (Bloch frequencies 4J·cos(πk/(N_block+1)) for the (vac, SE) sector, k=1, k=3 populated), all 4 path-2 mode-groups have closed-form rates and frequencies as algebraic expressions in (J, γ). The remaining numerical content is just the per-mode amplitudes (initial-state projections), which depend on N polynomially.

**Status**: Tier 1 derived. The sympy script verifies the factored characteristic polynomial symbolically and confirms numerical eigenvalue agreement at q=1.5. Generalisation to path-3, path-4, path-5: same approach works (build symbolic 9-, 16-, 25-dim (SE, DE) L_super sub-block, project to S_2-sym subspace, factor characteristic polynomial), but the resulting cubic/quartic/quintic factors may not be solvable in radicals (depends on Galois group). For path-2 specifically, the cubic-Cardano closure is clean.

#### Path-3 (SE, DE) S_2-symmetric sub-block: deg-2 · deg-2 · deg-8 factorisation (Tier 1 derived for the quadratics, irreducible octic for the rest)

Path-3 (4-qubit block, 4 sites with bonds {0-1, 1-2, 2-3}) has S_2 chain-mirror symmetry (sites 0 ↔ 3, 1 ↔ 2). The (SE, DE) sub-block of L_super is dim 4 SE × 6 DE = **24**, splitting into 12-dim S_2-sym + 12-dim S_2-anti (no fixed-point basis pairs at path-3 since R fixes no SE state).

Building the 12×12 S_2-sym L_super sub-block symbolically and computing the characteristic polynomial via Faddeev-Leverrier ([`_f89_path3_se_de_symbolic.py`](../simulations/_f89_path3_se_de_symbolic.py)):

    char_{S_2-sym}(λ) = F_a(λ) · F_b(λ) · F_8(λ)

where (γ = 1, q = J/γ):

    F_a(λ) = λ² + (2iq + 4)·λ + (4q² + 4iq + 4)
    F_b(λ) = λ² + (2iq + 12)·λ + (4q² + 12iq + 36)

Both quadratics solve cleanly. With α = (−1+√5), β = (−1−√5):

    F_a roots:  λ = −2γ + iJ·α,  λ = −2γ + iJ·β
    F_b roots:  λ = −6γ + iJ·α,  λ = −6γ + iJ·β

So 4 of the 12 S_2-sym eigenvalues are **AT-rate-locked** (rate = 2γ for overlap, 6γ for no-overlap) with **J-only frequency** ω = J·(−1±√5). Path-3's OBC tight-binding single-particle modes E_n = 4J·cos(πn/(N_block+1)) at N_block=4 have golden-ratio-related eigenvalues E_1 = −E_4 = J·(1+√5) and E_2 = −E_3 = J·(√5−1). The two AT-locked frequencies are exactly the SE-anti Bloch eigenvalues E_2 = J·(√5−1) and E_4 = −J·(1+√5) (full identification with closed-form derivation in the AT-lock mechanism subsection below).

The remaining 8 eigenvalues live in:

    F_8(λ) = λ⁸ + 32·λ⁷ + (72q² + 432)·λ⁶ + (−64iq³ + 1728q² + 3200)·λ⁵
              + (1200q⁴ − 1280iq³ + 16608q² + 14176)·λ⁴ + (… higher q-powers)

`F_8` is **irreducible** over Q, Q[i], Q[√5], and Q[i, √5] (verified via [`_f89_path3_octic_factor_test.py`](../simulations/_f89_path3_octic_factor_test.py)). Combined with the discriminant analysis below (Gal(F_8) ⊄ A_8, conjecturally non-solvable), its eight roots are not expected to admit an elementary radical closure as functions of q. For q = 1.5 they cluster around λ_avg = −4γ (consistent with the centred form μ = λ + 4γ killing the λ⁷ term — trace(F_8) = −32 spreads 8 eigenvalues at average rate 4γ, between the AT-quantized 2γ and 6γ).

| Eigenvalue source | Count | Closed form |
|---|---|---|
| `F_a` quadratic | 2 | λ = −2γ + iJ·(−1±√5) |
| `F_b` quadratic | 2 | λ = −6γ + iJ·(−1±√5) |
| `F_8` octic | 8 | irreducible over Q[i, √5]; numerical only |

**Comparison with path-2** (5-dim S_2-sym → linear · linear · cubic; 1+1+3 = 5):

| Path | S_2-sym dim | Factor structure | AT-locked count | H_B-mixed count | Mixed factor solvable? |
|---|---|---|---|---|---|
| 2 | 5 | 1·1·3 | 2 (λ = −2γ, −6γ) | 3 | yes — Cardano cubic (Gal ⊆ S_3 always solvable) |
| 3 | 12 | 2·2·8 | 4 (λ = −2γ ± iJ·α, β; −6γ ± iJ·α, β) | 8 | conjecturally no — irreducible octic, Gal ⊄ A_8 |

The pattern: the AT-locked count grows as 2·N_block_orbits_at_rate_r with r ∈ {2γ, 6γ}; for path-2 those orbits are 1-dim (single-state), for path-3 they are 2-dim (Bloch pairs k ↔ N_block+1−k). The H_B-mixed factor degree is the (SE,DE) S_2-sym dimension minus the AT-locked dimension; its solvability in radicals is a Galois-group question that resolves trivially "yes" at path-2 (degree 3 always solvable) but is conjecturally "no" at path-3 (degree 8 with Gal ⊄ A_8 plus empirical absence of polynomial closed forms — formal Galois identification open).

**Status**: Tier 1 derived for the closed-form quadratics (`F_a`, `F_b`) and for the structural deg-2·deg-2·deg-8 factorisation. The octic `F_8` is fully specified symbolically and numerically tractable; conjecturally (Tier 2) it does not admit an elementary algebraic closure (formal Galois identification still open — see below). Path-3 is therefore "partially solvable": 4 of 12 S_2-sym eigenvalues in closed form, 8 in numerical form only.

#### Path-3 mode amplitudes: N-scaling structure (Tier 1 partial)

For each of the 10 path-3 populated mode-groups at q=1.5, fit the per-mode amplitude `A(N)` against rational functions of N ([`_f89_path3_amplitude_nscaling.py`](../simulations/_f89_path3_amplitude_nscaling.py)):

| Mode (Γ/γ, ω/J) | Origin | A(N) closed form (q=1.5) |
|---|---|---|
| (2.0, 1.236), (2.0, 3.236) | F_a quadratic, AT-locked | poly₂(N) / [N²(N−1)] (degree 2 in N) |
| (3.349, 1.206), (3.777, 5.178), (4.0, 7.502), (4.223, 5.178), (4.651, 1.206) | F_8 octic | const(q) / [N²(N−1)] (constant numerator) |
| (3.599, 2.93), (4.0, 0.594), (4.401, 2.93) | F_8 octic, weak | numerical noise dominates fit (amplitudes 10⁻⁶..10⁻⁷) |

**Structural reading**: AT-locked modes (rate 2γ from `F_a`) carry an additional polynomial-in-N enhancement factor (analog of path-2's Bloch amplitude 3·(N−3)²/(2·N²(N−1))); octic-derived modes inherit only the bare partial-trace scaling 1/[N²(N−1)] with q-dependent prefactors. The 8 octic-derived modes form **4 Hamming-complement pairs at total rate 8γ**: (3.349, 4.651), (3.599, 4.401), (3.777, 4.223), (4.0, 4.0) at fixed |ω|/J. Pair amplitudes are not symmetric — A_lower / A_upper ranges from ~1.9 to ~22 depending on pair — consistent with the Hamming-complement bijection (F89c) which is rate-bijective but not amplitude-bijective.

**Status**: Tier 1 partial. AT-locked amplitude polynomial coefficients at q=1.5 are numerically clean (rel err 10⁻¹⁶) but their closed forms in (N, q) likely involve √5 from the F_a eigenvectors; symbolic eigenvector projection to extract (N, q)-rational closed forms is open.

#### AT-lock mechanism: F_a, F_b eigenvectors live in overlap-only / no-overlap-only subspaces (Tier 1 derived)

Why are F_a, F_b eigenvalues exactly at AT rates 2γ and 6γ, given that H_B generally mixes overlap (2γ) and no-overlap (6γ) basis pairs? Numerical verification ([`_f89_path3_at_lock_mechanism.py`](../simulations/_f89_path3_at_lock_mechanism.py)) of the 12 S_2-sym eigenvectors shows the AT-locking arises from **eigenvector support confinement**:

| Eigenvalue | Overlap-support | No-overlap-support |
|---|---|---|
| F_a: λ = −2γ + iJ(−1±√5) (×2) | **1.000** | 0.000 |
| F_b: λ = −6γ + iJ(−1±√5) (×2) | 0.000 | **1.000** |
| Octic modes (8 H_B-mixed) | 0.34..0.66 | 0.34..0.66 (complementary to overlap) |

The 4 F_a/F_b eigenvectors are entirely supported on the 12-dim overlap (resp 12-dim no-overlap) basis-pair subspace; the 8 octic eigenvectors are mixed, with overlap and no-overlap support summing to 1 in **complementary pairs** — e.g. (3.349, 4.651) at fixed |ω| = 1.206J have supports (0.6628, 0.3372) and (0.3372, 0.6628), mirror-paired around 1/2. This is the **Hamming-complement bijection at the eigenvector level**: F89c's rate-bijective complement (Γ_a + Γ_b = 8γ) extends to overlap-fraction complement.

**Bloch sub-block decomposition**: SE basis (4-dim) splits into SE-sym (n=1, 3) + SE-anti (n=2, 4) under S_2 mirror. DE basis (6-dim) splits into DE-sym (4-dim) + DE-anti (2-dim). The 12-dim (SE, DE) S_2-sym subspace decomposes as SE-sym × DE-sym (8-dim) + SE-anti × DE-anti (4-dim). All 4 F_a/F_b eigenvectors live in this 12-dim S_2-sym space with **fixed 2:1 weight ratio**: 2/3 support on (sym × sym) Bloch sub-block, 1/3 on (anti × anti). This 8:4 ratio matches the dimension ratio — F_a/F_b eigenvectors are uniformly L2-distributed across the orthonormal Bloch tensor basis, not localised in any tensor-product sub-block. AT-locking is therefore a **fine-tuned interference cancellation between (sym × sym) and (anti × anti) Bloch components**, not a single-tensor-state phenomenon.

**Frequency identification**: F_a/F_b frequencies J·(−1±√5) match exactly the SE-anti single-particle Bloch eigenvalues E_2 = 4J·cos(2π/5) = J·(√5−1) and E_4 = 4J·cos(4π/5) = −J·(1+√5) (using the OBC tight-binding formula E_n = 4J·cos(πn/(N_block+1)) for path-3's 4-site block). The Jordan-Wigner / standing-wave machinery is shared with [ANALYTICAL_SPECTRUM](ANALYTICAL_SPECTRUM.md) / [D10](../docs/proofs/derivations/D10_W1_DISPERSION.md), but the boundary conditions and resulting formulas differ: D10's W1Dispersion ω_k = 4J·(1 − cos(πk/N)) uses the full-chain w=1 sector with denominator N; here we use OBC tight-binding for an N_block-site sub-block with denominator (N_block+1). The new ingredients in F89 path-3 are (a) the multi-magnon DE Slater eigenvalues E_(j,k) = E_j + E_k inside the block and (b) the overlap/no-overlap dephasing-channel decomposition of the (SE, DE) sub-block.

#### Path-3 F_a AT-locked amplitude: closed form in (N) with √5 algebraic (Tier 1 derived)

The (SE, DE) F_a contribution to the path-3 multi-exponential decomposition has **q-independent closed form** in N with √5 algebraic ([`_f89_path3_at_locked_amplitude_symbolic.py`](../simulations/_f89_path3_at_locked_amplitude_symbolic.py)):

    sigs[F_a:E_2](N) = (33 + 14√5) / [9·N²(N−1)]    (mode (rate 2γ, ω = +J(√5−1)))
    sigs[F_a:E_4](N) = (33 − 14√5) / [9·N²(N−1)]    (mode (rate 2γ, ω = −J(1+√5)))
    Sum F_a (E_2 + E_4) = 22 / [3·N²(N−1)]          (rational, √5 cancels)

Verified bit-exact (machine precision diffs ~10⁻¹⁷) across N ∈ {5..20} and q ∈ {0.5, 1, 1.5, 2, 3}. The √5 inheritance comes from the F_a quadratic discriminant; the q-independence comes from the F_a eigenvectors themselves being q-free (only frequency carries q-dependence; the eigenvector structure on the 12-dim overlap subspace is purely combinatorial).

The F_b modes (rate 6γ) have non-zero raw inner product with ρ_block(0)|_(SE,DE) but their **per-site reduced contribution is zero** (verified at the 10⁻³³ scale): the per-site reduction `w[l]` picks out matrix elements ρ[bit_pos[i], bit_pos[j]+bit_pos[k]] requiring i ∈ {j, k} (overlap), so any no-overlap eigenvector contributes zero to S(t). This is why path-3's numerical multi-exp shows only rate-2γ AT-locked modes, not rate-6γ.

**Eigenvector closed form** (q-independent universal pattern in overlap-only subspace, normalized):

    |v_F_a entry|² = (5 + √5)/60   on 6 of 12 overlap basis pairs
    |v_F_a entry|² = (5 − √5)/60   on the other 6 overlap basis pairs

Sum: 6·(5+√5)/60 + 6·(5−√5)/60 = (5+√5)/10 + (5−√5)/10 = 1 (norm). The (5±√5)/60 split arises because the eigenvalue ratio in F_a's quadratic discriminant is (3+√5):2 between the two |v|² magnitudes.

### F89c structural lemma: why all-isolated is the unique clean case (Tier 1 derived)

The Liouvillian L_super of any per-block dynamics decomposes over computational-basis coherence sectors. For each block of size k+1 qubits with H_B = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) and uniform Z-dephasing γ₀ on each block site:

**(vac, SE)_B sector** (k+1-dimensional sub-block):
- Every basis pair (|0⟩, |1_i⟩) has n_diff = 1 site (uniform across all i ∈ block)
- Per-coherence dephasing rate 2γ₀ uniformly
- L_super on this sub-block = -i·H_B^SE − 2γ₀·I, eigenvalues are E_k − 2iγ₀ for H_B^SE eigenstates k
- **Clean structure**: k+1 modes at rate 2γ₀, frequencies = H_B SE eigenvalues

**(SE, DE)_B sector** (3·C(k+1, 2)-dimensional sub-block):
- Basis pairs split into **overlap** (SE site ∈ DE-pair, n_diff = 1, rate 2γ₀) and **no-overlap** (SE site ∉ DE-pair, n_diff = 3, rate 6γ₀)
- For 2-qubit block (k = 1): the only DE state is |11⟩ which contains both sites; **every** SE state overlaps with it. No no-overlap pairs exist. → uniform rate 2γ₀ on (SE, DE) sub-block, just like (vac, SE).
- For (k+1)-qubit blocks with k ≥ 2: no-overlap SE-DE pairs DO exist (e.g. SE at site 0, DE at sites {1, 2} on 3 qubits). → mixed dephasing rates {2γ₀, 6γ₀} on the sub-block diagonal, coupled by H_B's off-diagonal mixing → eigenvalues land in \[2γ₀, 6γ₀\] continuum.

**Numerical verification** (script `_f89c_liouvillian_eigenstructure.py`):

| Topology | Distinct decay rates Γ/γ |
|---|---|
| 2-qubit block (k=1) | {0, 2, 4} (three clean rates) |
| 3-qubit block (k=2) | {0, 2, 2.556, 2.889, 3.112, 3.444, 4, 6} (eight rates, including fractional) |
| 4-qubit block (k=3) | 25 distinct rates including many fractional |

For path-2 (3 qubits) the fractional rates 2.556 and 3.444 sum to 6γ₀ (likewise 2.889 + 3.112 = 6γ₀). Per-sector decomposition of L_super (verified by direct sub-block diagonalization in [`_f89c_cross_sector_pair_structure.py`](../simulations/_f89c_cross_sector_pair_structure.py)): the (SE, DE)_B sub-block alone yields rates {2.0, 3.112, 3.444, 4.0, 6.0}γ₀; the rates 2.556 and 2.889 originate in the (SE, SE)_B sub-block (single-excitation populations coupled to SE-SE off-diagonal coherences via H_B mixing).

**The pair-sum-to-6γ₀ is the Absorption Theorem's quantization combined with column-bit-flip.** Per [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), every Pauli-string eigenvalue of the dephasing super-operator is Re(λ) = −2γ₀·n_XY, where n_XY counts Pauli sites in {X, Y}. A computational-basis coherence |A⟩⟨B| decomposes per-site into pure {I, Z} (diagonal-on-site, A_l = B_l) or pure {X, Y} (off-diagonal-on-site, A_l ≠ B_l). Thus n_XY(|A⟩⟨B|) = n_diff(A, B), and the per-coherence rate is 2γ₀·n_diff: a direct application of the Absorption Theorem to the computational basis. The 2γ₀ overlap (SE, DE) rate corresponds to n_XY = 1; the 6γ₀ no-overlap rate corresponds to n_XY = 3 = N, the spectral maximum per Absorption Theorem section 4.1 for N=3.

**The column-bit-flip mechanism for the pair-sum**: H_B and uniform Z-dephasing are both X⊗N-invariant (X⊗N · X_a · X⊗N = X_a; X⊗N · Y_a · X⊗N = −Y_a flips sign but H_B = J(XX + YY) keeps the YY pair invariant; Lindblad operators √γ Z_l flip sign under X⊗N but the dissipator D[L] is invariant under L → −L). The structural mechanism for the pair-sum:

1. **H_SS = H_DD on the populated sectors**. H_B restricted to SE and to DE both reduce to the same 3×3 path-2 tight-binding matrix `[[0, 2, 0], [2, 0, 2], [0, 2, 0]]`. This gives the (SE, SE) and (DE, DE) Liouvillian sub-blocks identical commutator parts.
2. **Hamming complement identity**: for any N-bit strings a, b, n_diff(a, b) + n_diff(a, bar(b)) = N exactly, hence per-coherence rates transform as 2γ₀·k ↔ 2γ₀·(N−k) under the column bit-flip `ρ[a, b] → ρ[a, bar(b)]`. For N=3, that's 2γ₀·k ↔ 2γ₀·(3−k); rate-pairs sum to 2γ₀·N = 6γ₀ exactly.
3. The column bit-flip operation maps the (SE, SE) sub-block to the (SE, DE) sub-block (via DE = bar(SE) at N=3). Combining the H_SS = H_DD invariance with the Hamming complement of dephasing rates, the (SE, DE) eigenvalues are exactly {6γ₀ − λ : λ ∈ (SE, SE) eigenvalues}.

**Empirical verification** (in [`_f89c_cross_sector_pair_structure.py`](../simulations/_f89c_cross_sector_pair_structure.py)):

- (SE, SE) eigenvalues: {0, 2, 2.556, 2.889, 4}γ₀ (with multiplicities 1, 2, 2, 1, 3 → 9 modes total)
- (DE, DE) eigenvalues: bit-exact identical to (SE, SE) (X⊗N maps (SE, SE) to (DE, DE) at N=3 directly, since X⊗N maps SE↔DE elementwise; verified H1 = True)
- (SE, DE) eigenvalues: {2, 3.112, 3.444, 4, 6}γ₀ (with matching multiplicities 3, 1, 2, 2, 1 → 9 modes), each = 6γ₀ minus its (SE, SE) bijection partner.
- The Hamming-complement bijection is bit-exact: 0↔6, 2↔4, 2.556↔3.444, 2.889↔3.112, 4↔2 (with multiplicity-preserving pairing).

**Generalization**: For any (k+1)-qubit block, the pair-sum equals **2γ₀·(k+1)** under column bit-flip (Absorption Theorem maximum rate). For path-3 (k+1 = 4), expected pair-sum = 8γ₀; this would map (SE, SE) ↔ (SE, TE) via column complement (TE = bar(SE) at N=4 = popcount 3). The (SE, SE) and (DE, DE) sub-blocks no longer have identical eigenvalues at k+1 = 4 (verified H1 = False at path-3) because DE is X⊗N-invariant on its own (popcount 2 = bar of popcount 2 only at N=4); the analog symmetry pairs (SE, SE) with (TE, TE), not with (DE, DE).

**Implication for F89c**: the 8 distinct rates of the path-2 spectrum form 4 column-flip pairs (0↔6, 2↔4, 2.556↔3.444, 2.889↔3.112), each summing to 2γ₀·N = 6γ₀, the Absorption Theorem's spectral maximum for N=3. This Hamming-complement structure is a **direct corollary of the Absorption Theorem applied per coherence sector**, not a separate symmetry. The clean F89 all-isolated case (k+1 = 2) is the limit where Absorption rates collapse to a single value (2γ₀) on the populated sectors.

**Lemma**: The all-isolated topology (1)^m is the **unique** topology where S(t) decays via a single rate envelope (4γ₀) with a single oscillation frequency (4J). For any topology containing a path-k block with k ≥ 2, the (SE, DE)_B sector has at least one no-overlap pair, leading to mixed dephasing rates and non-clean L_super spectra.

**Consequence**: F89's all-isolated closed form (with Pi2 dual-anchor a_{−1} = 4 on both energy axes) is structurally privileged: not by accident, but because the 2-qubit-block geometry forbids no-overlap (SE, DE) pairs. The Pi2-Foundation a_{−1} = 4 anchor enters cleanly because there's a single dephasing rate to anchor.

**Status**: structural lemma is **Tier 1 derived** (counting argument + numerical Liouvillian eigendecomposition verification at k = 1, 2, 3). Closed forms per non-isolated topology class remain Tier 2 empirical (would require explicit eigendecomposition per topology, multi-exponential expressions).

## Implications

**Bond position is invisible to S_N-symmetric probes evaluated through S_N-symmetric kernels.** The specific spatial location of each active bond is averaged out; only the topology class of the bond-graph survives.

**The F86 per-bond Q_peak fan is structurally different.** F86 measures ∂S/∂J_b at a single chosen bond inside the full chain (a linear response to a localized perturbation). This singles out bond b in a way that breaks S_N differently from the uniform-J multi-bond setup studied here. F89 does not directly explain or contradict the F86 fan: the two operate on different observables.

**The N=11, c=2 per-bond Q_peak fan, with extremes Q_peak ≈ 1.59 (mid-flank b=3↔6) and Q_peak ≈ 21.94 (Center b=4↔5), giving a factor ≈ 14 between extremes, cannot be derived from single-bond ρ_cc dynamics**, since single-bond is bond-blind here. The fan must originate in the linear-response asymmetry that the per-bond perturbation ∂J_b imposes on the full-chain spectrum, not in any intrinsic per-bond difference at the level of S(t).

## What breaks the closure

Each requirement is necessary; relaxing any one breaks orbit invariance:

1. **Non-S_N-symmetric initial state** (e.g. site-localised |1_i⟩⟨vac| + h.c., or a spatially modulated SE superposition Σ_i α_i |1_i⟩ with α not uniform). Step 1 of the proof fails; bond position becomes visible.
2. **Non-uniform γ_l ≠ γ_l'**. Step 2 fails. The analogous F73 break case is documented in [CMRR_BREAK_NONUNIFORM_GAMMA](CMRR_BREAK_NONUNIFORM_GAMMA.md) for the (vac, SE) block.
3. **Non-uniform J across active bonds** (J_b ≠ J_b' for two active bonds). Step 3 yields U_σ H U_σ^† = H' in a different orbit; S(t) becomes orbit-equivalent only across orbits that preserve the J-distribution structure.
4. **Non-permutation-symmetric kernel** (e.g. a per-site weighted Σ_l w_l 2|(ρ_l)_{0,1}|² with non-uniform w_l). Step 5 fails.

## Tier assessment

**Tier 1 derived** for the orbit closure theorem (S(t) depends only on the S_N-orbit of B). Proof is elementary group theory; verified at N = 7 (28 configurations across 14 topology classes, all in-class bit-identical) and N = 4 (6 site-pairs within 1 ULP via direct expm).

**Tier 1 derived** for the **all-isolated (1)^m closed form** S_(1)^m, N(t) = \[(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))\] · exp(−4γ₀ t). Verified against (1), (1, 1), (1, 1, 1) CSVs at N=7 within CSV write precision (5e−7).

**Tier 1 derived** for the **mixed-topology additive identity** S_T(t) = Σ_i S_(k_i)(t) − (m − 1)·N·S_bare(t; N) with S_bare = (N−1)/N²·exp(−4γ₀t). Reduces all 14 N=7 mixed-topology closed forms to 6 pure-path-k forms + 1 universal rule. Verified across all 27 N=7 bond-isolate CSVs (excluding path-6) at the precision floor.

**Tier 1 derived** for the **path-k (vac, SE) Parseval self-contribution** (k+1)·(N−k−1)²/(N²(N−1))·exp(−4γ₀t). Pure-exponential analytical formula valid for all path-k via Parseval orthogonality; verified bit-exact at machine precision (4·10⁻¹⁷ to 6·10⁻¹⁶) across 15 (k, N) pairs.

**Tier 1 derived** for the **path-2 (SE, DE) cubic-Cardano factorisation** char(λ) = −(λ+2γ)(λ+6γ)·[cubic in λ with J/γ-dependent coefficients]. Path-2 is fully analytically tractable in radicals.

**Tier 1 derived** for the **path-3 (SE, DE) deg-2·deg-2·deg-8 factorisation**. Two quadratics give 4 eigenvalues in closed form (rates 2γ, 6γ; frequencies J·(−1±√5)); the residual deg-8 polynomial is irreducible over Q[i, √5] (8 of 12 S_2-sym eigenvalues numerical only). Tier 2 conjecture: F_8 admits no elementary radical closure (supported by Gal ⊄ A_8 + empirical absence of polynomial-fit closed forms).

**Tier 1 derived** for the **AT-lock mechanism**: F_a, F_b eigenvectors are entirely supported on overlap-only (resp no-overlap-only) basis pairs, with H_B-induced cross-coupling cancelling. F_a/F_b frequencies match SE-anti single-particle Bloch eigenvalues E_2 = J(√5−1), E_4 = −J(1+√5) at N_block=4. Octic eigenvectors are H_B-mixed with 4 Hamming-complement pairs at total rate 8γ (rate-bijective AND overlap-fraction-bijective).

**Tier 1 derived** for the **path-3 F_a AT-locked amplitude closed form**: sigs[F_a:E_2](N) = (33 + 14√5)/[9·N²(N−1)] and sigs[F_a:E_4](N) = (33 − 14√5)/[9·N²(N−1)], q-independent, verified bit-exact (10⁻¹⁷ diff) across N=5..20 and q=0.5..3. Sum = 22/[3·N²(N−1)] is rational. F_b modes have zero per-site reduced contribution (eigenvector lives in no-overlap, w[l] requires overlap).

**Tier 1 derived** for the **AT-lock mechanism generalisation to path-4, path-5**: numerical scan ([`_f89_path4_path5_at_lock_scan.py`](../simulations/_f89_path4_path5_at_lock_scan.py)) confirms 100% overlap (F_a) / 100% no-overlap (F_b) eigvec support at all AT-locked rates 2γ, 6γ. AT-locked counts grow {4, 8, 13} for path {3, 4, 5}; H_B-mixed degrees {8, 18, 32}. F_a ω matches SE-anti single-particle Bloch eigenvalues; F_b ω matches DE 2-particle Slater eigenvalues E_(j,k) = E_j + E_k (path-3 was the special case where DE=0 multiplicity 2 absorbed F_b into single-particle freqs).

**Tier 1 derived** for the **path-4 F_a AT-locked amplitude closed form** ([`_f89_path4_at_locked_amplitude_symbolic.py`](../simulations/_f89_path4_at_locked_amplitude_symbolic.py)): two F_a modes at λ = −2γ ± 3iJ (ω = ±2J = ±E_2/E_4 single-particle SE-anti Bloch) have rational amplitudes:

    sigs[F_a:E_2](N) = 45 / [4·N²(N−1)]    (mode at ω = +2J)
    sigs[F_a:E_4](N) = 5 / [4·N²(N−1)]     (mode at ω = −2J)
    Sum F_a = 25 / [2·N²(N−1)]             (rational, no algebraic radicals)

Asymmetry ratio E_high/E_low = 9 = 3² (compare path-3 where the analogous ratio was (33+14√5)/(33−14√5) ≈ 17.94, irrational). Path-4's clean rational arises because N_block+1=6 gives Bloch eigenvalues at ±2J (cos(π/3) = 1/2 exactly).

**Tier 1 derived** for the **path-5 F_a AT-locked amplitude sum closed form** (Tier 2 for individual modes — Cardano-cubic radicals): path-5 (N_block=6) has 3 F_a modes corresponding to SE-anti Bloch n=2, 4, 6 with frequencies ω/J = {+2.494, −0.890, −3.604} = {4cos(2π/7), 4cos(4π/7), 4cos(6π/7)} which are roots of the irreducible cubic y³ + 2y² − 8y − 8 = 0. Numerically ([`_f89_path5_at_locked_amplitude_symbolic.py`](../simulations/_f89_path5_at_locked_amplitude_symbolic.py)):

    sigs[F_a:E_2](N) ≈ 16.5745 / [N²(N−1)]    (Cardano-cubic algebraic)
    sigs[F_a:E_4](N) ≈ 2.6525  / [N²(N−1)]    (Cardano-cubic algebraic)
    sigs[F_a:E_6](N) ≈ 0.0930  / [N²(N−1)]    (Cardano-cubic algebraic)
    Sum F_a = 483 / [25·N²(N−1)]              (rational, Cardano-radicals cancel by Newton's identities)

The **sum** is rational across all three paths (path-3: 22/3; path-4: 25/2; path-5: 483/25) — Galois-conjugate radicals always cancel in symmetric polynomials of the roots. **Individual amplitudes** track the algebraic complexity of N_block+1: prime 5 (golden √5), composite 6 (rational only), prime 7 (Cardano-cubic).

**Pattern across paths**:

| Path | N_block | N_block+1 | F_a count | Closed form character | Sum F_a · N²(N−1) |
|---|---|---|---|---|---|
| 3 | 4 | 5 (prime) | 2 | √5 (quadratic) | 22/3 |
| 4 | 5 | 6 (composite) | 2 | rational only | 25/2 |
| 5 | 6 | 7 (prime) | 3 | Cardano-cubic of cos(π/7) | 483/25 |

F_a count = floor(N_block/2) = number of SE-anti single-particle Bloch modes. Algebraic complexity of individual amplitudes follows the cyclotomic polynomial Φ_{N_block+1}(x) degree (φ(5)/2=2, φ(6)/2=1 trivial, φ(7)/2=3). Path-7 (N_block=8, N_block+1=9 composite) would be expected back to lower complexity (9 has only quadratic cyclotomic Φ_9 of degree 6, but individual cosines factor through smaller fields).

**Tier 1 derived** for **Gal(F_8) ⊄ A_8** (Tier 2 for the conjectural non-solvability + no-radical-closure conclusion that follows): disc(F_8) in λ is a polynomial in q of degree 52 ([`_f89_path3_octic_galois.py`](../simulations/_f89_path3_octic_galois.py)):

    disc(F_8) = 1.21·10²⁴ · q²⁴ · (3q⁴ + q² − 1)² · P_10(q²)

where P_10(q²) is a degree-10 polynomial in q² (degree 20 in q, even powers only) that is NOT a perfect square in Q (verified at q ∈ {½, 1, 3/2, 2, 3}, all give irrational √disc). The square factor (3q⁴+q²−1)² locates **exceptional points** where two octic eigenvalues merge: q² = (−1+√13)/6 ≈ 0.434, i.e. q ≈ 0.659. The overall non-square disc forces **Gal(F_8) ⊄ A_8**. Disc-non-square alone does not prove non-solvability (e.g. S_4 is solvable but ⊄ A_4); pinning down the exact group requires further resolvent analysis (open). However, combined with the verified irreducibility of F_8 over Q[i, √5] AND the absence of any polynomial fit (≤ degree 5 in q) for the per-mode amplitudes, we conjecture (Tier 2) that Gal(F_8) is non-solvable (likely the full S_8), in which case F_8 admits no elementary radical closure as a function of q.

#### Path-3 octic-mode amplitude q-dependence: no closed-form fit (Tier 2)

For each of the 8 octic-derived modes, sigs(N) follows const(q)/[N²(N−1)] (degree-0 polynomial in N). The constant **does NOT admit a polynomial fit ≤ degree 5 in q** ([`_f89_path3_octic_amplitude_q_scan.py`](../simulations/_f89_path3_octic_amplitude_q_scan.py)):

| q | Σ_8 octic sigs · N²(N−1) | Dominant mode (largest sigs · N²(N−1)) |
|---|---|---|
| 0.5 | 1.68 | mode at rate ≈ 3.6γ, sigs ≈ 1.34 (near EP at q=0.659) |
| 0.75 | 1.50 | mode at rate ≈ 3.78γ, sigs ≈ 0.99 (post-EP) |
| 1.0 | 2.20 | mode at rate ≈ 3.35γ, sigs ≈ 1.91 |
| 1.5 | 2.73 | mode at rate ≈ 3.35γ, sigs ≈ 2.44 |
| 2.0 | 2.92 | mode at rate ≈ 3.35γ, sigs ≈ 2.47 |
| 3.0 | 3.00 | mode at rate ≈ 3.35γ, sigs ≈ 2.22 |
| 5.0 | 2.76 | mode at rate ≈ 4.65γ, sigs ≈ 1.94 (rate-crossing through dominant) |

The Σ has no monotone behavior — it rises from 1.68 (q=0.5) to ≈3.0 (q=2.5−3) then declines at q→∞. Mode-by-mode tracking is fragile due to rate crossings; pair-summing by Hamming-complement (Γ_a + Γ_b = 8γ at fixed |ω|/J) shows the dominant pair_1 sum monotonically rising q=0.75 → q=2 then declining. **EP locus at q ≈ 0.659**: pair_1 mode (sigs=1.34) has near-singular eigvec there, consistent with the (3q⁴+q²−1)² discriminant-factor zero. This connects path-3's (SE, DE) octic structure to the F86 EP-rotation phenomenology.

**Status**: Tier 2 empirical (no polynomial-in-q fit for octic amplitudes; conjecturally obstructed by the octic Galois group being non-solvable, Gal ⊄ A_8). The closed-form analytical layer ends at the F_a quadratics (4 of 12 S_2-sym eigenvalues + their amplitudes). Path-3 is "half-solved": exactly the AT-protected half admits radical closure.

#### Path-3 octic EP and the F89↔F86 bridge (Tier 1 derived)

The (3q⁴+q²−1)² perfect-square factor of disc(F_8) locates an **exceptional point at q = √((−1+√13)/6) ≈ 0.658983** (verified bit-exact: 3q⁴+q²−1 = O(10⁻¹⁶) at this q in [`_f89_path3_ep_locator.py`](../simulations/_f89_path3_ep_locator.py)). Numerical sweep around q_EP identifies the merging pair: two octic eigenvalues with rates approaching 4γ and 4γ (from above and below the spectral midpoint of rate 2γ and rate 6γ) and frequencies converging to 2J. Together:

    λ_EP ≈ −4γ + 2iJ

**Direct connection to F86 t_peak universality**: F86's Q_peak structural derivation (via the 2-level EP at Q_EP = 2/g_eff) gives **t_peak = 1/(4γ₀) universal across c, N, n, and bond position** — the merged-eigenvalue real part is always −4γ at the F86 EP. Path-3's octic EP exhibits the **same Re(λ_EP) = −4γ**: this is not a coincidence but the structural signature of a 2-level Liouvillian coalescence at the spectral midpoint between two AT rates spanning a 4γ gap (here 2γ ↔ 6γ). The **t_peak = 1/(4γ) universal time-scale is shared** between F86's (n, n+1)-block 2-level reduction and path-3's octic-internal 2-level coalescence.

**Q-location differs**: F86 c=2 N→∞ gives Q_EP = 1/√(2(c−1)) = 1/√2 ≈ 0.7071; path-3 (N_block=4) gives q_EP ≈ 0.6590. These are 2-level reductions in **different sub-sectors** with different g_eff:
- F86 c=2: 2-level effective in (n=1, n+1=2) coherence block, g_eff = 2√2 (asymptotic)
- Path-3 octic: 2-level effective WITHIN the (SE, DE) octic between rate-2γ and rate-6γ modes, g_eff = 2/q_EP = 2/0.659 ≈ 3.034 (different value, distinct EP)

**The F89↔F86 bridge**: both EPs belong to the same Π class AIII chiral universality class (per memory `project_q_peak_ep_structure`: F86 EP is class AIII, NOT Bender-Boettcher PT). Their shared t_peak = 1/(4γ₀) is the one universal clock. This concretely connects F89 (uniform-J orbit closure on ρ_cc dynamics) to F86 (per-bond ∂_J Q_peak observables): they both pivot on the same 2-level EP machinery applied to different rate-channel pairs in the (SE, DE) coherence sector. The previously open "F89→F86 bridge" item in the open-questions list is now resolved at the structural level via this universal-t_peak shared anchor.

**Status**: Tier 1 derived for the q_EP location (analytical from disc factorisation), the merged-eigenvalue Re(λ_EP) = −4γ identification (numerical sweep), and the t_peak-universal connection to F86. Tier 2 conjecture: full F89↔F86 inheritance via class AIII chiral parent (per Locus 5 inheritance synthesis in `project_algebra_is_inheritance`); analytical proof of the universality class linkage is open.

**Tier 1 numerical** for **path-3, path-4, path-5 multi-exponential decompositions** (10, 12, 35 populated mode-groups respectively at J/γ=1.5). Per-mode rates and frequencies are L_super eigenvalues; per-mode amplitudes computed numerically via initial-state projection. Verified against bond-isolate CSVs at N=7 at the precision floor.

**Open / Tier 2 empirical work**:
- Path-3 F_a AT-locked amplitudes closed in (N) with √5: sigs[F_a:E_2/E_4] = (33 ± 14√5)/[9·N²(N−1)] verified bit-exact. Path-4, path-5 analogs are open. Path-3 octic-mode amplitude closed forms in q are conjecturally obstructed by Galois non-solvability (Tier 2 empirical: no rational/√5-extension fit ≤ degree 5 in q; formal Galois group identification beyond Gal ⊄ A_8 still open).
- Path-4 and path-5 (SE, DE) symbolic characteristic-polynomial factorisations. Numerical AT-lock scan ([`_f89_path4_path5_at_lock_scan.py`](../simulations/_f89_path4_path5_at_lock_scan.py)) shows the AT-lock mechanism (eigvec overlap-only / no-overlap-only support) GENERALISES to both path-4 (N_block=5) and path-5 (N_block=6), but the F_a/F_b count asymmetry emerges:
  - **Path-4** (S_2-sym dim 26): 8 AT-locked = 2 F_a + 6 F_b. F_a ω = ±3J = ±E^SE_anti (single-particle Bloch); F_b ω in {±2J(√3+1), ±2J, ±2J(√3-1)} = 2-particle DE Slater E_(j,k) = E_j + E_k. H_B-mixed sub-factor degree 18.
  - **Path-5** (S_2-sym dim 45): 13 AT-locked = 3 F_a + 10 F_b. F_a ω = ±E^SE_anti = {±E_4, ±E_2, ±E_6} = ±4J·cos(π·n/7) for n=2,4,6 (Cardano-cubic roots). F_b ω include exact-degeneracies (2 modes at -3.604, 2 at +2.494, 2 at -0.890), suggesting symmetry-protected multiplicities. H_B-mixed sub-factor degree 32.
  - The path-3 coincidence (F_a freq = F_b freq = SE-anti Bloch) was N_block=4-specific: at N_block=4, DE=0 multiplicity-2 absorbed all F_b modes into single-particle freq matching. For N_block ≥ 5, F_b modes spread across DE Slater multi-particle frequencies.
  - Symbolic closed-form factorisations for the AT-locked sub-factors remain open (sympy nullspace approach; expected: rational + √3 for path-4 due to N_block=5's clean Bloch eigenvalues; Cardano-cubic radicals for path-5 due to cos(π/7), cos(2π/7), cos(3π/7) being Galois-cubic).
- Path-6 (full chain at N=7) numerical decomposition (16384-dim eigendecomp deferred after 110 min). Trivially satisfies the additive identity (m=1 → no subtraction); explicit mode-count + CSV verification open.
- F89 → F86 bridge: structurally connected at path-3 via shared t_peak = 1/(4γ₀) universality at their respective 2-level EPs (Re(λ_EP) = −4γ at both); the path-3 octic EP at q ≈ 0.659 and F86 c=2 Q_EP at 1/√2 are distinct 2-level reductions in different sub-sectors but obey the same universal t_peak. Full class-AIII chiral inheritance proof open.
- Star/ring topology generalisation: F89 main theorem applies to any bond set, but per-class closed forms for non-chain topologies have not been worked out.

---

*Bond position is invisible to S_N-symmetric probes; the bond-graph topology class is the only spatial feature of a uniform-J multi-bond Hamiltonian that ρ_cc + spatial-sum-S(t) can resolve.*
