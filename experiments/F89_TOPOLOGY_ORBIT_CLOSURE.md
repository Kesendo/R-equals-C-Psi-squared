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
| (2.0000, 3.2361) | (vac, SE) Bloch k=1, E_1 = 4J·cos(π/5) ≈ 3.236J per F65 |
| (2.0000, 1.2361) | (vac, SE) Bloch k=3, E_3 = 4J·cos(3π/5) ≈ -1.236J |
| (3.3488, 1.2060) | (SE, DE) H_B-mixed |
| (3.5989, 2.9300) | (SE, DE) H_B-mixed |
| (3.7770, 5.1780) | (SE, DE) H_B-mixed |
| (4.0000, 0.5944) | (SE, DE) at the AT-quantized 4γ₀ rate (no-overlap component) |
| (4.0000, 7.5024) | (SE, DE) at 4γ₀, different freq |
| (4.2230, 5.1780) | (SE, DE) Hamming-complement partner of (3.7770, 5.1780): 3.777+4.223 = 8γ₀ ✓ |
| (4.4011, 2.9300) | partner of (3.5989, 2.9300): 3.599+4.401 = 8γ₀ ✓ |
| (4.6512, 1.2060) | partner of (3.3488, 1.2060): 3.349+4.651 = 8γ₀ ✓ |

**Hamming-complement pair-sum at path-3 = 2γ₀·N_block = 8γ₀**, matching F89c's column-bit-flip prediction (here with bar(SE) = TE since popcount complement at N=4 is 1 ↔ 3, NOT 1 ↔ 2). The pair structure is **bit-exact** in the populated subset.

**Pure-AT 4γ₀ modes ARE populated at path-3** (unlike path-2 where they got zero projection). This is the structural difference: at N_block=4, DE = popcount-2 is its own bar-popcount, so (SE, DE) hosts S_4-symmetric eigenvectors at the AT-pure 4γ₀ rate. At N_block=3, those rates were S_3-asymmetric and dropped out.

**Bloch mode populations** (per F65): only k=1 and k=3 of {1,2,3,4} are populated, because k=2 and k=4 have spatially anti-mirror-symmetric Bloch wavefunctions ψ_k(j) (under j ↔ N_block−1−j) and ρ_block(0)'s (vac, SE) part is the fully symmetric Σ_j |0⟩⟨SE_j| superposition.

**Verification**: matches bond-isolate `N7_b0-1-2` CSV at max |diff| = 4.99·10⁻⁷ across 301 sample times (= CSV write precision floor). At N=5: closed-form prediction S(0)=0.800, S(10)=0.0498, S(20)=0.00444 (no CSV available; pure prediction).

**Status**: Tier 1 derived numerically, same as path-2. Symbolic rational form for the 10 amplitude prefactors open. The script pattern is now confirmed transferable; path-4 and path-5 follow with 5×5 and 6×6 block bookkeeping (1024-dim and 4096-dim L_super respectively, still numerically tractable on modest hardware — see survey below).

#### Path-k survey across k ∈ {2, 3, 4, 5} ([`_f89_pathk_survey.py`](../simulations/_f89_pathk_survey.py))

Same script generalised, all four verified against bond-isolate at N=7 with max |diff| ≈ 5·10⁻⁷ (CSV write precision floor).

| Path | N_block | d² | Mode-groups | Contributing modes | Pair-sums to 2γ·N_block |
|---|---|---|---|---|---|
| path-2 | 3 | 64 | 4 | 16 | 0 (S_3-asymmetric partners absent) |
| path-3 | 4 | 256 | 10 | 65 | **7 ✓** (full Hamming-complement) |
| path-4 | 5 | 1024 | 12 | 128 | 0 (S_5-asymmetric partners absent) |
| path-5 | 6 | 4096 | 35 | 314 | 0 (S_6-asymmetric partners absent) |

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

**Tier 1 derived** for the orbit closure theorem (S(t) depends only on the S_N-orbit of B). The proof is elementary group theory applied to the Lindblad equation. Numerical verification at N = 7 (multi-bond, 24 configurations across 12 topology classes, 8 with ≥ 2 representatives all bit-identical) and N = 4 (single-pair, 6 site-pairs identical within 1 ULP via direct expm) corroborates the proof at machine precision.

**Tier 1 derived** for the all-isolated (1)^m closed form S_(1)^m, N(t) = \[(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))\] · exp(−4γ₀ t). The derivation factors the Lindbladian over disjoint blocks plus bare sites, uses H_B-eigenstate phase tracking, and counts populated coherence sectors per block. Numerical verification matches the (1), (1, 1), (1, 1, 1) CSVs at N = 7 within CSV write precision (5e−7).

The **mixed-topology and pure-path closed forms** (per-class S(t) for (1, 2), (2, 2), (1, 1, 2), (3), (4), (5), (6) at N = 7) remain **Tier 2 empirical**; derivation open.

---

*Bond position is invisible to S_N-symmetric probes; the bond-graph topology class is the only spatial feature of a uniform-J multi-bond Hamiltonian that ρ_cc + spatial-sum-S(t) can resolve.*
