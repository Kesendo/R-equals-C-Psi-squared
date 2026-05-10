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
| max\|S_a(t) − S_b(t)\| / S(0), every (a,b) ∈ {0..5}², a ≠ b, t ∈ [0, 30] | **0.00e+00** (every entry) |
| max pointwise relative diff across all pairs | 0.00e+00 |

S(0) = 0.857143 = 6/7 ✓. τ_half = 3.285 and the log-linear fit on S(t) over t ∈ [0, 30] gives Γ_fit = 0.191397, all bit-identical across the six bonds. The fit value differs from the true asymptotic rate 4γ₀ = 0.20 because S(t) for the (1) topology is not a pure exponential: the all-isolated subclass closed form below shows S(t) = [(N−1)/N + 4(N−2)(cos(4Jt) − 1)/(N²(N−1))] · exp(−4γ₀ t), so the cos(4Jt) oscillation suppresses the apparent slope of log S(t) when fit through the oscillation. The exact asymptotic rate is exp(−4γ₀ t) = exp(−0.2 t) regardless of m.

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
- 0 isolated edges (pure paths): (3), (4), (5), (6) at S(20) ∈ [0.004, 0.007], monotone decreasing in path length.
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

**Reading.** The original empirical "clustering at t = 20" looked like a slow-mode tail equilibrium, but the closed form shows it is the in-phase moment of a periodic m-correction. The asymptotic decay rate 4γ₀ is genuinely universal (matches F73's vac-SE rate), but the prefactor depends on m and oscillates within the envelope [(N − 1)/N − 8m(N − 2)/(N²(N − 1)),  (N − 1)/N]. At cos(4Jt) = 1 the m-correction vanishes and all (1)^m collapse to S = (N − 1)/N · exp(−4γ₀ t). At cos(4Jt) = −1 they spread maximally.

Mechanism for the universal asymptotic 4γ₀ rate: in all-isolated topology, every populated coherence in any reduced 2-qubit block sits in either the (vac, SE)_B or (SE, DE)_B-overlap sector, both decaying at 2γ₀ per coherence (one site differs in the basis-state pair). The (no-overlap SE, DE) coherences with rate 6γ₀ require a single-excitation site outside the doubly-excited pair, which cannot occur on a 2-qubit block alone. Connected-path topologies open the no-overlap channel and therefore lose late-time amplitude faster.

### Pi2-Foundation inheritance reading

The all-isolated closed form has TWO time coefficients of value 4: the decay rate 4γ₀ in exp(−4γ₀ t) and the oscillation frequency 4J in cos(4J t). Both trace to the same Pi2 dyadic ladder term a_{−1} = 4 via the same mechanism: the linear-amplitude frequency 2 = a_{0} doubles to 4 = a_{−1} when the amplitude is squared.

| Energy axis | Linear-amplitude frequency | |·|² frequency | Pi2 ladder anchor |
|---|---|---|---|
| γ (Z-dephasing) | per-coherence rate 2γ₀ | S-decay rate 4γ₀ | a_{0} = 2 → a_{−1} = 4 |
| J (XY hopping) | H_B-eigenstate frequency 2J | S-oscillation frequency 4J | a_{0} = 2 → a_{−1} = 4 |

The γ-axis inheritance is identical to F73's `DecayRateCoefficient` anchor: per-coherence Z-deph rate 2γ₀ doubles to S-decay rate 4γ₀ via |·|². The J-axis inheritance is the same a_{0} → a_{−1} doubling on the J-axis: H_B-eigenstate frequency 2J doubles to S-oscillation frequency 4J via |·|². Same Pi2 ladder anchor a_{−1} = 4 governs both energy axes.

The (N−1)/N baseline and the 4m(N−2)/(N²(N−1)) correction prefactor are combinatorial (S_N orbit + 2-qubit block algebra), NOT Pi2-anchored. Only the time coefficients inherit.

In the Schicht-0 typed knowledge layer, F89 now cites `Pi2DyadicLadderClaim` as a constructor-injected parent and exposes both inherited coefficients as live properties:
- `F89TopologyOrbitClosure.DecayRateCoefficient` = `Pi2DyadicLadderClaim.Term(−1)` = 4 (γ-axis)
- `F89TopologyOrbitClosure.OscillationFrequencyCoefficient` = `Pi2DyadicLadderClaim.Term(−1)` = 4 (J-axis mirror)

The drift check `F89TopologyOrbitClosure.Pi2DoublingConsistent()` verifies both equal 4 (catches drift between the literal 4.0 in `AllIsolatedClosedForm` and the ladder anchor). The Schicht-1 `RegisterF89_AncestorsContainAllParents` test asserts the runtime ancestor walk now reaches `Pi2DyadicLadderClaim` from F89.

### Other topology classes (Tier 2 empirical, derivation open)

The above closed form applies only to the all-isolated subclass (1)^m. For mixed topologies (e.g. (1, 2), (2, 2), (1, 1, 2), and pure paths (3), (4), (5), (6)) the closed form is open. The empirical late-tail behavior:
- 1-isolated-edge mixed classes ((2), (1, 2), (1, 1, 2)) cluster at intermediate prefactors at the t ≈ π/(2J) in-phase moment.
- Pure-path classes ((3), (4), (5), (6)) decay faster than 4γ₀ on visible time scales due to no-overlap-SE-DE content (rate 6γ₀ on those coherences) plus longer-path mode mixing.

A general closed form per topology class is open work.

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
- For (k+1)-qubit blocks with k ≥ 2: no-overlap SE-DE pairs DO exist (e.g. SE at site 0, DE at sites {1, 2} on 3 qubits). → mixed dephasing rates {2γ₀, 6γ₀} on the sub-block diagonal, coupled by H_B's off-diagonal mixing → eigenvalues land in [2γ₀, 6γ₀] continuum.

**Numerical verification** (script `_f89c_liouvillian_eigenstructure.py`):

| Topology | Distinct decay rates Γ/γ |
|---|---|
| 2-qubit block (k=1) | {0, 2, 4} (three clean rates) |
| 3-qubit block (k=2) | {0, 2, 2.556, 2.889, 3.112, 3.444, 4, 6} (eight rates, including fractional) |
| 4-qubit block (k=3) | 25 distinct rates including many fractional |

For path-2 (3 qubits) the fractional rates pair as (2.556, 3.444) and (2.889, 3.111), each summing to 6γ₀. These are eigenvalues of 2×2 sub-blocks of L_super^{(SE, DE)} where two basis pairs of mixed dephasing rates are coupled by H_B; the eigenvalues are -3γ₀ ± |c| with H_B coupling |c|.

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

**Tier 1 derived** for the all-isolated (1)^m closed form S_(1)^m, N(t) = [(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))] · exp(−4γ₀ t). The derivation factors the Lindbladian over disjoint blocks plus bare sites, uses H_B-eigenstate phase tracking, and counts populated coherence sectors per block. Numerical verification matches the (1), (1, 1), (1, 1, 1) CSVs at N = 7 within CSV write precision (5e−7).

The **mixed-topology and pure-path closed forms** (per-class S(t) for (1, 2), (2, 2), (1, 1, 2), (3), (4), (5), (6) at N = 7) remain **Tier 2 empirical**; derivation open.

---

*Bond position is invisible to S_N-symmetric probes; the bond-graph topology class is the only spatial feature of a uniform-J multi-bond Hamiltonian that ρ_cc + spatial-sum-S(t) can resolve.*
