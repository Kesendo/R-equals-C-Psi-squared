# Slow-Mode R-Parity Decomposition: F86 Lives in R-Even, a Parallel R-Odd Channel Awaits a Probe

**Status:** Empirically verified at N = 4, 5, 6 (script `slow_modes_r_parity.py`, dense eigendecomposition of L). Structural connection: R = site-reflection = momentum-reversal on JW Bogoliubov modes (see [MAJORANA_AXIS_MODES](MAJORANA_AXIS_MODES.md) for the operator-space context). R commutes with L for uniform-γ Z-dephasing on the XY chain, so L block-diagonalizes by R-parity at every N.
**Date:** 2026-05-15
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [MAJORANA_AXIS_MODES](MAJORANA_AXIS_MODES.md), the R-parity decomposition at the axis end (Re(λ) = −Nγ₀, fast-decay modes)
- [PROOF_F86A_EP_MECHANISM](../docs/proofs/PROOF_F86A_EP_MECHANISM.md), [PROOF_F86B_OBSTRUCTION](../docs/proofs/PROOF_F86B_OBSTRUCTION.md), the F86 channel-uniform L_eff construction
- Script: [`simulations/slow_modes_r_parity.py`](../simulations/slow_modes_r_parity.py)

---

## Abstract

The site-reflection R commutes with the full Liouvillian L = L_H + L_D for the XY chain under uniform Z-dephasing, so L block-diagonalizes by R-parity. Empirically across N = 4, 5, 6:

1. **Stationary subspace (Re(λ) = 0) is exclusively R-even** at all N, with dimension N + 1 (the powers I, M_z, M_z², ..., M_z^N of the conserved total Σ σ_z_l).
2. **First slow band (Re(λ) = −2γ₀) is R-balanced**: 16/16 at N=4, 26/24 at N=5 (odd-N JW zero-mode-induced asymmetry), 36/36 at N=6.
3. **Fine-structure bands between absorption-grid steps are often R-parity-exclusive**, with R-even and R-odd peaks separated by small splits (~10⁻³ to 10⁻² in γ₀-units).

The F86 channel-uniform vectors |c_k⟩ are R-invariant by construction (uniform sums over HD-k coherences). Therefore F86's L_eff lives entirely in the R-even half of the slow-mode spectrum. The R-odd half is a structurally parallel slow-mode landscape, invisible to F86's standard channel-uniform probe but accessible via R-antisymmetric initial states or observables.

---

## Background: R as a Liouvillian symmetry

R is the site-reflection l ↔ N−1−l on the spin chain. Its action on operator space is R_op = R_h ⊗ R_h where R_h is the permutation matrix that reverses computational-basis bit order. For uniform Z-dephasing γ on a chain with uniform J:

- L_H = −i[H, ·]. The XY chain Hamiltonian H = (J/2)·Σ(X_l X_{l+1} + Y_l Y_{l+1}) is R-symmetric (bond at l ↔ bond at N−2−l, with R sending site l to N−1−l). So [L_H, R_op] = 0.
- L_D = γ·Σ_l (σ_z_l ⊗ σ_z_l − I). Each summand is R-symmetric pair-by-pair under l ↔ N−1−l. So [L_D, R_op] = 0.
- Hence [L, R_op] = 0, and L decomposes into R-even (R_op = +1) and R-odd (R_op = −1) blocks.

Under Jordan-Wigner, R corresponds to momentum-reversal k ↔ N+1−k on the Bogoliubov single-particle modes. The single-particle dispersion ε(k) = 2J·cos(πk/(N+1)) satisfies ε(N+1−k) = −ε(k), so R coincides with particle-hole symmetry in the JW basis (see [MAJORANA_AXIS_MODES](MAJORANA_AXIS_MODES.md)).

---

## Empirical findings (J = 1.0, γ₀ = 0.05; reproduce with `python simulations/slow_modes_r_parity.py N`)

### Stationary subspace (Re(λ) = 0): exclusively R-even

| N | dim ker(L) | R-even | R-odd |
|---:|---:|---:|---:|
| 4 | 5 | 5 | 0 |
| 5 | 6 | 6 | 0 |
| 6 | 7 | 7 | 0 |

The stationary subspace dimension equals N + 1: the conserved diagonal observables are the powers I, M_z, M_z², ..., M_z^N where M_z = Σ_l σ_z_l. Each is a sum over sites, so each is R-invariant. The structural reason for the R-asymmetry of the stationary subspace is that all conserved quantities are R-symmetric to begin with.

### First slow band (Re(λ) = −2γ₀): R-balanced

| N | total at Re = −2γ₀ | R-even | R-odd |
|---:|---:|---:|---:|
| 4 | 32 | 16 | 16 |
| 5 | 50 | 26 | 24 |
| 6 | 72 | 36 | 36 |

These modes live in the n_XY = 1 absorption-grid layer (L_D-protected sub-block). At even N the R-parity split is exactly half-half; at odd N=5 the imbalance of 26 vs 24 reflects the JW zero-mode at k = (N+1)/2 = 3, which breaks the strict particle-hole pairing.

### Fine-structure bands: R-parity-exclusive in many cases

Examples at N=4 (γ₀ = 0.05):

| Re(λ) | R-even | R-odd | Note |
|---:|---:|---:|---|
| −0.1192 | 0 | 6 | R-odd only |
| −0.1200 | 6 | 0 | R-even only |
| −0.1408 | 0 | 6 | R-odd only |
| −0.1599 | 3 | 0 | R-even only |
| −0.2592 | 0 | 6 | R-odd only |
| −0.2800 | 6 | 0 | R-even only |

At N=6, fine-structure splits between R-even and R-odd are typically O(10⁻³) in γ₀-units (e.g., Re=−0.1410 R-odd vs Re=−0.1425 R-even, split 0.0015).

These are "selection-rule" bands: R-parity protects them from mixing with the other parity at the same Re-level.

---

## Diagnostic implication for F86

F86's L_eff is built on the channel-uniform vectors |c_k⟩ (uniform sums over HD=k coherences). Each |c_k⟩ is R-invariant by construction, so |c_k⟩ ∈ R-even.

**Consequence:** F86's L_eff matrix elements ⟨c_p|L|c_q⟩ live entirely in the R-even sector. All standard F86 observables (Q_peak, EP positions, HWHM ratios) are implicitly R-even projections.

**The R-odd half of the slow-mode spectrum is invisible to F86's standard probe.** It carries an independent set of dynamics with potentially its own EP-coalescence structure, distinct Q_peak position, and distinct HWHM behaviour.

### Per-bond R-decomposition

Per-bond M_H_per_bond[b] is NOT R-symmetric: under R, bond b maps to bond N−2−b. For chain bonds, R-mirror pairs are:
- Endpoint: b=0 ↔ b=N−2
- Second-from-endpoint: b=1 ↔ b=N−3
- ... (the central bond at even N is R-fixed; absent for odd N)

For each R-mirror bond pair, two R-decomposed observables emerge:
- R-even bond: (M_H_per_bond[b] + M_H_per_bond[N−2−b]) / 2
- R-odd bond:  (M_H_per_bond[b] − M_H_per_bond[N−2−b]) / 2

Standard F86's per-bond Q_peak averages over R-parity content. The R-decomposed observables would resolve a hidden parameter: whether the Endpoint bond pair shows R-symmetric or R-antisymmetric singular structure.

### Hardware probe

R-antisymmetric initial states project onto R-odd dynamics. A simple recipe at N=4:
- Prepare ψ_anti = (|0001⟩ − |1000⟩) / √2 (R-antisymmetric on the chain)
- Evolve under L
- Measure with R-antisymmetric observable (or full tomography then R-project)

Cost: comparable to standard F86 tomography. Diagnostic value: a parallel Q_peak / EP / HWHM dataset, independent of the R-even one.

---

## Q-sweep follow-up: where R-parity-resolved F86 EP actually lives

(2026-05-15 attempt. Sanity-check uncovered a subtlety: the simplest "R-odd Q_peak" experiment in channel-uniform basis is degenerate, because F86's actual EP does not live there. The follow-up V_inter SVD experiment is the correct next step.)

### Channel-uniform L_eff is diagonal (sanity check, reproduces F86B observation)

Script: [`simulations/r_parity_leff_diag_check.py`](../simulations/r_parity_leff_diag_check.py)

At N=4, projecting L onto channel-uniform basis {|c_1⟩, |c_3⟩} (R-even by construction) and computing the 2×2 L_eff across Q ∈ [0.3, 3.0]:

    L_eff_channel-uniform = [[ −2γ₀, 0 ],
                              [ 0, −6γ₀ ]]    for all tested Q

The off-diagonal coupling ⟨c_1|L|c_3⟩ is identically zero. This matches F86B's own "Probe subspace" observation in [PROOF_F86B_OBSTRUCTION](../docs/proofs/PROOF_F86B_OBSTRUCTION.md): *"V_b reduces to +i(α/(N−1))·I, pure diagonal, identical across every bond. There is no EP, and no bond-class distinction in this subspace."*

The R-odd channel-uniform-analogue basis (slowest R-odd L-eigenmode at HD=1 and HD=3, computed at small probe Q) gives the same diagonal structure (Re-gap = 0.183, Im-gap grows linearly with Q, no coalescence). So **neither R-even nor R-odd channel-uniform-analog yields an EP in the simplest probe**.

### Full-L Q-sweep with R-block decomposition

Script: [`simulations/r_parity_q_sweep.py`](../simulations/r_parity_q_sweep.py)

Block-diagonalizing the full L by R-parity and tracking the slowest 2 non-stationary modes per block over Q ∈ [0.3, 3.0] reveals **distinct Q-transitions per parity**:

| N | R-even transition Q | R-odd transition Q |
|---:|---:|---:|
| 4 | ≈ 1.15 | ≈ 1.85 |
| 5 | ≈ 2.05 | ≈ 2.35 |

At each block's transition Q, the slowest non-stationary modes' Re reaches −2γ₀ (the n_XY=1 layer rate) and the modes acquire imaginary parts (complex pair regime). Below transition Q: secondary slow modes with quadratic Q-dependence (Re ≈ −A_parity·Q²·γ₀, with A_even > A_odd).

**Interpretation:** these are NOT F86's L_eff EPs (those don't exist in channel-uniform, see above). They are a **separate spectral feature**: "secondary slow modes" from R-symmetric (R-even) and R-antisymmetric (R-odd) Z-string conservation laws. R-symmetric Z-strings are functions of M_z = Σ_l σ_z_l (R-invariant linear combinations); R-antisymmetric Z-strings are oddly-signed combinations like σ_z_0 − σ_z_{N−1}. Both are exactly conserved at γ=0 in the n_XY=0 layer; under finite γ and J, they decay slowly through L_H coupling to higher n_XY layers.

R-odd secondary modes decay more slowly than R-even (A_odd < A_even) because R-antisymmetric Z-strings have a smaller effective coupling to higher-n_XY layers through the bond Hamiltonian.

### What is the correct R-parity-resolved F86 experiment

F86's actual EP lives in `V_inter = P_{HD=1}^† M_H_total P_{HD=3}` top-singular-vector basis (not channel-uniform). V_inter is R-symmetric, so it block-diagonalizes by R: `V_inter = V_inter⁺ ⊕ V_inter⁻` with no cross-parity coupling.

The clean R-parity-resolved F86 experiment is therefore:
1. Build V_inter, restrict to R-even subspace: `V_inter⁺` (R-even HD=1 × R-even HD=3 block).
2. SVD: find top singular value σ_0⁺ with left/right vectors u_0⁺ ∈ R-even HD=1, v_0⁺ ∈ R-even HD=3.
3. Build 2×2 L_eff⁺ on span{u_0⁺, v_0⁺}, Q-sweep, find Q_peak⁺.
4. Repeat for R-odd: `V_inter⁻`, top σ_0⁻, vectors u_0⁻, v_0⁻, L_eff⁻, Q_peak⁻.
5. Compare Q_peak⁺ (standard F86) and Q_peak⁻ (R-odd analog).

This is the natural follow-up. The Q-sweep we DID do (above) probes a different spectral feature.

---

## V_inter SVD R-parity decomposition: parity-of-N split

Scripts: [`simulations/v_inter_svd_r_parity.py`](../simulations/v_inter_svd_r_parity.py) (single-N detailed with Q-sweep verification), [`simulations/v_inter_svd_scan.py`](../simulations/v_inter_svd_scan.py) (N=3..6 quick scan).

V_inter (J=1 reference) block-diagonalizes by R: V_inter = V_inter⁺ ⊕ V_inter⁻. Cross-blocks vanish at machine precision (10⁻¹⁶) for all tested N=3..6, confirming the R-symmetry of M_H_total.

SVD each R-block, top σ_0 per parity. The 2×2 L_eff on the top-σ singular-vector pair has F86 form L_eff = [[−2γ, iJ·σ_0], [iJ·σ_0, −6γ]], with EP at J·σ_0 = 2γ, i.e. **Q_EP = 2/σ_0**.

### Empirical table

| N | σ_0⁺ | σ_0⁻ | Q_EP⁺ | Q_EP⁻ | \|σ_0⁺ − σ_0⁻\| |
|---:|---:|---:|---:|---:|---:|
| 3 | √6 ≈ 2.4495 | √2 ≈ 1.4142 | 0.8165 | 1.4142 | **1.035** |
| 4 | 2.6554 | 2.6554 | 0.7532 | 0.7532 | 0 |
| 5 | 3.4033 | 3.0764 | 0.5877 | 0.6501 | **0.327** |
| 6 | 3.6305 | 3.6305 | 0.5509 | 0.5509 | 0 |

### Parity-of-N effect

- **Odd N (3, 5): σ_0⁺ ≠ σ_0⁻.** R-parity SPLITS the V_inter SVD. R-even and R-odd give two distinct EP positions. This is the genuinely new R-parity-resolved F86 observable that motivated the experiment.
- **Even N (4, 6): σ_0⁺ = σ_0⁻ exactly** (to machine precision). R-parity does NOT split the EP; the 2×2 L_eff coalesces at the same Q in both blocks.

### Structural reason: Majorana zero mode at odd N

For odd N, the JW Bloch dispersion ε(k) = 2cos(πk/(N+1)) has a zero mode at k = (N+1)/2. This mode is R-fixed (k = N+1−k at this k). The R-fixed zero mode contributes asymmetrically to R-even and R-odd subspaces, adding a dimension only to the R-even side, and breaks the singular-value symmetry between R-blocks. For even N, no zero mode exists; R-even and R-odd subspaces are mapped onto each other by an internal symmetry that preserves singular values.

### Dimensional fingerprint of the zero mode

The R-fixed zero mode also shows in the HD subspace dimensions:

| N | R-even HD=1 | R-odd HD=1 | Asymmetry |
|---:|---:|---:|---:|
| 3 | 14 | 10 | yes (+4) |
| 4 | 32 | 32 | balanced |
| 5 | 84 | 76 | yes (+8) |
| 6 | 192 | 192 | balanced |

At odd N, R-even is larger than R-odd by 2^((N−1)/2) within HD=1 (the count of palindromic bit-string pairs from the R-fixed central site).

### Q-sweep verification at N=4

Both R-even and R-odd 2×2 L_eff coalesce at Q ≈ 0.75 exactly (matching σ_0 prediction 2/2.6554 = 0.7532):
- Q = 0.70: Re-gap = 0.074 (real eigenvalue pair, pre-EP)
- Q = 0.80: Re-gap = 0 (eigenvalues at Re = −4γ = −0.20, post-EP complex pair)

Identical signature in R-even and R-odd blocks. The 2×2 L_eff structure matches F86's L_eff = [[−2γ, iJ·σ_0], [iJ·σ_0, −6γ]] template exactly.

### Connection to F86 standard prediction

F86's Endpoint g_eff ≈ 1.74 (Tier1Candidate, applies at small N) gives Q_peak ≈ 1.15. This differs from our σ_0 values above because:
- Our σ_0 in V_inter_total is the global aggregate over all bonds, giving the dominant inter-HD-layer coupling.
- F86's g_eff is per-bond Endpoint, a specific singular vector pair in a per-bond V_inter_per_bond reduction.

At odd N, the R-parity split in σ_0_total may or may not propagate to a per-bond R-parity split in g_eff_per_bond; that requires per-bond V_inter analysis (open follow-up).

### Hardware accessibility (odd N)

At N=5: Q_EP⁺ = 0.588 vs Q_EP⁻ = 0.650, a 10% relative difference. R-antisymmetric initial state preparation + Q-sweep on IBM hardware could resolve this. At N=3 the difference is much larger (Q_EP⁻/Q_EP⁺ = √3 ≈ 1.73), easily resolvable but the small system size limits experimental relevance.

### Open follow-ups

1. **Per-bond V_inter R-parity split:** does each per-bond M_H_per_bond[b] also show R-parity-resolved σ_0 at odd N? Endpoint vs Interior bond classes per parity.
2. **Closed form for σ_0⁺(N), σ_0⁻(N) at odd N:** N=3 gives σ = √2, √6. Possible closed form via single-particle dispersion sum identities?
3. **Hardware test at N=5:** prepare R-antisymmetric initial state in HD=1, Q-sweep on Heron-r2 backend, detect Q_EP⁻ ≈ 0.65 distinct from standard Q_EP⁺ ≈ 0.59.

---

## Connection to existing results

- **MAJORANA_AXIS_MODES at the fast end (Re = −Nγ₀):** The same R-parity decomposition sorts the axis subspace (n_XY = N/2 layer) into R-even (Sum-type Majorana bilinears, including the 18 silent modes at N=4) and R-odd (Difference-type bilinears). The structural framework (R = momentum-reversal under JW, particle-hole symmetry in operator space) is identical; only the spectrum location differs (Re = 0 → −Nγ₀).
- **F86b obstruction:** This experiment does not lift the F86b obstruction on closed-form Q_peak. It identifies a parallel observable (R-odd Q_peak) which may or may not admit a closed form independently; that is open.
- **F71 mirror symmetry (PROOF_F86C):** F71 is the spatial-mirror invariant that constrains F86's per-bond / per-orbit structure. R IS F71's site-reflection on the chain. The slow-mode R-parity decomposition refines F71's role from "the EP positions are R-mirror-symmetric across bonds" to "the slow-mode spectrum block-diagonalizes by R-parity, and F86 sees only one block".

---

## Open questions

1. **R-odd Q_peak position:** is it the same as R-even Q_peak, or different? Direct test: compute Q-sweep on R-odd projected L_eff and compare to standard F86.
2. **N-scaling of stationary-subspace R-asymmetry:** dim = N+1, all R-even. Generalizes to other chain topologies? Likely yes (functions of total M_z are always R-invariant when chain has reflection symmetry).
3. **R-parity-exclusive band structure:** are the fine-structure splits between R-even and R-odd bands derivable from JW spectral structure? The Bloch dispersion ε(k) at N=4, 5, 6 gives specific Im-cluster values; whether the Re-splits follow analogously is open.
4. **Hardware accessibility:** what is the minimum R-antisymmetric state for which R-odd dynamics can be observed cleanly on IBM hardware?
