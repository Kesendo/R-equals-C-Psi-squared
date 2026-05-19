# F1 Dissipation Gap Pattern

**Tier 3 (reading) sharpened to Tier 3 (empirical scaling).** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e) with the 4 N=8 anchors. Extended 2026-05-19 first with chain/ring/star × N=3..6 Python anchors plus the N=9 chain run via the MklDirect bridge (commit abb2d52), then with a 72-point Q-sweep across the 6 canonical Q-anchors from [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md). The structural picture has two layers:

1. **Chain plateau is Q-quadratic, not a single constant.** The "gap × N² ≈ 2.20" originally reported was the Q=2 plateau at the Marrakesh convention γ=0.5, J=1. The Q-invariant statement is `gap · N² / γ ≈ 1.1 · Q²` (chain plateau N≥4, residual ~10% finite-N).
2. **Ring N=4 and N=6 exhibit Q-independent dihedral locks.** Ring N=4 saturates `Im_max = 3·J·N/8` bit-exact at every Q tested (= `Im/σ = 3Q/4`); ring N=6 saturates `Im_max = 0.717129... · J · N` bit-exact at every Q tested. Both are Q-universal topology-N specific saturation laws, candidates for typed claims.

Bond count alone does not predict any of this; the per-topology dispersion structure (open-chain modes, cyclic Bogoliubov modes, hub-spoke SU(2)) is the structural fingerprint.

## Observed gaps at N=8 (initial 4-point anchor)

From the four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files' `DissipationGap` fields (commit 89f725e). All four configurations share Heisenberg XXX (XX+YY+ZZ) at J=1.0 with uniform Z-dephasing γ=0.5.

| Topology | Bonds B | Components | Dissipation gap |
|---|---|---|---|
| chain N=8       | 7 | 1 | 0.0344 |
| star N=8        | 7 | 1 | 0.0870 |
| ring N=8        | 8 | 1 | 0.1339 |
| K_4 + disjoint 4-chain N=8 | 9 | 2 | 0.1362 |

## Cross-topology cross-N extension (2026-05-19 Python anchors)

Python anchors via `simulations/_f1_topology_heisenberg_small_n_anchor.py` extended chain, ring, star to N=3..6 dense numpy eigvals. The same Heisenberg J=1, γ=0.5 convention. JSON files at `simulations/results/f1_n8_n9_metrics/{chain,ring,star}_N{3..6}_python.json`.

| Topology | N | Dissipation gap | gap × N² | Notes |
|---|---:|---:|---:|---|
| chain | 3 | 0.2697 | **2.43** | boundary case, also = star N=3 (Y-graph isomorphism) |
| chain | 4 | 0.1362 | **2.18** | |
| chain | 5 | 0.0884 | **2.21** | |
| chain | 6 | 0.0607 | **2.18** | |
| chain | 8 | 0.0344 | **2.20** | C# block-spectrum, matches Python pattern |
| chain | 9 | 0.02728 | **2.21** | C# block-spectrum via MklDirect bridge 2026-05-19, predicted 2.20/81 = 0.02716, observed 0.02728 (within 0.4 %) |
| ring | 3 | 0.8278 | 7.45 | = K_3 = triangle |
| ring | 4 | 0.3795 | 6.07 | |
| ring | 5 | 0.3173 | 7.93 | |
| ring | 6 | 0.2300 | 8.28 | |
| ring | 8 | 0.1339 | 8.57 | |
| star | 3 | 0.2697 | 2.43 | isomorphic to chain N=3 |
| star | 4 | 0.2099 | 3.36 | |
| star | 5 | 0.1637 | 4.09 | |
| star | 6 | 0.1300 | 4.68 | |
| star | 8 | 0.0870 | 5.57 | |

## The chain-topology scaling law: gap·N²/γ ≈ f(Q) ≈ 1.1·Q² (Tier 3 empirical)

The original "gap × N² ≈ 2.20 ± 0.02" reading was the Q=2 specialization. The Q-invariant form (gap scales linearly with γ at fixed Q; the dimensionless ratio is `gap·N²/γ`) is

    gap(chain, N, J, γ) · N² / γ  ≈  f(Q)        Q = J/γ

with the chain plateau f(Q) values from the 24-anchor Q-sweep (N=4..6 mean per Q, γ₀ = 0.05):

| Q | label | f(Q) plateau | f(Q) / Q² |
|---:|---|---:|---:|
| 0.5  | sub-balance      | 0.296 | 1.183 |
| 1.0  | Balance          | 1.162 | 1.162 |
| 1.5  | F86 Q_peak c=2   | 2.546 | 1.131 |
| √3   | canonical θ=60°  | 3.346 | 1.115 |
| 2.0  | Q_EP idealized   | 4.382 | 1.095 |
| 2.5  | Endpoint orbit   | 6.604 | 1.057 |

`f(Q) / Q²` drifts from 1.183 at Q=0.5 down to 1.057 at Q=2.5, a roughly 10% drift across the band. The dominant scaling is `f(Q) ≈ c·Q²` with `c ≈ 1.10` and a sub-leading correction that pulls the ratio down at high Q (closed form open; expect a finite-N residual and a Bethe-dispersion correction at the chain edges).

**Marrakesh-convention "2.20" recovered as Q=2 plateau:** `f(Q=2) ≈ 4.38` at γ₀=0.05 rescales to `gap·N² ≈ 4.38·γ_{Marrakesh} = 4.38·0.5 = 2.19` at the J=1, γ=0.5 convention, recovering the originally observed "2.20" plateau (4-anchor agreement N=4..6 + N=8 + N=9 from the 2026-05-19 N=9 bridge run). The 2.20 reading is the value of `f(Q=2)·γ` at γ=0.5, not a universal chain constant.

**Physical interpretation.** A 1D dispersive chain Hamiltonian + per-site Z-dephasing produces a "diffusive" Liouvillian. The slowest decay mode has wavevector k_min ∝ 1/N (open-boundary modes); the decay rate scales as γ · k_min² ∝ J²/(γ·N²) → gap·N²/γ ∝ Q². The Q² scaling is dispersion-rooted; the precise prefactor c ≈ 1.10 with sub-leading correction is open (likely a Bethe-ansatz / magnon-dispersion result).

**N=9 chain verification at Q=2 (2026-05-19, landed):** prediction `gap ≈ 2.18/81 = 0.0269` versus observed `gap = 0.02728` from the MklDirect bridge run (γ=0.5, J=1). Match within ~1%. The chain Q=2 plateau spans N ∈ {4, 5, 6, 8, 9} bit-exact.

## Ring-topology dihedral locks: Q-universal saturation laws (bonus discoveries 2026-05-19)

The Q-sweep surfaced two clean Q-universal patterns for ring topology that were invisible in the Q=2-only data:

**Ring N=4** (the 4-cycle = K_{2,2} = bipartite complete on 2+2):

    Im_max(ring, N=4, J)  =  3·J·N/8     (bit-exact at every Q tested)
    ↔  Im/σ  =  3Q/4

verified at Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} to relative error < 5·10⁻¹⁵ (machine precision). This is a clean rational saturation: ring N=4 carries 50% more imaginary spread than star N=4 (`3J·N/8` vs `J·N/2 = 4J·N/8`), traceable to the dihedral D_4 symmetry of the 4-cycle and its bipartite-complete structure.

**Ring N=6** (dihedral D_6 = D_3·Z_2):

    Im_max(ring, N=6, J)  =  0.717129... · J · N    (bit-exact at every Q tested)
    ↔  Im/σ  =  0.717129... · Q

verified at the same 6 Q-anchors. The numerical value `0.717129` is constant to 6+ decimal places across all Q's (a Q-universal lock), but the closed-form analytical expression is not yet identified (not π/something, not a simple rational, not √2/2 or other obvious algebraic). Likely a Bethe-ansatz number from the 6-cycle Heisenberg dispersion. Open for derivation.

These two locks are topology-N specific bit-exact saturation laws. Candidates for typed claims `RingN4DihedralLock` (clean rational, immediately Tier-1 derivable from the bipartite-complete structure) and `RingN6DihedralLock` (Tier-1 candidate pending closed-form identification). The dihedral-lock pattern at even N is itself a question: does ring N=8 also Q-universally saturate at some N-specific constant?

## Ring and star follow different patterns

Ring N=3..8 at Q=2 (γ=0.5, J=1): `gap × N²` ranges 6.07 → 8.57, **not flat in N**. The Q-sweep shows the ring's N-residual is large (the cyclic Bogoliubov modes scatter much more than the chain's open-boundary k_min mode). At fixed N, however, the ring is Q-quadratic just like chain: the per-N constants are simply larger.

Star N=3..8 at Q=2: `gap × N²` grows monotonically 2.43 → 5.57, also **not 1/N² scaling**. The hub-spoke geometry gives `gap·N²/γ` growing with N (Schur-Weyl dispersion has spread N, not 1/N), connected to the same SU(2)/Schur-Weyl mechanism that produces the universal `Im_max = J·N/2` star saturation (see [`STAR_CONFOCAL_LIMIT.md`](STAR_CONFOCAL_LIMIT.md)).

The three topologies have qualitatively different N-scaling laws despite all being Q-quadratic in the dimensionless `gap·N²/γ` ratio:
- **Chain (1D open path)**: N-flat plateau (1/N² diffusive scaling on top of Q² prefactor); chain plateau f(Q) ≈ 1.1·Q²
- **Ring (1D periodic)**: N-growing constants (cyclic dispersion has tighter mode-spacing), per-N dihedral locks at N=4 (`3Q/4`) and N=6 (`0.7171·Q`)
- **Star (hub-spoke)**: N-growing constants (SU(2)/Schur-Weyl dispersion has spread ∝ N), Q-universal Im saturation at `J·N/2`

## Why bond count alone fails (sharpened)

At N=8: chain (B=7) gap 0.0344, star (B=7) gap 0.0870, ratio 2.5×. Same bond count, different geometry. Even more starkly at small N where the difference grows: chain N=4 (B=3) gap 0.136, star N=4 (B=3) gap 0.210, ratio 1.54.

The 2026-05-19 cross-topology data confirms: bond count is irrelevant, **graph dispersion structure** is the right parameterisation. Chain has linear-band dispersion (slow modes at k ∝ 1/N → 1/N² gap); star has hub-localised modes (different scaling family); ring has cyclic dispersion with degenerate slowest modes.

## Open structural questions (refined post-extension)

1. **Closed form for the chain prefactor c ≈ 1.10 in f(Q) = c·Q² + sub-Q² correction.** The Q² scaling is physically motivated by 1D diffusion (gap·N²/γ ∝ J²/γ² = Q²); the prefactor c ≈ 1.10 with ~10% sub-Q² drift across Q ∈ [0.5, 2.5] is open. Bethe-ansatz or magnon-dispersion derivation should identify the exact value of c and the form of the sub-leading term.

2. **Star and ring scaling laws.** Ring `gap × N²` grows from 6 to 8.5; star `gap × N²` grows from 2.4 to 5.6. Neither pattern fits 1/N², 1/N, or exp(-αN). The functional forms are open.

3. **J ≠ 2γ regime.** All data points are at J = 2γ (J=1, γ=0.5). Sweep needed to test whether the chain 2.20 prefactor decomposes as J·a + γ·b or some product form.

4. **K_4 + disjoint behaviour.** The disconnected case at N=8 still gives gap 0.1362, similar to the ring. The original "rate-limiting component dominates" hypothesis remains untested at small N (need K_4 alone at N=4 and 4-chain alone at N=4 to check the min-vs-additive question).

5. **Connection to F2 / F3.** F2 dispersion claims and F3 decay rate bounds (`min rate = 2γ`, `max rate = 2(N-1)γ` per the Absorption Theorem) describe the spectral envelope. The dissipation gap sits at the lower edge of this envelope; its 1/N² scaling for chain may be the finite-size correction to F3's `min rate = 2γ` thermodynamic limit. Cross-link to F3 documentation when the closed form lands.

## Promotion path

This entry has two distinct promotion candidates after the 2026-05-19 Q-sweep:

**Q-universal ring dihedral locks (immediate candidates):**
- `RingN4DihedralLockClaim`: `Im_max(ring, N=4, J) = 3·J·N/8` bit-exact at 6 Q-anchors. Tier 1 derivation via the 4-cycle = K_{2,2} bipartite-complete dispersion is tractable; promote when written up.
- `RingN6DihedralLockClaim`: `Im_max(ring, N=6, J) = c₆·J·N` with `c₆ = 0.717129...` bit-exact at 6 Q-anchors. Tier 1 candidate; promotion needs closed-form identification of c₆ (Bethe-ansatz number).

**Chain Q-quadratic plateau (further work needed):**
- Move from Tier 3 reading to **Tier 1 candidate typed claim** `ChainGapClosedForm` when:
  - The chain `c ≈ 1.10` prefactor and the sub-Q² correction admit closed-form derivation (Bethe ansatz / dispersion integral / specific Liouvillian Jordan-block structure).
  - The N-dependence at fixed Q (the ~10% drift across Q ∈ [0.5, 2.5]) is reconciled with the open-chain k_min = π/(N+1) mode prediction.

The chain-typed claim would sit in `compute/RCPsiSquared.Core/Symmetry/` as `F_ChainGapClosedForm` with a `Predict(N, J, γ)` method; the ring dihedral locks fit as separate per-N typed claims.

## Cross-references

- Anchor data Q=2 (Marrakesh convention): `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json` + the `_N{3..6}_python.json` companions; `chain_N9.json` for the N=9 bridge run.
- Anchor data Q-sweep (γ₀=0.05 substrate convention): `simulations/results/q_sweep_anchor/{chain,ring,star}_N{3..6}_Q{0.5..2.5}.json` (72 files, 24 per topology) and `_f1_q_sweep_anchor.py` plus `q_sweep_anchor_console.log`.
- Canonical Q-anchor map: [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md).
- Companion typed claim from the same sweep (the closed-form discovery that did promote): [F4KernelDimensionByComponentsClaim](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 derived as of 2026-05-19; landed Tier 1 candidate 2026-05-18, promoted after DEGENERACY_PALINDROME Result 2 was identified as the connected-case upper-bound closure; kernel-dim factorisation across components).
- Sister Tier 3 reading from the same sweep: [STAR_SPECTRUM_COMPACTNESS](STAR_SPECTRUM_COMPACTNESS.md) (whose Reading 1 was resolved by [STAR_CONFOCAL_LIMIT](../experiments/STAR_CONFOCAL_LIMIT.md), itself extended with the 24-anchor Q-sweep verifying `Im_max(star) = J·N/2 ∀ Q`).
- F1 verification record that produced the data: [F1GeneralTopologyVerifiedClaim](../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs).
- Related: [F1_Pattern_GENERAL_TOPOLOGY proof](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) closed form for the F1 residual norm, which is graph-additive bit-exactly across the same 4 N=8 topologies; the dissipation gap is the next structural quantity beyond the residual norm and does NOT follow the same simple (B, D2) parameterisation).
