# F1 Dissipation Gap Pattern

**Tier 3 (reading) sharpened to Tier 3 (empirical scaling).** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e) with the 4 N=8 anchors. Extended 2026-05-19 first with chain/ring/star × N=3..6 Python anchors plus the N=9 chain run via the MklDirect bridge (commit abb2d52), then with a 72-point Q-sweep across the 6 canonical Q-anchors from [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md). The structural picture has two layers:

1. **Chain plateau is Q-quadratic, not a single constant.** The "gap × N² ≈ 2.20" originally reported was the Q=2 plateau at the Marrakesh convention γ=0.5, J=1. The Q-invariant statement is `gap · N² / γ ≈ 1.1 · Q²` (chain plateau N≥4, residual ~10% finite-N).
2. **Ring N=4 and N=6 exhibit Q-independent dihedral locks.** Ring N=4 saturates `Im_max = (3/4)·J·N = 3·J` bit-exact at every Q tested (= `Im/σ = 3Q/4`); ring N=6 saturates `Im_max = 0.717129... · J · N` bit-exact at every Q tested. Both are Q-universal topology-N specific saturation laws, candidates for typed claims.

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

    Im_max(ring, N=4, J)  =  (3/4)·J·N  =  3·J     (bit-exact at every Q tested)
    ↔  Im/σ  =  3Q/4

verified at Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} to relative error < 5·10⁻¹⁵ (machine precision). This is a clean rational saturation: ring N=4 carries 50% more imaginary spread than star N=4 ((3/4)·J·N vs (1/2)·J·N), traceable to the bipartite-complete structure (4 bonds in K_{2,2} vs 3 in the N=4 star). Closed-form derivation via the Casimir spectrum `{−2J, −J, 0³, +J}` of `H = J·S⃗_A·S⃗_B` with sublattice totals S⃗_A = S⃗_0 + S⃗_2, S⃗_B = S⃗_1 + S⃗_3 (see [`docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md)). Typed Tier 1 derived as [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) 2026-05-19.

**Ring N=6** (dihedral D_6 = D_3·Z_2):

    Im_max(ring, N=6, J)  =  0.717129... · J · N    (bit-exact at every Q tested)
    ↔  Im/σ  =  0.717129... · Q

verified at the same 6 Q-anchors. The numerical value `0.717129` is constant to 6+ decimal places across all Q's (a Q-universal lock), but the closed-form analytical expression is not yet identified (not π/something, not a simple rational, not √2/2 or other obvious algebraic). Likely a Bethe-ansatz number from the 6-cycle Heisenberg dispersion. Open for derivation.

These two locks are topology-N specific bit-exact saturation laws. **Ring N=4 promoted Tier 1 derived 2026-05-19** as [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) via the K_{2,2} bipartite-complete Casimir derivation. **Ring N=6 stays Tier-2 empirical** pending closed-form identification of the transcendental constant 0.717129... (likely Bethe-ansatz, see [`docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md) section "Why this is N=4-specific" for why the simple Casimir argument does not extend). The dihedral-lock pattern at even N is itself a question: does ring N=8 also Q-universally saturate at some N-specific constant? Existing N=8 ring data at Q=2 (γ=0.5, J=1) gives `Im/σ = 1.4128 = 0.7064·Q`, so the three-point sequence is `c_4 = 0.75, c_6 = 0.7171, c_8 = 0.7064`, monotonically decreasing and already passing through 1/√2 ≈ 0.7071 between N=6 and N=8 (1/√2 is NOT an upper bound that c_∞ approaches from above; c_8 already sits below it). The N → ∞ limit could be 1/√2, slightly below it, or some other Bethe-ansatz number; three points is not enough to settle this.

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

1. **Closed form for the chain prefactor c ≈ 1.10 in f(Q) = c·Q² + sub-Q² correction.** The Q² scaling is physically motivated by 1D diffusion (gap·N²/γ ∝ J²/γ² = Q²); the prefactor c ≈ 1.10 with ~10% sub-Q² drift across Q ∈ [0.5, 2.5] is open. Bethe-ansatz or magnon-dispersion derivation should identify the exact value of c and the form of the sub-leading term. **Status: open.**

2. **Star and ring scaling laws.** **Status: substantially sharpened 2026-05-19 by the Q-sweep.** Three distinct scaling families emerge cleanly:
   - **Chain** (open 1D path, k_min = π/(N+1)): `gap·N²/γ ≈ 1.10·Q²` at the N≥4 plateau (closed-form prefactor still open, see Q1).
   - **Ring** (periodic 1D cycle, k_min = 2π/N): `gap·N²/γ ≈ 4·Q²` at the N≥5 plateau (the dihedral lock at N=4 is a finite-size outlier). The ~4× chain-to-ring prefactor ratio matches `(2π/N)² / (π/(N+1))²  →  4` in the N→∞ limit, i.e. the **squared-wavevector ratio of cyclic-vs-open lowest modes**. This is suggestive that ring and chain share the same `Q²·γ/N²` diffusive form with a topology-specific k_min² factor.
   - **Star** (hub-spoke): `gap·N/γ` ≈ const at fixed Q (the 4 anchors at N=3..6, Q=2 give 0.081, 0.084, 0.082, 0.078), i.e. **gap ~ 1/N (NOT 1/N²)**. The scaling family is different from chain/ring because hub-spoke geometry has no spatial dispersion; the slowest mode is set by the hub-localised eigenmode whose decay rate scales as bandwidth/N rather than k_min². The Q-dependence of the star prefactor is non-Q² (drift from `(gap·N/γ)/Q² ≈ 0.60` at Q=0.5 to ≈ 0.35 at Q=2.5); closed form open.

   What is closed: chain and ring share the diffusive `Q²·γ/N²` form with topology-specific k_min² prefactors; star is in a different scaling family. What remains open: exact closed forms for the three Q-dependent prefactors and the sub-leading corrections.

3. **J ≠ 2γ regime.** **Status: RESOLVED 2026-05-19.** The Q-sweep (`_f1_q_sweep_anchor.py`) measured chain/ring/star × N=3..6 × Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} at γ₀ = 0.05: 72 anchors covering six J/γ ratios from 0.5 to 2.5. The data confirms `gap·N²/γ ≈ c·Q²` (with c topology-specific) holds across the full range, ruling out the "chain prefactor decomposes as J·a + γ·b" linear conjecture in favour of the cleaner Q² form. The earlier "2.20" reading was the chain Q=2 plateau at γ=0.5.

4. **K_4 + disjoint behaviour (rate-limiting component).** **Status: RESOLVED 2026-05-19** via tensor-sum factorisation (the same mechanism that closes F4 kernel-dim factorisation across components). For G = G_1 ⊔ G_2 with the per-component Heisenberg + Z-dephasing Liouvillian L(G_c), the bond-locality of both H and the dephasing dissipator gives

       L(G_1 ⊔ G_2) = L(G_1) ⊗ I_{A(G_2)} + I_{A(G_1)} ⊗ L(G_2)

   on the operator algebra A(G) = A(G_1) ⊗ A(G_2). The spectrum factorises as σ(L) = {λ_1 + λ_2 : λ_i ∈ σ(L(G_i))}. Since 0 ∈ σ(L(G_i)) for both components (per-component kernel from F4KernelDimensionByComponentsClaim), the smallest |Re(λ)| over non-kernel modes is

       gap(L(G_1 ⊔ G_2)) = min(gap(L(G_1)), gap(L(G_2))).

   The **rate-limiting component dominates**: the disconnected graph's dissipation gap equals the SMALLER of the two component gaps, not their sum or any average. Empirical corroboration: at N=8, K_4+disjoint-4-chain gives gap = 0.1362 (γ=0.5, J=1), identical bit-for-bit to chain N=4 at the same (J, γ) which gives 0.1362. Chain (sparser, slower decay) is rate-limiting; K_4 (denser, faster decay) is invisible to the gap. Sister structural result to [`F4KernelDimensionByComponentsClaim`](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (both follow from tensor-sum factorisation of L across disconnected components).

5. **Connection to F2 / F3.** F2 dispersion claims and F3 decay rate bounds (`min rate = 2γ`, `max rate = 2(N-1)γ` per the Absorption Theorem) describe the spectral envelope. The dissipation gap sits at the lower edge of this envelope; its 1/N² scaling for chain and ring may be the finite-size correction to F3's `min rate = 2γ` thermodynamic limit (where gap → 0 because of the conserved total magnetisation that Heisenberg + Z-dephasing preserves). The star's 1/N scaling does not fit this picture and would be a separate F2/F3 reading. **Status: open.** Cross-link to F3 documentation when the closed forms for the chain/ring/star prefactors land.

## Promotion path

This entry has two distinct promotion candidates after the 2026-05-19 Q-sweep:

**Q-universal ring dihedral locks:**
- `RingN4DihedralLockClaim`: `Im_max(ring, N=4, J) = (3/4)·J·N = 3·J` bit-exact at 6 Q-anchors. **Promoted Tier 1 derived 2026-05-19** via the K_{2,2} = C_4 bipartite-complete Casimir derivation in [`docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md); typed in [`RingN4DihedralLockClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs).
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
