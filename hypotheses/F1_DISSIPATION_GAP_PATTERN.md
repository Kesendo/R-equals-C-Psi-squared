# F1 Dissipation Gap Pattern

**Tier 3 (reading) sharpened to Tier 3 (empirical scaling).** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e) with the 4 N=8 anchors. Extended 2026-05-19 first with chain/ring/star × N=3..6 Python anchors plus the N=9 chain run via the MklDirect bridge (commit abb2d52), then with a 72-point Q-sweep across the 6 canonical Q-anchors from [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md). The structural picture has two layers:

1. **Chain plateau is Q-quadratic, not a single constant.** The "gap × N² ≈ 2.20" originally reported was the Q=2 plateau at γ=0.5, J=1. The Q-invariant statement is `gap · N² / γ ≈ 1.1 · Q²` (chain plateau N≥4; the ~10% is the drift of f(Q)/Q² ACROSS Q, 1.183 at Q=0.5 down to 1.057 at Q=2.5, not a finite-N residual, which runs about 1.4% over N = 4..9).
2. **Ring N=4 and N=6 exhibit Q-independent dihedral locks.** Ring N=4 saturates `Im_max = (3/4)·J·N = 3·J` at every Q tested, to machine precision (= `Im/σ = 3Q/4`); ring N=6 saturates `Im_max = 0.717129... · J · N` the same way. Both are Q-universal topology-N specific saturation laws; ring N=4 is typed Tier 1 as `RingN4DihedralLockClaim`.

Bond count alone does not predict any of this; the per-topology dispersion structure (open-chain modes, cyclic Bogoliubov modes, hub-spoke SU(2)) is the structural fingerprint.

## Observed gaps at N=8 (initial 4-point anchor)

From the four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files' `DissipationGap` fields (commit 89f725e). All four configurations share Heisenberg XXX, H = (J/4)·Σ(XX+YY+ZZ), at J=1.0 with uniform Z-dephasing γ=0.5.

| Topology | Bonds B | Components | Dissipation gap |
|---|---|---|---|
| chain N=8       | 7 | 1 | 0.0344 |
| star N=8        | 7 | 1 | 0.0870 |
| ring N=8        | 8 | 1 | 0.1339 |
| K_4 + disjoint 4-chain N=8 | 9 | 2 | 0.1362 |

## Cross-topology cross-N extension (2026-05-19 Python anchors)

Python anchors via `simulations/f1_topology_heisenberg_small_n_anchor.py` extended chain, ring, star to N=3..6 dense numpy eigvals. The same Heisenberg J=1, γ=0.5 convention. JSON files at `simulations/results/f1_n8_n9_metrics/{chain,ring,star}_N{3..6}_python.json`.

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

`f(Q) / Q²` drifts from 1.183 at Q=0.5 down to 1.057 at Q=2.5, a roughly 10% drift across the band. The dominant scaling is `f(Q) ≈ c·Q²` with `c ≈ 1.10` and a sub-leading correction that pulls the ratio down at high Q (its Q → 0 value is closed, see Q1; the form of the drift is open).

**The γ=0.5, J=1 "2.20" recovered as Q=2 plateau:** `f(Q=2) ≈ 4.38` at γ₀=0.05 rescales to `gap·N² ≈ 4.38·γ = 4.38·0.5 = 2.19` at the J=1, γ=0.5 convention, recovering the originally observed "2.20" plateau (4-anchor agreement N=4..6 + N=8 + N=9 from the 2026-05-19 N=9 bridge run). The 2.20 reading is the value of `f(Q=2)·γ` at γ=0.5, not a universal chain constant.

**Physical interpretation.** A 1D dispersive chain Hamiltonian + per-site Z-dephasing produces a "diffusive" Liouvillian. The slowest decay mode has wavevector k_min ∝ 1/N (open-boundary modes); the decay rate scales as (J²/γ)·k_min² ∝ J²/(γ·N²) → gap·N²/γ ∝ Q². The Q² scaling is dispersion-rooted; the prefactor's Q → 0 value is closed (Q1) and its sub-leading correction is open.

**N=9 chain verification at Q=2 (2026-05-19, landed):** prediction `gap ≈ 2.18/81 = 0.0269` versus observed `gap = 0.02728` from the MklDirect bridge run (γ=0.5, J=1). Match within ~1%. The chain Q=2 plateau spans N ∈ {4, 5, 6, 8, 9}, holding to about 1.5% across that range (2.18 to 2.21 in the table above); it is an empirical plateau, not an exact identity.

## Ring-topology dihedral locks: Q-universal saturation laws (bonus discoveries 2026-05-19)

The Q-sweep surfaced two clean Q-universal patterns for ring topology that were invisible in the Q=2-only data:

**Ring N=4** (the 4-cycle = K_{2,2} = bipartite complete on 2+2):

    Im_max(ring, N=4, J)  =  (3/4)·J·N  =  3·J     (machine precision at every Q tested)
    ↔  Im/σ  =  3Q/4

verified at Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} to relative error < 1·10⁻¹⁴ (machine precision; the largest is 5.1·10⁻¹⁵). This is a clean rational saturation: ring N=4 carries 50% more imaginary spread than star N=4 ((3/4)·J·N vs (1/2)·J·N). It is not traceable to the bipartite-complete structure and not to bond count: K_4, which is not bipartite and has 6 bonds, reaches the same 3J, and so does the diamond with 5, ten of the 38 connected labelled graphs in all. The star's J·N/2 is the minimum of ΔE_max over connected graphs (exhaustive at N ≤ 6) and 3J is the maximum at N=4; the ratio is that of the two ends. Closed-form derivation via the Casimir spectrum `{−2J, −J³, 0⁷, J⁵}` (multiplicities summing to 16 = 2⁴) of `H = J·S_A·S_B` with sublattice totals S_A = S_0 + S_2, S_B = S_1 + S_3 (see [`docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md)). Typed Tier 1 derived as [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) 2026-05-19.

**Ring N=6** (dihedral D_6 = D_3·Z_2):

    Im_max(ring, N=6, J)  =  0.717129... · J · N    (machine precision at every Q tested)
    ↔  Im/σ  =  0.717129... · Q

verified at the same 6 Q-anchors. The numerical value `0.717129` is constant to 6+ decimal places across all Q's (a Q-universal lock). **Closed form identified 2026-07-28:** `c₆ = (5 + √13)/12`, a quadratic surd. It comes out of the exact characteristic polynomial of the S_z = 0 sector, whose ground-state factor is `λ² + 8λ − 36` in units of 4H, giving `E₀ = −(2 + √13)/2` and `ΔE_max = (5 + √13)/2` ([`simulations/ring_n4_lock_gate.py`](../simulations/ring_n4_lock_gate.py) §7 checks the 4H-scaled factor itself; [`simulations/ring_dihedral_lock_limit.py`](../simulations/ring_dihedral_lock_limit.py) STAGE 2 checks the same root as `12c² − 10c + 1` in c). No Bethe ansatz is needed at this N; what needs it is only the N → ∞ limit.

These two locks are topology-N specific closed forms for ΔE_max, saturated to machine precision. **Ring N=4 promoted Tier 1 derived 2026-05-19** as [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) via the K_{2,2} bipartite-complete Casimir derivation. **Ring N=6** carries the quadratic surd c₆ = (5+√13)/12 (the factorisation above), so no Bethe ansatz is needed at that N (see [`docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md) section "Why this is N=4-specific" for why the simple Casimir argument does not extend). The dihedral-lock pattern at even N is itself a question: does ring N=8 also Q-universally saturate at some N-specific constant? Existing N=8 ring data at Q=2 (γ=0.5, J=1) gives `Im/σ = 1.4128 = 0.7064·Q`, so the three-point sequence is `c_4 = 0.75, c_6 = 0.7171, c_8 = 0.7064`, decreasing and already passing through 1/√2 ≈ 0.7071 between N=6 and N=8 (1/√2 is NOT an upper bound that c_∞ approaches from above; c_8 already sits below it). The N → ∞ limit could be 1/√2, slightly below it, or some other Bethe-ansatz number; three points is not enough to settle this. **RESOLVED 2026-06-04: it is ln 2.** The bound of [`PROOF_RING_N4_DIHEDRAL_LOCK`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md) §5 rests on D being Hilbert-Schmidt self-adjoint, which the Hermiticity of the Z_l buys; the saturation needs uniform γ over the full single-site set plus a fully polarised H extreme in addition, and the ring meets both at every N. That is [F148](../docs/ANALYTICAL_FORMULAS.md#f148-the-imaginary-reach-is-the-hamiltonian-spread-on-every-graph-minted-2026-07-31), Tier 1 and general, so the lock is purely the Hamiltonian's spectral spread at every N, not only where it has been checked: `c_N = 1/4 − E₀(N)/(J·N)` with E₀ the AFM Heisenberg-ring ground state. The per-bond ground energy has the Hulthén limit `E₀/(J·N) → 1/4 − ln 2`, so `c_∞ = ln 2 = 0.69315` (NOT 1/√2, which the even branch merely crosses between N=6 and N=8). Confirmed N=4..16 with a clean 1/N² approach ([`simulations/ring_dihedral_lock_limit.py`](../simulations/ring_dihedral_lock_limit.py)); c₆ = 0.71713 reproduces the empirical 0.717129, validating the reduction. The per-N value is just the ring ground state (no simpler closed form beyond the low-N minimal polynomials); only the limit is closed. **The even N are one branch of two:** the frustrated odd rings sit BELOW ln 2 and rise toward it (c₅ = 0.623607, c₇ = 0.657883, c₉ = 0.671922, c₁₁ = 0.678994), so the sequence over all N is not monotone.

## Ring and star follow different patterns

Ring N=3..8 at Q=2 (γ=0.5, J=1): `gap × N²` ranges 6.07 → 8.57, **not flat in N**. The Q-sweep shows the ring's N-residual is large (the cyclic Bogoliubov modes scatter much more than the chain's open-boundary k_min mode). At fixed N, however, the ring is Q-quadratic just like chain: the per-N constants are simply larger.

Star N=3..8 at Q=2: `gap × N²` grows monotonically 2.43 → 5.57, also **not 1/N² scaling**. The hub-spoke geometry gives `gap·N²/γ` growing with N (Schur-Weyl dispersion has spread N, not 1/N), connected to the same SU(2)/Schur-Weyl mechanism that produces the star's closed form `Im_max = J·N/2` (see [`STAR_CONFOCAL_LIMIT.md`](../experiments/STAR_CONFOCAL_LIMIT.md)).

The three topologies have qualitatively different N-scaling laws in the dimensionless `gap·N²/γ` ratio; all three are Q-quadratic as Q → 0, where the Zeno end gives gap = Q²γλ₁/8 in this entry's book, and only roughly so at the Q of the sweep, the star least:
- **Chain (1D open path)**: N-flat plateau (1/N² diffusive scaling on top of Q² prefactor); chain plateau f(Q) ≈ 1.1·Q²
- **Ring (1D periodic)**: N-growing constants (cyclic dispersion has tighter mode-spacing), per-N dihedral locks at N=4 (`3Q/4`) and N=6 (`0.7171·Q`)
- **Star (hub-spoke)**: N-growing constants (SU(2)/Schur-Weyl dispersion has spread ∝ N), Q-universal Im saturation at `J·N/2`

## Why bond count alone fails (sharpened)

At N=8: chain (B=7) gap 0.0344, star (B=7) gap 0.0870, ratio 2.5×. Same bond count, different geometry. Even more starkly at small N where the difference grows: chain N=4 (B=3) gap 0.136, star N=4 (B=3) gap 0.210, ratio 1.54.

The 2026-05-19 cross-topology data confirms: bond count is irrelevant, **graph dispersion structure** is the right parameterisation. Chain has linear-band dispersion (slow modes at k ∝ 1/N → 1/N² gap); star has, at Q = 2, hub-localised modes (a different scaling family there); ring has cyclic dispersion with degenerate slowest modes.

## Structural questions: status

Five structural questions were posed when this entry first landed (2026-05-18, after the SLOW_N8 sweep). Their status:

- **Q1 chain prefactor closed form:** closed at Q → 0, c(0) = N²(1 − cos(π/N))/4 (detail 1); the drift with Q OPEN
- **Q2 star and ring scaling laws:** SHARPENED (three distinct scaling families identified at finite Q; closed forms open beyond their Q → 0 values, which the Zeno end closes for the chain and the star, and for the ring by the outside route it recognizes)
- **Q3 J ≠ 2γ regime:** **RESOLVED 2026-05-19** via the Q-sweep covering Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5}
- **Q4 K_4 + disjoint rate-limiting component:** **RESOLVED 2026-05-19** via tensor-sum factorisation
- **Q5 F2 / F3 connection:** **RESOLVED 2026-05-19** via the per-joint-popcount-sector diagnostic ([`experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md`](../experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md)): slow mode lives in the central diagonal popcount block `(⌈N/2⌉, ⌈N/2⌉)`, is 93-97% n_XY=0 (near-stationary I/Z-only) with a small n_XY=2 magnon admixture, and the Absorption Theorem `gap = 2γ·⟨n_XY⟩_slow` is confirmed to a relative deviation of 1.8e-15 / 6.8e-14 / 1.0e-13 at N = 4 / 5 / 6. The empirical form 0.55·Q²/N² of ⟨n_XY⟩ at Q = 2 matches the measurements to ~1%. F1² + F2 + F3 + F50 + Absorption Theorem compose at finite Q: X⊗N, the square of F1, pairs the diagonal sectors around the center, F50 puts a 2γ floor under the off-diagonal popcount sectors (a floor their slowest mode attains, not a rate the sector carries as a whole), F2 supplies the magnon's frequency scale from the `(0,1)` coherence block, F3/Absorption gives Re(λ) from ⟨n_XY⟩ in the diagonal sector where the gap actually lives.

Detail per question:

1. **Closed form for the chain prefactor c ≈ 1.10 in f(Q) = c·Q² + sub-Q² correction.** The Q² scaling is physically motivated by 1D diffusion (gap·N²/γ ∝ J²/γ² = Q²); the prefactor c ≈ 1.10 with ~10% sub-Q² drift across Q ∈ [0.5, 2.5] is open in its drift. The Q → 0 end is closed: at strong dephasing the populations relax under the exclusion graph's Laplacian, blind to the anisotropy, whose gap in every block of the chain is the path's ([PROOF_WEIGHT1_DEGENERACY](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end), the Zeno end), so D06's asymptote in this entry's book H = (J/4)·Σ σ_i·σ_j gives c(0) = N²(1 − cos(π/N))/4 = 1.1716, 1.1936, 1.2058 at N = 4, 5, 6 (mean 1.190), and π²/8 = 1.2337 as N → ∞; the full Liouvillian at Q = 0.5 gives the mean 1.1829 over N = 4..6, the 1.183 of the table (gate row Z8 of [`f50_zeno_end_ferromagnet.py`](../simulations/f50_zeno_end_ferromagnet.py) prints both). **Status: the sub-Q² term is open.**

2. **Star and ring scaling laws.** **Status: substantially sharpened 2026-05-19 by the Q-sweep.** Three distinct scaling families emerge cleanly:
   - **Chain** (open 1D path, k_min = π/N, the lowest Neumann cosine of the population diffusion): `gap·N²/γ ≈ 1.10·Q²` at the N≥4 plateau (the prefactor closed at Q → 0, its drift open, see Q1).
   - **Ring** (periodic 1D cycle, k_min = 2π/N): `gap·N²/γ ≈ 4·Q²` at the N≥5 plateau (the dihedral lock at N=4 is a finite-size outlier). The ~4× chain-to-ring prefactor ratio matches `(2π/N)² / (π/N)² = 4`, at Q → 0 the ratio of the two Laplacian gaps `(1 − cos(2π/N))/(1 − cos(π/N)) = 2(1 + cos(π/N)) → 4` (the ring's gap equality holding by the outside route the Zeno end recognizes, and read in gate row Z7 of [`f50_zeno_end_ferromagnet.py`](../simulations/f50_zeno_end_ferromagnet.py)), i.e. the **squared-wavevector ratio of cyclic-vs-open lowest modes**. This is suggestive that ring and chain share the same `Q²·γ/N²` diffusive form with a topology-specific k_min² factor.
   - **Star** (hub-spoke): `gap·N/γ` ≈ const at Q = 2 (the 4 anchors at N=3..6 give 1.62, 1.68, 1.64, 1.56), i.e. **gap ~ 1/N (NOT 1/N²)** there; toward Q → 0 the gap turns N-independent instead, the star being a tree whose Zeno gap is Q²γ/8 in this book, the same in every block (λ₁ = 1 of the star), and at Q = 0.5 it lies within 10% of that for N = 3 to 6 (gap/γ = 0.0314 to 0.0293 against 0.03125, `simulations/results/q_sweep_anchor/`), so the 1/N reading belongs to Q = 2 and the N-independence to Q → 0 at fixed N. The scaling family is different from chain/ring because hub-spoke geometry has no spatial dispersion; at Q = 2 the slowest mode is set by the hub-localised eigenmode whose decay rate scales as bandwidth/N rather than k_min². The Q-dependence of the star prefactor is non-Q² (at N = 5 a drift from `(gap·N/γ)/Q² ≈ 0.60` at Q=0.5 to ≈ 0.35 at Q=2.5); closed form open.

   What is closed: chain and ring share the diffusive `Q²·γ/N²` form with topology-specific k_min² prefactors; star is in a different scaling family at finite Q. What remains open: the three prefactors beyond their Q → 0 values (the chain's and the star's closed by the tree theorem, the ring's by the outside route the Zeno end recognizes, read in gate row Z7), and the sub-leading corrections.

3. **J ≠ 2γ regime.** **Status: RESOLVED 2026-05-19.** The Q-sweep (`f1_q_sweep_anchor.py`) measured chain/ring/star × N=3..6 × Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} at γ₀ = 0.05: 72 anchors covering six J/γ ratios from 0.5 to 2.5. The data confirms `gap·N²/γ ≈ c·Q²` (with c topology-specific) holds across the full range, ruling out the "chain prefactor decomposes as J·a + γ·b" linear conjecture in favour of the cleaner Q² form. The earlier "2.20" reading was the chain Q=2 plateau at γ=0.5.

4. **K_4 + disjoint behaviour (rate-limiting component).** **Status: RESOLVED 2026-05-19** via tensor-sum factorisation (the same mechanism that closes F4 kernel-dim factorisation across components). For G = G_1 ⊔ G_2 with the per-component Heisenberg + Z-dephasing Liouvillian L(G_c), the bond-locality of both H and the dephasing dissipator gives

       L(G_1 ⊔ G_2) = L(G_1) ⊗ I_{A(G_2)} + I_{A(G_1)} ⊗ L(G_2)

   on the operator algebra A(G) = A(G_1) ⊗ A(G_2). The spectrum factorises as σ(L) = {λ_1 + λ_2 : λ_i ∈ σ(L(G_i))}. Since 0 ∈ σ(L(G_i)) for both components (per-component kernel from F4KernelDimensionByComponentsClaim), the smallest |Re(λ)| over non-kernel modes is

       gap(L(G_1 ⊔ G_2)) = min(gap(L(G_1)), gap(L(G_2))).

   The **rate-limiting component dominates**: the disconnected graph's dissipation gap equals the SMALLER of the two component gaps, not their sum or any average. Empirical corroboration: at N=8, K_4+disjoint-4-chain gives gap = 0.1362 (γ=0.5, J=1), agreeing with chain N=4 at the same (J, γ) to the four decimals both were quoted at. The EQUALITY is exact, but as a consequence of the factorisation above; the two numbers are separate eigensolver outputs and were never compared below that precision. Chain (sparser, slower decay) is rate-limiting; K_4 (denser, faster decay) is invisible to the gap. Sister structural result to [`F4KernelDimensionByComponentsClaim`](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (both follow from tensor-sum factorisation of L across disconnected components).

5. **Connection to F2 / F3.** **Status: RESOLVED 2026-05-19** via the per-joint-popcount-sector diagnostic (see [`experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md`](../experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md) for the sector data + Pauli-weight-distribution evidence). The slow mode lives in the **central diagonal popcount sector** `(⌈N/2⌉, ⌈N/2⌉)` (the largest block, of dimension C(N, ⌈N/2⌉)²), is **93-97% n_XY=0** (near-stationary I/Z-only Pauli content), with a small **n_XY=2 magnon admixture** of weight `w_2 ≈ 0.275·Q²/N²`. The Absorption Theorem `Re(λ_slow) = -2γ·⟨n_XY⟩_slow` is exact as a theorem, and the two independent float routes to it (the eigensolver's Re λ and the Pauli-basis projection's ⟨n_XY⟩) agree to a relative deviation of 1.8e-15 / 6.8e-14 / 1.0e-13 at N = 4 / 5 / 6, γ=0.5, J=1, which in absolute terms is 2.5e-16 / 6.0e-15 / 6.3e-15. That is a float floor between the two routes, and its N-dependence is NOT established: the relative column mostly tracks the 1/N² gap in its denominator, and the absolute one is flat from N=5 to N=6. With `⟨n_XY⟩_slow = 2·w_2 ≈ Q²/N² × c` where c ≈ 0.55 is the empirical magnon-mixing coefficient (closed at Q → 0 by the Zeno end, see Q1; the drift with Q open). F1² + F2 + F3 + F50 + Absorption Theorem compose cleanly: X⊗N, the square of F1, pairs the diagonal sectors around the center, F50 puts a 2γ floor under the off-diagonal popcount sectors (their slowest mode sits on it; the sectors themselves spread above it, and only the two END sectors `(0,1)` and `(N−1,N)` sit at −2γ as a whole), F2 supplies the magnon's frequency scale from its own object the `(0,1)` coherence block, and F3/Absorption Theorem reads Re(λ) from `⟨n_XY⟩` in the diagonal popcount sector where the gap actually lives.

   **Three-line summary of the precise statements (full data in the CHAIN_GAP_SECTOR_DIAGNOSTIC experiment doc):**
   - **F2** is the oscillation-frequency formula `ω_k = 4J·(1 − cos(π·k/N))` of the `(0,1)` coherence block: `Im(λ)` only, and only for that block, not for the off-diagonal popcount sectors as a class. At finite Q it does not predict the dissipation gap; at Q → 0 it does, the Zeno gap being its slowest frequency times J/γ in the Pauli book, one Laplacian running both (F152).
   - **F3 / Absorption Theorem** is `Re(λ) = -2γ·⟨n_XY⟩` for any Lindblad eigenmode. Its "min rate = 2γ" applies only to pure-w=1 modes; the diagonal-popcount slow mode has fractional ⟨n_XY⟩ ≈ 0.55·Q²/N², far below 1.
   - **The dissipation gap is the magnon-admixture weight × 4γ.** The slow mode is a 93-97% I/Z near-stationary operator dressed with a small n_XY=2 single-magnon (XX or YY pair) admixture; the empirical `w_2 ≈ 0.275·Q²/N²` (perturbation-theory order of magnitude `Q²/N²` with empirical 0.275 prefactor) and the Absorption Theorem turns `4γ·w_2 = 1.10·γ·Q²/N²` into the decay rate.

   **External literature anchors (2026-05-19 web search):**

   - **Medvedyeva, Essler, Prosen (PRL 2016, [arXiv:1606.09122](https://arxiv.org/abs/1606.09122))** derives the exact Bethe-ansatz spectrum of the periodic XX chain with on-site Z-dephasing: `gap = 2π²·J²/(γ·L²)`, i.e. `gap·N²/γ ≈ 19.74·Q²` for periodic XX. Divided by the factor 4 between the two books that CHAIN_GAP_SECTOR_DIAGNOSTIC applies, 19.74 is π²/2 = 4.93, the ring's Zeno-end coefficient N²(1 − cos(2π/N))/4 → π²/2 in this entry's book, which the Zeno generator gives for XX and XXX alike. Our chain 1.10·Q² sits below it by the ring-to-chain factor 4, which does not see the ZZ term, by finite N, and by the finite-Q drift, which at Q = 2 the ZZ term turns downward: in the gap's units the XX chain lies there at 1.265, 1.248, 1.241 for N = 4, 5, 6, above its Zeno values 1.172, 1.194, 1.206, and the Heisenberg chain at 1.089, 1.105, 1.092, below them (gate row Z8). The 4× ring/chain ratio in our data matches MEP's `(2π/N)² / (π/N)²` wavevector-squared ratio structurally.
   - **Bortz, Stolze (PRB 2008, [arXiv:cond-mat/0612382](https://arxiv.org/abs/cond-mat/0612382))** for the inhomogeneous central-spin model: "oscillation frequency is proportional to the number of spins, whereas the amplitude behaves like 1/N". Independent literature confirmation that our star scaling family is a documented physical regime (`StarImMaxBoundClaim` Im_max = J·N/2 frequency ∝ N + gap·N/γ const at Q = 2, amplitude ∝ 1/N).
   - **Žnidarič ([arXiv:2311.07375](https://arxiv.org/abs/2311.07375))** confirms that local Z-dephasing on XX gives diffusive transport with N⁻² gap scaling; our chain/ring sit in this regime.

   **Sub-leading work needed (not blocking the resolution):**
   - **Closed-form derivation of the 0.55 magnon-mixing coefficient's drift.** Its Q → 0 value is second-order perturbation theory on the kernel of `L_D`, N²(1 − cos(π/N))/8, half of Q1's gap coefficient, with no ZZ factor at that order; the drift with Q, which the ZZ term turns downward at Q = 2, is what remains.

## Promotion path

**Already promoted (typed Tier-1-derived elsewhere):**

- [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs): the Im_max bound on the 4-cycle, Tier 1 derived 2026-05-19 via K_{2,2} = C_4 Casimir derivation in [`PROOF_RING_N4_DIHEDRAL_LOCK.md`](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md). This is an Im_max bound, not a gap closed form, but it surfaced from the same Q-sweep.
- [`StarImMaxBoundClaim`](../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs): the universal Im_max = J·N/2 saturation on any star, Tier 1 derived 2026-05-19 via SU(2)/Schur-Weyl hub-leaf Casimir in [`PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md`](../docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md). Same caveat: Im_max bound, not gap.

**Structurally resolved in this doc, not separately typed:**

- **Disconnected-graph gap factorisation:** `gap(L(G_1 ⊔ G_2)) = min(gap(L(G_1)), gap(L(G_2)))` via tensor-sum factorisation of the bond-local Liouvillian (Q4 above). This is the gap-side companion to [`F4KernelDimensionByComponentsClaim`](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (kernel-dim product factorisation across components). Could be typed as a `DisconnectedGapMinClaim` Tier 1 derived if a downstream consumer needs `.Predict(componentGaps)` as a typed API; for now the structural argument lives in Q4's resolution text and the F4 claim documents the underlying tensor-sum mechanism.

**Still candidate, not yet promoted:**

- `RingN6DihedralLockClaim`: `Im_max(ring, N=6, J) = c₆·J·N` with `c₆ = 0.717129...` at 6 Q-anchors, machine precision. **Reframed 2026-06-04:** c₆ is the N=6 AFM Heisenberg-ring ground state via `c_N = 1/4 − E₀(N)/(J·N)`, which by F148 is the Liouvillian's imaginary reach at every N. **Closed form 2026-07-28:** that ground state is exactly `−(2+√13)/2`, so `c₆ = (5+√13)/12`; the neighbours stay algebraic too (c₈ the largest root of `512c³ − 640c² + 232c − 25`, c₁₀ of a sextic). What IS closed is the limit: `c_∞ = ln 2` (Hulthén), the sequence c_4=0.75, c_6=0.7171, c_8=0.7064 converging to ln2 ≈ 0.6931 with a 1/N² approach along the EVEN N, merely crossing 1/√2 between N=6 and N=8; the frustrated odd rings rise to the same limit from below (c_5=0.623607, c_7=0.657883, c_9=0.671922). See §"Ring N=6" above + `PROOF_RING_N4_DIHEDRAL_LOCK` §"The N → ∞ limit".

- `ChainGapClosedForm`: gap_chain(N, J, γ) = c(Q) · γ/N² with c(Q) ≈ 1.10·Q² at the N≥4 plateau (Q1 above). Promote to Tier 1 candidate when:
  - The sub-Q² correction to the chain prefactor admits closed-form derivation (its Q → 0 value is closed, Q1).
  - The N-dependence at fixed finite Q is reconciled with its Q → 0 form N²(1 − cos(π/N))/4 (Q1).

- `RingGapClosedForm`: gap_ring(N, J, γ) = c(Q) · γ/N² with c(Q) ≈ 4·Q² at the N≥5 plateau (the dihedral lock at N=4 is a finite-size outlier; Q2 above). Same diffusive-form premise as chain with a 4× prefactor that matches the (2π/N)² / (π/N)² wavevector-squared ratio. Promotion needs the same closed-form derivation as chain.

The chain and ring gap closed forms would naturally share a parent abstraction `DiffusiveGapClosedForm` since both fit `gap = (k_min² · Q² + sub-leading) · γ/N²`. Star is in a separate scaling family at finite Q (gap ~ 1/N) and would not fit that abstraction; at Q → 0 the tree theorem gives chain and star one form, the Laplacian gap 2J²λ₁/γ in the Pauli book.

## Cross-references

- Anchor data Q=2 (γ=0.5, J=1): `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json` + the `_N{3..6}_python.json` companions; `chain_N9.json` for the N=9 bridge run.
- Anchor data Q-sweep (γ₀=0.05 substrate convention): `simulations/results/q_sweep_anchor/{chain,ring,star}_N{3..6}_Q{0.5..2.5}.json` (72 files, 24 per topology) and `f1_q_sweep_anchor.py` plus `q_sweep_anchor_console.log`.
- Canonical Q-anchor map: [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md).
- Companion typed claim from the same sweep (the closed-form discovery that did promote): [F4KernelDimensionByComponentsClaim](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 derived as of 2026-05-19; landed Tier 1 candidate 2026-05-18, promoted after DEGENERACY_PALINDROME Result 2 was identified as the connected-case upper-bound closure; kernel-dim factorisation across components).
- Sister Tier 3 reading from the same sweep: [Star Spectrum Compactness](STAR_SPECTRUM_COMPACTNESS.md) (whose Reading 1 was resolved by [Star Confocal Limit](../experiments/STAR_CONFOCAL_LIMIT.md), itself extended with the 24-anchor Q-sweep verifying `Im_max(star) = J·N/2 ∀ Q`).
- F1 verification record that produced the data: [F1GeneralTopologyVerifiedClaim](../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs).
- Related: [the F1 general-topology proof](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) closed form for the F1 residual norm, whose graph-additivity is exact by an algebraic argument, the Frobenius orthogonality of the bond-resolved residuals, and is verified numerically through N = 6, while the same 4 N=8 topologies check the F1 palindromic pairing; the dissipation gap is the next structural quantity beyond the residual norm and does NOT follow the same simple (B, D2) parameterisation).
