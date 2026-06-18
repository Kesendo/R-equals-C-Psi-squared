# The Structural Ceiling: closed forms for the topology gap rate

**Status:** Tier 1 derived (analytical proof + gate-exact numerical verification to machine precision, N=4..8). Resolves the `topology_band_edge` arc NextStep items 1, 2, 3.
**Date:** 2026-06-16
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Statement:** `g2(K_N) = 4/N` (N≥5), `g2(star_N) = 4/(N−1)` (N≥6), `g2(K_4) = 2 − 2/√3`, where `g2 = strict_gap / 2γ` is the high-Q structural ceiling.
**Typed claim:** [`StructuralCeilingClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/StructuralCeilingClaim.cs) (Tier 1 derived)
**Verifier:** [`topology_ceiling_rep_derivation.py`](../../simulations/topology_ceiling_rep_derivation.py) (gate-first, all stages exact)
**Builds on:** [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) (the floor `Re = −2γ⟨n_XY⟩`); the [topology band edge](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (the regime where the band edge loses gap-dominance).
**Formula registry:** [F122 in `ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).

---

## What this is about

A small quantum network, a handful of spins wired together, is being watched. The watching is not passive: light from outside reads the system out, and the reading makes its internal rhythms fade. Different patterns of coherence fade at different speeds, and the one that lingers longest is the system's memory, the last rhythm still ringing after the others have gone quiet.

On a plain chain the longest-ringing pattern is a single gentle wave, the most spread-out way one excitation can sit along the line. It rings at a fixed pace that the watching light cannot hurry: the chain keeps its clock, and the memory is that clock.

Wire the network more tightly and something new appears. A dimmer pattern, one that barely shows itself to the light, outlives even that ringing wave, and it does not ring at all; it just sits and fades, slowly. Because it lasts longer than the ringing pattern, it puts a lid on how long the audible memory can survive. That lid is the structural ceiling, and this document works out, for each shape of network, exactly how low the lid sits and why.

The answer is a story about symmetry. The more interchangeable a network's parts are, the more ways these dim, nearly hidden patterns have to fold themselves away from the watching light, and the lower the ceiling drops. The fully-connected network, where every part touches every other, drops it furthest; the star, one hub with many spokes, less so; the plain chain, with no symmetry to exploit, never grows such a pattern at all and keeps ringing freely.

Two surprises sit inside the story. At one particular small size, two networks that look nothing alike, the fully-connected one and the ring, go strange in the very same way and for the very same reason: a pattern built from two excitations, rather than one, takes the lead. And a tidy single rule that seems to predict every ceiling from how symmetric a network is almost works, then quietly fails for the ring, so we keep the honest shape-by-shape answers rather than the clean law that is not quite true.

## Abstract

Under uniform Z-dephasing, the slowest non-steady Liouvillian mode of an XY spin network sits at a decay rate `2γ·g2`, where `g2 = ⟨n_XY⟩` of that mode (the Absorption Theorem). For the chain `g2 = 1`: the single-excitation band edge `|vac⟩⟨ψ_k|` (`n_XY = 1`) is the slowest, and the clock reads it. On more connected graphs a darker mode appears below the band edge, capping the gap at `g2 < 1`, a structural ceiling. This document derives the closed forms.

The mechanism is degenerate perturbation theory at `γ ≪ J`. The decay rates are the eigenvalues of `N_XY` (the X/Y-counting operator, diagonal in the coherence basis with entry `hamming(a,b)`) restricted to each eigenspace of the Hamiltonian super-operator `ad_H`. Two kinds of slow mode result: the **band edge**, an `ad_H`-eigenblock inside the `(0,1)` excitation sector where every coherence has `hamming = 1`, so `N_XY = I` there and the rate is exactly `2γ` (`g2 = 1`) at all Q; and the **ceiling**, a coherence in the `ad_H` kernel (the H-commutant, `[H,A]=0`) with `⟨n_XY⟩ < 1`. The ceiling lives in the largest *degenerate single-particle level*: its coherences are the darkest because they can be made almost diagonal.

For the complete graph `K_N` that level is the `(N−1)`-fold `−J` band (the `S_N` standard representation), and the darkest coherence has `⟨n_XY⟩ = 2(1 − λ₂) = 4/N`, where `λ₂ = (N−2)/N` is the second principal-angle overlap between the commutant and the population (diagonal) operators. For the star `K_{1,N−1}` the level is the `(N−2)`-fold `0`-eigenvalue leaf manifold, giving `4/(N−1)`. Both are exact to machine precision (`K_5,6,7 = 4/5, 2/3, 4/7`; `star_6,7,8 = 4/5, 4/6, 4/7`). The connectivity ordering of the ceiling onset, chain never `<` star `N≥6` `<` complete `N≥5`, is the growth of that degeneracy with edge count.

`N = 4` is the shared outlier of `K_4` and the ring: the `4/N` ladder reaches `1` exactly at `N = 4`, so the single-excitation ladder no longer dips below the band edge, and the ceiling is set instead by the `(2,2)` half-filling two-excitation sector, the same sector that makes ring-4 special. There `K_4 = 2 − 2/√3 ≈ 0.8453` (below the floor) while ring-4 `= 1` (co-occupying the floor). The tempting universal law `4/(m+1)` in the degeneracy `m` fits complete and star but is **not** general: the ring's Fourier-degenerate manifold breaks it.

## 1. Setup and the high-Q mechanism

Take the open XY network `H = (J/2) Σ_{(i,j)∈E} (X_iX_j + Y_iY_j)` under uniform Z-dephasing at rate `γ`, with Liouvillian `L = L_H + L_D`, `L_H = −i[H,·]`, `L_D(·) = γ Σ_l (Z_l · Z_l − ·)`. Define the **gap rate**

    g2 = (smallest nonzero decay rate of L) / (2γ).

**Two diagonal facts.** In the computational coherence basis `|a⟩⟨b|`:

- `L_D(|a⟩⟨b|) = −2γ · hamming(a,b) · |a⟩⟨b|`, since `Z_l|a⟩⟨b|Z_l = (−1)^{a_l+b_l}|a⟩⟨b|`. So `N_XY := −L_D/2γ` is diagonal with eigenvalue `hamming(a,b) = ` the X/Y site count `n_XY` of that coherence.
- `ad_H = [H,·]` preserves the bigrading `(n_left, n_right) = (popcount a, popcount b)`, because `H` conserves the total excitation number. The whole problem is block-diagonal in `(p,q)`.

**The high-Q limit (`γ ≪ J`).** Treat `L_D` as a perturbation of `L_H`. To leading order the decay rates are the eigenvalues of the part of `N_XY` that commutes with `ad_H`, i.e. for each `ad_H`-eigenspace `Ω` the rates are `2γ ·` the eigenvalues of `P_Ω N_XY P_Ω`. The high-Q `g2` is the smallest nonzero such eigenvalue over all blocks. This reproduces the full Liouvillian `g2` at `Q = J/γ = 1000` to `O(1/Q) ≈ 10⁻⁶` for every topology and `N` tested (verifier Stage 0).

**The band-edge floor.** The `(0,1)` sector (and its conjugate `(1,0)`) consists of the coherences `|vac⟩⟨1_i|`, each with `hamming = 1`. So `N_XY = I` on this whole sector and `L = L_H − 2γ·I` there *exactly* (all Q): the eigenvalues are `iE_k − 2γ`, all at `Re = −2γ`. A band-edge mode at `g2 = 1` is therefore always present, so

    g2 ≤ 1   for every graph and N.

**The ceiling.** A structural ceiling (`g2 < 1`) occurs iff some `ad_H`-kernel mode (an operator with `[H,A]=0`) has `⟨n_XY⟩ < 1`. The kernel is spanned by coherences `|ψ_a⟩⟨ψ_b|` between degenerate Hamiltonian eigenstates. The darkest such coherence lives among the eigenstates of the largest degenerate level, because a larger level allows an operator that is more nearly diagonal in the site basis (lower X/Y content).

## 2. Theorem 1: the complete graph, `g2(K_N) = 4/N`

In the single-excitation sector the Hamiltonian acts as `H_SE = J·A`, `A` the adjacency matrix. For `K_N`, `A = 𝟙 − I` where `𝟙 = N|w⟩⟨w|` is the all-ones matrix and `|w⟩ = N^{−1/2} Σ_i |e_i⟩` is the uniform single-excitation state. The spectrum of `A` is `N−1` (once, eigenvector `|w⟩`) and `−1` (multiplicity `N−1`, eigenspace `w^⊥`).

The H-commutant inside the `(1,1)` sector is `C = ℂ|w⟩⟨w| ⊕ B(w^⊥)` (operators preserving the two `A`-eigenspaces). The only `n_XY = 0` (steady) element of `C` is the identity `I_SE` (a diagonal site operator commutes with `A` only if it is a multiple of `I`). We need the smallest nonzero eigenvalue of `N_XY` restricted to `C`.

On the `(1,1)` sector `N_XY` acts as `(N_XY M)_{ij} = 2(1−δ_{ij}) M_{ij}`, i.e. `N_XY = 2(Id − Diag)` where `Diag` keeps the site-diagonal of `M`. Hence on `C`

    min nonzero eig of N_XY|_C = 2(1 − λ₂),

with `λ₂` the second-largest eigenvalue of `P_C Diag P_C` (the largest, `λ₁ = 1`, is the steady `I_SE`). By Schur's lemma the `S_N` symmetry makes every traceless diagonal direction (the standard rep) have the same overlap, so we compute one. Take `d = |e_1⟩⟨e_1| − |e_2⟩⟨e_2|` (traceless, `‖d‖² = 2`). Its component in `C` is `P_⊥ d P_⊥` with `P_⊥ = I − |w⟩⟨w|` (the `|w⟩⟨w|` component vanishes since `⟨|w⟩⟨w|, d⟩ = 0`). Writing `u_i = P_⊥|e_i⟩`,

    ‖u_i‖² = 1 − 1/N,     ⟨u_1,u_2⟩ = −1/N,
    ‖P_⊥ d P_⊥‖² = ‖u_1‖⁴ + ‖u_2‖⁴ − 2|⟨u_1,u_2⟩|²
                 = 2((N−1)/N)² − 2/N² = 2·N(N−2)/N² = 2(N−2)/N.

Therefore `λ₂ = ‖P_C d‖²/‖d‖² = (N−2)/N`, and

    g2(K_N) = 2(1 − λ₂) = 2·(2/N) = 4/N.    ∎

This is `< 1` exactly when `N ≥ 5`. At `N = 4` it equals `1` (the band edge), so the `(1,1)` ladder produces no ceiling there; see §4. Verifier Stage 1 confirms `K_5,6,7 = 4/5, 2/3, 4/7` to `< 10⁻⁹`.

## 3. Theorem 2: the star, `g2(star_N) = 4/(N−1)`

For `star_N = K_{1,N−1}` (center site `0`, leaves `1..N−1`) the adjacency spectrum is `±√(N−1)` (once each) and `0` with multiplicity `N−2`. The largest degenerate level is the `0`-eigenvalue **leaf manifold** `S₀ = { v : v_0 = 0, Σ_{i≥1} v_i = 0 }`, the sum-zero leaf vectors, of dimension `N−2`. It sits inside the `(N−1)`-dimensional leaf space exactly as `w^⊥` sits inside `ℂ^N` in Theorem 1, now with `N → N−1`. The darkest coherence in `B(S₀)` (whose elements are supported entirely on leaf sites, so their site-diagonal is the leaf-diagonal) has, by the identical principal-angle computation,

    λ₂ = ((N−1) − 2)/(N−1) = (N−3)/(N−1),
    g2(star_N) = 2(1 − λ₂) = 4/(N−1).    ∎

This is `< 1` exactly when `N ≥ 6` (`star_5 = 4/4 = 1`, the band edge). Verifier Stage 3 confirms `star_6,7,8 = 4/5, 4/6, 4/7` to `< 10⁻⁹`. **This corrects** the band-edge arc's tentative reading that the star ceiling "saturates at 0.80, N-independent": it is `4/(N−1)`, and `0.80 = 4/5` was the `N = 6` value alone.

## 4. Theorem 3: the `N = 4` outlier is the half-filling sector

The `(1,1)` ladders of both `K_4` (`4/4`) and ring-4 (`4/3 > 1`) fail to produce a ceiling at `N = 4`: `4/N = 1` exactly meets the band-edge floor. The actual `N = 4` ceiling comes from the **`(2,2)` half-filling two-excitation sector**, the same sector that makes ring-4 special:

- **`K_4`**: the darkest `(2,2)` commutant coherence has `g2 = 2 − 2/√3 = 2(1 − 1/√3) ≈ 0.845299`, strictly below the floor. (The diagonal weight of that mode is `1/√3`.)
- **ring-4**: the darkest `(2,2)` coherence has `g2 = 1` exactly, co-occupying the floor (the known ring-4 "co-occupied floor": a `(2,2)` mode at `Re = −2γ` with `|Im| = 2√2·J >` band edge `2J`, so the clock reads `2√2·J`, not the band edge).

**Where `2√2·J` comes from, and why it is a mirror pair.** The ring-4 `(2,2)` two-excitation energies are `{−2√2, 0, 0, 0, 0, +2√2}·J`, two readings of the same fact:

- *The value* `2√2·J` is the top of the **anti-periodic** free-fermion band of the even-parity sector. Under Jordan-Wigner the parity string wraps the ring, so the even (two-excitation) sector carries anti-periodic boundary conditions and the single-fermion momenta shift to `k = π(2m+1)/N`. For `N = 4` this gives single-fermion energies `±√2·J` (`= 2J cos(π/4)`) and a two-fermion top `2√2·J`, which **exceeds the periodic single-excitation band edge** `2J` (`= 2J cos 0`, read in the odd sector). That overshoot is why the even-sector mode, not the band edge, holds the floor's clock hand.
- *The symmetry* is the palindrome about `0`: the **chiral K-mirror** `K·H·K = −H ⟹ E ↔ −E` ([`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs), `PROOF_K_PARTNERSHIP.md` — the chain's "second mirror", a sibling of F1 acting on `H`, not the Liouvillian), which holds because the even cycle `C_4` is **bipartite**. So `±2√2·J` are a chiral mirror pair, both at `Re = −2γ`. The odd rings `C_3, C_5` are *not* bipartite, their spectra are not `±`-symmetric, and this is the bipartite root of the even/odd ring dichotomy (`PROOF_K_PARTNERSHIP.md` §XXZ, `PROOF_RING_N4_DIHEDRAL_LOCK.md`).

So the `N = 4` anomaly is one mechanism in two topologies: the single-excitation ladder vacates the sub-floor region at `N = 4`, leaving the half-filling sector to set the gap. Verifier Stage 2 gates `K_4 (1,1) = 1.0`, the global minimum in `(2,2)`, and `= 2 − 2/√3` to `< 10⁻⁷`; Stage 2b gates the ring-4 `(2,2)` spectrum `= {±2√2, 0⁴}·J` against the anti-periodic two-fermion sums (the chiral palindrome).

## 5. The connectivity ordering, and what is not universal

The ceiling onset is the growth of the largest single-particle degeneracy `m` with edge count:

| topology | `m` (max adjacency degeneracy) | `(1,1)` ladder | ceiling onset |
|---|---|---|---|
| chain (path) | `1` (no degeneracy) | `> 1`, never dips | never |
| ring (cycle) | `2` | `4/3 > 1`, never dips | only `N=4` via `(2,2)` |
| star `K_{1,N−1}` | `N−2` | `4/(N−1)` | `N ≥ 6` |
| complete `K_N` | `N−1` | `4/N` | `N ≥ 5` (and `N=4` via `(2,2)`) |

More edges → larger degenerate manifold → darker coherence → lower ceiling → earlier onset. It is tempting to read off a universal `g2 = 4/(m+1)` (it gives `4/N` for `K_N` with `m=N−1`, and `4/(N−1)` for the star with `m=N−2`). **This is not a general law.** The ring's degenerate manifold is built from Fourier modes `e^{2πikj/N}`, whose embedding in the site basis differs from the complete/star case, and the principal-angle overlap is different: ring-5's `(1,1)` commutant is `1.6 ≠ 4/3`, ring-7's is `12/7 ≠ 4/3`. The per-family closed forms `4/N` and `4/(N−1)` are the real results; the degeneracy gives the *intuition* (more degeneracy → darker) but not the *value*. (Verifier final diagnostic prints the ring values to make the break explicit.)

## 6. Relation to the band edge and the Absorption Theorem

The structural ceiling is the Re-side companion of the topology band edge. The band edge `= J·ρ` (`ρ` the adjacency spectral radius) is the Im side: when it is the gap mode the clock reads `ω_mem = J·ρ` at `Re = −2γ` (`g2 = 1`). The ceiling is the regime where a darker `[H,A]=0` coherence undercuts it, dragging the strict gap to `g2 < 1` and the clock to `ω_mem = 0` (overdamped). Both rest on the Absorption Theorem `Re = −2γ⟨n_XY⟩`: the band edge is its `⟨n_XY⟩ = 1` line, the ceiling its darkest commutant value. Nothing here uses the chain gap-dominance proof, so the closed forms were Tier 1 derived independently of it (that proof has since landed, [PROOF_CHAIN_GAP_DOMINANCE](PROOF_CHAIN_GAP_DOMINANCE.md), lifting `TopologyBandEdgeClaim` to Tier 1 derived too).

## 7. Corollary: the star has no coherence horizon (the low-Q regime)

The chain has a **coherence horizon** `Q*(N)`: below it the band edge stops being the slowest oscillating mode and an overdamped mode holds the gap (the clock reads `ω_mem → 0`). That horizon is a *dispersion* effect, a square-root exceptional point where the `{0,2}`-coherence pair of the dispersive single-particle band coalesces; it reduces to the single-excitation (Haken-Strobl) Liouvillian and scales as `Q*(N) ~ 2N/π` ([`CoherenceHorizonClaim`](../../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs), `PROOF_COHERENCE_HORIZON_SLOPE.md`).

**The star has no such horizon.** Its single-particle band is *flat*: the star adjacency spectrum is `±√(N−1)` (once each) and `0` with multiplicity `N−2`. No dispersion, no coalescence, no EP. The chain's Haken-Strobl SE-EP does not port: applied to the star it predicts a spurious horizon (e.g. `N=4` at `Q ≈ 261`) that the full Liouvillian flatly contradicts (the star is already protected far below it, at `Q = 20`).

Instead the *same* `(1,1)`-commutant value `4/(N−1)` that sets the high-Q ceiling governs the star at **every** `Q`:

- `4/(N−1) > 1` (`N = 4`): the sub-band mode is above the `2γ` floor, so the band edge protects down to a low-Q *crossing* near `Q ≈ 1.9` (a real, overdamped mode crosses the floor, not an EP).
- `4/(N−1) = 1` (`N = 5`, marginal): the sub-band mode sits *exactly* on the floor and is pulled marginally below at finite `Q` as `g2 = 1 − 1/Q²`. There is no horizon and no oscillating memory at any finite `Q`, only asymptotic protection. (An apparent `Q* ≈ 316` read off a finite tolerance is an artifact of this `1/Q²` approach.)
- `4/(N−1) < 1` (`N ≥ 6`): the structural ceiling `g2 = 4/(N−1)` (§2–§3).
- `N = 3` is the path `P_3`, i.e. a chain, the lone exception, with the genuine `{0,2}`-EP at `√2`.

So the low-Q question is **subsumed** by the high-Q ceiling: one number, `4/(N−1)`, decides protection at all `Q`, and there is no separate star `Q*(N)` closed form. The chain's coherence horizon is the special case of a *dispersive* band; flat-band graphs do not have one. Verifier: [`star_no_coherence_horizon.py`](../../simulations/star_no_coherence_horizon.py) (gate-first: chain SE-EP fidelity, the spurious star SE-EP, the `4/(N−1)` governance, the real-mode dichotomy).

Numerical corollary (Tier 2): it is the *hub*, not the flatness, that is decisive. A wheel graph (the same hub with a fully dispersive leaf-ring of tunable strength `ε`) still has no coherence horizon at any `ε`, robust across `N = 5, 6, 7`; the hub always hosts a real, zero-frequency survivor that outlives the dispersive leaf modes. So the no-horizon extends from the flat-band star to *any* hub-graph (the natural generalisation, not proven analytically here). See [`THE_HUB_KILLS_THE_HORIZON.md`](../../experiments/THE_HUB_KILLS_THE_HORIZON.md), verifier [`wheel_qstar_bandwidth.py`](../../simulations/wheel_qstar_bandwidth.py).
