# PROOF F86a: EP Mechanism (Q_EP, t_peak)

**Status:** Tier 1 derived. Q_EP = 2/g_eff, t_peak = 1/(4γ₀), bit-exact verified against full block-L numerics (universal across all tested c, N, n, bond position). The 2-level reduction itself is **heuristic** (see Tier label below). Local-vs-global EP connection (2026-05-06): Tier 2 verified at c=2 via Petermann-K sweep N=5..8 (max K = 2384.7 at N=7, ≈ 6× above FRAGILE_BRIDGE's K=403); analytic continuation along complex γ remains open. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs`.
**Date:** 2026-05-02 (Statement 1); 2026-05-06 (local-vs-global EP connection, Tier 2 verified at c=2 via Petermann-K sweep N=5..8).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86a, the EP mechanism**: the 2-level rate-channel exceptional point behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14.
**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and shared references.
**F-entry:** [F86a in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity); siblings [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md).

---

## Statement 1 (EP mechanism). [Tier 1 for EP location and t_peak; 2-level reduction heuristic]

For a uniform N-qubit XY chain with Z-dephasing γ₀ at each site, the (n, n+1) coherence block contains c = min(n, N−1−n) + 1 pure dephasing rates 2γ₀·HD with HD ∈ {1, 3, ..., 2c−1} (F74). For each pair of adjacent channels at HD = 2k−1 and HD = 2k+1 (k = 1, ..., c−1), the two-level effective Liouvillian with inter-channel coupling J·g_eff has eigenvalues

    λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The discriminant vanishes at the **exceptional point** (EP)

    J·g_eff = 2γ₀     ⟺     Q_EP = 2 / g_eff

At the EP, λ_±(k) = −4γ₀·k. The slowest mode (k = 1) gives the universal e-folding time

    t_peak = 1 / (4γ₀)

independent of c, N, n, and bond position.

**Tier label.** The EP location Q_EP = 2/g_eff and t_peak = 1/(4γ₀) are bit-exact verified against full block-L numerics (universal across all tested c, N, n, bond position), Tier 1. The 2-level reduction itself is **heuristic**: in the channel-uniform basis M_H_eff is diagonal (off-diagonals exactly zero, established empirically across N=3..9 in `_eq022_b1_channel_projection.py`), so the EP physics lives in the orthogonal complement. The effective 2-level form reproduces the EP location and degenerate-eigenvalue structure, but the explicit basis change from full block-L to the reduced model is not derived here.

---

## Proof of Statement 1 (EP mechanism)

This is elementary 2×2 matrix algebra applied to a Liouvillian sub-block.

### Setup

Let n be fixed, c = c(n, N), and consider the (n, n+1) coherence block of L = L_H + L_D with H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) and L_D the uniform Z-dephasing dissipator at rate γ₀.

By F74 and F2b, at J = 0 the block-restricted L is diagonal with eigenvalues 2γ₀·HD for HD ∈ {1, 3, ..., 2c−1}, each with multiplicity ≥ 1. As J grows, H couples the rate channels to one another. Adjacent rate channels, those differing by Δ(HD) = 2, couple at first order in J (the bond flips two adjacent bits, changing HD by 0 or ±2). Non-adjacent channels (Δ(HD) ≥ 4) couple only at higher order in J.

### Two-level effective

Restrict to the subspace spanned by the two "slowest" relevant rate channels: HD = 2k−1 and HD = 2k+1 for some k ∈ {1, ..., c−1}. Within this subspace the effective Liouvillian, in a suitably chosen rate-channel basis, takes the form

    L_eff(k) = [ −2γ₀(2k−1)     iJ·g_eff   ]
               [  iJ·g_eff      −2γ₀(2k+1) ]

The diagonal entries are real-negative (rate channels at J = 0). Both off-diagonal entries carry the same imaginary sign +iJ·g_eff (not the anti-Hermitian opposite-sign pattern). This same-sign-imaginary structure is what produces an EP at finite coupling. Verified numerically: the opposite-sign pattern (+iJg, −iJg) gives discriminant 4γ₀² + J²·g_eff² with no EP; only the same-sign pattern produces an EP at J·g_eff = 2γ₀.

This 2-level form is "PT-phenomenology-like" (EP at finite coupling, spectral flow) but algebraically inside **class AIII chiral** per [`experiments/PT_SYMMETRY_ANALYSIS.md`](../../experiments/PT_SYMMETRY_ANALYSIS.md), distinct from Bender-Boettcher PT (Π is linear; classical PT requires anti-linear operators). The local EP at Q_EP = 2/g_eff is the 2-level rate-channel instance of the chiral classification established for the full Liouvillian; the Hopf bifurcation in [`hypotheses/FRAGILE_BRIDGE.md`](../../hypotheses/FRAGILE_BRIDGE.md) is the global instance, with Petermann factor K = 403 signaling an EP in the complex γ plane.

**2026-05-06 local-vs-global EP connection (Tier 2 verified at c=2).** The local-2-level-EP and global-complex-γ-EP are the same algebraic object: same-sign-imaginary 2×2 form, AIII chiral classification, read at two residuals of the F1 palindrome `Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ for the local instance; Σγ = 0 for the global instance, where the gain side cancels the loss side). A Petermann-K sweep on the real Q axis at c=2 N=5..8 (data-pinning probe `compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`) records max K = 1333.6 (N=5, odd), 337.9 (N=6, even), 2384.7 (N=7, odd), 795.4 (N=8, even); by N=7 the spike sits ≈ 6× above FRAGILE_BRIDGE's K = 403 ballpark, with within-parity monotonic growth (odd 1.79× per step; even 2.36× per step) and a 2-4× odd/even asymmetry empirically confirming A3's σ_0 R-even/R-odd-degeneracy prediction (chain-mirror R splitting of σ_0 at even N; see `compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs`). Reading: the local F86 EP is a real-axis hit of the same EP whose near-singularity FRAGILE_BRIDGE detects in the complex-γ direction. The structural connection is an analytic continuation along complex γ, specified but not yet executed in code; explicit modulated gain-loss in `LindbladPropagator` or a closed-form K(N) at the EP would promote the Tier 2 to Tier 1. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (Tier2Verified Claim, four PetermannSpikeWitness entries pinning the table; PendingDerivationNote names the analytic-continuation gap).

(The heuristic-vs-Tier-1 split for this statement is in the Tier label of the Statement above.)

### Eigenvalues

L_eff(k) has

    Trace = −2γ₀(2k−1) + (−2γ₀(2k+1)) = −8γ₀·k
    det   = (−2γ₀(2k−1))·(−2γ₀(2k+1)) − (iJg_eff)·(iJg_eff)
          = 4γ₀²·(4k² − 1) − (i²J²g_eff²)
          = 4γ₀²·(4k² − 1) + J²·g_eff²

Eigenvalues λ_± = (Trace/2) ± √((Trace/2)² − det):

    (Trace/2)² − det = 16γ₀²·k² − 4γ₀²·(4k² − 1) − J²·g_eff²
                     = 4γ₀² − J²·g_eff²

Therefore

    λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The discriminant vanishes when J·g_eff = 2γ₀, giving the EP at Q_EP = 2/g_eff and degenerate λ_± = −4γ₀·k. Beyond the EP the eigenvalues form a complex-conjugate pair around the real centre −4γ₀·k.

Numerical verification of this 2-level form: at α = J·g_eff = 2γ₀ the eigenvalues of `[[−2γ₀, +i·2γ₀], [+i·2γ₀, −6γ₀]]` coalesce at λ = −4γ₀ (verified by direct diagonalisation). For α > 2γ₀ they become λ = −4γ₀ ± i·√(α² − 4γ₀²), the post-EP oscillation observed in the (n, n+1) block dynamics.

### t_peak universality

The slowest mode is k = 1 with λ_±(1) = −4γ₀ at the EP. Its e-folding time is 1/(4γ₀). Higher k decay faster (1/(8γ₀), 1/(12γ₀), ...) and are masked by k=1 in the long-time observable. Therefore

    t_peak = 1/(4γ₀)

for any (c, N, n, bond position), the slowest k=1 EP universally sets the J-derivative peak time. ∎

---

## Pointers

**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and the shared reference list.
**Sibling theorems:** [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md) (F86b), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md) (the g_eff obstruction proof), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md) (F86c).
**Chiral classification anchor:** [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear).
**Global EP instance:** [FRAGILE_BRIDGE](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in complex γ plane).
**Empirical anchor:** [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/` carries `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `LocalGlobalEpLink` (Tier2Verified, four PetermannSpikeWitness entries). CLI: `rcpsi inspect --root f86 --with-measured`.
