# PROOF F86a: EP Mechanism (Q_EP, t_peak)

**Status:** Tier 1 derived. Q_EP = 2/g_eff (definitional, g_eff := σ_0), t_peak = 1/(4γ₀), bit-exact verified against full block-L numerics (the decay timescale and the value, universal across all tested c, N, n, bond position). The 2-level reduction itself is **heuristic** (see Tier label below). The genuine EP is the toy 2×2 reduction; the local-vs-global EP connection (2026-05-06) was **retracted 2026-06-21** (F86a-retraction) and is now an **OpenQuestion**: the full block-L is genuinely non-normal near Q_peak but has NO defective EP on the real Q axis (eigenvalues simple), and the Petermann magnitudes are grid-sensitive. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (OpenQuestion).
**Date:** 2026-05-02 (Statement 1); 2026-05-06 (local-vs-global EP connection, since retracted 2026-06-21 → OpenQuestion).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86a, the EP mechanism**: the 2-level rate-channel exceptional point behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14.
**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md): three-theorem overview and shared references.
**F-entry:** [F86a in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity); siblings [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md).

---

## Abstract

The Q_peak phenomenon has a mechanism, and this proof writes it down. Inside each coherence block of an XY chain under Z-dephasing, adjacent dephasing channels couple through the Hamiltonian by an effective coupling g_eff. The two-channel effective Liouvillian has a discriminant that vanishes at a specific Q value: an exceptional point where the two eigenvalues coalesce, both their real and imaginary parts colliding at once. The Q at which this happens is exactly Q_EP = 2/g_eff. The time at which the corresponding coherent response peaks (the EP-time) is t_peak = 1/(4γ₀), universal across chromaticity, chain length, and bond position.

Both numbers are bit-exact across every (c, N, n, bond) configuration we have tested. The proof itself is a clean two-level reduction inside each (n, n+1) coherence block, with the discriminant analysis dropping out the EP location and t_peak as algebraic consequences of the rate-channel pairing structure.

A subtlety the proof tracks honestly: the 2-level reduction is heuristic when projected onto a higher-dimensional block. The Q_EP value (definitional, g_eff := σ_0) and t_peak survive that projection bit-exactly against the full-block decay numerics. But the local-vs-global EP connection (whether the toy 2×2 EP is a faithful defective-EP of the full-block dynamics) was **retracted 2026-06-21** and is now OPEN: an artifact-free re-verification (Riesz spectral-projector norm) found the full block-L has NO eigenvalue coalescence on the real Q axis, its eigenvalues stay simple, so there is no real-axis defective EP to "hit". The block is genuinely NON-NORMAL there (large but FINITE Petermann factor), and the earlier peak magnitudes ("6× above FRAGILE_BRIDGE") were grid artifacts (K swings 2–4× over ΔQ = 1e-3). The genuine EPs are only the toy 2×2 reduction and the SEPARATE Σγ=0 gain-loss system (FRAGILE_BRIDGE); whether the full block shares their defective-EP structure off the real axis is open.

The diagnostic upshot is that Q_EP = 2/g_eff is the structural anchor of the whole F86 family. F86b inherits the universal shape from the EP rotation; F86c inherits the spatial-mirror invariance from the symmetry of the coupling structure; F86_block confirms that the dressed Q_peak values cannot be closed-form even though the underlying mechanism can. F86a is the mechanism; everything downstream is the consequence.

---

## Retraction (2026-06-21)

The "exceptional point **on the real Q axis**" reading of this proof is **retracted**. An artifact-free re-verification (Riesz spectral-projector norm) found the full (n, n+1)-coherence block has **no eigenvalue coalescence on the real Q axis**, its eigenvalues stay simple (nearest-neighbour gap ~0.25–0.35), so there is no real-axis defective EP. The block IS genuinely strongly non-normal there (a near-EP shadow; cond(V) up to ~268; ‖P‖ = √375 on a simple eigenvalue at N=5), with the Petermann factor large but **finite**, exactly `experiments/PT_SYMMETRY_ANALYSIS.md`'s reading ("no real-axis EP; large-but-finite Petermann signals a nearby EP in the complex plane"). The earlier real-axis Petermann-K evidence, "6× above FRAGILE_BRIDGE", K = 2384.7, the within-parity growth law, the odd/even parity asymmetry, is dropped as **grid artifacts** (K swings 2–4× over ΔQ = 1e-3).

**Survives unchanged:** `t_peak = 1/(4γ₀)`; `Q_EP = 2/g_eff` as a **definition** (g_eff := σ₀); the **toy 2×2** reduction's genuine EP (the Abstract's "exceptional point where the two eigenvalues coalesce" refers to this heuristic toy, not the full block); the **separate Σγ = 0 gain-loss EP** (`FRAGILE_BRIDGE`, K = 403); the centered-L_c **class-AIII-chiral** label.

**Open:** whether the full Σγ = N·γ₀ block has a defective EP off the real axis at all (the nearest complex-Q coalescences found 2026-06-21 are themselves *diabolic*). Full detail in the "local-vs-global EP connection: RETRACTED" note in §Proof below; typed as `LocalGlobalEpLink` (OpenQuestion).

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

This 2-level form is "PT-phenomenology-like" (EP at finite coupling, spectral flow); the centered Liouvillian L_c sits in **class AIII chiral** per [`experiments/PT_SYMMETRY_ANALYSIS.md`](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π linear, Π⁴=I, {Π, L_c}=0), distinct from Bender-Boettcher PT (classical PT requires anti-linear operators). The EP at Q_EP = 2/g_eff is the genuine EP of this toy 2×2 rate-channel reduction, a SEPARATE object, not "the instance of" the full-block classification (in isolation its eigenvalues coalesce at centre −4γ₀·k and carry no λ↔−λ pairing). The Hopf bifurcation in [`hypotheses/FRAGILE_BRIDGE.md`](../../hypotheses/FRAGILE_BRIDGE.md) is a DISTINCT genuine EP, the SEPARATE Σγ=0 gain-loss system (centre 0, exact λ↔−λ pairing), Petermann factor K = 403 in the complex γ plane. Whether the full Σγ=N·γ₀ block-L shares a defective-EP structure with these two is OPEN (see the local-vs-global note below).

**2026-05-06 local-vs-global EP connection: RETRACTED 2026-06-21 (now OpenQuestion).** The shared algebraic substrate survives: the same-sign-imaginary 2×2 form and the AIII chiral classification of the centered L_c, read at two residuals of the F1 palindrome `Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ vs Σγ = 0, where the gain side cancels the loss side). What was retracted is the reading built on the real-Q Petermann-K sweep. The original sweep at c=2 N=5..8 (data-pinning probe `compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`) recorded max K = 1333.6 / 337.9 / 2384.7 / 795.4 across N = 5 / 6 / 7 / 8 and was read as "a real-axis hit of the same EP FRAGILE_BRIDGE detects in complex γ, ≈ 6× above its K = 403, with an odd/even parity asymmetry confirming a σ_0 R-degeneracy prediction". An artifact-free re-verification (Riesz spectral-projector norm ‖P‖) refuted this: the full block-L has **NO eigenvalue coalescence on the real Q axis**, its eigenvalues stay SIMPLE there (nearest-neighbour gap ~0.25–0.35), so there is no real-axis defective EP to hit. The block IS genuinely strongly **NON-NORMAL** on the real axis (cond(V) = 48.7 / 50.9 / 268.5 at N = 5 / 6 / 7; ‖P‖ = 19.4 = √375 on a simple eigenvalue at N=5), with the Petermann factor large but FINITE. RETRACTED as grid artifacts: the peak magnitudes (K swings 2–4× over ΔQ = 1e-3), the "6× above FRAGILE_BRIDGE", the within-parity growth rates, and the odd/even parity asymmetry. The genuine EPs are ONLY the toy 2×2 reduction and the SEPARATE Σγ = 0 gain-loss system (FRAGILE_BRIDGE, K=403). OPEN QUESTION: whether the full Σγ = N·γ₀ block has an off-axis defective EP at all: a 2026-06-21 search found the nearest complex-Q coalescences are themselves DIABOLIC (‖P‖ = 1), not defective. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (OpenQuestion Claim; the four PetermannSpikeWitness rows are retained only as a cautionary non-normality record, not EP evidence).

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

**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md): three-theorem overview and the shared reference list.
**Sibling theorems:** [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md) (F86b), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md) (the g_eff obstruction proof), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md) (F86c).
**Chiral classification anchor:** [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md) (centered L_c is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear; "no real-axis EP, Petermann large but finite, nearby EP in the complex plane").
**Separate genuine EP (Σγ=0):** [FRAGILE_BRIDGE](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in complex γ plane), a DISTINCT gain-loss system, not a "global instance" of the full block-L.
**Empirical anchor:** [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/` carries `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `LocalGlobalEpLink` (OpenQuestion since the F86a-retraction 2026-06-21; four PetermannSpikeWitness entries retained as a cautionary non-normality record). CLI: `rcpsi inspect --root f86 --with-measured`.
