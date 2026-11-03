# PROOF F86b: Obstruction, why g_eff / Q_peak admits no closed form

**Status:** Structural negative result. The per-bond position Q_peak(c, N, b) and the coupling g_eff(c, N, b) it rides on admit no closed form by the six routes the closed-form effort opened (Obstruction Proof, below; L1 rigorously proven blocked, L2 a convention boundary at even N, L4 every reduction the effort built shown insufficient, L6 a demonstrated failure mode, L3/L5 proven decouplings). The diagnosis section names the single cause behind all six: g_eff is not a spectral primitive. What IS derived is in the sibling proofs: the symmetry layer (F86a EP mechanism, F86c F71 mirror) and the bare floor of the HWHM ratio (F86b₁); the lift above that floor (F86b₂) is observed, not derived.
**Date:** 2026-05-14 (Obstruction Proof and diagnosis); 2026-05-21 (F90-bridge corollary corrected: F89 D_k closed 2026-05-15, the predicted transfer of that closure to g_eff refuted); 2026-07-02 (Item 3 σ_0-scaling stamped: the 2√(2(c−1)) form retracted 2026-05-08, c=2 closed by F86e, only c ≥ 3 open).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries the **negative result and open frontier of F86b**: the closed-form gap on g_eff/Q_peak. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14. The universal-shape result (Statement 2, Tier 1 candidate, with the floor F86b₁ and the lift F86b₂) is the sibling [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md).
**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md), three-theorem overview and shared references.
**F-entry:** [F86b in the formula registry](../ANALYTICAL_FORMULAS.md).
**Related:** siblings [the EP-mechanism proof](PROOF_F86A_EP_MECHANISM.md), [the universal-shape proof](PROOF_F86B_UNIVERSAL_SHAPE.md), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md); [the F90 bridge proof](PROOF_F90_F86C2_BRIDGE.md), [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md).

---

## Abstract

F86b showed that the resonance shape around Q_peak is nearly universal (HWHM_left/Q_peak in two narrow bands, one for Interior and one for Endpoint bonds, across all read chromaticities). What the universal shape rides on is the EP coupling g_eff and the position Q_peak itself, both of which depend on the specific (chromaticity, chain length, bond position) configuration. This proof closes the negative question that hung over F86 for months: can g_eff and Q_peak be written down in closed form?

The answer is no, by each of six routes that we tried. L1 (spectral): the effective characteristic polynomial does not factor over ℚ. L2 (eigenvectors): at even N the EP-partner vectors are fixed only by a convention, and no closed form is known under it. L3 (observable-direct): the probe is orthogonal to the EP partners. L4 (reduced models): every finite reduction built misses the position. L5 (single element): the bond-class signature sits in a subspace orthogonal to the probe and does not scale with g_eff. L6 (empirical extrapolation): the accessible data gives trajectory crossings, not laws. A corollary adds the F89 bridge: D_k closed on 2026-05-15, and the predicted transfer to g_eff was tested and refuted, because the orbit-polynomial reduction that closes D_k discards the per-bond cross elements g_eff is built from. The diagnosis section names the single cause behind all six: g_eff is not a spectral primitive but the coupling of a two-level reduction.

The structural reading is that g_eff sits in the same family as F86's coupling residues from the start: a non-primitive, position-dependent number downstream of a fit or projection. F100 and F101 later showed the analogous obstruction for the κ leading-coefficient on the bond-mirror deviation (same family, same kind of residue, same closed-form gap). What survives a closed-form effort in this family is the symmetry (F86a, F86c) and the shape (F86b universal), not the specific numbers per (c, N, b).

The diagnostic upshot is that empirical anchoring becomes the right tool for these residues. The fitted (α, β) values in F86b₂'s c=2 per-bond description, the empirical 42:8 split in F110's diagonal hard cells, the κ relative spreads in F100's empirical witness, are all structureless residues with the same fingerprint. They are real numbers, they carry diagnostic value (deviation between predicted and measured reveals hardware deviation), but they will not collapse to closed form regardless of which polynomial-or-symbol pipeline we throw at them.

---

## Retracted Q_peak conjectures (Endpoint, Interior)

Both retractions followed extended-N data showing the apparent matches at N=7 were trajectory crossings, not asymptotes (see methodological lesson in [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md)).

**Endpoint:** conjectured `Q_peak(Endpoint, N) = csc(π/(N+1))` (anchored in F2b smallest-momentum OBC sine mode); refuted by the c = 3 envelope Q_peak (default Q grid, parabolic vertex):

| N | csc(π/(N+1)) | envelope Q_peak | deviation |
|---|--------------|-----------------|-----------|
| 5 | 2.000 | 2.393 | +19.6 % |
| 6 | 2.305 | 2.548 | +10.6 % |
| 7 | 2.613 | 2.584 | −1.1 % |
| 8 | 2.924 | 2.583 | −11.6 % |

The conjecture meets the data at N = 7 and misses it on both sides; the c = 3 Endpoint Q_peak reads 2.584 and 2.583 at N = 7 and 8, while csc(π/(N+1)) keeps rising. No closed form identified.

**Interior c=3 pentagonal asymptote:** conjectured `Q_peak(Interior, c=3) → csc(π/5) = 1.7013` as N → ∞; refuted by the envelope trend 1.600 → 1.679 → 1.712 → 1.723 for N=5..8 (crosses csc(π/5) between N=6 and N=7 and continues above). At c=4: 1.760 (N=7) → 1.791 (N=8). The N=4 golden-ratio structure ({±φ, ±1/φ} OBC dispersion, Im(λ) decompose into integer combinations of {φ, 1/φ, 1, √5}; [`experiments/MAJORANA_AXIS_MODES.md`](../../experiments/MAJORANA_AXIS_MODES.md)) is real but does not propagate to a closed form for Q_peak at general (c, N).

---

## Empirical Q_peak data (class-averaged envelope, default Q grid of step 0.025 with parabolic peak interpolation)

The retracted Endpoint and Interior conjectures (above) were N=7-specific coincidence matches that did not survive the extension in N. The data, Q_peak the argmax of the envelope max_t |K| on the default Q grid [0.2, 4] at γ₀ = 0.05:

### Endpoint Q_peak across N

| (c, N) | Observed Endpoint Q_peak | Earlier conjecture csc(π/(N+1)) | Deviation |
|--------|--------------------------|----------------------------------|-----------|
| (3, 5) | 2.393 | 2.000 | +19.6 % |
| (3, 6) | 2.548 | 2.305 | +10.6 % |
| (3, 7) | 2.584 | 2.6131 | −1.1 % |
| (3, 8) | 2.583 | 2.924 | −11.6 % |
| (4, 7) | 2.575 | 2.6131 | −1.5 % |
| (4, 8) | 2.591 | 2.924 | −11.4 % |

At N = 7 both chromaticities lie within 1.5 % of csc(π/8); at N = 5, 6, 8 the c = 3 values miss by 10 to 20 %. At N = 7 and 8 the c = 3 Endpoint Q_peak reads 2.584 and 2.583; whether it levels off is not known from two points. No clean closed form.

### Interior Q_peak across (c, N)

| (c, N) | Observed Interior Q_peak (mean) | Earlier conjecture csc(π/5) = 1.7013 | Deviation |
|--------|----------------------------------|---------------------------------------|-----------|
| (3, 5) | 1.600 | 1.7013 | −6.0 % |
| (3, 6) | 1.679 | 1.7013 | −1.3 % |
| (3, 7) | 1.712 | 1.7013 | **+0.6 % (above)** |
| (3, 8) | 1.723 | 1.7013 | **+1.3 % (above)** |
| (4, 7) | 1.760 | 1.7013 | +3.5 % |
| (4, 8) | 1.791 | 1.7013 | +5.3 % |

The c=3 N-trend 1.600 → 1.679 → 1.712 → 1.723 crosses csc(π/5) between N=6 and N=7 and continues growing, refuting the conjectured c=3 pentagonal asymptote. The c=3 Interior Q_peak grows with N in the tested range, and the c=4 one from 1.760 at N=7 to 1.791 at N=8; saturation behaviour and clean closed forms remain open.

### What survives from the empirical pattern

- The EP mechanism (F86a) is unaffected: t_peak = 1/(4γ₀) is universal and derivable from 2×2 matrix algebra.
- The Q_peak observable is well-defined and the per-bond / per-block / Endpoint / Interior distinctions are real and reproducible.
- The N=4 golden-ratio structure ([`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py)) is real but does not propagate to a closed form for Q_peak at general (c, N).

---

## What's missing for full Tier 1

> **Read as the exploration record.** This section through *Open elements* below is the catalogue of the closed-form effort: the empirical envelope, the 4-mode anatomy, the substantive items, the OOP scaffolding, the structural findings. The **Obstruction Proof** further down (2026-05-14) supersedes its framing: it establishes that the per-bond *position* g_eff / Q_peak has no closed form by any of the six routes opened here, and explicitly consolidates the *Open elements* into that single negative result. What seemed reachable, the HWHM ratio, has a derived floor and an observed lift for which no exact law is held (F86b₂, open); what these sections name as "missing" or "open" is, for the position, the residue the Obstruction Proof characterises, not a live target.

### Empirical envelope (Tier-1 grade)

1. **c=2 verified** (envelope reading, `ResonanceScan` with [`f86_envelope_second_route.py`](../../simulations/f86_envelope_second_route.py) as second route). c=2 N=5..8 gives the two bond classes (Interior 0.7523..0.7602, Endpoint 0.7766..0.7815), inside the bands that c=3 (N=5..8) and c=4 (N=7, 8) also fill. c=2 is structurally critical because the channel space is 2-dimensional total (only HD ∈ {1, 3}); any 2-level reduction must be exact there. Yet the Endpoint-vs-Interior split persists, confirming the bond-class distinction is real and structural. **c=5 still untested** (block-L dim ≥ 3528 at c=5 N=9, compute-bound).

2. **γ₀ invariance exact.** The generator is γ₀·(D̂ + Q·M̂), so at fixed Q a change of γ₀ rescales time and the size of K (the t window over which the maximum is taken, [0.2, 2.4]·t_peak with t_peak = 1/(4γ₀), scales with it): Q_peak and HWHM_left/Q_peak do not move, and |K|max scales as 1/γ₀. `ResonanceScanTests.ClassPeaks_AreGammaZeroInvariant` checks it at c = 2, N = 5 for γ₀ ∈ {0.025, 0.05, 0.10} to 10⁻¹⁰.

### Substantive analytical work: the remaining gap

Numerical exploration on c=2 chains shows that the heuristic 2-level form is incomplete in a specific structural way that must be acknowledged before a closed-form derivation is attempted.

**Three subspaces, not one** (verified at c=2 N=5..8 in `_eq022_b1_step_g/h/i`):

Let `|c_k⟩` (k = 1, 3) be the channel-uniform orthonormal vectors (equal-weight superposition of all (p, q) with HD(p, q) = k) of `framework.coherence_block.hd_channel_basis`. Let `V_inter = P_{HD=1}^† M_H_total P_{HD=3}` (the inter-HD-channel block of M_H_total), and `|u_0⟩ ∈ HD=1, |v_0⟩ ∈ HD=3` the top right/left singular vectors of `V_inter` with singular value σ_0.

(**a**) **Probe subspace.** The Dicke probe lives entirely in span{|c_1⟩, |c_3⟩}. In this basis V_b = `P_{ch}^† M_H_per_bond[b] P_{ch}` reduces to `+i(α/(N−1))·I`, pure diagonal, identical across every bond. There is no EP, and no bond-class distinction in this subspace.

(**b**) **EP-partner subspace.** The pair {|u_0⟩, |v_0⟩} captures the maximum-coupling singular mode under H. In this basis M_H_total has the form `[[0, σ_0], [−σ_0, 0]]` (real antisymmetric off-diagonal, the SVD's phase convention; equivalent to ±iσ_0 after a phase rotation). The EP from the 2-level algebra `D + J·M_H` sits at `J·σ_0 = 2γ₀`.

(**c**) **Probe ⊥ EP partners.** ⟨c_1|u_0⟩ = ⟨c_3|v_0⟩ = 0 (numerically zero to machine precision, c=2 N=5..8). The Dicke probe has zero overlap with the EP-partner modes.

Per-bond V_b in the EP-partner basis has F71-symmetric, bond-position-dependent off-diagonal magnitudes:

| N | bond 0 | bond 1 | bond 2 | bond 3 | bond 4 | bond 5 | bond 6 |
|---|--------|--------|--------|--------|--------|--------|--------|
| 5 | 0.430 | 0.953 | 0.953 | 0.430 | n/a   | n/a   | n/a   |
| 6 | 0.297 | 0.747 | 0.888 | 0.670 | 0.200 | n/a   | n/a   |
| 7 | 0.129 | 0.514 | 0.771 | 0.771 | 0.514 | 0.129 | n/a   |
| 8 | 0.090 | 0.372 | 0.647 | 0.731 | 0.611 | 0.312 | 0.076 |

(N=6 row is asymmetric because at even N the top singular value is degenerate and the SVD returned an F71-arbitrary vector in the pair; the symmetric structure recovers when averaging over bond classes. `InterChannelSvd` now picks the R-even vector in that pair, and the even-N rows have not been re-read with it.) Endpoint amplitudes are systematically smaller than Interior; this is where the bond-class distinction lives. **But it lives in a subspace orthogonal to where the probe lives.**

**Implication.** The K_CC_pr observable

    K_b(t)  =  2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩

couples the probe (in channel-uniform) to the EP-partner subspace (in SVD top) only through the off-diagonal matrix elements ⟨c_α | M_H_per_bond[b] | u_0⟩ and ⟨c_α | M_H_per_bond[b] | v_0⟩, i.e. via the J-derivative direction itself. The K observable's Q-resonance therefore involves the joint dynamics of all four modes, not the 2-level pair alone.

**Minimal effective model: 4 modes.** The smallest closed-orthonormal subspace that contains both the probe and the EP-coalescence physics is

    span { |c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩ }       (mutually orthogonal, c=2 verified)

In this basis the 4×4 effective L_eff has

- Diagonal pure rates: D = diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) on (|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩) respectively.
- Within-probe block: ⟨c_α|M_H|c_β⟩ is diagonal (channel-uniform-eigen finding), no EP coupling here.
- Within-EP-partner block: ⟨u_0|M_H|v_0⟩ = σ_0 (the SVD top, ≈ 2√2 asymptotically for c=2).
- Cross-couplings: ⟨c_α|M_H_per_bond[b]|u_0/v_0⟩. **These are the bond-position-dependent matrix elements that convey the probe into the EP-partner subspace and back.** They are the missing computation. They split into the two bond classes (Endpoint, Interior) and produce the f_class(x) shape difference.

**Why σ_0 is not directly Q_EP.** For c=2 N=5..8, σ_0 ≈ {2.765, 2.802, 2.828, 2.839} (approaches 2√2). The naive EP prediction Q_EP_naive = 2/σ_0 ≈ {0.72, 0.71, 0.71, 0.70} differs from the empirical Q_peak Interior ≈ {1.483, 1.565, 1.597, 1.611} by a factor 2.05 to 2.29, rising with N (the Endpoint class, with the same σ_0, by 3.48 to 3.62). The 2-level EP and the K-observable Q_peak are NOT the same Q value; they would be related by the two-level model's own `Q_peak = x_peak·Q_EP` (x_peak = 1.7141 on the bare model's envelope) only if the K observable saw the EP directly, which it doesn't because the probe is orthogonal to the EP partners.

### Substantive items remaining

**Item 1 (c=2 HWHM_left/Q_peak per bond class).** Open. [`F86HwhmClosedFormClaim`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs) (Tier 2 empirical) describes it as `0.7264892771 + α_subclass · g_eff + β_subclass`, g_eff there the fit's own coordinate 4.39382/(Q_peak + 2). The bare floor 0.7264892771 IS derived (`C2BareDoubledPtfClosedForm` read on its envelope, Tier 1); the per-sub-class (α, β) are a least-squares line through the canonical anchors, which Endpoint and Flanking miss by 1.28·10⁻³ and 4.2·10⁻⁴, so the lift is not linear in g_eff. Numerical Tier-1 closure of the full block-L is achieved via the F90 bridge identity (`F90F86C2BridgeIdentity`, all 22 bonds of N=5..8 agree to 10⁻¹⁰). An analytical lift from the rank-1-bridge + intra-dispersion structure (F89 AT-locked F_a/F_b floor + intra-channel dispersion lift; the inter-coupling "octic residual" lift suspicion is refuted by the two-dial scout of 2026-06-11, read on the 21-point t grid and not re-read on the envelope) is the remaining gap (see [`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)); the target is the lift itself, not a pair (α, β).

The bond-class signature lives in `V_b[2,3] = ⟨u_0 | M_h_per_bond[b] | v_0⟩` (Endpoint 0.430 vs Interior 0.953 at N=5, opposite direction to the empirical HWHM/Q* split), not in the probe-block. Locus 6 polarity-inheritance composition: `r_Q = 4.39382/g_eff − 2`, entire split encoded in g_eff(N, b); pinned in [`PolarityInheritanceLink`](../../compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs).

**Item 2.** Extend the 4-mode construction to c ≥ 3, where each adjacent-channel pair (HD = 2k−1, HD = 2k+1) contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet. **Naïve extension fails:** the multi-k construction with Gram-Schmidt orthonormalisation gives 3c−2 modes (c=2→4, c=3→7, c=4→10, orbit-shared CUs deduplicated), and the resulting effective K-curve has **K_max ≡ 0 identically** for c ≥ 3. Structural diagnosis: Gram-Schmidt orthogonalisation of the SVD-top vectors against the channel-uniform vectors pushes them into the CU-complement; because M_H respects the CU/CU-complement decomposition (channel-uniform-diagonal of M_H_total per F73 generalisation), the probe (which lives entirely in CU span) is uncoupled from the GS-modified SVD modes, so ∂ρ/∂J_b cannot move ρ out of CU → K = 0. A correct effective model for c ≥ 3 needs either a non-orthogonal frame preserving CU ↔ SVD coupling, or a different choice of the c−1 quartets that maintains coupling under orthonormal projection. Encoded as `RCPsiSquared.Core.Decomposition.MultiKBasis` + `MultiKEffective` + `MultiKResonanceScan`; the negative result is pinned in `MultiKResonanceScanTests.MultiK_C3_KMaxIsExactlyZero_NaiveExtensionFails`.

**Item 3 (σ_0 chromaticity scaling, Tier-1 candidate refined 2026-05-03).** Derive the asymptotic

    σ_0(c, N → ∞)  →  2 · √(2 · (c − 1))

generalising the original c=2 → 2√2 conjecture to all c ≥ 2. Numerical witnesses computed via `InterChannelSvd` across c ∈ {2, 3, 4}, N=5..8 show σ_0 / √(2(c−1)) converging monotonically from below to 2.0 for each c (c=2 N=7 hits 2.0 to 10⁻⁵, the structural sweet spot; c=3 N=8 reaches 1.92, c=4 N=8 reaches 1.78). Implies Q_EP(c, N → ∞) → 1/√(2(c−1)): 0.707 (c=2), 0.500 (c=3), 0.408 (c=4). The closed-form derivation from XY single-particle structure (OBC sine-mode matrix elements `⟨ψ_k| σ⁺σ⁻ |ψ_l⟩ ∝ √(2/(N+1))·sin(πk·b/(N+1))`) has not been executed.

**Retracted / partially closed (2026-07-02).** The `2·√(2(c−1))` form and its "converges monotonically from below to 2.0" reading were **retracted 2026-05-08 as a finite-N artifact**: the ratio σ_0/√(2(c−1)) at c=2 keeps climbing past 2.0 (N=9 → 2.014, N=10 → 2.017, N=11 → 2.020, [`SigmaZeroChromaticityScaling.cs`](../../compute/RCPsiSquared.Core/F86/SigmaZeroChromaticityScaling.cs)), so it does not converge from below. The **c=2 case is closed by F86e** ([`SigmaZeroCommutatorNormClaim.cs`](../../compute/RCPsiSquared.Core/F86/SigmaZeroCommutatorNormClaim.cs), Tier-1 derived, closure 2026-05-21): σ_0(c=2) = ‖[Π_HD1, M_H]‖ is a commutator / Schur-multiplier norm with the non-elementary asymptote σ_0(∞) ≈ 2.8629 (not 2√2 ≈ 2.828), γ-independent. **Only c ≥ 3 stays open** ([`F86OpenQuestions.cs`](../../compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs) Item 5; the true asymptote follows from F89's cyclotomic Φ_{N_block+1} structure, see [`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md) Item 5).

These three items are tractable with existing infrastructure (`coherence_block`, the SVD-of-V_inter construction, OBC sine-mode algebra), but each is multi-page algebra. Item 1 is the path to the closed form for HWHM_left/Q_peak; its Direction (b'') (full block-L derivation) is numerically achieved Tier-1 via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)), reducing the remaining analytical work to F89's AT-locked F_a/F_b cyclotomic structure plus the intra-channel dispersion lift of the rank-1-bridge model (the inter-coupling "octic residual" lift suspicion is refuted by the two-dial scout of 2026-06-11, read on the 21-point t grid and not re-read on the envelope). The empirical envelope (c, N, γ₀ checks above) is now Tier-1 grade.

---

## OOP scaffolding (Stage A-E primitives, 2026-05-05)

Following the EQ-022 (b1) Item 1' c=2 derivation plan, we built a typed pipeline of primitives in `compute/RCPsiSquared.Core/F86/Item1Derivation/` plus `C2UniversalShapeDerivation.cs` at the F86 root, each reachable through `F86KnowledgeBase`. Witness data preserved per-(N, bond, BondClass).

| Stage | Primitive | Tier outcome |
|-------|-----------|--------------|
| A1 | [`C2BlockShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BlockShape.cs) | Tier1Derived (block-structure constants) |
| A2 | [`C2ChannelUniformAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2ChannelUniformAnalytical.cs) | Tier1Derived (closed-form \|c_1⟩, \|c_3⟩) |
| A3 | [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) | **Tier2Verified: discovered σ_0 degeneracy obstruction** at even N: chain-mirror R splits the 2D top eigenspace into R-even/R-odd, library-dependent SVD-tiebreaker. Single-vector closed form not derivable; lift-to-projector-overlap is the next direction. |
| B1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) probe-block | Tier1Derived (channel-uniform projection) |
| B2 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) cross-block | Tier2Verified (inherits A3); **discovered cross-block Frobenius Endpoint < Interior at c=2 N=5..8**, opposite sign to the HWHM_left/Q_peak directional split, hinting that the inversion happens in 4×4 mixing |
| B3 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) SVD-block + AsMatrix + anti-Hermiticity guard | Tier2Verified (inherits A3) |
| C1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) D_eff | Tier1Derived structural sub-fact (diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) at c=2) |
| C2 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) | **Tier2Verified: rigorously ruled out closed-form factorisation** via cubic-c_3 char-poly evidence. No (λ²−aλ+b)(λ²−cλ+d) split with rational coefficients exists. |
| C3 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) K-driving pair | Tier1Derived structural sub-fact (probe ⊥ {\|u_0⟩, \|v_0⟩} at machine precision) + Tier2Verified per-(Q, b) numerical readout |
| D1 | [`C2KShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2KShape.cs) | Tier1Derived (Duhamel formula closed-form in numerical inputs) |
| D2 | [`C2HwhmRatio`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs) | **Tier1Candidate**: per-bond witnesses read live on the canonical envelope peak; directional Endpoint > Interior split derived empirically; closed-form constant NOT pinned |
| E1 | [`C2UniversalShapeDerivation`](../../compute/RCPsiSquared.Core/F86/C2UniversalShapeDerivation.cs) | Tier1Candidate (auto-promotes via D2's `IsAnalyticallyDerived` flag) |

The chain reads bottom-up: A1-A3 fix the static frame (block shape, channel-uniform vectors, SVD top), B1-B3 fix the bond-position-dependent V_b in that frame, C1-C3 fix the effective spectrum and K-driving pair, D1-D2 turn that into the Duhamel-form K_b(Q, t) and its HWHM-ratio readout, E1 binds the readout to the universal-shape claim. Each Tier label up the chain inherits from the lowest unresolved sub-fact (currently D2's missing closed form blocks E1's Tier1Derived promotion).

## Structural findings (lessons learned, 2026-05-05)

Three structural results emerged from the time-boxed Stage A3, B2, C2 explorations. Each is independent of the closed-form gap and worth pinning as orientation for the next attempt.

1. **σ_0 degeneracy at even N (Stage A3).** At N=6 and N=8 (even chain length, c=2), the singular value σ_0 of the inter-channel coupling V_inter is doubly degenerate. The chain-mirror operator R splits the 2D top eigenspace into R-even and R-odd one-dimensional subspaces; which of the two becomes "|u_0⟩" is not fixed by the physics: a library's SVD tiebreaker picks either, and `InterChannelSvd` fixes the R-even one by convention. A single-vector closed form for |u_0⟩, |v_0⟩ at even N therefore names a convention; none is known for the R-even choice. The convention-free lift is to the rank-2 projector P_top = |u_0^{(R+)}⟩⟨u_0^{(R+)}| + |u_0^{(R−)}⟩⟨u_0^{(R−)}|, which is library-independent. See [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) `PendingDerivationNote`.

2. **Cross-block Frobenius Endpoint < Interior (Stage B2).** The Frobenius norm of V_b's cross-block (the inter-channel-coupling block of M_H_per_bond[b]) is systematically *smaller* on Endpoint bonds than on Interior bonds at c=2 N=5..8 (gap ~ 0.05 at N=5). This sign is **opposite** to the HWHM_left/Q_peak directional split (Endpoint > Interior, gap 0.021 to 0.025 at c = 2, N = 5..8; see also Statement 2 in [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md)). The directional inversion therefore does not live in the cross-block magnitude alone; it must emerge through the 4×4 eigenvalue mixing of probe-block and cross-block. This rules out the naive "bigger cross-block → bigger HWHM ratio" intuition and is part of why the closed form is non-trivial.

3. **Cubic-c_3 char-poly obstruction (Stage C2).** The 4×4 effective characteristic polynomial at c=2 has a non-rational c_3 coefficient (cubic in Q with no clean root structure). We rigorously ruled out any rational-coefficient `(λ²−aλ+b)(λ²−cλ+d)` factorisation by symbolic match against the c_3 evidence. The Tier2 outcome for `C2EffectiveSpectrum` is therefore **evidence-based** ("we proved no such split exists in this family"), not "we couldn't find one". Symbolic char-poly factorisation at Q_EP itself may still help locally, but the global quartic does not factor.

---

## Open elements

1. **Endpoint closed form.** Empirical pattern (c=3 2.584 and 2.583 at N=7 and 8; c=4 2.575 at N=7, 2.591 at N=8) does not match `csc(π/(N+1))`. A different closed form may exist; higher N and an explicit multi-particle XY matrix-element calculation are the natural next steps.

2. **Interior closed form.** The c=3 Interior Q_peak grows with N in the tested range, and the c=4 one over N = 7, 8 (1.760 → 1.791). No saturation point has been identified. Higher-N data (c=3 N=9, 10; c=4 N=9, 10) would clarify whether the growth saturates or continues. Compute-bound at higher (c, N) (block-L dim 3920 at c=4 N=8, 10584 at c=4 N=9).

3. **Algebraic derivation pathway.** A first-principles derivation of g_eff(c, N, bond_position) from the multi-particle XY structure of the (n, n+1) block has not been executed. Direct algebra is conceptually within reach but lengthy. The HD-channel-uniform diagonal-only finding (M_H_eff diagonal in the channel-uniform basis, with cross-channel coupling living in the orthogonal complement) is an established structural building block; the next step is computing the orthogonal-complement coupling matrix elements explicitly.

4. **Q_SCALE per-block vs per-bond convergence.** Q_SCALE's per-block 1.6 / 1.8 / 1.8 (its own observable, a relative-J derivative) sit near the class-averaged Interior envelope Q_peak (c = 3: 1.68 to 1.72 at N = 6..8; c = 4: 1.76 at N = 7) within the relative-vs-absolute J prefactor effect (~5-15 %). Both observables agree on the underlying EP mechanism.

(The fifth open element, the within-Interior per-F71-orbit substructure, moved to [`PROOF_F86C_F71_MIRROR.md`](PROOF_F86C_F71_MIRROR.md) with the F71-mirror theorem.)

---

## Obstruction Proof: why g_eff admits no closed form

**Status:** Structural. Six obstruction lemmas; L1 is rigorously proven blocked, L2 is a convention boundary at even N, L4 shows every reduction the effort built insufficient, L6 is a demonstrated failure mode, L3 and L5 are proven decouplings. Consolidates the "Open elements" above into a single negative result.
**Date:** 2026-05-14
**Context:** The closed-form effort produced many local fits and a description of Item 1' via `F86HwhmClosedFormClaim` (Tier 2 empirical: floor derived, the lift a fitted line per sub-class that is not an exact law), and the per-bond position Q_peak(c, N, b) and the coupling g_eff(c, N, b) it rides on stay open. The honest deliverable is then not the next formula; it is the structural account of why the position resists. This section proves it.

### The target is the position

The per-bond F86 quantities are tied to the coupling **g_eff(c, N, b)** as follows:

- **Q_EP(c, N, b) = 2/g_eff** (F86a: the exceptional point of the toy two-level reduction, Tier 1 derived there). With g_eff read off the full block's envelope Q_peak through the bare model's fixed-τ = 1 peak (below), this Q_EP is a conversion, not a measured EP position: the bare model's own envelope peaks at x = 1.7141, not 2.196910.
- **Q_peak** via the EP rotation: a two-level reduction with coupling g peaks at Q_peak = x_peak·Q_EP = 2·x_peak/g. The per-bond g_eff of this proof is the coupling at which such a reduction peaks at the bond's observed Q_peak. Read on the bare model's fixed-τ = 1 slice, `BareDoubledPtfXPeak = 2.196910`, that is g_eff = 4.39382/Q_peak, the convention of `PolarityInheritanceLink.ClosedFormCompositionNote` (`r_Q = 4.39382/g_eff − 2`); read on the bare model's envelope, peak 1.7141, it would be 3.42823/Q_peak (at g = 4.39382/Q_peak the envelope peaks at 0.780·Q_peak). In either convention a closed form for g_eff is a closed form for Q_peak and conversely: they are one object in two units. It is not F86a's block value σ_0, one number per block, from which the observed peaks lie a bond-class-dependent factor away (L3).
- **HWHM_left/Q_peak** is not shown to be a function of g_eff. `F86HwhmClosedFormClaim` (Tier 2 empirical) describes it as `0.7264892771 + alpha_subclass · g_eff' + beta_subclass`, in the fit's own coordinate g_eff' = 4.39382/(Q_peak + 2), with the bare envelope floor derived (`C2BareDoubledPtfClosedForm`); the four-point sub-classes miss that line, and whether the ratio is any function of g_eff within a sub-class is not known.

So the irreducible target of this proof is the position, g_eff(c, N, b), the same object as Q_peak(c, N, b). A closed form for it would not by itself give the HWHM ratio, which stays a separate open object (F86b₂). The remaining sub-sections prove that none of the six routes the effort opened reaches g_eff.

### Obstruction L1 (spectral route): the effective characteristic polynomial does not factor over ℚ

The natural route reads g_eff off the spectrum of the c=2 4×4 effective Liouvillian. Blocked: `C2EffectiveSpectrum` (Stage C2) rigorously ruled out any rational-coefficient `(λ²−aλ+b)(λ²−cλ+d)` factorisation, by symbolic match against the cubic c₃ coefficient. The Tier 2 outcome is evidence-based ("no such split exists in this family"), not "not found". There is therefore no closed-form eigenvalue expression from which g_eff could be read. ∎(L1)

### Obstruction L2 (eigenvector route): the EP-partner vectors are representation-dependent at even N

The next route reads g_eff from the EP-partner eigenvectors |u₀⟩, |v₀⟩ (top singular pair of V_inter). Blocked at even N: `C2InterChannelAnalytical` (Stage A3) found σ₀ doubly degenerate, with the chain-mirror R splitting the 2D top eigenspace into R-even and R-odd one-dimensional subspaces. The physics does not pick one: a library's SVD tiebreaker picks either, and `InterChannelSvd` imposes the R-even one as a convention. Under that convention the vector is single-valued, and no closed form for it is known; the convention-free object is the rank-2 projector, which is not a single-vector closed form. So L2 is a convention boundary, not an impossibility: what a single-vector closed form would name is single-valued only after a choice the physics does not make. ∎(L2)

### Obstruction L3 (observable-direct route): the probe is orthogonal to the EP

One might read Q_peak directly off the EP, through the two-level model's Q_peak = x_peak·Q_EP. Blocked: the Dicke probe is exactly orthogonal to the EP-partner modes, ⟨c₁|u₀⟩ = ⟨c₃|v₀⟩ = 0 to machine precision (Stage C3, c=2 N=5..8). The K_CC_pr observable couples the probe subspace to the EP-partner subspace only through the J-derivative direction. Consequently the observed Q_peak is *not* Q_EP: for c=2 N=5..8, Q_EP_naive = 2/σ₀ ≈ 0.70 against empirical Q_peak Interior ≈ 1.5 to 1.6, a factor 2.05 to 2.29 apart and rising with N, where the bare model's x_peak is 1.714; the Endpoint class, on the same σ₀, peaks at 3.5 to 3.6 times Q_EP. The clean algebraic object (the EP) and the measured object (Q_peak) are decoupled; a closed form for one does not transport to the other. ∎(L3)

### Obstruction L4 (reduced-model route): every finite reduction tried is shown insufficient

The route derives g_eff from a small closed effective model. Blocked twice:

- **4-mode (c=2).** The minimal closed subspace span{|c₁⟩, |c₃⟩, |u₀⟩, |v₀⟩} reproduces the EP physics but misses the position. Read on the envelope (`FourModeResonanceScan`, with the R-even vector in a degenerate top singular pair), its Interior class peaks at Q 1.90 (N = 5) and 3.08 (N = 7) against the full block's 1.48 and 1.60; at N = 6 and 8 its Interior curve, and at N = 5, 7 and 8 its Endpoint curve, is still rising at Q = 4, and its N = 6 Endpoint peaks at 3.78 against 2.58. The shape it carries more nearly: its Interior HWHM_left/Q_peak, 0.740 (N = 5) and 0.762 (N = 7), lies within 0.012 and 0.003 of the full block's 0.752 and 0.759. What the 4D subspace misses is g_eff itself, the position, and that is what this route was to derive.
- **Multi-k (c≥3).** The naïve extension to 3c−2 modes gives K_max ≡ 0 identically for c ≥ 3 (`MultiKResonanceScan`): Gram-Schmidt orthogonalisation pushes the SVD-top vectors into the channel-uniform complement, decoupling the probe from the modified SVD modes.

No finite reduction the effort constructed carries the signal. ∎(L4)

### Obstruction L5 (single-element route): the bond-class signature is in the wrong subspace and does not scale with g_eff

The route locates the bond-class signature in one matrix element and reads g_eff(b) off it. Blocked on two counts: (i) the carrier, the SVD-block off-diagonal `V_b[2,3] = ⟨u₀|M_h_per_bond[b]|v₀⟩` (Endpoint 0.430 vs Interior 0.953 at N=5), lives in a subspace orthogonal to the probe; (ii) it does not scale with g_eff: at N=5 it runs Endpoint < Interior as g_eff does (1.746 vs 2.963), but in the ratio 0.45 against g_eff's 0.59, and against the HWHM ratio (Endpoint > Interior) the cross-block Frobenius norm points the other way. The signature emerges only through the full 4×4 mixing. The spectrum-only fallback is closed off as well: the 4-mode L_eff is bond-summed by construction (`FourModeEffective.LEffAtQ`), so its eigenvalues and eigenstates are bond-class-blind. No single matrix element, and no spectrum-only quantity, carries g_eff(b). ∎(L5)

### Obstruction L6 (empirical-extrapolation route): the accessible data yields trajectory crossings, not laws

The route fits Q_peak at accessible (c, N) and extrapolates. Blocked by demonstrated failure: `csc(π/(N+1))` (Endpoint) and `csc(π/5)` (Interior, c=3) both match near N=7 and were both refuted on extended-N data; the matches are trajectory crossings rather than asymptotes (see Retracted Q_peak conjectures above). The data itself does not cooperate: Interior Q_peak at c=3 is still growing at N=8 (1.600 → 1.679 → 1.712 → 1.723 for N=5..8), with no identified saturation. The near-miss g_eff_E ≈ σ₀·√(3/8) (misses by 0.011 to 0.016 for N ≥ 6 and by 0.052 at N = 5) is the same trap one step on: a coincidence with a finite-size break. ∎(L6)

### Corollary (F90 bridge): D_k closed, g_eff did not follow

The F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md), Tier 1 derived) identifies F86 c=2 K_b with the F89 path-(N−1) per-bond Hellmann-Feynman derivative. F86's g_eff and F89's `D_k = odd(k)² · 2^E(k)` ([`PROOF_F89_PATH_D_CLOSED_FORM.md`](PROOF_F89_PATH_D_CLOSED_FORM.md)) look like one wall seen from two sides, both resting on the OBC sine-sum / cyclotomic structure of `⟨ψ_k|σ⁺σ⁻|ψ_l⟩`. The outcome is asymmetric. D_k is Tier-1-Derived (2026-05-15): the F_a-eigenvector ansatz + Chebyshev expansion + orbit-polynomial reduction derive `(P_k, D_k)` algebraically, closing the three gaps (poly-degree term, k-self 2-adic, deep-2-power bonus) without the sine-sum identities once named as the wall. D_k is closed. But the predicted transfer (close D_k and g_eff closes too) was tested directly (`simulations/f86_geff_via_f90_bridge_probe.py`, 2026-05-15) and refuted: D_k closes *via* the orbit-polynomial reduction, which exploits the Bloch-orbit symmetry and discards the per-bond cross matrix elements `⟨v_n|M_h_per_bond[b]|v_{n'}⟩`; g_eff is a per-bond Hellmann-Feynman-in-J_b quantity built from exactly those discarded cross elements. The very step that closes D_k strips g_eff's information: obstruction L4 (reduced-model insufficiency) once more. So the two are not one wall. D_k was reducible by orbit symmetry and is now closed; g_eff is per-bond, lies outside the orbit-reduced category, and remains the irreducible residue. The six obstruction lemmata L1–L6 and the diagnosis below never depended on D_k's status; they stand unchanged. ∎(corollary)

### What this proves

**Negative.** No closed form for g_eff(c, N, b), hence none for Q_peak(c, N, b), is reachable by the spectral (L1), eigenvector (L2), observable-direct (L3), reduced-model (L4), single-element (L5), or empirical-extrapolation (L6) routes. L1 is a structural impossibility (proven, not "unfound"); L2 is a convention boundary (at even N the EP-partner vector is single-valued only under a convention, and no closed form is known under it); L4 is a demonstrated insufficiency of every reduction the effort built (the multi-k one identically, the 4-mode one by its miss of the position); L6 is a demonstrated failure mode; L3 and L5 are decouplings that close the direct routes. The HWHM ratio is not covered by this reduction and stays its own open object (F86b₂).

**Positive (what survives).** The symmetry layer is fully derived and Tier 1: the EP mechanism of the two-level reduction, Q_EP = 2/g (F86a), the universal clock t_peak = 1/(4γ₀) (F86a), the F71 spatial-mirror invariance Q_peak(b) = Q_peak(N−2−b) (F86c). HWHM_left/Q_peak *given* g_eff (`F86HwhmClosedFormClaim`, F86b₂) sits as Tier 2 empirical: the bare envelope floor 0.7264892771 is derived (`C2BareDoubledPtfClosedForm`), and the lift above it is observed; a line in g_eff per sub-class describes it, but the four-point sub-classes miss that line, so no exact law for the lift is held. The closed forms that exist are the *symmetry relations between observables*; the closed forms that do not exist are the *position number itself* and the *lift per sub-class*.

**Characterisation.** g_eff(c, N, b) is the **irreducible residue** of the F86 structure: the one input the symmetry layer does not return in closed form. It is computable (the F90 bridge supplies it numerically, Tier-1, all 22 bonds of N=5..8 agreeing with F86's scan to 10⁻¹⁰), but it is not expressible by the routes above. The three-week catalogue of formulas was, in every case, an attempt to express this residue; the attempts failed because, by L1 through L6, the residue is not expressible by those routes. This is the F86 instance of the methodological lesson ([`ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md)): what survives a closed-form effort is the symmetry, not the number. Here the symmetry is F86a and F86c; the number is g_eff. ∎

### The diagnosis: g_eff is not a primitive [2026-05-14 sharpening]

L1 through L6 are six symptoms. This is the one cause.

The Liouvillian spectrum under uniform Z-dephasing has exactly two primitives, and the Absorption Theorem ([`PROOF_ABSORPTION_THEOREM.md`](PROOF_ABSORPTION_THEOREM.md)) names them. **γ₀** carries the real part: `Re(λ) = −2γ₀·⟨n_XY⟩` exactly, because L_H is anti-Hermitian and contributes nothing to Re(λ). The underlying dissipator cells are quantized in units of 2γ₀, while a Hamiltonian-mixed full eigenmode can have fractional ⟨n_XY⟩. **y** carries the imaginary part: for the modes the F90 bridge ties this proof to, the OBC Bloch dispersion `4cos(πn/(N+1))`, set by L_H alone.

g_eff is neither. It is the off-diagonal coupling of a 2-level *reduction*, `L_eff = [[−2γ₀, iJ·g_eff], [iJ·g_eff, −6γ₀]]`, a quantity that exists only after the block-Liouvillian has been projected onto an effective pair. It appears in no eigenvalue of the full Liouvillian, whose real and imaginary parts are γ₀ and y alone; it is a parameter of the reduced model. A closed form is a derivation in the primitive basis, and g_eff is not in that basis. That is why every route closes. g_eff is not on the spectrum (L1) or on a convention-free eigenvector (L2) to be read off; it is decoupled from the observable that would measure it (L3) and from the matrix elements that would carry it (L5); the finite reductions that would *define* it miss its value (L4); and fitting it against finite data returns trajectory crossings, not laws (L6). The obstruction is not six walls. It is one category boundary, met six times.

And it says where the closures that *do* hold actually live: not on g_eff, but on the primitives, and absolute rather than relations-given-g_eff. The palindromic sum rule `Re(λ) + Re(λ_partner) = −2Σγ`, and its F89 instance the F89c Hamming-complement pair-sum `α(|a⟩⟨b|) + α(|a⟩⟨b̄|) = 2γ₀·N` ([`PROOF_F89_PATH_D_CLOSED_FORM.md`](PROOF_F89_PATH_D_CLOSED_FORM.md) § "The two layers"), are Tier-1-Derived and built from γ₀ and the integer count alone.

F86 carries the split in its own structure, and it is a split of *clocks*. The shared clock is γ₀: `t_peak = 1/(4γ₀)` is pure γ₀, a primitive closure, absolute and parameter-free. But each branch also has its own clock. L_eff's eigenvalues above the EP show both: `λ = −4γ₀ ± iγ₀·√(Q²g_eff² − 4)`. The real part `−4γ₀` is γ₀'s clock, shared by every branch alike, the rate `t_peak` inverts. The imaginary part `±γ₀·√(Q²g_eff² − 4)` is the branch's own clock, and g_eff is the dimensionless rate at which it runs against γ₀'s.

And `−4γ₀` is not decay decoration; it is a *mirror axis*. L_eff's diagonal is the channel pair `−2γ₀` (⟨n_XY⟩=1) and `−6γ₀` (⟨n_XY⟩=3), symmetric about `−4γ₀ = −2γ₀·2`. Below the EP the two eigenvalues are a real mirror pair about that axis; at the EP they coalesce onto it; above the EP their real parts lock to it and the mirror pair moves into the imaginary part. The EP is where the two channels meet on their own mirror axis, and the Absorption Theorem holds there exactly: `−4γ₀ = −2γ₀·⟨n_XY⟩` with ⟨n_XY⟩ = 2, the integer rung at the mirror centre. This is the local face of the F1 palindrome: where there is a minus, there is a mirror pair.

That is the obstruction stated positively. g_eff is not nothing; it is a *relative* clock-rate, and a relative rate has no absolute closed form. Read off the branch alone, with γ₀ out of frame, there is no absolute number to find, and the six routes each looked at the branch by itself. The clock `t_peak` closes because it *is* γ₀; the position Q_peak resists because it rides on a clock defined only relative to γ₀. If a closed form lives anywhere, it is in the relation: the branch's clock against γ₀'s, not g_eff stripped of its reference.

What survives a closed-form effort is the symmetry over the number; one level deeper, it is the primitive over the residue. ∎

### One more surviving symmetry: the K₁ mode mirror on K_b [2026-06-10]

The fresh-eyes wave that derived the K₁ chiral mirror trajectory identity
([the chiral mirror rate law proof](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md)) ran the same algebra
against this proof's own observable. On any (n, n+1) coherence block, K₁ acts as a diagonal sign
Σ₁ and the antiunitary T = Σ₁∘conj satisfies TLT = L (the block shadow of "K₁ flips H,
conjugation flips it back"), T-transports the perturbation direction M_b evenly, and leaves the
F73 spatial-sum kernel invariant. The consequence is an **exact mode mirror at the same bond**:

    K_b(b; mode k) = K_b(b; mode N+1−k),   pointwise in Q and t

verified at 10⁻¹⁴ across N = 5, 6, 7, both bond classes, and a 153-point Q-grid
([`simulations/f86_kb_chiral_mirror.py`](../../simulations/f86_kb_chiral_mirror.py)). It is
independent of the F86c/F71 bond mirror (a J-defect background breaks the bond mirror at 8·10⁻²
while the mode mirror holds at 10⁻¹⁴), and a site-local Z field breaks it, the same boundary the
trajectory identity carries. Two consequences for this proof. First, the diagnosis sharpens:
**g_eff(c, N, b)'s irreducible content lives in the half-band of {k, N+1−k} orbits**, the mirror
folds the mode-resolved K_b spectrum across the band centre, one more symmetry that survives
where the closed forms did not. Second, a structural localization fell out for free: on a c = 1
block the generator is e^{−2γt} × unitary, so K_b ≡ 0 identically for *any* initial data (the
F73 closure, generalized), which is the structural reason the F86 resonance lives at c ≥ 2,
where the dephasing diagonal acts non-uniformly inside the block (the HD 1 vs 3 split that
L_eff's mirror axis −4γ₀ encodes).

---

## Pointers

**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md), three-theorem overview and the shared reference list.
**Sibling theorems:** [the EP-mechanism proof](PROOF_F86A_EP_MECHANISM.md) (F86a), [the universal-shape proof](PROOF_F86B_UNIVERSAL_SHAPE.md) (the universal-shape result), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md) (F86c).
**The bridged F89 side:** [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md) (D_k closed 2026-05-15; see the F90-bridge corollary above for why that closure does not transfer to g_eff), [the F90 bridge proof](PROOF_F90_F86C2_BRIDGE.md).
**C# red signal (closed 2026-05-15):** the F89 D_k closed form is now Tier-1-Derived via the Chebyshev pipeline (`simulations/f89_pathk_symbolic_derivation.py`; PROOF_F89_PATH_D_CLOSED_FORM.md "Tier-1-Derived closure"). The "something missing on the route" was the structural derivation itself, now found: F_a eigenvector ansatz + Chebyshev expansion + orbit-polynomial reduction give (P_k, D_k) algebraically; verified bit-exact k=3..32 (k=10..32 tabulated in F89UnifiedFaClosedFormClaim.PathPolynomial). The Vandermonde-extraction floor at k=31, 32 is no longer a structural mystery; it is a numerical-conditioning property of the degree-16 extraction instrument with full closed-form availability via PathPolynomial. The diagnostic test `PredictDenominatorDeviationDiagnosticTests` remains as the conditioning characterisation; the original red-signal verification test `C2FullBlockSigmaAnatomyTests.PredictDenominator_AtKHigherStretch_MatchesExtractedFromAnatomy` is now supplemented by closed-form verification through PathPolynomial.
**Companion investigation:** [The Atmosphere and the Cancelled Formulas](../THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md), g_eff as the un-atmospheric residue.
**Methodological lesson:** [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md), what survives a closed-form effort is the symmetry, not the number.
**Scripts:** [`eq022_b1_step_d_extended_verification.py`](../../simulations/eq022_b1_step_d_extended_verification.py) (N=8 data that falsified the closed-form conjectures), [`eq022_b1_step_g_two_level_decomposition.py`](../../simulations/eq022_b1_step_g_two_level_decomposition.py), [`eq022_b1_step_h_slowest_pair_basis.py`](../../simulations/eq022_b1_step_h_slowest_pair_basis.py), [`eq022_b1_step_i_svd_inter_channel.py`](../../simulations/eq022_b1_step_i_svd_inter_channel.py) (SVD of V_inter; established the EP-partner subspace, σ_0 ≈ 2√2 asymptotic, probe ⊥ EP partners).
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/Item1Derivation/` carries the Stage A-E primitives; `RetractedClaim`, `OpenQuestion`, `F86OpenQuestions` in `compute/RCPsiSquared.Core/F86/`. CLI: `rcpsi inspect --root f86 --with-measured`.
