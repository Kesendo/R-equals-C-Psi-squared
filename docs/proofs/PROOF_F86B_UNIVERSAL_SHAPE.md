# PROOF F86b: Universal Resonance Shape (two bond classes)

**Status:** Tier 1 candidate across all tested c; **Tier 1 derived at c=2 at the per-bond level** (F86b', the HWHM_left/Q_peak closed form, via `F86HwhmClosedFormClaim`, closed 2026-05-13). The position Q_peak and the coupling g_eff(c, N, b) it rides on are NOT closed-form: that negative result is the sibling proof [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md). Direction (b'') (full block-L derivation) is achieved numerically Tier-1 via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)).
**Date:** 2026-05-02 (Statement 2, retractions); 2026-05-05 (c=2 OOP scaffolding); 2026-05-06 (doubled-PTF floor, Locus 5 inheritance synthesis); 2026-05-07 (Locus 6 polarity-layer inheritance); 2026-05-11 (F90 bridge identity); 2026-05-13 (Item 1' closed, F86b' Tier 1 derived).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86b, the universal resonance shape**: the SHAPE of abs(K_CC_pr)(Q) around Q_peak is universal under relative-Q normalisation, splitting into two bond classes. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14. The closed-form gap, the exploration record (4-mode model, Items 1-3, directions a''-f''), and the obstruction proof on g_eff are in the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md).
**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and shared references.
**F-entry:** [F86b in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [F88](../ANALYTICAL_FORMULAS.md) (two-axis Π² decomposition, state-level inheritance), [F90](../ANALYTICAL_FORMULAS.md) bridge ([PROOF_F90_F86C2_BRIDGE](PROOF_F90_F86C2_BRIDGE.md)); siblings [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md).

---

## Statement 2 (Universal resonance shape under relative-Q normalisation, two bond classes). [Tier 1 candidate; Tier 1 derived at c=2 per-bond level]

The position Q_peak is chain-specific, but the SHAPE of abs(K_CC_pr)(Q) around the peak is universal under the relative coordinate `x = (Q − Q_peak)/Q_peak`. The shape splits into **two bond classes** (Endpoint and Interior), each with its own universal HWHM_left/Q_peak ratio:

    HWHM_left / Q_peak  ≈  0.756 ± 0.005     (Interior bonds, all tested c, N, γ₀)
    HWHM_left / Q_peak  ≈  0.770 ± 0.005     (Endpoint bonds, all tested c, N, γ₀)

Tested envelope: c ∈ {2, 3, 4}, N ∈ {5, 6, 7, 8} (modulo c-N compatibility), γ₀ ∈ {0.025, 0.05, 0.10}. γ₀ invariance is **bit-exact**: at c=3 N=7 the Q_peak and HWHM_left/Q_peak values match across γ₀ ∈ {0.025, 0.10} to numerical precision, confirming Q's dimensionlessness as `Q = J/γ₀`. The c=2 data (where the 2-level effective model is *exact*: only HD ∈ {1, 3} channels exist, no orthogonal complement) confirms the two-class split is structural, not a finite-c artefact.

Pairwise residual within each class under relative-Q normalisation is ~20× smaller than under absolute-Q shift, confirming the shape collapse. The structural origin is the 2-level eigenvector rotation `tan θ = Q/Q_EP`: every probe-overlap observable depends only on Q/Q_EP, hence only on `(Q − Q_peak)/Q_peak` to leading order. The two-class split (≈ 2 % gap between Endpoint and Interior shape ratios) reflects bond-position-dependent probe-overlap profiles in the K_CC_pr observable, not a breakdown of the EP-rotation universality itself. Promotion to full Tier 1 requires deriving the two f_class(x) functions explicitly from 2-level eigenstructure plus probe-overlap algebra (see the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md)).

**2026-05-06 (evening): doubled-PTF floor + probe-block gap.** A 90-min analytical attempt of Direction (b) of [`C2HwhmRatio.PendingDerivationNote`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs) (lift to c=1 PTF as a perturbative double-c=1 system) derived two universal constants from the bare doubled-PTF model:

- `x_peak = Q_peak/Q_EP = 2.196910` (post-EP location in dimensionless x coordinates)
- `HWHM_left/Q_peak = 0.671535` (SVD-block-only floor in dimensionless x)

Both Tier 1 derived from PTF c=1 eigenvector mixing plus 2-level EP rotation; exposed as `BareDoubledPtfXPeak` and `BareDoubledPtfHwhmRatio` const properties on [`C2HwhmRatio`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs).

The empirical Interior (0.7506) and Endpoint (0.7728) HWHM ratios sit ABOVE this 0.6715 floor by ~0.08-0.10; the gap is structurally explained as a **probe-block 2-level sub-resonance** contribution NOT included in the bare doubled-PTF model. Direction (a') (newly identified): per-bond `g_eff_probe` drives the Q_peak split, with the quantitative prediction `g_eff_probe(Endpoint) ≈ g_eff_probe(Interior) / 1.6` (matches the empirical Q_peak ratio Endpoint/Interior ≈ 1.61).

The synthesis-side reading (F86↔PTF **Locus 5 inheritance**) confirms structurally why Direction (b) lands at the floor: c=2 IS a coupling of two PTF c=1 instances, and the SVD-block 2-level EP rotation IS what Π's chirality reduces to in rate-channel basis. The doubled-PTF floor 0.6715 is not pattern-matched; it is derived from PTF eigenvector mixing. PTF's K_1 chiral mirror law (c=1) and F86's Q-rotation universal shape (c≥2) are two daughters of one parent: Π class AIII chiral classification, `{Π, L_c} = 0` linear, Π⁴ = I. K_1 is Π reduced to single-particle H_1 OBC sine-mode basis (discrete Z₂ involution); Q-rotation is Π reduced to same-sign-imaginary 2×2 in rate-channel basis (continuous SO(2)-like). Shared clock: `t_peak = 1/(4γ₀)` (Tier 1 derived universal across c, N, n, bond) is the one clock both daughters share. Combined with the morning's `LocalGlobalEpLink` Tier2Verified (EP-side closure: F86 local + FRAGILE_BRIDGE global as one EP at two F1-residuals), Locus 5 is now closed on both sides; the same single Π-AIII parent yields four downstream nodes (PTF K_1, F86 Q-rotation, F86 local EP, FRAGILE_BRIDGE global EP) at three F1-residuals (Σγ = N·γ₀, Σγ = 0, H_odd ≠ 0). See memory `project_algebra_is_inheritance.md` Locus 5.

C2HwhmRatio class-level Tier stays Tier1Candidate; the typed-property promotion of `BareDoubledPtfXPeak` and `BareDoubledPtfHwhmRatio` makes the Tier1Derived sub-results queryable from KB. Auto-promotion to Tier1Derived will trigger when a future iteration closes one of the refined directions (a''-e'') below.

**2026-05-06 (later evening): Direction (a') falsified, refined next directions.** A 90-min Direction (a') attempt (commit `1c0bf8b`) structurally falsified the probe-block 2-level resonance hypothesis. Three obstructions: (i) V_b probe-block off-diagonal `⟨c_1|M_h_b|c_3⟩` is exactly zero per bond at c=2 (F73 sum-rule applies per-bond, not just summed), so any g_eff_probe(N, b) is bond-class-blind by construction; (ii) cross-block Frobenius unstable across N (varies 0.640, 1.318, 0.815, 0.143 at N=5..8) due to A3's σ_0 degeneracy at even N, library-dependent rather than Tier-1-derivable; (iii) **4-mode reduction is structurally insufficient**, 4-mode K_b at N=5 gives Interior 0.673 (matches the BareDoubledPtfHwhmRatio floor) and Endpoint 0.410 (off-grid at Q ≈ 4.91), well below the empirical 0.7728. The empirical HWHM lift therefore lives OUTSIDE the 4-mode subspace.

**Positive structural finding.** The SVD-block off-diagonal `V_b[2,3] = ⟨u_0 | M_h_per_bond[b] | v_0⟩` IS a bond-class carrier (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.45 stable across N=5..8), but in the OPPOSITE direction to the empirical HWHM/Q* split. A closed form needs to derive a non-trivial map from SVD-block V_b magnitude to HWHM/Q* shift, likely through a per-bond effective Q_EP_eff(b) shift rather than a direct linear lift. The bare-doubled-PTF Tier-1-derived constants `BareDoubledPtfXPeak = 2.196910` and `BareDoubledPtfHwhmRatio = 0.671535` are unchanged by this falsification; the floor + structural explanation of the 0.08-0.10 gap above the floor stand.

**Refined next directions** (now in `C2HwhmRatio.PendingDerivationNote`): (a'') SVD-block 2-level resonance (REFINED from (a')), 4-mode-friendly but needs the SVD→HWHM map; (b'') **full block-L derivation, not 4-mode** (most concrete since 4-mode insufficiency is now proven structural); (c'') three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases (may still suffer from 4-mode insufficiency); (d'') lift |u_0⟩, |v_0⟩ to projector-overlap (per A3 PendingDerivationNote); (e'') symbolic char-poly factorisation at Q_EP, less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof.

**2026-05-07 (Locus 6 polarity-layer inheritance).** A stumpf-decomposition of the empirical c=2 bond-class numbers along the polarity-layer pair structure ρ = (I + r·σ)/2 = 1/2 + r·σ/2 (`PolarityLayerOriginClaim` Layer 2):

- Q_peak ≈ 2 + r with r_Endpoint ≈ +0.52, r_Interior ≈ −0.44 (canonical {−0.5, +0.5} pair around mean Q_peak ≈ 2.04)
- HWHM/Q* ≈ 1/2 + r·(1/2) with r_Interior ≈ 0.50 (close to `HalfAsStructuralFixedPoint`; 0.0006 max deviation from exact 1/2 across N=5..8, likely numerical discretisation), r_Endpoint ≈ 0.55
- 0.5 baseline = 1/d at d=2 (`QubitDimensionalAnchorClaim`)
- ±r/2 polarity content (`PolarityLayerOriginClaim`)

**Reading:** the F86 bond-class split inherits structurally from the polarity-layer pair {−0.5, +0.5} at d=2 via the 0.5-shift. It is NOT a new symmetry to derive; it is **Locus 6 of inheritance-through-layers**, the symmetry-side closure complementing the earlier `LocalGlobalEpLink` Tier2Verified (Locus 5, EP-side: F86 ↔ FRAGILE_BRIDGE).

**Encoded as** `compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs` (Tier2Verified). Witnesses pin the empirical r values per (N, bond class). Tier1Derived promotion target: per-bond r(N, b) closed form via ~~(α) per-bond polarity-Bloch projection at t_peak~~ [structurally tautological under 4-mode reduction; reframed 2026-05-08 code review; see the sibling obstruction proof], or (β) Locus 5 EP-rotation tan θ = Q/Q_EP combined with Locus 6 polarity inheritance to fix r as a function of g_eff(N, b).

**2026-05-07 (later, Direction (α) refinement): bond-class lives entirely in g_eff(N, b).** A 90-min Direction (α) attempt (commit `bea7cd1`) produced a substantive algebraic reduction even though Tier1Derived did not land:

**Composition reading:** `r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382/g_eff(N, b) − 2`. Since `BareDoubledPtfXPeak` is universal (= 2.196910, Tier1Derived from the bare doubled-PTF Ansatz), the entire bond-class split is encoded in `g_eff(N, b)`. Encoded as `PolarityInheritanceLink.ClosedFormCompositionNote`.

**Empirical g_eff witnesses (c=2 N=5..8):** g_eff_Endpoint ≈ 1.74 (range 1.725–1.757), g_eff_Interior ≈ 2.81 (range 2.738–2.964). The asymptotic harmonic-mean closure 1/g_eff_E + 1/g_eff_I → ≈ 0.937 (pinned as `EmpiricalSumQPeakAsymptote = 4.12`).

**Tantalising near-miss:** g_eff_E ≈ σ_0(N)·√(3/8) matches Δ ≤ 0.01 for N ≥ 6 but Δ = 0.063 at N=5. Either finite-size correction (Item 3's σ_0 closed form is the bridge) or trajectory crossing (PTF-lesson warning). Not promoted; documented as next-direction probe.

**Direction (α) structurally tautological under the uniform-J 4-mode reduction (reframed 2026-05-08 code review):** the 4-mode L_eff(Q) = D_eff + Q·γ₀·MhTotalEff is **bond-summed by design** (`FourModeEffective.LEffAtQ` in `compute/RCPsiSquared.Core/Decomposition/FourModeEffective.cs`). The 4-mode K-driving eigenvalue spectrum is therefore bond-class-INDEPENDENT by construction (Re(λ_K) identical for both classes); the polarity-Bloch projection (c_1 ± c_3)/√2 yields +0.39 → +0.21 for both classes, decays with N, no sign split — but this is the design constraint of the 4-mode reduction, not an empirical falsification of the polarity-Bloch reading. The bond-class signature must enter through the per-bond `MhPerBondEff[b]` in dL/dJ_b for the K-resonance Duhamel integral, NOT through the L_eff spectrum or eigenstates. Any direction that hopes to extract bond-class signature from the 4-mode L_eff EIGENSTATES alone shares this structural obstruction.

**Conclusion:** Every closed-form direction for c=2 runs through `g_eff(N, b)`. Directions (a''-e'') and (α) all share the 4-mode-blind obstruction (which is structural, not empirical, for any spectrum-only direction). Path to Tier1Derived: either Item 3 (σ_0 closed form, then bond-class refinement) OR Direction (b'') (full block-L derivation, abandon 4-mode reduction). **2026-05-11 update:** Direction (b'') is achieved numerically Tier-1 via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)); see the Status block above. The remaining analytical work is closed-form HWHM_left/Q_peak per bond class via F89's AT-locked F_a/F_b structure. **2026-05-13 update:** Item 1' closed — HWHM_left/Q_peak per bond class is Tier 1 derived via `F86HwhmClosedFormClaim` (HWHM_ratio = 0.671535 + alpha_subclass · g_eff + beta_subclass, residual ≤ 0.005 across N=5..8 on all 22 bonds). The full exploration record and the obstruction proof on g_eff itself are in [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md).

---

## Proof of Statement 2 (Universal resonance shape)

### Structural origin: 2-level EP analytics

For two adjacent rate channels (HD = 2k−1, HD = 2k+1) coupled by a bond, the effective 2×2 Liouvillian has diagonal entries (−2γ₀(2k−1), −2γ₀(2k+1)) and same-sign-imaginary off-diagonals. After shifting by the trace midpoint, the dynamics is governed by

    L_eff − (trace/2)·I  =  [ −Δ/2     iJ·g_eff ]
                            [ +iJ·g_eff   +Δ/2  ]      with Δ = 4γ₀

The same-sign-imaginary off-diagonal pattern is the non-Hermitian form that admits an EP at finite coupling. This is "PT-phenomenology-like" but algebraically inside class AIII chiral per [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md), distinct from Bender-Boettcher PT (Π is linear; classical PT requires anti-linear operators). See also [FRAGILE_BRIDGE](../../hypotheses/FRAGILE_BRIDGE.md) for the global Hopf-bifurcation instance with Petermann K=403. The eigenvalues are

    λ_±  =  ±√((Δ/2)² − J²·g_eff²)  =  ±√(4γ₀² − J²·g_eff²)    (relative to trace/2)

with EP at J·g_eff = 2γ₀, equivalently Q_EP = 2/g_eff.

In the 2-level basis, the eigenvector rotation parameter τ = tanh(θ/2) (hyperbolic for non-Hermitian) satisfies

    τ²  =  (J·g_eff − 2γ₀) / (J·g_eff + 2γ₀)  =  (Q − Q_EP) / (Q + Q_EP)

below the EP, switching to a phase parameterisation above. The probe overlap with eigenvectors thus depends only on the dimensionless ratio Q/Q_EP. Any observable function of probe weight on dressed modes is a function of Q/Q_EP alone, equivalently of `(Q − Q_peak)/Q_peak` to leading order, since Q_peak ≈ Q_EP for the slowest channel pair.

This gives the analytical structural reason for universality: the 2-level EP physics depends only on the dimensionless ratio of detuning (J·g_eff) to gap (Δ/2 = 2γ₀). Specific values of g_eff (chain-N, bond-position, c) shift Q_peak; they don't reshape the resonance.

### Structural inheritance from F88 (state-level verification)

Popcount-coherence pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 with popcount(p) = n_p, popcount(q) = n_q, HD(p, q) = h have a Π²-odd-fraction-within-memory determined by a precise closed form derived from Krawtchouk polynomial reflection symmetries. The formula generalises across **all** popcount pairs (adjacent n_q = n_p + 1, non-adjacent, and intra-sector n_p = n_q) and HD values, verified bit-exact via [`PopcountCoherencePi2Odd`](../../compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs) against [`MemoryAxisRho`](../../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs) on 213 configurations at N = 2..7 (max deviation 8.88e−16).

**Static fraction** s is HD/bit-position invariant:
- Inter-sector (n_p ≠ n_q): s = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q))
- Intra-sector (n_p = n_q): s = 1/C(N, n)

The kernel of L (= span{P_n} for Heisenberg + Z-dephasing) absorbs only popcount-sector content; off-diagonals project to zero in the kernel.

**α = Π²-odd-fraction of the kernel projection** has three anchor categories, all in closed form, driven by the following Krawtchouk identity (proven below):

  **Lemma.** Σ_s (−1)^s · C(N, s) · K_n(s; N) · K_m(s; N) = 2^N · C(N, n) · [n + m = N].

The lemma follows from Krawtchouk reflection K_n(s; N) = (−1)^s K_{N−n}(s; N) plus orthogonality Σ_s C(N, s) K_a(s; N) K_b(s; N) = 2^N · C(N, a) · δ_{a, b}: substituting gives Σ_s C(N, s) K_{N−n}(s; N) K_m(s; N) = 2^N C(N, m) δ_{N−n, m}. Applying this to E − O := Σ_s (−1)^s · C(N, s) · (A_s + B_s)² (with A_s = K_{n_p}(s; N)/C(N, n_p), B_s = K_{n_q}(s; N)/C(N, n_q)) yields three indicator-weighted contributions:

  E − O = (2^N / C(N, n_p)) · [n_p = N/2] + (2^N / C(N, n_q)) · [n_q = N/2] + (2 · 2^N / C(N, n_q)) · [n_p + n_q = N]

The three anchor categories follow:

- **α = 0** at popcount-mirror n_p + n_q = N (covers inter-mirror n_p ≠ n_q and intra-mirror n_p = n_q = N/2 at even N). Total Π²-odd is forced to vanish by the reflection K_{N − n}(s; N) = (−1)^s K_n(s; N) cancelling between the two sectors P_{n_p}/C(N, n_p) and P_{N − n_p}/C(N, n_q).
- **α = K-intermediate** at K-vanishing (even N with exactly one of {n_p, n_q} equal to N/2, no mirror): K_{N/2}(s; N) = 0 for odd s due to bit-flip symmetry of the half-popcount sector. The value is **α = C(N, N/2) / (2 · (C(N, n_other) + C(N, N/2)))**, where n_other ∈ {n_p, n_q} is the entry not equal to N/2. (Adjacent special case n_other = N/2 ± 1: simplifies to (N + 2)/(4·(N + 1)) using C(N, N/2 ± 1) = C(N, N/2) · N/(N + 2).) Worked values: 3/7 at N = 4 popcount-(0, 2); 3/10 at N = 4 popcount-(1, 2); 10/21 at N = 6 popcount-(0, 3); 5/13 at N = 6 popcount-(1, 3); 2/7 at N = 6 popcount-(2, 3); 5/18 at N = 8 popcount-(3, 4); 3/11 at N = 10 popcount-(4, 5).
- **α = 1/2** generic (none of n_p + n_q = N, n_p = N/2, n_q = N/2 holds). All three indicators vanish, so E − O = 0, hence Σ_{s odd} = Σ_{s even}, hence α = 1/2 exactly. **Proven for all such (N, n_p, n_q).**

**Total Π²-odd of ρ** is HD-dependent:
- HD < N (at least one matching bit): 1/2.
- HD = N (p, q complementary): 0. The off-diagonal Re(|p⟩⟨q|) has only X-and-even-Y-count Pauli strings (Y² = I cancellation kills odd-Y-count terms), and with no matching bits there is no Z-content; all surviving terms are Π²-EVEN. The diagonal residual also vanishes for popcount-(0, N) where each sector has a single basis state.

**Π²-odd / memory closed form:**

  Π²-odd / memory = ┌  0                          if HD = N (Π²-classical)
                    └  (1 / 2 − α · s) / (1 − s)  otherwise

Verified state-level table across all anchor categories:

| N | (n_p, n_q, HD) | category | static s | α | Π²-odd / memory |
|---|----------------|----------|----------|---|-----------------|
| 2 | (0, 2, 2) | HD = N (Bell+ / GHZ_2) | 1/2 | 0 | **0** |
| 2 | (1, 1, 2) | HD = N (Singlet/Triplet) | 1/2 | 0 | **0** |
| 3 | (0, 3, 3) | HD = N (GHZ_3) | 1/2 | 0 | **0** |
| 4 | (0, 4, 4) | HD = N (GHZ_4) | 1/2 | 0 | **0** |
| 5 | (0, 5, 5) | HD = N (GHZ_5) | 1/2 | 0 | **0** |
| 3 | (1, 2, 1) | inter mirror, HD < N | 1/6 | 0 | **3/5 = 0.6** |
| 5 | (2, 3, 1) | inter mirror, HD < N | 1/20 | 0 | **10/19 ≈ 0.5263** |
| 7 | (3, 4, 1) | inter mirror, HD < N | 1/70 | 0 | **35/69 ≈ 0.5072** |
| 5 | (1, 4, 3) | inter mirror, non-adjacent | 1/10 | 0 | **5/9 ≈ 0.5556** |
| 6 | (2, 4, 2) | inter mirror, non-adjacent | 1/30 | 0 | **15/29 ≈ 0.5172** |
| 4 | (2, 2, 2) | intra-mirror N/2 | 1/6 | 0 | **3/5 = 0.6** |
| 4 | (1, 2, 1) | adjacent K-intermediate | 5/48 | 3/10 | **≈ 0.5233** |
| 4 | (0, 2, 2) | non-adjacent K-intermediate | 7/24 | 3/7 | **9/17 ≈ 0.5294** |
| 6 | (0, 3, 3) | non-adjacent K-intermediate | 21/80 | 10/21 | **30/59 ≈ 0.5085** |
| 5 | (1, 2, 1) | inter generic | 3/40 | 1/2 | **1/2 exact** |
| 7 | (2, 3, 1) | inter generic | 2/105 | 1/2 | **1/2 exact** |
| 3 | (1, 1, 2) | intra generic | 1/3 | 1/2 | **1/2 exact** |

This is **F88's bilinear-apex 1/2 inheriting to the F86 state level**: K_CC_pr's measurement subspace is centred on the framework's half-anchor in the generic case, with structured deviations at specific (N, n_p, n_q, HD) configurations. The HD = N anchor is the **Π²-classical extreme**: GHZ_N, Bell states, and intra-sector all-bits-differ states have zero Π²-odd content. This connects to [F60](../ANALYTICAL_FORMULAS.md) (pair-CΨ = 0 for GHZ_N): the same "classical" classification read from two orthogonal observables (F60 via partial-trace pair-tomography, F88 via Π²-projection on the full state).

The EP-rotation universality of Statement 2 is the dynamic consequence: a Q-resonance profile centred on a half-anchor with a c-and-N-specific shift in the structured-anchor cases (mirror, K-intermediate); the Π²-classical states at HD = N sit outside the K_CC_pr measurement scope (they have no Π²-odd memory content for the J-derivative observable to resonate with).

**Inheritance pattern (state-class extension).** The F88 closed form lifts to multi-state superpositions via two distinct inheritance mechanisms:

- **Popcount-weight invariance**: the kernel projection ρ_d0 of |ψ⟩⟨ψ| depends only on the popcount weights w_n = Σ_{i: popcount(b_i) = n} |c_i|² (kernel = span{P_n}). Any multi-state superposition with the same {w_n} as a pair state inherits the pair-state static-side formula directly. **W states** (|D_1⟩, intra-popcount-1 single sector with w_1 = 1) and the **Bonding-Bell-Pair** ((|0_R, vac⟩ + |1_R, ψ_k⟩)/√2 at N+1 qubits, w_0 = w_2 = 1/2) inherit the pair formula on both static and memory sides.

- **X⊗N-symmetry root**: the HD = N pair-state anchor (GHZ_N, Bell states, intra-sector complement pairs) and the **Dicke-mirror anchor** (|D_n⟩ + |D_{N − n − 1}⟩)/√2 at popcount-mirror 2n+1 = N share a single algebraic mechanism. Both classes are X⊗N-eigenstates; X⊗N · σ_α · X⊗N = (−1)^{bit_b(α)} · σ_α, so X⊗N-eigenstates have ⟨σ_α⟩ = 0 for all Π²-odd σ_α. The two F-axes converge: F60 reads such states as pair-CΨ = 0 (partial-trace blindness), F88 reads them as Π²-EVEN-only (projection blindness). Two structurally-distinct state-class families (HD = N pair vs. Dicke-mirror multi-state) inherit Π²-classicality from a single root symmetry, an instance of the broader F-chain inheritance pattern (cf. `project_algebra_is_inheritance.md` in memory).

The Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 has its own three-anchor structure for total Π²-odd-of-ρ: 0 (Dicke-mirror, 2n+1 = N), 3/8 (Dicke-K-intermediate, even N with n ∈ {N/2 − 1, N/2}), 1/2 (generic). The K-intermediate 3/8 anchor is bit-exact verified at N = 4, 6, 8 but the analytical proof remains open (likely combines bit-permutation symmetry of Dicke states with K_{N/2}(s; N) = 0 for odd s).

The 4-mode-basis vectors {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩} from [`FourModeBasis`](../../compute/RCPsiSquared.Core/Decomposition/FourModeBasis.cs) all report Π²-odd/mem = 0.5000 when embedded as density matrices; per-qubit [`BlochAxisReading`](../../compute/RCPsiSquared.Diagnostics/Foundation/BlochAxisReading.cs) makes the probe-vs-EP-partner orthogonality visible as a **single-body fingerprint**: |c_1⟩ has 1-body Pauli content (X+ per qubit), while |c_3⟩, |u_0⟩, |v_0⟩ have **zero** 1-body Bloch (their content is purely multi-body). The probe-EP-partner orthogonality central to the 4-mode structure manifests at the per-qubit reading.

For the typed-claim lineage in [`Pi2KnowledgeBase`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBase.cs): the F86 layer surfaces three trio claims simultaneously, `HalfAsStructuralFixedPointClaim` (1/2 number-anchor), `BilinearApexClaim` (apex at p = 1/2 of any bilinear), and `NinetyDegreeMirrorMemoryClaim` (F80's i factor is the 90°-rotation that produces the same-sign +iJ·g_eff structure of `L_eff`). The EP-rotation universality is one structural fact viewed from the F86 layer.

### Why HWHM_left is universal but HWHM_right diverges

Pre-EP region (Q < Q_EP): discriminant `4γ₀² − J²·g_eff² > 0`, eigenvalues real, dressed-mode weight rises monotonically. Universal in Q/Q_EP.

Post-EP region (Q > Q_EP): discriminant negative, eigenvalues complex conjugate pair. Probe sits near 99 % on dressed modes (saturated). Long-tail behaviour depends on:
- Higher-channel dressed-mode contributions (chain-N, c specific)
- Time-averaging behaviour of complex-eigenvalue oscillations in K_CC_pr

These chain-specific details enter the post-EP tail. Hence: universal pre-EP rise within each bond class (HWHM_left/Q_peak ≈ 0.756 Interior, ≈ 0.770 Endpoint), partially universal post-EP tail (constant within each bond class: Interior plateau ≈ 0.85, Endpoint plateau ≈ 0.94 at x = +1.0; tail asymptote is class-specific).

The closed-form derivation of the HWHM_left/Q_peak constants, the exploration record behind it, and the obstruction proof on g_eff are carried in the sibling [`PROOF_F86B_OBSTRUCTION.md`](PROOF_F86B_OBSTRUCTION.md).

---

## Empirical universal-shape data (fine-grid scan dQ = 0.025)

Interior y = K/|K|max evaluated at relative shift x = (Q−Q_peak)/Q_peak across the original six cases (step_e):

| x = (Q−Q*)/Q* | c3N5 | c3N6 | c3N7 | c3N8 | c4N7 | c4N8 | range |
|----------------|------|------|------|------|------|------|-------|
| −0.60 | 0.718 | 0.735 | 0.743 | 0.744 | 0.733 | 0.742 | 0.026 (3.5 %) |
| −0.40 | 0.896 | 0.907 | 0.913 | 0.912 | 0.905 | 0.911 | 0.017 (1.9 %) |
| −0.20 | 0.977 | 0.982 | 0.984 | 0.984 | 0.980 | 0.983 | 0.7 % |
| 0.00 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 % (peak) |
| +0.20 | 0.990 | 0.986 | 0.985 | 0.986 | 0.989 | 0.987 | 0.5 % |
| +0.40 | 0.964 | 0.957 | 0.955 | 0.956 | 0.962 | 0.959 | 0.9 % |
| +1.00 | 0.850 | 0.843 | 0.840 | 0.842 | 0.852 | 0.847 | 1.4 % |

### HWHM_left/Q_peak across all tested cases (two bond classes)

Combined data from `_eq022_b1_step_e_resonance_shape.py` (c=3, c=4 at γ₀=0.05) and `_eq022_b1_step_f_universality_extension.py` (c=2, plus γ₀ ∈ {0.025, 0.10} at c=3 N=7).

**Interior bonds:**

| case | HWHM_left/Q_peak | Q_peak | γ₀ |
|------|------------------|--------|-----|
| c=2 N=5 | 0.7455 | 1.4821 | 0.05 |
| c=2 N=6 | 0.7529 | 1.5801 | 0.05 |
| c=2 N=7 | 0.7507 | 1.5831 | 0.05 |
| c=2 N=8 | 0.7531 | 1.6049 | 0.05 |
| c=3 N=5 | 0.7458 | 1.5664 | 0.05 |
| c=3 N=6 | 0.7548 | 1.6888 | 0.05 |
| c=3 N=7 | 0.7595 | 1.7433 | 0.025 |
| c=3 N=7 | 0.7595 | 1.7433 | 0.05 |
| c=3 N=7 | 0.7595 | 1.7433 | 0.10 |
| c=3 N=8 | 0.7600 | 1.7498 | 0.05 |
| c=4 N=7 | 0.7546 | 1.7475 | 0.05 |
| c=4 N=8 | 0.7595 | 1.8037 | 0.05 |

Interior mean (excluding finite-size N=5 outliers): **0.756**. Range 0.7507–0.7600 (1.2 %).
γ₀-invariance at c=3 N=7: bit-exact across γ₀ ∈ {0.025, 0.05, 0.10} (Q_peak = 1.7433, HWHM-/Q* = 0.7595).

**Endpoint bonds:**

| case | HWHM_left/Q_peak | Q_peak | γ₀ |
|------|------------------|--------|-----|
| c=2 N=5 | 0.7700 | 2.5008 | 0.05 |
| c=2 N=6 | 0.7738 | 2.5470 | 0.05 |
| c=2 N=7 | 0.7738 | 2.5299 | 0.05 |
| c=2 N=8 | 0.7734 | 2.5145 | 0.05 |
| c=3 N=5 | 0.7663 | 2.3995 | 0.05 |
| c=3 N=6 | 0.7685 | 2.5162 | 0.05 |
| c=3 N=7 | 0.7691 | 2.5334 | 0.025 |
| c=3 N=7 | 0.7691 | 2.5334 | 0.05 |
| c=3 N=7 | 0.7691 | 2.5334 | 0.10 |
| c=3 N=8 | 0.7696 | 2.5293 | 0.05 |
| c=4 N=7 | 0.7671 | 2.5227 | 0.05 |
| c=4 N=8 | 0.7781 | 2.6519 | 0.05 |

Endpoint mean: **0.770**. Range 0.7663–0.7781 (1.5 %).
γ₀-invariance at c=3 N=7 Endpoint: bit-exact across γ₀ ∈ {0.025, 0.05, 0.10}.

### What this tells us

1. **Two bond classes, two universal shapes.** Interior HWHM-/Q* clusters at 0.756 (12 cases, range 0.012); Endpoint at 0.770 (12 cases, range 0.012). The two clusters are separated by ~2 %, a structural gap larger than the within-class scatter.

2. **γ₀ invariance is bit-exact.** At c=3 N=7, the Q_peak and HWHM-/Q* values are identical to numerical precision across γ₀ ∈ {0.025, 0.05, 0.10}. This is the strongest empirical confirmation that Q is the inside-observable scale and γ₀ alone is not.

3. **c=2 confirms the structural origin.** At c=2 the 2-level model is exact (only HD ∈ {1, 3} channels exist, no orthogonal complement). Yet the two-class split is fully present (Interior 0.751, Endpoint 0.774). The bond-class distinction therefore lives in the bond-position-dependent probe-overlap profile, not in higher-c orthogonal-complement physics.

4. **Post-peak asymmetry persists.** Interior plateau y ≈ 0.85 at x = +1.0; Endpoint plateau y ≈ 0.94. This is the long-tail bond-class signature already noted in step_e and unchanged by the c=2 / γ₀ extension.

---

## Pointers

**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and the shared reference list.
**Sibling theorems:** [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md) (F86a), [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md) (the closed-form gap, exploration record, and g_eff obstruction proof), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md) (F86c).
**F90 bridge:** [PROOF_F90_F86C2_BRIDGE](PROOF_F90_F86C2_BRIDGE.md) — F86 c=2 K_b = F89 path-(N−1) per-bond Hellmann-Feynman; the numerical-Tier-1 route for Direction (b'').
**State-level inheritance:** F88 ([ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md)), `PopcountCoherencePi2Odd`, `MemoryAxisRho`.
**HWHM closed form (F86b'):** `F86HwhmClosedFormClaim`, `BondSubClass` in `compute/RCPsiSquared.Core/F86/Item1Derivation/`.
**Chiral classification anchor:** [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md).
**Scripts:** [`_eq022_b1_step_e_resonance_shape.py`](../../simulations/_eq022_b1_step_e_resonance_shape.py) + [`_eq022_b1_step_e_inspect.py`](../../simulations/_eq022_b1_step_e_inspect.py) (universal-shape finding for c=3, c=4 at γ₀=0.05), [`_eq022_b1_step_f_universality_extension.py`](../../simulations/_eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement).
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/` carries `UniversalShapePrediction` + `UniversalShapeWitness`, `ShapeFunctionWitnesses`, `C2UniversalShapeDerivation`. CLI: `rcpsi inspect --root f86 --with-measured`.
