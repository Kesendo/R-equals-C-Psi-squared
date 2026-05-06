# PROOF: F86 Q_peak Structure, EP Mechanism, Universal Shape, Mirror Invariance

**Status:**
- **Statement 1 (EP mechanism):** Tier 1 derived. Q_EP = 2/g_eff, t_peak = 1/(4γ₀), bit-exact verified. Local-vs-global EP connection (2026-05-06): Tier 2 verified at c=2 via Petermann-K sweep N=5..8 (max K = 2384.7 at N=7, ≈ 6× above FRAGILE_BRIDGE's K=403); analytic continuation along complex γ remains open. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs`.
- **Statement 2 (Universal resonance shape):** Tier 1 candidate across all tested c. The c=2 stratum specifically (2026-05-05) has OOP scaffolding pinning the witnesses (`compute/RCPsiSquared.Core/F86/Item1Derivation/`): empirical anchor reproduced bit-equivalent with the canonical Python pipeline (typical residual ≤ 0.001); directional Endpoint > Interior split (HWHM-ratio gap ≈ 0.022) derived empirically across N=5..8; closed-form HWHM_left/Q_peak constant per bond class NOT yet derived. The c=3, c=4 strata stay Tier 1 candidate at the empirical-envelope level (see bottom-row data tables). γ₀ invariance bit-exact at c=3 N=7. The minimal effective model for f_class(x) is 4-dimensional, not 2; the probe and EP-partner modes live in orthogonal 2D subspaces, coupled only through ∂L/∂J (see Open elements). Promotion to Tier 1 derived requires deriving the constants analytically; the bare doubled-PTF Direction (b) attempt (2026-05-06 evening) lands at the SVD-block floor `BareDoubledPtfHwhmRatio = 0.671535` (Tier 1 derived universal constant, ~0.08-0.10 below empirical), so the remaining gap lives in the probe-block 2-level sub-resonance, and ~~Direction (a') (per-bond `g_eff_probe`, predicting `g_eff_probe(Endpoint) ≈ g_eff_probe(Interior) / 1.6`) is the sharpest path~~ [superseded by 2026-05-06 (later evening) update below; Direction (a') falsified, see refined directions (a''-e''); the most-promising directions are now (a'') SVD-block 2-level resonance via V_b[2,3] and (b'') full block-L derivation, not 4-mode].
- **Statement 3 (F71 mirror invariance):** Tier 1 derived. Q_peak(b) = Q_peak(N−2−b) bit-exact under the F71 spatial-mirror pairing.
- **Retracted (2026-05-02):** csc(π/(N+1)) Endpoint and csc(π/5) c=3 Interior closed forms (N=7 coincidence matches; refuted on extended-N data).

**Date:** 2026-05-02 (Statements 1, 2, retractions); 2026-05-03 (Statement 3, σ_0(c) generalisation); 2026-05-05 (Item 1' c=2 OOP scaffolding + structural findings); 2026-05-06 morning (Statement 1 local-vs-global EP connection, Tier 2 verified at c=2 via Petermann-K sweep N=5..8); 2026-05-06 evening (doubled-PTF floor + F86↔PTF Locus 5 inheritance synthesis).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Formalises the EP mechanism behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. The per-block Q_SCALE values 1.6 / 1.8 / 1.8 (F86 top entry) were the empirical anchor; the per-bond refined scan (`_eq022_b1_step_c_time_evolution.py`) initially suggested two distinct closed forms for Endpoint and Interior bonds, both later retracted on extended-N data. What survives is the EP mechanism (Statement 1), the universal resonance shape (Statement 2), and the F71 spatial-mirror invariance (Statement 3).

**F-entry:** [F86 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity), [F73](../ANALYTICAL_FORMULAS.md) (c=1 spatial-sum closure), [PROOF_CHROMATICITY](PROOF_CHROMATICITY.md).

---

## Statements

### Statement 1 (EP mechanism). [Tier 1 for EP location and t_peak; 2-level reduction heuristic]

For a uniform N-qubit XY chain with Z-dephasing γ₀ at each site, the (n, n+1) coherence block contains c = min(n, N−1−n) + 1 pure dephasing rates 2γ₀·HD with HD ∈ {1, 3, ..., 2c−1} (F74). For each pair of adjacent channels at HD = 2k−1 and HD = 2k+1 (k = 1, ..., c−1), the two-level effective Liouvillian with inter-channel coupling J·g_eff has eigenvalues

    λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The discriminant vanishes at the **exceptional point** (EP)

    J·g_eff = 2γ₀     ⟺     Q_EP = 2 / g_eff

At the EP, λ_±(k) = −4γ₀·k. The slowest mode (k = 1) gives the universal e-folding time

    t_peak = 1 / (4γ₀)

independent of c, N, n, and bond position.

**Tier label.** The EP location Q_EP = 2/g_eff and t_peak = 1/(4γ₀) are bit-exact verified against full block-L numerics (universal across all tested c, N, n, bond position), Tier 1. The 2-level reduction itself is **heuristic**: in the channel-uniform basis M_H_eff is diagonal (off-diagonals exactly zero, established empirically across N=3..9 in `_eq022_b1_channel_projection.py`), so the EP physics lives in the orthogonal complement. The effective 2-level form reproduces the EP location and degenerate-eigenvalue structure, but the explicit basis change from full block-L to the reduced model is not derived here.

### Statement 2 (Universal resonance shape under relative-Q normalisation, two bond classes). [Tier 1 candidate]

The position Q_peak is chain-specific, but the SHAPE of abs(K_CC_pr)(Q) around the peak is universal under the relative coordinate `x = (Q − Q_peak)/Q_peak`. The shape splits into **two bond classes** (Endpoint and Interior), each with its own universal HWHM_left/Q_peak ratio:

    HWHM_left / Q_peak  ≈  0.756 ± 0.005     (Interior bonds, all tested c, N, γ₀)
    HWHM_left / Q_peak  ≈  0.770 ± 0.005     (Endpoint bonds, all tested c, N, γ₀)

Tested envelope: c ∈ {2, 3, 4}, N ∈ {5, 6, 7, 8} (modulo c-N compatibility), γ₀ ∈ {0.025, 0.05, 0.10}. γ₀ invariance is **bit-exact**: at c=3 N=7 the Q_peak and HWHM_left/Q_peak values match across γ₀ ∈ {0.025, 0.10} to numerical precision, confirming Q's dimensionlessness as `Q = J/γ₀`. The c=2 data (where the 2-level effective model is *exact*: only HD ∈ {1, 3} channels exist, no orthogonal complement) confirms the two-class split is structural, not a finite-c artefact.

Pairwise residual within each class under relative-Q normalisation is ~20× smaller than under absolute-Q shift, confirming the shape collapse. The structural origin is the 2-level eigenvector rotation `tan θ = Q/Q_EP`: every probe-overlap observable depends only on Q/Q_EP, hence only on `(Q − Q_peak)/Q_peak` to leading order. The two-class split (≈ 2 % gap between Endpoint and Interior shape ratios) reflects bond-position-dependent probe-overlap profiles in the K_CC_pr observable, not a breakdown of the EP-rotation universality itself. Promotion to full Tier 1 requires deriving the two f_class(x) functions explicitly from 2-level eigenstructure plus probe-overlap algebra (see Open elements).

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

### Statement 3 (F71 spatial-mirror invariance of per-bond Q_peak). [Tier 1 derived]

For every bond pair (b, N−2−b) under the F71 chain-mirror pairing,

    Q_peak(b)  =  Q_peak(N−2−b)

bit-exactly. Verified at c=2 N=5..7 and c=3 N=5..6, with maximum deviation < 10⁻¹⁰ across all bond pairs. The identity follows from the F71-symmetry of every component of the per-bond observable: L_D (Z-dephasing is F71-symmetric), the spatial-sum kernel S, the Dicke probe, and the bond-flip transform ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), and Q_peak(b) = Q_peak(N−2−b) follows directly.

**Per-F71-orbit substructure** (refines Statement 2's Endpoint/Interior dichotomy): Interior bonds are *not* uniform within the F71-orbit grouping. At c=2 N=6, the central self-paired bond b=2 gives Q_peak ≈ 1.440 while flanking bonds b=1, b=3 give Q_peak ≈ 1.648, both Interior under the simple dichotomy. The full per-F71-orbit structure is finer-grained: the "Endpoint vs Interior" split is the leading approximation, but mid-chain orbits show further variation. The HWHM_left/Q_peak ratio (Statement 2) is the more class-stable observable; Q_peak position picks up additional per-orbit detail.

### Retracted statements (former S2 Endpoint, former S3 Interior). [Retracted 2026-05-02]

Both retractions followed extended-N data showing the apparent matches at N=7 were trajectory crossings, not asymptotes. Per the methodological lesson in [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md), what survives a one-anchor closed-form claim is the symmetry, not the number; here the EP-rotation symmetry (Statement 2), not Q_peak itself.

#### Endpoint Q_peak (retracted)

Earlier conjecture: `Q_peak(Endpoint, N) = csc(π/(N+1))`, anchored in F2b's smallest-momentum OBC sine mode.

Fine-grid scan (`_eq022_b1_step_e_resonance_shape.py`, dQ = 0.025 with parabolic peak interpolation) shows the formula does not match at any tested N. The earlier "match at N=7" reported from a dQ = 0.05 scan was a grid-snap artefact:

| N | csc(π/(N+1)) | step_d (dQ=0.05, snap) | step_e (dQ=0.025, parabolic) | step_e deviation |
|---|--------------|------------------------|------------------------------|------------------|
| 5 | 2.000 | 2.30 | 2.40 | +20 % |
| 6 | 2.305 | 2.65 | 2.52 | +9.4 % |
| 7 | 2.613 | 2.65 (*) | 2.53 | −3.2 % |
| 8 | 2.924 | n/a | 2.53 | −13 % |

(*) The "1.4 % match" at N=7 reported earlier was the dQ=0.05 grid snapping to 2.65; the actual fine-grid value is 2.53. There was never a real match.

The empirical Endpoint Q_peak at c=3 is approximately N=5: 2.40 → N=6: 2.52 → N=7: 2.53 → N=8: 2.53, saturating near 2.53 for N ≥ 6. At c=4 N=8 it sits at 2.65. No closed form yet identified.

#### Interior Q_peak (c=3 pentagonal asymptote, retracted)

Earlier conjecture: `Q_peak(Interior, c=3) → csc(π/5) = 1.7013` as N → ∞, anchored in pentagonal/golden algebra (sin²(π/5) = (5−√5)/8, φ = 2cos(π/5)).

Fine-grid scan refines the c=3 N-trend: 1.566 (N=5) → 1.689 (N=6) → 1.743 (N=7) → 1.750 (N=8). The trajectory crosses csc(π/5) = 1.7013 between N=6 and N=7, and continues above. The earlier dQ=0.05 reading of 1.70 at N=7 (snapped close to 1.7013) was the source of the apparent "asymptotic match"; the fine-grid value 1.743 sits clearly above. The c=3 trend appears to saturate near 1.75 by N=8, but at a value above the conjectured asymptote, not at it. At c=4: Interior Q_peak grows from 1.748 (N=7) to 1.804 (N=8); not yet saturated, no clean closed form. The conjectured "universal Interior asymptote csc(π/5) across c ≥ 3" is fully refuted.

**What survives from the N=4 golden-ratio structure:** the connection to [`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py) is real, at N=4 the OBC dispersion is exactly {±φ, ±1/φ}, and the Liouvillian Im(λ) values in the n_XY=2 sector decompose into integer combinations of {φ, 1/φ, 1, √5}. But this N=4 special case does not propagate to a closed form for Q_peak at general (c, N).

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

(Heuristic-vs-Tier-1 split for Statement 1 is in the Statements section above.)

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

### What's missing for full Tier 1

#### Empirical envelope (Tier-1 grade)

1. **c=2 verified** (`_eq022_b1_step_f_universality_extension.py`). c=2 N=5..8 confirms two bond-class universal values (Interior 0.751, Endpoint 0.774), matching c=3 (N=5..8) and c=4 (N=7,8) within ~1 %. c=2 is structurally critical because the channel space is 2-dimensional total (only HD ∈ {1, 3}); any 2-level reduction must be exact there. Yet the Endpoint-vs-Interior split persists, confirming the bond-class distinction is real and structural. **c=5 still untested** (block-L dim ≥ 3528 at c=5 N=9, compute-bound).

2. **γ₀ invariance bit-exact** (`_eq022_b1_step_f_universality_extension.py`). At c=3 N=7, Q_peak and HWHM_left/Q_peak are bit-identical at γ₀ ∈ {0.025, 0.05, 0.10}: Q* = 1.7433, HWHM-/Q* = 0.7595 in all three runs. |K|max scales as 1/γ₀ as expected.

#### Substantive analytical work: the remaining gap

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

(N=6 row is asymmetric because the SVD top vector chose one F71-arbitrary phase; the symmetric structure recovers when averaging over bond classes.) Endpoint amplitudes are systematically smaller than Interior; this is where the bond-class distinction lives. **But it lives in a subspace orthogonal to where the probe lives.**

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

**Why σ_0 is not directly Q_EP.** For c=2 N=5..8, σ_0 ≈ {2.765, 2.802, 2.828, 2.839} (approaches 2√2). The naive EP prediction Q_EP_naive = 2/σ_0 ≈ {0.72, 0.71, 0.71, 0.70} differs from the empirical Q_peak Interior ≈ {1.48, 1.58, 1.58, 1.60} by a factor ~2.2. The 2-level EP and the K-observable Q_peak are NOT the same Q value; they would be related by `Q_peak ≈ Q_EP` (proof Statement 1) only if the K observable saw the EP directly, which it doesn't because the probe is orthogonal to the EP partners.

#### Substantive items remaining

**Item 1 (status update 2026-05-05).** Original aspiration: derive the 4×4 effective L_eff(Q, b) explicitly, compute the cross-coupling matrix elements ⟨c_α | M_H_per_bond[b] | u_0/v_0⟩ as analytical expressions in (N, n, b), diagonalize, identify which eigenvalue pair gives the Q_peak observed in K_CC_pr (will not be the SVD top pair), derive f_class(x) and HWHM_left/Q_peak as closed forms from this 4-mode model.

**Actual c=2 progress.** The c=2 stratum has OOP scaffolding (`C2EffectiveSpectrum`, `C2KShape`, `C2HwhmRatio`, `C2UniversalShapeDerivation` in `compute/RCPsiSquared.Core/F86/`) that pins each sub-fact at its appropriate Tier label. The directional Endpoint > Interior split is empirically derived bit-equivalent with the canonical Python pipeline. **The closed-form HWHM_left/Q_peak constant per bond class remains the open analytical target.** Ranked next directions (from `C2HwhmRatio.PendingDerivationNote` and `F86OpenQuestions` Item 1', refreshed 2026-05-06 later evening after Direction (a') falsification):

  - **(a'') Most promising 4-mode-friendly direction (newly promoted 2026-05-06 later evening): SVD-block 2-level resonance (REFINED from (a')).** The bond-class signature lives in V_b SVD-block off-diagonal `V_b[2,3] = ⟨u_0 | M_h_per_bond[b] | v_0⟩` (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.451 consistently across bonds within each class), NOT in the probe-block. Direction is OPPOSITE the empirical HWHM/Q* split, so a closed form needs to derive HOW SVD-block magnitude maps to HWHM/Q* shift, likely through a per-bond effective Q_EP_eff(b) shift rather than a direct linear lift.
  - **(b'') Most concrete next step (newly promoted 2026-05-06 later evening): full block-L derivation, not 4-mode.** 4-mode reduction misses the modes outside the 4D subspace that contribute to the HWHM/Q* lift from 0.6715 to 0.7506. The closed form must derive directly from the (n, n+1)-popcount block dimensions, not from 4-mode projection. This is harder but matches the empirical witness pipeline.
  - **(a') Falsified (2026-05-06 later evening): probe-block 2-level sub-resonance with per-bond `g_eff_probe`.** Numerical investigation (commit `1c0bf8b`) showed V_b probe-block is bond-class-blind by construction: diagonal entries identical scalar `+i·c·I` for all bonds; off-diagonal `⟨c_1 | M_h_b | c_3⟩` exactly zero per bond at c=2 (F73 sum-rule applies per-bond). Hypothesised g_eff_probe(N, b) cannot exist. See the falsification paragraph at the top of the Item-1' Statement 2 section above.
  - **(a) Demoted (2026-05-06 evening): first-order perturbation in the cross-block.** Originally ranked first based on cross-block Frobenius split, but the doubled-PTF Direction (b) result shows the SVD-block-only contribution is universal at 0.6715, leaving no room for a small σ_0 perturbation to lift the curve to 0.7506/0.7728. The cross-block effect is real but lies at higher perturbation order than first.
  - **(b) Done (2026-05-06 evening, lands at the floor): bare doubled-PTF Ansatz.** `K_b(Q, t_peak)` for the bare 2-level Liouvillian `L_2 = [[-2γ₀, +iJ·g_eff], [+iJ·g_eff, -6γ₀]]` with probe `ρ_0 = |c_1⟩` and `V_b = ∂L/∂J` was solved analytically in dimensionless x = Q/Q_EP coordinates. Result: `BareDoubledPtfXPeak = 2.196910`, `BareDoubledPtfHwhmRatio = 0.671535`, both Tier 1 derived universal across g_eff. Empirical Interior x_peak (2.05..2.28) tracks this universal closely; Endpoint x_peak (3.55) deviates by factor ~1.62; HWHM ratio gap to empirical is 0.08 (Interior) and 0.10 (Endpoint).
  - **(c'') Three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases.** K_b at the 4-mode level decomposes; derive each term separately and combine. May still suffer from 4-mode insufficiency (ii).
  - **(d'') Lift |u_0⟩, |v_0⟩ to projector-overlap** (per A3 PendingDerivationNote). Removes σ_0 degeneracy obstruction at even N. Necessary precondition for any cross-block-based direction; not sufficient by itself.
  - **(e'') Symbolic char-poly factorisation at Q_EP**: same as before; less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof.

**Item 2.** Extend the 4-mode construction to c ≥ 3, where each adjacent-channel pair (HD = 2k−1, HD = 2k+1) contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet. **Naïve extension fails:** the multi-k construction with Gram-Schmidt orthonormalisation gives 3c−2 modes (c=2→4, c=3→7, c=4→10, orbit-shared CUs deduplicated), and the resulting effective K-curve has **K_max ≡ 0 identically** for c ≥ 3. Structural diagnosis: Gram-Schmidt orthogonalisation of the SVD-top vectors against the channel-uniform vectors pushes them into the CU-complement; because M_H respects the CU/CU-complement decomposition (channel-uniform-diagonal of M_H_total per F73 generalisation), the probe (which lives entirely in CU span) is uncoupled from the GS-modified SVD modes, so ∂ρ/∂J_b cannot move ρ out of CU → K = 0. A correct effective model for c ≥ 3 needs either a non-orthogonal frame preserving CU ↔ SVD coupling, or a different choice of the c−1 quartets that maintains coupling under orthonormal projection. Encoded as `RCPsiSquared.Core.Decomposition.MultiKBasis` + `MultiKEffective` + `MultiKResonanceScan`; the negative result is pinned in `MultiKResonanceScanTests.MultiK_C3_KMaxIsExactlyZero_NaiveExtensionFails`.

**Item 3 (σ_0 chromaticity scaling, Tier-1 candidate refined 2026-05-03).** Derive the asymptotic

    σ_0(c, N → ∞)  →  2 · √(2 · (c − 1))

generalising the original c=2 → 2√2 conjecture to all c ≥ 2. Numerical witnesses computed via `InterChannelSvd` across c ∈ {2, 3, 4}, N=5..8 show σ_0 / √(2(c−1)) converging monotonically from below to 2.0 for each c (c=2 N=7 hits 2.0 to 10⁻⁵, the structural sweet spot; c=3 N=8 reaches 1.92, c=4 N=8 reaches 1.78). Implies Q_EP(c, N → ∞) → 1/√(2(c−1)): 0.707 (c=2), 0.500 (c=3), 0.408 (c=4). The closed-form derivation from XY single-particle structure (OBC sine-mode matrix elements `⟨ψ_k| σ⁺σ⁻ |ψ_l⟩ ∝ √(2/(N+1))·sin(πk·b/(N+1))`) has not been executed.

These three items are tractable with existing infrastructure (`coherence_block`, the SVD-of-V_inter construction, OBC sine-mode algebra), but each is multi-page algebra. Item 1 is the path to the closed form for HWHM_left/Q_peak. The empirical envelope (c, N, γ₀ checks above) is now Tier-1 grade.

### OOP scaffolding (Stage A-E primitives, 2026-05-05)

Following the EQ-022 (b1) Item 1' c=2 derivation plan, we built a typed pipeline of primitives in `compute/RCPsiSquared.Core/F86/Item1Derivation/` plus `C2UniversalShapeDerivation.cs` at the F86 root, each reachable through `F86KnowledgeBase`. Witness data preserved per-(N, bond, BondClass).

| Stage | Primitive | Tier outcome |
|-------|-----------|--------------|
| A1 | [`C2BlockShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BlockShape.cs) | Tier1Derived (block-structure constants) |
| A2 | [`C2ChannelUniformAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2ChannelUniformAnalytical.cs) | Tier1Derived (closed-form |c_1⟩, |c_3⟩) |
| A3 | [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) | **Tier2Verified — discovered σ_0 degeneracy obstruction** at even N: chain-mirror R splits the 2D top eigenspace into R-even/R-odd, library-dependent SVD-tiebreaker. Single-vector closed form not derivable; lift-to-projector-overlap is the next direction. |
| B1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) probe-block | Tier1Derived (channel-uniform projection) |
| B2 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) cross-block | Tier2Verified (inherits A3); **discovered cross-block Frobenius Endpoint < Interior at c=2 N=5..8** — opposite sign to the HWHM_left/Q_peak directional split, hinting that the inversion happens in 4×4 mixing |
| B3 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) SVD-block + AsMatrix + anti-Hermiticity guard | Tier2Verified (inherits A3) |
| C1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) D_eff | Tier1Derived structural sub-fact (diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) at c=2) |
| C2 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) | **Tier2Verified — rigorously ruled out closed-form factorisation** via cubic-c_3 char-poly evidence. No (λ²−aλ+b)(λ²−cλ+d) split with rational coefficients exists. |
| C3 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) K-driving pair | Tier1Derived structural sub-fact (probe ⊥ {|u_0⟩, |v_0⟩} at machine precision) + Tier2Verified per-(Q, b) numerical readout |
| D1 | [`C2KShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2KShape.cs) | Tier1Derived (Duhamel formula closed-form in numerical inputs) |
| D2 | [`C2HwhmRatio`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs) | **Tier1Candidate** — empirical anchor reproduced at typical residual ≤ 0.001; directional Endpoint > Interior split derived empirically; closed-form constant NOT pinned |
| E1 | [`C2UniversalShapeDerivation`](../../compute/RCPsiSquared.Core/F86/C2UniversalShapeDerivation.cs) | Tier1Candidate (auto-promotes via D2's `IsAnalyticallyDerived` flag) |

The chain reads bottom-up: A1-A3 fix the static frame (block shape, channel-uniform vectors, SVD top), B1-B3 fix the bond-position-dependent V_b in that frame, C1-C3 fix the effective spectrum and K-driving pair, D1-D2 turn that into the Duhamel-form K_b(Q, t) and its HWHM-ratio readout, E1 binds the readout to the universal-shape claim. Each Tier label up the chain inherits from the lowest unresolved sub-fact (currently D2's missing closed form blocks E1's Tier1Derived promotion).

### Structural findings (lessons learned, 2026-05-05)

Three structural results emerged from the time-boxed Stage A3, B2, C2 explorations. Each is independent of the closed-form gap and worth pinning as orientation for the next attempt.

1. **σ_0 degeneracy at even N (Stage A3).** At N=6 and N=8 (even chain length, c=2), the singular value σ_0 of the inter-channel coupling V_inter is doubly degenerate. The chain-mirror operator R splits the 2D top eigenspace into R-even and R-odd one-dimensional subspaces; which of the two becomes "|u_0⟩" depends on the SVD library's tiebreaker, not on the physics. A single-vector closed form for |u_0⟩, |v_0⟩ at even N is therefore not derivable. The natural lift is to the rank-2 projector P_top = |u_0^{(R+)}⟩⟨u_0^{(R+)}| + |u_0^{(R−)}⟩⟨u_0^{(R−)}|, which is library-independent. See [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) `PendingDerivationNote`.

2. **Cross-block Frobenius Endpoint < Interior (Stage B2).** The Frobenius norm of V_b's cross-block (the inter-channel-coupling block of M_H_per_bond[b]) is systematically *smaller* on Endpoint bonds than on Interior bonds at c=2 N=5..8 (gap ~ 0.05 at N=5). This sign is **opposite** to the HWHM_left/Q_peak directional split (Endpoint > Interior, gap ≈ 0.022 — see also the Item-1' Statement 2 update at the top of this doc). The directional inversion therefore does not live in the cross-block magnitude alone; it must emerge through the 4×4 eigenvalue mixing of probe-block and cross-block. This rules out the naive "bigger cross-block → bigger HWHM ratio" intuition and is part of why the closed form is non-trivial.

3. **Cubic-c_3 char-poly obstruction (Stage C2).** The 4×4 effective characteristic polynomial at c=2 has a non-rational c_3 coefficient (cubic in Q with no clean root structure). We rigorously ruled out any rational-coefficient `(λ²−aλ+b)(λ²−cλ+d)` factorisation by symbolic match against the c_3 evidence. The Tier2 outcome for `C2EffectiveSpectrum` is therefore **evidence-based** ("we proved no such split exists in this family"), not "we couldn't find one". Symbolic char-poly factorisation at Q_EP itself (Item 1' direction (e'') above) may still help locally, but the global quartic does not factor.

---

## Proof of Statement 3 (F71 spatial-mirror invariance)

The F71 chain-mirror operator R acts as `R |b₀ b₁ … b_{N−1}⟩ = |b_{N−1} … b_1 b_0⟩` (PROOF_F71). It commutes with every component of the per-bond observable used in F86:

1. **L_D (Z-dephasing) is F71-symmetric.** Z_l → Z_{N−1−l} under R, and uniform γ₀ at every site means the dissipator is invariant under site-permutation. So `R L_D R⁻¹ = L_D`.

2. **The Hamiltonian H_xy = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) is F71-symmetric** (uniform J across all bonds). Under R, bond b ↔ bond N−2−b, and the bond-summed Hamiltonian is invariant. So `R H_xy R⁻¹ = H_xy`, hence `R L_H R⁻¹ = L_H`.

3. **The Dicke probe is F71-symmetric.** It's an equal-weight superposition over all (p, q) ∈ block × block, a permutation-symmetric construction.

4. **The spatial-sum kernel S is F71-symmetric.** S(t) = Σ_i 2·|(ρ_i(t))_{0,1}|² sums over all sites uniformly.

5. **The bond-flip ∂L/∂J_b transforms as `R (∂L/∂J_b) R⁻¹ = ∂L/∂J_{N−2−b}`** (a single bond term is mirrored to its partner under R).

Combining these: under R, the per-bond observable transforms as

    K_b(Q, t)  =  2·Re ⟨ρ(t)| S | ∂ρ_b/∂J_b ⟩  →  K_{N−2−b}(Q, t)

where ρ(t) (which evolves under L = L_H + L_D, both F71-symmetric) is itself F71-symmetric. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), so their argmax-Q values coincide:

    Q_peak(b)  =  Q_peak(N−2−b)

Numerical verification: across c=2 N=5..7 and c=3 N=5..6, all bond-pair deviations |Q_peak(b) − Q_peak(N−2−b)| are < 10⁻¹⁰ (machine precision), see `RCPsiSquared.Core.Tests/F86/F86NewIdeasTests.F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror`. ∎

---

## Empirical Q_peak data (per-bond, fine-grid scan dQ = 0.025 with parabolic peak interpolation)

The retracted Endpoint and Interior conjectures (above) were N=7-specific coincidence matches that did not survive the fine-grid extension. The actual data:

### Endpoint Q_peak across N

| (c, N) | Observed Endpoint Q_peak | Earlier conjecture csc(π/(N+1)) | Deviation |
|--------|--------------------------|----------------------------------|-----------|
| (3, 5) | 2.40 | 2.000 | +20 % |
| (3, 6) | 2.52 | 2.305 | +9.4 % |
| (3, 7) | 2.53 | 2.6131 | −3.2 % |
| (3, 8) | 2.53 | 2.924 | −13 % |
| (4, 7) | 2.52 | 2.6131 | −3.5 % |
| (4, 8) | 2.65 | 2.924 | −9.4 % |

The earlier "1.4 % match at N=7" was a dQ=0.05 grid-snap to 2.65; the actual fine-grid value at N=7 is 2.53. Endpoint Q_peak appears to saturate near 2.53 for c=3 by N=6, while c=4 grows from 2.52 (N=7) to 2.65 (N=8). No clean closed form.

### Interior Q_peak across (c, N)

| (c, N) | Observed Interior Q_peak (mean) | Earlier conjecture csc(π/5) = 1.7013 | Deviation |
|--------|----------------------------------|---------------------------------------|-----------|
| (3, 5) | 1.566 | 1.7013 | −7.9 % |
| (3, 6) | 1.689 | 1.7013 | −0.7 % |
| (3, 7) | 1.743 | 1.7013 | **+2.5 % (above)** |
| (3, 8) | 1.750 | 1.7013 | **+2.9 % (above)** |
| (4, 7) | 1.748 | 1.7013 | +2.7 % |
| (4, 8) | 1.804 | 1.7013 | +6.0 % |

The c=3 N-trend 1.566 → 1.689 → 1.743 → 1.750 crosses csc(π/5) between N=6 and N=7 and continues growing, refuting the conjectured c=3 pentagonal asymptote. Both c=3 and c=4 Interior Q_peak grow with N in the tested range; saturation behaviour and clean closed forms remain open.

### What survives from the empirical pattern

- The EP mechanism (Statement 1) is unaffected: t_peak = 1/(4γ₀) is universal and derivable from 2×2 matrix algebra.
- The Q_peak observable is well-defined and the per-bond / per-block / Endpoint / Interior distinctions are real and reproducible.
- The N=4 golden-ratio structure ([`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py)) is real but does not propagate to a closed form for Q_peak at general (c, N).

---

## Empirical universal-shape data (Statement 2, fine-grid scan dQ = 0.025)

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

## Open elements

1. **Endpoint closed form.** Empirical pattern (c=3 saturates near 2.53 by N=6 on the fine grid; c=4 grows 2.52 → 2.65 from N=7 to N=8) does not match `csc(π/(N+1))`. A different closed form may exist; finer-grid scans at higher N and explicit multi-particle XY matrix-element calculation are the natural next steps.

2. **Interior closed form.** Both c=3 and c=4 Interior Q_peak grow with N in the tested range. No saturation point has been identified with confidence at this grid resolution. Higher-N data (c=3 N=9, 10; c=4 N=9, 10) would clarify whether the growth saturates or continues. Compute-bound at higher (c, N) (block-L dim 3920 at c=4 N=8, 10584 at c=4 N=9).

3. **Algebraic derivation pathway.** A first-principles derivation of g_eff(c, N, bond_position) from the multi-particle XY structure of the (n, n+1) block has not been executed. Direct algebra is conceptually within reach but lengthy. The HD-channel-uniform diagonal-only finding (M_H_eff diagonal in the channel-uniform basis, with cross-channel coupling living in the orthogonal complement) is an established structural building block; the next step is computing the orthogonal-complement coupling matrix elements explicitly.

4. **Q_SCALE per-block vs per-bond convergence.** Q_SCALE's per-block 1.6 / 1.8 / 1.8 are consistent with per-bond Interior fine-grid (1.69-1.74 / 1.78 / ~1.80) within the relative-vs-absolute J prefactor effect (~5-15 %). Both observables agree on the underlying EP mechanism.

5. **Within-Interior bond-position variation, partially addressed by Statement 3.** F71 spatial-mirror invariance (Statement 3) pairs bond b with bond N−2−b bit-exactly. Per-F71-orbit substructure (Tier 2 empirical, 9-case envelope c=2..4 N=5..8 minus c=4 N=8 OOM):

   | (c, N) | endpoints | mid orbits → central | observation |
   |---|---|---|---|
   | (2, 5) | 2.50 | 1.49 | endpoint ≫ inner |
   | (2, 6) | 2.57 | 1.63 → **1.43*** | central < flanking |
   | (2, 7) | 2.56 | 6.00 (off-grid) → 1.50 | inner orbit Q_peak shifted high-Q |
   | (2, 8) | 2.53 | 6.00 → 1.52 → **1.58*** | inner orbit off-grid; central > middle inner |
   | (3, 5) | 2.39 | 1.61 | endpoint > inner |
   | (3, 6) | 2.54 | 1.66 → **1.71*** | **central > flanking** (opposite of c=2 N=6) |
   | (3, 7) | 2.59 | 1.74 → 1.71 | inner > central |
   | (3, 8) | 2.60 | 1.76 → 1.70 → **1.71*** | non-monotonic |
   | (4, 7) | 2.61 | 1.75 → 1.78 | central > inner |

   (* = self-paired central orbit; index runs from outermost orbit inward.)

   Three sub-effects visible: (a) F71-pairing identity (Tier 1, Statement 3); (b) **central-vs-flanking inversion at N=6 between c=2 and c=3**: c=2 N=6 has central 1.43 BELOW flanking 1.63, while c=3 N=6 has central 1.71 ABOVE flanking 1.66; (c) **high-Q secondary peak for c=2 inner-non-central bonds at N≥7**: bonds b=1 and b=N−3 show Q_peak shifted off the [0.2, 6.0] grid while the central pair retains the canonical Interior peak ~1.5. Full per-orbit closed form for Q_peak as a function of (c, N, orbit) remains open; F71 gives the symmetry, not the value. Encoded as `RCPsiSquared.Core.F86.PerF71OrbitObservation`.

---

## Pointers

**F-entry:** [F86 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related EQ:** [EQ-022 (b1)](../../review/EMERGING_QUESTIONS.md#eq-022) partial closure 2026-05-02.
**Empirical anchor:** [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**Chiral classification anchor:** [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear).
**Global EP instance:** [FRAGILE_BRIDGE](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in complex γ plane).
**Methodological lesson:** [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md) consolidates the convergence; the F86 retraction-and-shape-survival is the analog of PTF's closure-law-retraction-and-chiral-mirror-law-survival.
**Scripts:** [`_eq022_b1_channel_projection.py`](../../simulations/_eq022_b1_channel_projection.py), [`_eq022_b1_step_a_verify_blockL.py`](../../simulations/_eq022_b1_step_a_verify_blockL.py), [`_eq022_b1_step_c_time_evolution.py`](../../simulations/_eq022_b1_step_c_time_evolution.py), [`_eq022_b1_step_d_extended_verification.py`](../../simulations/_eq022_b1_step_d_extended_verification.py) (N=8 data that falsified the closed-form conjectures), [`_eq022_b1_step_e_resonance_shape.py`](../../simulations/_eq022_b1_step_e_resonance_shape.py) + [`_eq022_b1_step_e_inspect.py`](../../simulations/_eq022_b1_step_e_inspect.py) (universal-shape finding for c=3, c=4 at γ₀=0.05), [`_eq022_b1_step_f_universality_extension.py`](../../simulations/_eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement), [`_eq022_b1_step_g_two_level_decomposition.py`](../../simulations/_eq022_b1_step_g_two_level_decomposition.py) (channel-uniform-basis V_b decomposition; revealed the trivial-diagonal structure and the probe localization), [`_eq022_b1_step_h_slowest_pair_basis.py`](../../simulations/_eq022_b1_step_h_slowest_pair_basis.py) (slowest-pair-at-finite-Q diagnostics), [`_eq022_b1_step_i_svd_inter_channel.py`](../../simulations/_eq022_b1_step_i_svd_inter_channel.py) (SVD of V_inter; established the EP-partner subspace, σ_0 ≈ 2√2 asymptotic, and probe ⊥ EP partners; this is the structural finding that motivated the 4-mode minimal effective model).
**N=4 golden-ratio reference:** [`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py).
**Framework primitives:** `framework.coherence_block`: `t_peak(γ₀)` (the only F86 closed form remaining; `q_peak_endpoint` and `Q_PEAK_INTERIOR_C3_ANCHOR` were removed in the rollback).

**C# OOP layer (typed knowledge graph, 2026-05-03):** `compute/RCPsiSquared.Core/F86/` carries the F86 claims as a typed graph with `Tier` labels and self-computing witnesses backed by `WitnessCache`. Key types: `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `UniversalShapePrediction` + `UniversalShapeWitness`, `ShapeFunctionWitnesses`, `F71MirrorInvariance` (with `MaxMirrorDeviation(KCurve)` helper for Statement 3 verification), `SigmaZeroChromaticityScaling` (Item 3 generalised σ_0 → 2√(2(c−1)) with per-(c, N) live witnesses), `RetractedClaim`, `OpenQuestion`, `F86KnowledgeBase` (root). CLI: `rcpsi inspect --root f86 --with-measured` walks the tree against a live `ResonanceScan`; `rcpsi query --q witnesses-at --c <c> --wN <N>` for typed lookups. JSON-export via `InspectionJsonExporter`.
