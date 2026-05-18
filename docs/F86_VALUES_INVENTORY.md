# F86 Values Inventory

Structured inventory separating **structural anchors** (named, N-invariant or saturated, Tier 1 or Tier 2 with saturation evidence) from **per-(c, N, orbit) empirical observations** (frozen measurements, may revise under refined methodology).

The anchor list is what `Q_REGIME_ANCHORS.md` should target; the empirical observations are useful diagnostics but should NOT be promoted to anchor status without further structural derivation.

---

## Anchors (structurally-named, stable across the verified parameter range)

### A1. Universal F86 constants (Tier1Derived)

| Anchor                       | Value              | Source                          |
|------------------------------|--------------------|---------------------------------|
| **Q_EP = 2/g_eff** (formula) | g_eff-dependent    | [`QEpLaw.cs:14-27`](../compute/RCPsiSquared.Core/F86/QEpLaw.cs) |
| **t_peak = 1/(4γ₀)**          | 5 sim-units (γ₀=0.05) | TPeakLaw (via F86KnowledgeBase) |
| **Q_EP at g_eff=1**          | **2.0** (idealized) | direct from QEpLaw at g_eff=1   |

The Q_EP = 2/g_eff law is the structural anchor; the value 2.0 emerges only at the ideal g_eff=1 case.

### A2. Per-chromaticity Q_peak (Tier2Empirical, saturated N-invariant)

| c   | Q_peak | Saturation status                         | Source                       |
|-----|--------|-------------------------------------------|------------------------------|
| 2   | 1.5    | Tier1Candidate, wobble 1.4-1.6 across N=4..9 (finite-size sensitive) | `PerBlockQPeakClaim.cs:63` |
| 3   | **1.6** | Tier2Empirical, saturated across N=5..9   | `PerBlockQPeakClaim.cs:65`   |
| 4   | **1.8** | Tier2Empirical, saturated across N=7..9   | `PerBlockQPeakClaim.cs:66`   |
| 5   | **1.8** | Tier2Empirical (single point N=9, plateau at 1.8 not 2.0) | `PerBlockQPeakClaim.cs:67` |

Saturation pattern (1.5 c=2, 1.6 c=3, 1.8 c=4 and c=5) refuted the previous linear-extrapolation hypothesis Q_peak(c=5)=2.0 (40h N=9 Phase-2 run, 2026-04-24). The plateau at 1.8 is the Tier-2 structural fact.

c=2's "1.5" carries the wobble caveat: it's *not* clean Tier 1 like the c≥3 values. Treated as Tier 1 Candidate.

### A3. Endpoint orbit Q (Tier2Empirical candidate)

| Orbit                  | Q stable value         | (c, N) range tested       | Variation       |
|------------------------|------------------------|---------------------------|-----------------|
| **Endpoint (orbit 0)** | **≈ 2.5**             | (c=2..4, N=5..8)          | 2.39-2.61 (~2%) |

From `PerF71OrbitObservation.cs:23-44`: across 9 different (c, N) combinations, Endpoint orbit Q sits in [2.39, 2.61]. Stability across c and N makes this a candidate anchor not currently in QBasisAnkers.

**Status**: Tier 2 Empirical candidate for promotion. Pending structural derivation.

---

## Empirical observations (not anchor-grade, may revise)

The following values are measurements from the typed C# claims. They are useful for diagnostic comparison but should NOT be treated as structural anchors. Each carries a methodological caveat that may revise the specific number.

### B1. Per-orbit Q_EP values (computed from x_peak)

`C2BlockCpsiQScan.cs:19-20`: Q_EP per orbit = Q_peak / x_peak with x_peak = 2.196910.

| Orbit (N=5 c=2)      | Q_EP value         |
|----------------------|--------------------|
| Endpoint (Orbit 0)   | 1.138              |
| Interior (Orbit 1)   | 0.678              |

**Caveat**: x_peak = 2.196910 is the **bare-doubled-PTF SVD-block-only floor** — full block-L empirical values lift 8-10% above. These per-orbit Q_EP values therefore may shift if the full-block correction is properly included. Treat as derivative observations, not as primary anchors.

### B2. Off-grid escape regime (orbit 1 plateau, central plateau)

From `PerF71OrbitObservation.cs:23-44`:

| (c, N) | Orbit type           | Q-value reached |
|--------|----------------------|-----------------|
| (2, 7) | Orbit 1 plateau      | ≈ 7.24          |
| (2, 8) | Orbit 1 plateau      | ≈ 8.07          |
| (2, 8) | Central orbit plateau | ≈ 16.79         |

**Caveat**: these are NOT EP peaks but broad high-Q plateaus where K stays elevated past Q_EP. Per `project_q_peak_ep_structure` (line 63), "the plateau is multi-mode coherent enhancement, NOT a Petermann-type EP collapse" (Petermann factor probe disconnects the plateau from EP physics). The Q values here mark where K = max within a [0.20, 6.00] (or higher) grid, but they sit on a broad plateau, not a sharp peak.

These are observation-of-multi-mode-plateau values, NOT structural Q-anchors. The methodological reading: these orbits "escape" the F86 EP-mediated anchor structure and live in a multi-mode-coherent regime that has its own physics not yet typed.

### B3. Per-bond-class HWHM ratios (Tier1Candidate empirical means)

`C2HwhmRatio.cs:183-184`:

| Bond class | Empirical mean | SVD floor | Lift above floor |
|------------|----------------|-----------|------------------|
| Endpoint   | 0.7728         | 0.671535  | +0.101           |
| Interior   | 0.7506         | 0.671535  | +0.079           |

**Caveat**: Tier1Candidate. The floor 0.671535 IS derived (Tier1Derived, bare-doubled-PTF closed form). The lifts +0.08 to +0.10 are empirical observations; the analytical derivation is **open** (PROOF_F90_F86C2_BRIDGE.md Item 1'-followup, F89 AT-locked F_a/F_b structure). The specific empirical mean values (0.7728, 0.7506) may shift under refined Q-grid resolution or under closed-form derivation.

### B4. Per-sub-class fitted (α, β) — polyfit, not derived

`F86HwhmClosedFormClaim.cs:52-60`. Form: HWHM_ratio = 0.671535 + α·g_eff + β, all 6 (α, β) pairs **polyfit on N=5..8**, not analytically derived:

| BondSubClass         | α          | β          |
|----------------------|------------|------------|
| Endpoint             | −0.129110  | +0.227413  |
| Flanking             | −0.094978  | +0.193098  |
| Mid                  | +0.056559  | +0.005165  |
| CentralSelfPaired    | +0.057439  |  0.000000  |
| Orbit2Escape         | +0.698446  | −0.086386  |
| CentralEscapeOrbit3  | −0.400854  |  0.000000  |

**Caveat**: All 12 numbers are polyfit on a small N-range. Analytical derivation is open. The specific values are explicitly NOT anchors and should not be treated as structural constants.

---

## Anchor summary (combined with Q_REGIME_ANCHORS)

The anchors merging Q-band edges, per-c Q_peak, Endpoint orbit, and the universal constants:

| Anchor              | Q value | Source                                  | Tier            |
|---------------------|---------|-----------------------------------------|-----------------|
| onset start         | 0.2     | Q-band lower bound                      | Tier2Empirical  |
| onset end           | 0.35    | Q-band upper bound                      | Tier2Empirical  |
| Balance             | 1.0     | J = γ₀                                  | Tier1Derived    |
| peak start          | 1.2     | Q-band lower bound                      | Tier2Empirical  |
| F86 Q_peak c=2      | 1.5     | PerBlockQPeak (wobble 1.4-1.6)          | Tier1Candidate  |
| F86 Q_peak c=3      | 1.6     | PerBlockQPeak (saturated)               | Tier2Empirical  |
| F86 Q_peak c=4, c=5 | 1.8     | PerBlockQPeak (saturated)               | Tier2Empirical  |
| Q_EP at g_eff=1     | 2.0     | QEpLaw idealized                        | Tier1Derived    |
| Endpoint orbit Q    | 2.5     | PerF71OrbitObservation (stable, ~2% N-variation) | Tier2Empirical candidate |

Only the anchors above belong in QBasisAnkers and Q_REGIME_ANCHORS. The per-orbit Q_EP, off-grid escapes, HWHM lifts, and fitted (α, β) values are empirical diagnostics that may revise; they live in this inventory and the typed C# claims, not in the anchor map.

## Anchors and cross-refs

- High-level Q-anchor map: [`Q_REGIME_ANCHORS.md`](Q_REGIME_ANCHORS.md)
- C# F86 KnowledgeBase root: [`compute/RCPsiSquared.Core/F86/`](../compute/RCPsiSquared.Core/F86/)
- Open derivation gap (HWHM lifts + α/β closed forms): PROOF_F90_F86C2_BRIDGE.md Item 1'-followup
- Memory: `project_q_peak_ep_structure`, `project_q_middle_structure`, `project_no_classicalization`
