# F112 Hardware Lens on Kingston: What a Magnitude Cannot See

**Status: the model family is wrong for this data, so this experiment identifies no noise channel.** What it does show is the F112 and F113 machinery running end to end on Liouvillians fitted to hardware trajectories. The channel question is open and needs a redo.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude
**Data:**
- `data/ibm_block_cpsi_saturation_may2026/...20260508T032749Z.json` (Kingston q13-q14, 5 t-points, idle)
- `data/ibm_cusp_slowing_april2026/...20260416_212042.json` (Kingston A_mid q124-q125 + B_high, 2 pair-runs × 6 t-points, idle)
- `data/ibm_chain_gamma0_april2026/...20260419_110200.json` (Kingston Q12-Q19 chain, 4 pair-runs × 9 t-points, (J/2)(XX+YY) bond Hamiltonian via Trotter, J = 11.55 rad/μs; the JSON field is named `J_rad_per_us` but holds 1.838 MHz, as its own README says)
- `data/ibm_f95_angle_steering_may2026/...omega0.130_...json` and `...omega0.250_...json` (Kingston q82-q83 + q13-q14, 2 omega × 2 pairs × 6 t-points, with applied RZ Z-drive on ONE qubit of each pair)

**Scripts:**
- [`simulations/f112_block_cpsi_analysis.py`](../simulations/f112_block_cpsi_analysis.py) (single-dataset block_cpsi)
- [`simulations/f112_hardware_lens_multi.py`](../simulations/f112_hardware_lens_multi.py) (multi-dataset extension)

**Lens:** F112 polarity-asymmetry diagnostic on fitted effective Liouvillians

## Setup

These are the Tier-A datasets that combine full ρ tomography on a qubit pair with a multi-time-point delay sweep, so F112's `polarity_coordinates_from_L` diagnostic can be evaluated on a fitted effective L without new QPU spend.

The block_cpsi protocol:
- Initial state: (|D_0⟩ + |D_1⟩) / √2 = (|00⟩ + (|01⟩+|10⟩)/√2) / √2 on Kingston qubits 13, 14
- Pure-decoherence idle (no applied Hamiltonian during the delay)
- 5 t-points: 0, 120, 240, 360, 480 μs
- 16 Pauli expectations per t-point (full 2-qubit tomography)
- Calibrated T2_min = 480 μs, so the D[Z] rate reproducing it is γ_eff = 1/(2·T2) ≈ 0.00104 per μs (this line read `1/T2 ≈ 0.00208` until 2026-08-05, which is twice the Lindblad rate; see `docs/GLOSSARY.md`, "The T2 → γ conversion"). The distinction matters below, because a ratio is taken against this number and the 1.72× imported on the next line was computed on 1/(2·T2): two rulers in one setup was exactly the defect.
- Documented anomaly (in [`IBM_BLOCK_CPSI_SATURATION.md`](IBM_BLOCK_CPSI_SATURATION.md); the data directory has no README): hardware C_block decays ~1.72× faster than pure-T2 predicts

Two questions are asked of this data, and they are independent. The **structural** one: does the fitted Liouvillian sit inside F112's typed Tier1Derived scope (Hermitian H + bit_b-homogeneous c), and does the polarity balance hold? The **fit-quality** one: which candidate channel actually explains the trajectory?

## Method

For each candidate noise model, fit the model's parameters to the ρ(t) trajectory via Frobenius² least-squares (Nelder-Mead), then compute the F112 polarity asymmetry on the fitted L:

| Model | Parameters | F112 scope |
|---|---|---|
| `pure_Z` | γ_Z per qubit | in scope (Z is bit_b-homogeneous) |
| `Z + T1` | + σ⁻ amplitude damping per qubit | outside scope (σ⁻ = (X + iY)/2 is bit_b-mixed) |
| `Z + ZZ` | + ZZ-crosstalk Hamiltonian | in scope (ZZ is bit_b-homogeneous, Hermitian) |
| `Z + T1 + ZZ` | combined | outside scope (T1 component) |
| `Z + h_y` | + single-site Y transverse field | in scope (Y is bit_b-homogeneous, Hermitian) |

11 pair-runs × 5 models = 55 fitted Liouvillians. For each fitted L, compute `polarity_coordinates_from_L(L_pauli, N=2, σ=Σγ_Z)` and read `asymmetry`.

## Why the fit-quality question is not answered here

Three defects, in rising order of how much they cost.

**The collapse operator was inverted.** Both scripts built `SIGMA_MINUS = [[0,0],[1,0]]`, commented `|0><1|`. That matrix is `|1⟩⟨0|`, the raising operator. The Kingston trajectories relax: populations on the block_cpsi run move monotonically into |00⟩ (0.485 → 0.865 over 480 μs), and Z-dephasing cannot move a diagonal at all, so amplitude damping was the only channel in the set that could carry that flow, pointed the wrong way. With rates clamped non-negative the optimiser drove γ_T1 to exactly zero, which the first reading took at face value and reported as amplitude damping being rejected by the data. Fixed.

**The superoperator and the state disagreed about stacking.** The generator was assembled row-stack (`H⊗I − I⊗Hᵀ`) while the state was propagated column-stack (`rho.flatten('F')`), so the commutator actually integrated was with −Hᵀ. Invisible for the Y-field and the free-sign ZZ coupling; it flips any fixed real-symmetric H, which includes both the f95 Z-drive and the chain_gamma0 bond Hamiltonian (J/2)(XX+YY). Fixed. Note the polarity diagnostic is deliberately NOT fed the corrected propagator directly: the sign of the F112 asymmetry is pinned to a row-stack L read against an `order='F'` Pauli transform (see `PauliBasis.cs`), so the scripts convert to that representative at the diagnostic call and keep the correct propagator for the fit.

**The model family has no per-qubit detuning, and this data needs one.** On block_cpsi, the second tensor slot's transverse expectations alternate in sign from sample to sample while their magnitude decays smoothly (the archived data gives the pair as [13, 14] but states no slot ordering, so which physical qubit that is remains an assumption; the structural argument below does not depend on it, but two of the numerical comparisons do, and they say so where they stand):

```
t (μs)      0      120     240     360     480
<I,X>    0.665  -0.352   0.149  -0.064   0.029
<I,Y>    0.010  -0.136   0.188  -0.080   0.024
```

That is a coherence turning about half a revolution per 120 μs sample, aliased by the sampling. Qubit 13 turns too, more slowly and without aliasing: its transverse angle walks −35.3°, −35.2°, −36.2° and then −48.9° across the four steps. Neither dephasing nor amplitude damping can turn a coherence at all, and a ZZ coupling turns both qubits at nearly the same rate, which these two do not do. So no member of the family above can represent what the data plainly does, and whatever the optimiser cannot represent lands in the parameters it has. Against a sign-alternating sequence the least-squares optimum for a purely decaying model is **zero**, so the optimiser drives that qubit's γ_Z upward until the model predicts no transverse coherence at all and then stops. It is not fitting a rate; it is switching a prediction off. What it gives up is not negligible: at the plateau the model predicts zero for *every* Pauli carrying X or Y in that slot, and their measured content is 0.1461. That is a fixed quantity, so its share of the residual depends on which model you divide by, and the honest way to quote it is both ends: **18.0%** of `pure_Z`'s loss 0.8101, the worst-fitting model, and **98.8%** of `Z_plus_T1_plus_ZZ`'s 0.1479, the best-fitting one. In the model this experiment would otherwise call its winner, essentially the entire residual is transverse content the fit refuses to explain. Within the discarded set the largest single term is ⟨Z,X⟩ at 0.0435, ahead of ⟨I,X⟩ at 0.0378, so the two rows tabulated above are not the whole of it.

**The runaway value must not be read as a rate, and a measured rate says why.** The aliased slot's transverse *magnitude* √(⟨I,X⟩² + ⟨I,Y⟩²) decays cleanly, 0.665 → 0.377 → 0.240 → 0.102 → 0.037, log-linear with **T₂ ≈ 170 μs**, i.e. γ_Z ≈ 0.0029 per μs. So that qubit's coherence is emphatically alive at the first sampled delay (57% of its t=0 value at 120 μs), and its dephasing is identified by the data at about 2.8× the calibrated rate. The fitted parameter, by throwing away the sign, throws away the magnitude with it and lands a factor of 54 away (0.1597 against 0.0029), or 85× further out again in the one model that runs to 13.59. Read as a bound it would imply T₂ < 10 μs on a qubit whose trajectory in this same file gives T₂ ≈ 170 μs, which is the misreading to avoid. (The 170 μs is this document's log-linear fit to that trajectory, not a number the file states; the file's stated calibration is `t2_min_us_calibration` = 480 μs, and it is a Z-only envelope, so the `Z_plus_T1` fit's σ⁻ = 0.00321 would move the dephasing-only part to ≈ 0.0021.)

**Where the retracted number came from.** This paragraph quoted "between 0.138 and 0.147 per μs, about 140×" until 2026-08-05, and before that "about 70×", the latter because it divided by 1/T2 rather than by the Lindblad rate. Neither was a fitted coupling, and the reason is arithmetic rather than physics: around γ_Z ≈ 0.15 the model's predicted coherence at the first delay, exp(−2·γ_Z·120), falls under the ULP of the objective, and from there upward the objective is bit-identical to its γ_Z → ∞ limit. So the simplex stops, and the retracted range sits at exactly that onset. **The published "140×" was a readout of machine epsilon.** No figure is quoted for the onset itself, because it carries no digits: it is model-dependent, and at fixed model it moves by 8% under a relative 10⁻⁸ jitter of the other rate. The precision-free statement is the flattening: at 0.05 per μs the objective is already within 1.96·10⁻⁶ relative of its limit. Multi-start returns 0.1597, 0.1483, 0.1515 and 13.59 per μs for the four models at RMS unchanged to six digits, which is the plateau showing itself.

Removing the ratio does not weaken the point being made, it sharpens it: a parameter that has stopped denoting a rate is a stronger statement than one landing at 140×. Two things corroborate the reading directly. `Z_plus_hy` has no plateau at all: its second Z rate sits at an interior minimum near 0.0152, with γ_Z → ∞ about 6% worse. **That escape has nothing to do with the transverse alternation**, and the decomposition says so rather than a story about it. The interior minimum's total advantage over the plateau, taken with the other two parameters held at their fitted values, is +0.0275, and it splits into **+0.0277 gained on the slot-1 *longitudinal* Paulis** (⟨I,Z⟩, ⟨X,Z⟩ and ⟨Z,Z⟩ lead, contributing +0.0116, +0.0113 and +0.0062) **minus 0.0002 LOST on the transverse block**, where the fit is actually worse than the plateau's exact zero. The part exceeds the whole because the two have opposite signs. ⟨I,X⟩ itself accounts for 3% of the gain. So the one model that escapes does not escape by explaining what the section is about.

`Z_plus_ZZ` does turn the transverse coherence, and still plateaus, for the reason given where the aliasing is described above: it turns both slots by the *same* angle, reaching −34.4° and −34.3° at the first sample. That tracks slot 0's measured angle at that sample (−34.2°) almost exactly, which is precisely why it cannot also carry slot 1, whose step is the near-half-revolution above. One knob, two different rates. So neither escape reproduces the differential alternation the data shows, which is what "no member of the family can represent it" means concretely.

*(Two tidy causal stories were written for the `Z_plus_hy` exception and both were withdrawn under review, on 2026-08-05: first that it is the only model that can turn a coherence, which is backwards, and then that its Y-field gives ⟨I,X⟩ a second way down, which the decomposition above refutes. The measurement was right both times and the mechanism laid over it was the defect. It is left as a decomposition on purpose.)* And in `Z_plus_T1` the σ⁻ rates are identified and land near calibration, 0.00256 and 0.00321 against 1/T₁ = 0.00299 and 0.00341 in slot order q13, q14: only the second Z rate is broken. `simulations/f112_block_cpsi_analysis.py` now prints a per-parameter identifiability profile beside the fit, so a flat direction cannot be read as a measurement again; every identified parameter in that table moves the objective by a relative 5·10⁻² to 4.5 under a factor-of-ten rescaling, against exactly 0 in the ×10 column for the runaway one in three models and −1.8·10⁻¹⁶ in the fourth, where the fit sits a hair below its own onset.

One ruler note, since this document's own lesson is which ruler, and it is the first of the two places that lean on the slot ordering. The calibrated 1/(2·T2) ≈ 0.00104 uses **T2_min = 480 μs, which is q13's**; q14, the aliased qubit the comparison is about, has its own calibrated T₂ = 511.5 μs (`simulations/results/block_cpsi_run_planner_2026_05_08.txt`). That does move a number above: against q14's own ruler the measured 2.8× becomes **3.01×**. The second place is the σ⁻ comparison above, which pairs the fitted rates with 1/T₁ in slot order q13, q14; read in the other order it would give 25% and 7% agreement instead of 14% and 6%.

A re-fit with the corrected scripts is therefore not enough. Every RMS ranking this experiment produced is a comparison within the wrong family, and no replacement channel identification is offered. **A redo needs:** one free Z coefficient per qubit, converged multistart, model comparison at matched parameter count (AIC/BIC), and uncertainties. Multistart landed for the single-dataset script on 2026-08-05 and changed no RMS there: that table was already nesting-clean under single-start, and multistart confirmed it rather than fixing it. The nesting violations, a superset fitting worse than its subset, are in `f112_hardware_lens_multi.py`, which is still single-start (on `chain_gamma0 Q13-Q14`, `Z_plus_T1` gives RMS 0.628598 against `pure_Z`'s 0.628012). Adding multistart there alone would make those rows quieter without making any ranking mean more, so it waits for the rest of the redo.

Two of the four datasets need more than that. The chain_gamma0 "pair-runs" are four ADJACENT PAIRS of one five-qubit chain carrying a propagating single excitation, so each is an open, overlapping marginal and no closed two-qubit Lindblad describes it, whatever channels are added. And on block_cpsi the per-qubit detuning is not identifiable from a uniform 120 μs grid: it aliases into equally good minima spaced π/Δt apart. That dataset needs non-uniform delays before the missing term can be fitted at all.

## What the structural thread shows

**In-scope balance is a theorem being evaluated, not a measurement.** F112 says the asymmetry is exactly 0 for any in-scope model at any parameter values. All 33 in-scope fits give 0 bit-exactly, which is one identity evaluated 33 times rather than 33 independent confirmations. It is still worth having: it is an end-to-end check that the pipeline, the Π construction and the Pauli transform agree with the closed form on Liouvillians nobody hand-built.

**The nonzero readings are where F113 says they must be, but the split is set by the model file, not the chip.** Every out-of-scope fit with nonzero asymmetry is an f95 run, the one dataset whose protocol applies a Z-drive at all, and it applies it to one qubit of the pair. That is what F113 predicts, since the asymmetry reads H only through Tr(Z_l H). But the single-site Z enters those models as a known drive handed to the fitter, and no other model in the set contains a single-site Z at all, so the split follows by construction. Once the missing detuning is added, every run carries a single-site Z and the split is expected to disappear. The individual magnitudes came from the withdrawn fits and go with them.

The part that is clean is the synthetic isolation below, because it evaluates the formula at chosen parameters instead of fitted ones (ω = 0.13, γ_T1 = 0.001, N = 2). It is reproduced by the F113 derivation script (`simulations/f113_break_formula_derivation.py`), not by the two fitters:

```
Case A: Z-drive H + Z-deph             rel asym = 0.000        (in F112 scope: BALANCED)
Case B: idle H + σ⁻ T1                 rel asym = 0.000        (out of scope, still BALANCED)
Case C: Z-drive H + σ⁻ T1              rel asym = 3.85e-03     (out of scope, BROKEN)
Case D: idle H + Z-deph + σ⁻ T1        rel asym = 0.000        (out of scope, still BALANCED)
Case E: Z-drive H + Z-deph + σ⁻ T1     rel asym = 3.85e-03     (out of scope, BROKEN)
```

C and E agree exactly, so the Z-dephasing term contributes nothing to the +i / −i content; the breaker is the Z-drive against the σ⁻ collapse operator, through [Z, σ⁻] = +2σ⁻.

## Why none of this showed

The quantity the structural thread reports is `rel = |asymmetry| / ‖M‖²`, an absolute value, identical for σ⁻ and σ⁺. A magnitude cannot see a direction. The RMS ranking did move when the operator was corrected, so the pipeline was not blind to it in principle; what was missing is any check that pointed OUTWARD, comparing a fitted rate against the device's own calibration. That is the transferable lesson, and it is worth more than the reading it cost.

## Connections

- **F112** (`docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md`, `compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs`): the typed Tier1Derived theorem, unaffected by anything here.
- **F113** (`compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBreakMagnitude.cs`, [`F113_BREAK_MAGNITUDE_FORMULA.md`](F113_BREAK_MAGNITUDE_FORMULA.md)): the closed form for the break magnitude, and the selector Tr(Z_l H) behind the split above.
- **[F113 T1 extraction on Kingston](F113_T1_EXTRACTION_KINGSTON.md)**: a γ_T1 extraction on the f95 data with the correct operator. Its rates survive a per-qubit detuning term (they move by under 3%), which is the check this survey has not done for its own.
- **Block-CΨ saturation lens** ([`IBM_BLOCK_CPSI_SATURATION.md`](IBM_BLOCK_CPSI_SATURATION.md), Confirmations entry `block_cpsi_saturation_kingston_may2026`): the original lens on this dataset. Its measured numbers do not depend on the noise-model question and are unaffected; its attribution of the 1.72× gap to gate noise and readout is its own. This experiment does not adjudicate it, but one constraint is worth carrying across: the 1.72× is a decay RATE against a free-intercept fit, so a time-independent MULTIPLICATIVE state-preparation or readout error cannot move it, being absorbed by the free amplitude. (A time-independent additive offset is not absorbed and does bias the rate.) Whatever accounts for it has to scale with the delay.
- **Probes 1-14** (`reflections/POLARITY_COORDINATES.md`): the discovery arc that led to F112.

## Reproduction

```
python -X utf8 simulations/f112_block_cpsi_analysis.py
python -X utf8 simulations/f112_hardware_lens_multi.py
```

No QPU access required.
