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

**The model family has no per-qubit detuning, and this data needs one.** On block_cpsi, the second tensor slot's transverse expectations alternate in sign from sample to sample while their magnitude decays smoothly (the archived data gives the pair as [13, 14] but states no slot ordering, so which physical qubit that is remains an assumption; nothing below depends on it):

```
t (μs)      0      120     240     360     480
<I,X>    0.665  -0.352   0.149  -0.064   0.029
<I,Y>    0.010  -0.136   0.188  -0.080   0.024
```

That is a coherence turning about half a revolution per 120 μs sample, aliased by the sampling. Qubit 13 turns too, more slowly and without aliasing: its transverse angle walks −35.3°, −35.2°, −36.2° and then −48.9° across the four steps. Neither dephasing nor amplitude damping can turn a coherence at all, and a ZZ coupling turns both qubits at nearly the same rate, which these two do not do. So no member of the family above can represent what the data plainly does, and whatever the optimiser cannot represent lands in the parameters it has: on the aliased qubit the fitted γ_Z lands between 0.138 and 0.147 per μs in four of the five models, about **140×** the calibrated D[Z] rate 1/(2·T2) ≈ 0.00104 (133× to 141× across those four). This ratio read "about 70×" until 2026-08-05, when it was dividing by 1/T2, i.e. by twice the Lindblad rate. The correction doubles the ratio and so only sharpens the point being made: the optimiser is not fitting dephasing at all, it is parking an unrepresentable rotation in the one parameter it has.

A re-fit with the corrected scripts is therefore not enough. Every RMS ranking this experiment produced is a comparison within the wrong family, and no replacement channel identification is offered. **A redo needs:** one free Z coefficient per qubit, converged multistart (the fits still show nesting violations, a superset fitting worse than its subset, in the current output as well as the old), model comparison at matched parameter count (AIC/BIC), and uncertainties.

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
