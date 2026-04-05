# Trajectory Dwell Time at the CΨ = 1/4 Crossing

**Status:** Verified (analytical + numerical, April 5, 2026)
**Source script:** [critical_slowing_trajectory_dwell.py](../simulations/critical_slowing_trajectory_dwell.py)
**Main experiment:** [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md)

---

## Result

For Bell+ under Z-dephasing (Formula 25), the dwell time near the
CΨ = 1/4 crossing is:

    t_dwell(δ) = 2δ / |dCΨ/dt|_{t_cross}

where

    dCΨ/dt = -2γ · f · (1 + 3f²) / 3,    f_cross = 0.8612241

This gives:

    t_dwell = 2δ / (1.851701 · γ) = 1.080088 · δ / γ

## K-invariance

In rescaled time K = γt, the dwell K-interval is:

    K_dwell = γ · t_dwell = 1.080088 · δ

This is independent of γ. Numerically verified across γ = 0.1 to 10.0
with standard deviation < 2 × 10⁻¹⁷ (machine precision).

## Crossing point

The crossing f*(1 + f*²) = 3/2 is solved by f_cross = 0.8612240997,
giving K_cross = γ · t_cross = -ln(f_cross)/4 = 0.03735.

## Physical interpretation

The dwell time measures how long the system remains near the saddle-node
bifurcation during decoherence. It scales as:

- **δ** (linearity): symmetric about the crossing to leading order
- **1/γ** (K-invariance): faster dephasing compresses the time window
  but the K-measure stays constant
- **|dCΨ/dt|⁻¹**: the derivative at the crossing point fully determines
  the dwell behavior; higher-order corrections appear only at δ ≥ 10⁻²

This answers [Open Question #2](../docs/WEAKNESSES_OPEN_QUESTIONS.md)
("Crossing speed dependence"): yes, d(CΨ)/dt at the crossing moment
fully determines the post-crossing convergence timescale via the dwell
time formula.

## Numerical verification

| γ   | t_dwell (δ=10⁻³) | K_dwell      | Ratio to prediction |
|-----|-------------------|--------------|---------------------|
| 0.1 | 0.01080097        | 0.00108010   | 1.0000              |
| 0.5 | 0.00216019        | 0.00108010   | 1.0000              |
| 1.0 | 0.00108010        | 0.00108010   | 1.0000              |
| 2.0 | 0.00054005        | 0.00108010   | 1.0000              |
| 10  | 0.00010801        | 0.00108010   | 1.0000              |

| δ      | t_dwell/t_pred | Accuracy |
|--------|----------------|----------|
| 10⁻²   | 1.0008         | 0.08%    |
| 10⁻³   | 1.0000         | <0.01%   |
| 10⁻⁴   | 1.0000         | <0.01%   |
