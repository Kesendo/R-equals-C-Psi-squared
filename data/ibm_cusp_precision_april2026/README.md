<!-- QUARTER-CURRENT -->
# IBM Kingston finite radial-precision run, April 26 2026

Current reading: the immutable rows sample a normalized-purity/coherence proxy
near a scalar radial value. The saved phase is unavailable; the saved payload has no phase observable, so the
run does not establish a real-axis path, zero imposed phase, head-on cusp approach,
or physical recurrence trajectory.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the original cusp-precision description follows with the
measured payload and acquisition caveats intact.

# IBM Kingston radial-scalar precision run, April 26 2026

A curated copy of one Bell⁺-pair run sampled densely around the scalar value
CΨ=1/4. It tests an F25 magnitude fit; it does not measure a recurrence cusp.

## What is in this directory

| File | Description |
|------|-------------|
| `cusp_precision_ibm_kingston_20260426_115939.json` | The 19-delay trajectory: per-delay scalar CΨ magnitude, the F25 fit, and the pointwise residuals (the headline run) |
| `cusp_precision_ibm_kingston_20260426_115053.json` | The same morning's second run, same pair (14,15), 17 delays over the WIDE window 2.26–39.6 μs (job `d7mtuoraq2pc73a24pjg`): covers the April-16 crossing region and finds the coherence already dead there (CΨ ≈ 0.05 at 22.6 μs); its F25 fit γ_fit = 0.015039 is the source of the headline run's γ_calib = 0.015 |
| `cusp_precision_ibm_kingston_20260426_113534.json` | The morning's first run, auto-selected pair (149,150), 15 delays around its echo-predicted 25.0 μs (job `d7mtnmlqrg3c738ln1c0`): the same collapse on a second pair, γ_fit = 0.0269 = 18.0× its echo-γ |

The three runs are one same-day chain: run 1 (pair 149,150, echo-predicted window)
found the coherence dead across its whole sampled window (CΨ ≈ 0.01–0.02
over 12.5–37.6 μs); run 2 (pair 14,15 forced, wide window)
measured the fast decay across 2.3–39.6 μs and fit γ = 15.04/ms; run 3 placed 19
dense delays from that fit (`--gamma-override 15.0`) and confirmed γ_fit = 14.98/ms
(RMS 0.0097). The γ_calib = 0.015 in the headline JSON is therefore a same-day
in-situ fit, not a calibration value; the pair's T2-echo γ that day was 0.00165/μs,
which is the 9.08× gap tracked in EQ-025.

## Why this run is here

This is the experiment with points sampled densely near and straddling the
selected radial line at `1/4`. None of the 19 saved scalar values is exactly
`1/4`. Where the April-16
[cusp-slowing run](../ibm_cusp_slowing_april2026/) used only six delays per pair to bracket the
crossing (but saved full density matrices, so it carries phase), this run sampled the crossing
**densely**: 19 delays at factors 0.05 to 3.0 of the model-selected `t_cross`,
with about eight of them straddling the radial line (factors 0.85 to 1.15,
CΨ approximately 0.22 to 0.26). It is a point-by-point finite F25 fit check,
not an exact-quarter or recurrence-cusp measurement.

**It carries no phase.** The JSON stores only a scalar `cpsi` magnitude per delay,
not the 4×4 density matrices. Therefore phase is unknown: the rows cannot place
the trajectory on any axis, infer Ω=0, or determine a crossing argument. The
phase-carrying April-16 run is a separate dataset; this one is a dense radial
magnitude trace.

## Experiment summary

- **Backend:** ibm_kingston (Heron r2)
- **Date:** 2026-04-26
- **Job ID:** `d7mu36lqrg3c738lnda0`
- **Pair:** qubits 14, 15 (T2_min ≈ 303 μs, γ ≈ 0.00165 1/μs); γ_calib 0.015, γ_fit 0.014977
- **Protocol:** Bell⁺ (H + CX), native delay, 9-Pauli 2-qubit tomography, 2048 shots per basis, 171 circuits total
- **Delays:** 19, at factors [0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.96, 0.99, 1.0, 1.01, 1.04, 1.08, 1.15, 1.3, 1.5, 1.8, 2.3, 3.0] of t_cross
- **Result:** F25 = f·(1+f²)/6, f = exp(−4γt), fits with γ the only free parameter; `fit_rms` = 0.0097; `F_CROSS` = 0.8612, `K_CROSS` = 0.037350 (γ·t at the crossing)

## Structure of the JSON

```
{
  "mode": "hardware", "backend": "ibm_kingston", "job_id": "d7mu36lqrg3c738lnda0",
  "pair": { "qubits": [14, 15], "gamma_per_us": ... },
  "gamma_calib": 0.015, "gamma_fit": 0.014977,
  "delay_factors": [ ... 19 ... ],
  "cpsi_data": [ { "factor": ..., "t_us": ..., "cpsi": ... }, ... 19 ... ],
  "pointwise_residual_calib": [ ... ], "pointwise_residual_fit": [ ... ],
  "fit_rms": 0.0097, "F_CROSS": 0.8612, "K_CROSS": 0.037350
}
```

## Provenance

Curated copy of the full run, which lives external to this repo in the IBM tomography directory
(`ibm_quantum_tomography/run_cusp_precision.py`, results folder). See
[CRITICAL_SLOWING_AT_THE_CUSP.md](../../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) (the April-26
precision update) for the full analysis. The JSON here is the authoritative record of the values
the repo's plots read.
