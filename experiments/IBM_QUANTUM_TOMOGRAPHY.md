# IBM Quantum Tomography: First Hardware Test of R = CΨ²

## Discovery Date
2026-02-09

## Summary

First empirical test of the C·Ψ = ¼ boundary on real quantum hardware. State tomography
on IBM Torino (Heron processor, 133 qubits) confirms that a single qubit's C·Ψ product
crosses the ¼ boundary during decoherence.

Two critical findings emerged:
1. **The ¼ boundary crossing is real.** Measured on a physical qubit, not just simulated.
2. **T₂ ≠ T₂*.** IBM calibration reports T₂ from Hahn echo. Our experiment is a free
   induction decay, governed by T₂* ≈ T₂/2.7. The correct comparison timescale matters.

Additionally, the generalized crossing equation was derived to account for T₁ relaxation,
extending the pure dephasing cubic x³ + x = ½ to arbitrary T₁/T₂ ratios.

## Hardware

| Parameter | Value |
|-----------|-------|
| Backend | ibm_torino (Heron r2, 133 qubits) |
| Qubit | #52 (selected for best T₂) |
| T₁ (calibration) | 221.2 μs |
| T₂ (Hahn echo, calibration) | 298.2 μs |
| T₂* (free induction decay, measured) | 110.1 μs |
| Shots per circuit | 8,192 |
| Delay points | 25 (0 to 895 μs) |
| Total runtime | ~8 minutes |

## Experiment Design

1. **Prepare** |+⟩ = (|0⟩+|1⟩)/√2 using a Hadamard gate
2. **Wait** for variable delay time t (25 logarithmically spaced points, 0 to 3×T₂)
3. **Tomograph** full state reconstruction via Qiskit `StateTomography` (X, Y, Z basis, maximum-likelihood estimation)
4. **Extract** density matrix ρ(t) at each delay
5. **Compute** C(t) = Tr(ρ²), Ψ(t) = 2|ρ₀₁|, product C·Ψ

This is a free induction decay (FID): no refocusing pulses, no dynamical decoupling.
The qubit simply decoheres under ambient noise while we watch.

## Results

### Raw Density Matrix Evolution

![IBM Torino full analysis](../visualizations/ibm_tomography/ibm_torino_full_analysis.png)

*Six-panel analysis of IBM Torino qubit 52. Top-left: density matrix elements showing
T₁ relaxation (ρ₁₁ → 0) and T₂* dephasing (|ρ₀₁| → 0). Top-center: coherence decay
on log scale, showing T₂* = 110 μs vs calibration T₂ = 298 μs. Top-right: C·Ψ trajectory
normalized to T₂*. Bottom-left: zoomed crossing region. Bottom-center: tomography fidelity.
Bottom-right: result summary.*

### The ¼ Crossing

| Quantity | Value |
|----------|-------|
| C·Ψ at t = 0 | 0.885 (ideal: 1.000) |
| Crossing time t* | 114.7 μs |
| t*/T₂* (measured) | 1.041 |
| t*/T₂* (generalized prediction, r = 0.456) | 0.936 |
| t*/T₂* (pure dephasing prediction) | 0.858 |
| Deviation from generalized | 11.3% |
| Asymptotic purity C∞ | 0.740 (ideal: 0.500) |

### Key Observations

**The crossing exists.** C·Ψ starts at 0.885 (imperfect |+⟩ preparation) and falls
monotonically through ¼, reaching thermal equilibrium below 0.05. The transition is smooth,
as Lindblad theory predicts.

**The initial state is imperfect.** Gate fidelity 0.97 means C₀ = 0.942, Ψ₀ = 0.940.
The product C·Ψ₀ = 0.885 instead of the ideal 0.500. This shifts the crossing time
but does not eliminate it.

**The asymptotic state is not maximally mixed.** Purity stabilizes at ~0.74, not 0.50.
This is caused by readout errors and thermal population imbalance. Qubit 52 relaxes toward
|0⟩ (ground state), creating a bias that inflates the measured purity at long times.

**Coherence decays 2.7× faster than T₂ suggests.** The calibration T₂ = 298 μs comes
from Hahn echo, which refocuses low-frequency (1/f) noise. Our FID experiment sees all
noise sources. The measured T₂* = 110 μs, with the difference explained by pure dephasing:
1/T₂* = 1/T₂ + 1/T_φ, giving T_φ ≈ 175 μs.

## The T₂ vs T₂* Distinction

This is the most important experimental lesson:

```
T₂ (Hahn echo):  Refocuses slow noise. Measures only fast dephasing.
T₂* (FID):       Sees ALL noise. The relevant timescale for unprotected qubits.

Relationship:  1/T₂* = 1/T₂ + 1/T_φ

For ibm_torino qubit 52:
  T₂  = 298 μs  (calibration, echo)
  T₂* = 110 μs  (measured, FID)
  T_φ = 175 μs  (derived: low-frequency noise contribution)
```

When comparing the C·Ψ = ¼ crossing time with theory, **always use T₂***. The calibration
T₂ is the wrong timescale for free decoherence experiments.

This distinction is well-known in NMR and quantum computing, but easy to miss when pulling
calibration data from IBM's backend properties.

## The Generalized Crossing Equation

The original cubic x³ + x = ½ assumes pure dephasing (T₁ → ∞). Real hardware has finite T₁.

For a single qubit |+⟩ under combined T₁ + T₂* decay:

```
ρ₀₀(t) = 1 - ½ e^{-t/T₁}
ρ₁₁(t) = ½ e^{-t/T₁}
|ρ₀₁(t)| = ½ e^{-t/T₂*}
```

With b = e^{-t/T₂*} and r = T₂*/T₁:

```
C = 1 - b^r + b^{2r}/2 + b²/2
Ψ = b

Crossing equation:  [1 - b^r + b^{2r}/2 + b²/2] · b = ¼
```

![Generalized crossing](../visualizations/ibm_tomography/generalized_crossing.png)

*Left: crossing time t*/T₂ as function of r = T₂/T₁. Center: C·Ψ trajectories for
different r values. Right: crossing points trace the C·Ψ = ¼ hyperbola in C-Ψ space.*

### Special cases

| r = T₂*/T₁ | Equation | t*/T₂* | Physical regime |
|-------------|----------|--------|-----------------|
| r → 0 | b³ + b = ½ | 0.858 | Pure dephasing (T₁ ≫ T₂) |
| r = 0.456 | Numerical | 0.936 | IBM Torino qubit 52 |
| r = 0.5 | Mixed | 0.950 | T₁ = 2T₂ |
| r = 1 | 4b³-4b²+4b = 1 | 1.141 | T₁ = T₂ |

### Physical interpretation

T₁ relaxation drives population toward |0⟩, creating asymmetry (ρ₀₀ > ρ₁₁). An asymmetric
state has higher purity than a symmetric one (|0⟩ has C = 1). This counteracts the purity
loss from dephasing, delaying the crossing. Larger r means stronger T₁ effect, later crossing.

The ¼ boundary itself is invariant. Only the arrival time depends on the decoherence channel.

### Polynomial approximation (max error < 0.001)

```
t*(r) ≈ 0.858 + 0.012r + 0.375r² − 0.019r³ − 0.084r⁴
```

## Error Budget

| Source | Effect | Magnitude |
|--------|--------|-----------|
| Gate errors | C₀ < 1, Ψ₀ < 1 | ~6% initial state error |
| Readout errors | Purity floor at ~0.74 | Shifts asymptotic C·Ψ |
| Non-exponential decay | 1/f noise → Gaussian envelope | ~5% on T₂* fit |
| Tomography overhead | Extra gates for X/Y/Z measurement | Additional decoherence during measurement |
| Delay quantization | dt rounding to hardware clock | < 0.1%, negligible |

The 11.3% deviation between measured (1.041) and predicted (0.936) crossing times is
consistent with the combined error budget.

## What This Proves

1. **The ¼ boundary crossing is a real physical phenomenon**, not an artifact of simulation.
2. **The generalized crossing equation extends the framework** to arbitrary decoherence channels.
3. **T₂* is the operationally relevant timescale** for free decoherence.
4. **The framework survives contact with real hardware**, with deviations attributable to
   known error sources.

## What This Does Not Prove

1. The 11% deviation has not been fully explained. It could indicate model limitations.
2. A single qubit on one backend is not statistical evidence. Multiple qubits/backends needed.
3. The imperfect initial state complicates direct comparison with the ideal-state prediction.
4. We have not tested the two-qubit case (entanglement dynamics) on hardware.

## Relation to Universal Quantum Lifetime

This experiment extends [UNIVERSAL_QUANTUM_LIFETIME.md](UNIVERSAL_QUANTUM_LIFETIME.md):

- The cubic x³ + x = ½ is the **r → 0 limit** of the generalized equation.
- The published T₁/T₂ data analysis remains valid for platforms where T₁ ≫ T₂ (trapped ions, NV centers).
- For superconducting qubits (r ≈ 0.3–0.8), the generalized equation gives the correct prediction.
- All previous simulations used pure dephasing and remain valid within that model.

## Next Step: Final Hardware Run (March 2026)

IBM Quantum free tier gives 10 min/month. After that: $96/min. This means one more
run and then real hardware experiments end unless funding changes.

The final run targets three questions with 75 circuit batches in ~8 minutes:

1. **Repeat crossing measurement** with points clustered near t* (12-point FID tomo)
2. **Independent T₂* via Ramsey** (15 points, single-basis, 3× cheaper than tomo)
3. **Second qubit with different r = T₂/T₁** (8-point crossing-focused tomo)

Two data points on the r-curve. One independent T₂* confirmation. No budget for
error mitigation, Hahn echo, or two-qubit entanglement.

---

*Back to [experiments overview](README.md) | Related: [Universal Quantum Lifetime](UNIVERSAL_QUANTUM_LIFETIME.md)*
