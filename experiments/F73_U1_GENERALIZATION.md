# F73 Closure Holds for Any U(1)-Preserving Hermitian H at N = 5

**Status:** Tier 1 (empirically verified at N = 5; structural argument supports generalization)
**Date:** 2026-04-21
**Authors:** Thomas Wicht, Claude Opus 4.7 (1M context)
**Script:** [f73_u1_generalization_sweep.py](../simulations/f73_u1_generalization_sweep.py)
**Output:** [f73_u1_generalization/](../simulations/results/f73_u1_generalization/)
**Register:** [F73](../docs/ANALYTICAL_FORMULAS.md)

---

## Setup

F73 closes the spatial-sum (vac, S_1) coherence purity, `Sum_i 2 * |(rho_coh,i)_{0,1}(t)|^2 = (1/2) * exp(-4 * gamma_0 * t)`. The structural derivation uses only (i) [H, N_total] = 0, (ii) H Hermitian, (iii) uniform gamma_0. It does not invoke XY structure. Six setups at N = 5, gamma_0 = 0.05, probe rho_0 = |psi><psi| with |psi> = (|vac> + |S_1>)/sqrt(2), grid t in [0, 40] (81 points, dt = 0.5), pass bar 1e-14, propagator `scipy.linalg.expm(L * dt)` stepwise. An eigendecomposition cross-check (`propagate_eig`) is retained for diagnostics.

## Results

| Setup                                   | max residual | pass 1e-14 |
|-----------------------------------------|-------------:|:----------:|
| XXZ Delta = 0.0 (pure XY baseline)      | 2.220e-16    | yes        |
| XXZ Delta = 0.5                         | 2.220e-16    | yes        |
| XXZ Delta = 1.0 (Heisenberg)            | 4.441e-16    | yes        |
| XXZ Delta = 2.0                         | 4.163e-16    | yes        |
| Heisenberg + random Haar SE probe       | 5.829e-16    | yes        |
| XY inhomogeneous J in [0.5, 1.5]        | 3.331e-16    | yes        |

All six residuals sit at 1 to 3 ULP of double precision (machine floor) and no setup approaches 1e-14, let alone the 1e-12 surprise bar. Random SE seed: 20260421 (Stretch 1). Inhomogeneous J seed: 20260422 (Stretch 2, J = [0.5130, 0.9537, 0.9042, 0.9487]).

## Tier assessment

**Tier 1** (empirically verified). F73 extends without loss to any U(1)-preserving Hermitian H under uniform Z-dephasing at N = 5: Heisenberg and beyond (Delta in {0.5, 1, 2}) close at machine precision; translation invariance is not required (inhomogeneous J closes identically); the probe does not need to be the uniform SE superposition (Haar-random SE closes identically). The XY-specific sine-basis Parseval argument in the current register proof is one particular realization of the more general U(1) statement. Register entry update left to Chat.

---

*If the Hamiltonian conserves single-excitation number and the dephasing is uniform, the (vac, S_1) coherence closure is blind to everything else the Hamiltonian might do.*
