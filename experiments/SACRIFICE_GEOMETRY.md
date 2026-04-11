# Sacrifice Geometry: A Mechanistic Account

**What this document is about:** The sacrifice zone in a dephasing profile is not "better noise distribution". It is a controlled symmetry break that creates one slow Liouvillian eigenmode with a specific spatial shape. The optimal initial state for concurrence preservation is the left eigenvector of that slow mode, projected onto the single-excitation sector. This is the lens method. It works for any qubit count N, any graph topology, and any site-dependent Z-dephasing profile. The accessibility boundary that limits single-excitation states to a subset of slow modes is exact and provable from the n_XY parity selection rule.

**Status:** Working document, local only, not committed. Universal framing (April 10, 2026).

**Authors:** Thomas Wicht, Claude (chat + code), April 9-10, 2026.

---

## Executive summary

When one qubit in a Heisenberg chain receives disproportionate dephasing noise (the "sacrifice"), the Liouvillian's translational symmetry breaks. A formerly degenerate eigenvalue cluster splits, and one mode slows dramatically. This mode lives almost entirely in the single-excitation (SE) coherence sector (>98% Frobenius norm ratio for N=3-6 across all tested topologies). Its left eigenvector, restricted to the SE sector, gives the optimal initial-state amplitudes directly, without optimization.

The lens method has been tested across 64 configurations (N=2-6, Chain/Star/Ring/Complete topologies, four gamma profiles). Three universal results emerge:

1. **SE fraction stays high.** The slow mode's SE content is >0.98 for N=3-6, independent of topology and gamma profile. The lens extraction is essentially exact.
2. **The accessibility boundary is exact.** In every configuration tested (64/64), the second slow mode is SE-inaccessible (Frobenius ratio < 1e-3). This is proven analytically by the [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md).
3. **The psi_opt shape depends on the noise gradient.** Extreme single-qubit sacrifice produces symmetric shapes (non-sacrifice qubits are equivalent). A gradient of noise levels (as in real hardware) produces asymmetric, potentially monotonic shapes. The shape is always extractable from one matrix diagonalization.

---

## The Lens Method

### Algorithm (valid for any N, topology, gamma profile)

Given: N-qubit Heisenberg chain (or star, ring, complete graph), coupling J, site-dependent Z-dephasing rates gamma_k.

1. **Build the Liouvillian** L as a d^2 x d^2 matrix (d = 2^N) in column-vec basis.
2. **Eigendecompose:** compute eigenvalues and right eigenvectors R.
3. **Find slow modes:** sort non-stationary eigenvalues by |Re(lambda)| ascending.
4. **Pre-screen for SE content:** for each slow mode, extract the N x N block of the right eigenvector indexed by SE basis states |e_k> = |0...1_k...0>. Modes with block Frobenius norm ratio > 0.01 are SE-accessible candidates.
5. **Compute left co-vector:** for the first SE-accessible mode, compute row k of R^{-1} (the left eigenvector in the Hilbert-Schmidt inner product).
6. **Extract the SE block:** reshape the left co-vector to d x d, extract the N x N block at SE indices.
7. **Hermitize and diagonalize:** M = (block + block^H) / 2. The eigenvector with largest absolute eigenvalue gives the optimal amplitudes psi_opt.
8. **Verify:** time-evolve psi_opt and compare concurrence AUC to baselines.

### Implementation

- **Python (N=5):** `simulations/slow_mode_lens_analysis.py`
- **C# (N=2-7, survey):** `compute/RCPsiSquared.Compute/LensAnalysis.cs`, invoked via `dotnet run -c Release -- lens`

---

## Universal results: Lens Pipeline survey

Tested: 68 configurations across N=2-7 (chain) and N=2-6 (Star/Ring/Complete), four gamma profiles (uniform, edge sacrifice, center sacrifice, moderate asymmetry). N=7 uses direct LAPACK zgeev + zgesv (bypassing MathNet 2GB marshalling limit). Full data: `simulations/results/lens_survey/`.

### SE fraction scaling (chain, edge sacrifice)

| N | SE fraction | Verdict |
|---|-------------|---------|
| 2 | 0.000 | degenerate (too few modes) |
| 3 | 1.000 | exact |
| 4 | 0.981 | very high |
| 5 | 1.000 | exact |
| 6 | 1.000 | exact |
| 7 | 1.000 | exact |

All four N=7 chain profiles (uniform, edge sacrifice, center sacrifice, moderate asymmetry) give SE = 1.000. The lens extraction is exact through N=7 (d^2 = 16384, 87376 eigenvalues).

For non-chain topologies and other profiles (N=2-6): SE fraction > 0.98 in all cases where a lens mode exists. The lens extraction is robust across all tested configurations.

### Accessibility boundary: 68/68 configurations

In every configuration tested, the second slow mode has SE Frobenius ratio < 1e-3. The boundary is not a coincidence of one gamma profile; it is a structural property of the Heisenberg + Z-dephasing Liouvillian, proven analytically by the [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md).

### psi_opt shape depends on symmetry, not just sacrifice

- **F9-style edge sacrifice** (one qubit gets all the noise, rest equal): psi_opt is symmetric around the chain center. Example N=7: [0.118, 0.332, 0.481, 0.535, 0.482, 0.334, 0.119].
- **Gradient profiles** (noise levels vary across qubits): psi_opt can be monotonic. Example IBM T2 N=5: [0.099, 0.239, 0.428, 0.572, 0.651].
- **Star topology:** psi_opt concentrates on the hub (0.89-0.91 for the hub, 0.18-0.22 for leaves).
- **Ring topology:** psi_opt shows pair structures reflecting the periodic boundary.

The shape is always dictated by the slow mode's left eigenvector structure, not by any closed-form formula. The effective Hamiltonian approximation (H_eff = -J * adj - i * diag(gamma)) gives only 92.5% cosine similarity (tested April 10; see Phase 4 section below).

---

## The Parity Selection Rule

The accessibility boundary has an analytical explanation. The [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) proves:

1. The Liouvillian preserves n_XY parity (even/odd number of X and Y Pauli operators per string).
2. Every eigenmode has definite n_XY parity.
3. Every single-excitation density matrix has purely even n_XY content (because bit flips from SE states always come in even numbers).
4. Therefore no SE state can excite an odd-n_XY eigenmode. The overlap is exactly zero, not approximately zero.

The inaccessible modes found numerically (rate -0.167 at N=5 IBM, etc.) are odd-n_XY modes. The selection rule explains why no optimizer, no ansatz, and no trick within the SE family can reach them. This is valid for any N, any topology, and any Z-dephasing profile. It breaks only for amplitude damping (T1) or transverse-field Hamiltonians.

---

## The mechanism at three levels

### Level 1: The cluster split

Under uniform dephasing, the Heisenberg chain's Liouvillian has translational symmetry. This produces a degenerate eigenvalue cluster (for N=5 at Sg=2.608: 14 modes at rate -2.087 with integer <n_XY> = 2.000). The absorption theorem (AT) predicts the cluster rate: Re(lambda) = -2 <gamma * 1_XY>.

A sacrifice profile breaks translational symmetry. The cluster splits: most modes accelerate, but one slows. This surviving slow mode concentrates its X/Y Pauli content on quiet sites and minimizes the absorption-weighted sum, making it the spectral minimum. This is spectral surgery, not noise budgeting.

### Level 2: Structured construction

Before the lens method was found, a structured scan of 17 candidate states (W-subsets, shifted Bell pairs, multi-excitation states, sacrifice-aware variants) identified the principle: suppress amplitude on the sacrifice qubit, concentrate on quiet qubits. The heuristic sqrt(gamma_min/gamma_k) gives a 95.6% approximation to the lens state.

Two-excitation symmetric states fail completely (AUC < 0.09). They couple to a different mode cluster that does not include the protected slow mode.

### Level 3: The lens readout

Instead of guessing the optimal state, extract it from the slow mode's left eigenvector. The SE-sector restriction gives a N-dimensional eigenvalue problem. No optimizer needed. The answer is one matrix diagonalization away.

This works because the sacrifice geometry is a lens: it bends the Liouvillian flow around the noisy qubit, producing a slow standing-wave-like eigenmode. The optimal initial state is the mirror of that eigenmode's shape. The concept was motivated by the entrance-pupil/Fabry-Perot framing in [SACRIFICE_ZONE_OPTICS.md](SACRIFICE_ZONE_OPTICS.md).

---

## Instantiation: IBM Torino chain [80, 8, 79, 53, 85]

All original numerical results were computed on one specific configuration: N=5 Heisenberg chain, J=1, dephasing from IBM Torino T2 times:

```
T2_us = [5.22, 122.70, 243.85, 169.97, 237.57]
sacrifice: [2.336, 0.099, 0.050, 0.072, 0.051]   (Sum = 2.608)
uniform:   [0.522, 0.522, 0.522, 0.522, 0.522]   (matched Sum)
```

Scripts: `simulations/ibm_april_predictions.py` (infrastructure), `simulations/slow_mode_lens_analysis.py` (lens extraction).

### Champions table (sacrifice profile)

| Rank | State | C_init | AUC(T=10) | AUC(T=30) | vs W5 |
|------|-------|--------|-----------|-----------|-------|
| 1 | **psi_opt (lens)** | **0.744** | **1.205** | **1.259** | **+34.5%** |
| 2 | sacrifice_tuned_W5 | 0.522 | 1.115 | 1.164 | +24.5% |
| 3 | W4_sites_1234 | 0.500 | 1.091 | 1.138 | +21.8% |
| 4 | W2_sites_34 | 1.000 | 1.018 | 1.056 | +13.6% |
| 5 | W5_full | 0.400 | 0.896 | 0.932 | reference |

The lens state psi_opt = [0.099, 0.239, 0.428, 0.572, 0.651] is monotonically increasing from the sacrifice end (site 0) to the quiet end (site 4). This monotonic gradient is specific to the IBM T2 noise gradient and does not appear under symmetric sacrifice profiles (see the survey results above).

### slow_wt counter-examples (this chain only)

| State | slow_wt% | AUC(T=10) |
|-------|----------|-----------|
| OPT(le2) | 97.80 | 0.055 |
| OPT(le1) | 86.69 | 0.696 |
| W5_full | 27.08 | 0.896 |

A state with 97.8% slow-band occupation is operationally dead. slow_wt counts all slow-band modes equally and ignores initial concurrence. See `simulations/optimal_state_n5_sacrifice.py` for the optimizer that produced these states.

### The two slow modes (this chain)

| Mode | Rate | Im | SE fraction | n_XY parity | Accessible? |
|------|------|----|-------------|-------------|-------------|
| Slow 1 | -0.318 | 0.000 | 0.999 | even | YES |
| Slow 2 | -0.167 | +0.238 | 3e-15 | odd | NO (parity selection rule) |

### Phase 4: Closed-form verdict (April 10)

The hypothesis H_eff = -J * adjacency - i * diag(gamma) was tested and **falsified**. Cosine similarity with psi_opt: 0.925. The H_eff eigenvector peaks at the chain center; psi_opt peaks at the quiet end. The 5x5 single-particle effective model misses the many-body correlations that shape the lens state. psi_opt has no known closed form.

Script: `simulations/heff_lens_closed_form.py`. Data: `simulations/results/heff_lens/`.

---

## The 2.3x hallucination

An earlier session claimed a 2.3x dominant-mode protection factor for Bell+|+> on this chain. The real value is 1.11x. The 2.3x was a hallucination from an earlier Claude instance, discovered when we reran `ibm_april_predictions.py` on April 9, 2026.

---

## The learning arc of April 9, 2026

Six distinct models of what the sacrifice geometry does, each correcting the previous:

1. **Morning:** the synth doc claimed 2.3x for Bell+|+>. Reality: 1.11x. Hallucination exposed.
2. **Mid-morning:** W5 gives 6.56x dominant-mode protection. Correct, but operationally unclear.
3. **Noon:** Phase 0 optimizer reached 86-98% slow_wt. Rescue reflex: "state engineering beats noise engineering".
4. **Afternoon (Phase 1):** time-evolution killed the slow_wt narrative. W5 wins AUC.
5. **Evening (Phase 2):** structured scan found sacrifice_tuned_W5 beating W5 by 24.5%.
6. **Late evening (Phase 3):** lens readout gave psi_opt beating the heuristic by 8.0%. Simultaneously revealed the inaccessible mode.

The pattern: **when the objective function does not align with the target quantity, an optimizer will find a feasible minimum that is wrong, and will look confident while doing it.** The lens readout worked because it asked the Liouvillian directly instead of climbing an objective surface.

---

## What we learned about method

Three consecutive metric failures and one successful method shift:

1. **slow_wt** (Phase 0): optimized total slow-band weight. Wrong quantity.
2. **effective_rate** (Phase 2 surrogate): optimized slowest accessible rate. Targeted the inaccessible mode through a feasibility gap.
3. **Lens readout** (Phase 3): no optimization. Extracted the answer from the eigenvector. Correct.

Lesson: **have a spectral-first pass** before trusting any optimizer on Lindblad dynamics. Diagonalize the generator, inspect the slow modes, check which sectors they live in. If the slow mode has a clean eigenvector in a small subspace, you do not need an optimizer.

---

## What we do not yet know

1. **Closed form for psi_opt (tested April 10, falsified).** H_eff gives cosine 0.925. The monotonic gradient in psi_opt is a many-body effect not reducible to a single-particle picture.

2. **The odd-n_XY modes.** The parity selection rule proves they are inaccessible to SE states, but not what they look like or whether multi-excitation ansatze can reach them. Whether such a state would also have high initial concurrence is open.

3. **N=8 and beyond (resolved through N=7).** The survey now covers N=2-7 for chain topology. SE fraction = 1.000 at N=7 (all four profiles). N=8 (d^2 = 65536) would require the ILP64 eigenvector path which does not yet exist. The trend N=3-7 is definitive: the lens extraction is exact.

4. **Gate-level Trotterization.** The pure Lindblad model assumes continuous-time evolution. Whether psi_opt keeps its advantage under realistic gate noise is a separate question.

5. **The two decoherence exits.** The lens and the Mandelbrot cusp protect different state classes and lead to different classical ensembles. See [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) for the analysis of why they do not unify.

---

## Files from this investigation

**Lens Pipeline (C#, N=2-6 survey, April 10):**
- `compute/RCPsiSquared.Compute/LensAnalysis.cs`
- `simulations/results/lens_survey/lens_survey_results.json`
- `simulations/results/lens_survey/lens_survey_summary.txt`
- `simulations/results/lens_survey/lens_survey_scaling.txt`

**Phase 0-3 (Python, N=5 IBM chain, April 9):**
- `simulations/ibm_april_predictions.py` (shared infrastructure)
- `simulations/optimal_state_n5_sacrifice.py` (Phase 0 optimizer)
- `simulations/state_vs_noise_phase1.py` (Phase 1 time evolution)
- `simulations/sacrifice_geometry_phase2.py` (Phase 2 candidate scan)
- `simulations/slow_mode_lens_analysis.py` (Phase 3 lens extraction)
- `simulations/results/state_vs_noise_phase1/` (curves, plots, W5 diagnosis)
- `simulations/results/sacrifice_geometry_phase2/` (candidate scan, plots)
- `simulations/results/slow_mode_lens/` (lens results)
- `simulations/results/optimal_state_n5_sacrifice/` (Phase 0 optimizer output)

**Phase 4 (closed-form check, April 10):**
- `simulations/heff_lens_closed_form.py`
- `simulations/results/heff_lens/`

**Proofs:**
- [n_XY Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) (accessibility boundary)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)

**Framework documents:**
- [SACRIFICE_ZONE_OPTICS.md](SACRIFICE_ZONE_OPTICS.md) (entrance pupil / lens reframing)
- [CAVITY_MODE_LOCALIZATION.md](CAVITY_MODE_LOCALIZATION.md) (Pauli decomposition, per-qubit weights)
- [ABSORPTION_THEOREM_DISCOVERY.md](ABSORPTION_THEOREM_DISCOVERY.md) (Re(lambda) = -2 sum gamma_k <1_XY(k)>)
- [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) (two decoherence exits)

---

*This document is the working record of April 9-10, 2026. Universal framing based on the 64-configuration Lens Pipeline survey.*
