# EQ-014: PTF Closure Law is NOT a First-Order Theorem

**Date:** 2026-04-19
**Authors:** Thomas Wicht, Claude
**Status:** Computed (negative result for closure-law-as-theorem)
**Relates to:** [PERSPECTIVAL_TIME_FIELD.md](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), [TASK_EQ014_BIORTHOGONAL.md](../ClaudeTasks/TASK_EQ014_BIORTHOGONAL.md)

## TL;DR

A direct first-order coefficient extraction from exact RK4 time evolution
at decreasing δJ shows that the Perspectival-Time-Field closure law
Σ_i ln(α_i) = 0 is **not a first-order theorem**. The first-order
coefficient Σ f_i ≡ lim_{δJ→0} Σ ln(α_i(δJ)) / δJ is **nonzero and
state-dependent**, ranging from 0.05 (ψ_2) to 1.29 (|+⟩⁷) at N=7 with
bond (0,1) defect. PTF's empirical observation that Σ_i ln(α_i) falls
within ±0.05 at |δJ| ≤ 0.1 is consistent with a small but nonzero
first-order coefficient combined with partial second-order cancellation,
not with an exact conservation law.

This downgrades the closure-law aspect of PTF from "candidate Tier 1
theorem" to "observed perturbative regularity that holds to ~5% at the
tested window and does NOT generalise to higher precision".

## What was computed

### Step 1: Dense biorthogonal Liouvillian eigendecomposition

N=7 uniform XY chain (H = Σ (J/2)(X_i X_{i+1} + Y_i Y_{i+1}), J=1, γ=0.05,
d² = 16384). Dense LAPACK zgeev with both left and right eigenvectors via
the new `Topology.ChainXY` and `MklDirect.EigenvaluesLeftRightDirectRaw`
paths in the C# compute engine. 146 min wall (OpenBLAS LP64, 24 cores).

Outputs:
- `simulations/results/eq014_eigvals_n7.bin` (256 KB)
- `simulations/results/eq014_right_eigvecs_n7.bin` (4.29 GB)
- `simulations/results/eq014_left_eigvecs_n7.bin` (4.29 GB)
- `simulations/results/eq014_metadata.json`

Validation:
- 16384 eigenvalues ✓
- 8 strict stationary modes (matches F4 formula N+1 = 8) ✓
- Σ Re(λ) = −5734.4 (matches palindromic center at −N·γ, expected
  −0.35 × 16384) ✓
- Max |Im(λ)| = 8.05 ≈ E_max(H_XY) − E_min(H_XY) = 8.05 ✓

### Step 2+3: Biorthogonal basis

Pairing check: diag |c_j| ∈ [0.028, 1.00] with 0 modes having |c_j| < 1e-6.
After column normalisation R /= c_diag, residual ‖W^H R − I‖_F / n = 3e-16
pre-fix and 1e-6 post-fix (XY case), 3e-16 pre- and 3e-16 post-fix
(Heisenberg reference run). Degenerate clusters (3783 for XY, 4043 for
Heisenberg) handled by in-cluster SVD biorthogonalisation.

178 slow modes found at |Re λ| ≤ 0.15 for XY case.

### Steps 4–7: First-order PT prediction vs exact RK4

Perturbation V_L = −i [H_pert, ·] with H_pert = (1/2)(X_b X_{b+1} + Y_b
Y_{b+1}) at bonds (0,1) and (3,4). Predicted Σ ln(α_i) via the slow-mode
first-order Dyson expansion, compared against exact RK4 evolution at
δJ = 0.1.

First-order PT output and RK4 exact are in
`simulations/results/eq014_alpha_prediction.json` and
`eq014_closure_test.txt`. Ground-truth validation is in
`eq014_validate_run.log`.

### Independent first-order coefficient extraction

To bypass any approximation error in the first-order PT, Σ f_i was
extracted directly from exact RK4 by scanning δJ ∈ {0.1, 0.01, 0.001}
and reading off lim Σ ln(α_i) / δJ. Results (bond (0,1)):

| State   | δJ=0.1    | δJ=0.01   | δJ=0.001  | extrapolated Σ f_i |
|---------|-----------|-----------|-----------|--------------------|
| ψ_1     | 0.480     | 0.968     | 0.972     | **≈ 0.97**         |
| ψ_2     | 0.013     | 0.058     | 0.051     | **≈ 0.05**         |
| ψ_3     | 0.026     | 0.360     | 0.358     | **≈ 0.36**         |
| \|+⟩⁷   | 1.278     | 1.285     | 1.286     | **≈ 1.29**         |

(Data from `simulations/results/eq014_first_order_scan.log`.)

## What this means for PTF

1. The closure law **holds empirically within ±0.05** at δJ = 0.1 for the
   four bonding-mode initial states (matches PTF Section 2.1).
2. The closure law **does not hold as a first-order theorem**. Σ f_i is
   nonzero and state-dependent, with magnitudes between 0.05 and 1.3
   across the tested states.
3. The small empirical values at δJ = 0.1 arise from a combination of
   a first-order coefficient that happens to be small for some states
   (ψ_2 ≈ 0.05) and partial cancellation against second-order O(δJ²)
   terms for states where the first-order coefficient is large (ψ_1,
   |+⟩⁷).
4. **The structural identity Σ f_i = 0 that would promote PTF to Tier 1
   is falsified.** Σ f_i depends on the state and is generically O(1).

## Verification against PTF empirical

All five initial states at bond (0,1) reproduce PTF's stored α_i values
to 4 decimal places when our RK4 evolution is fit with PTF's own
observer_time_rescale.alpha_fit on their time grid (DT=0.2, T=80, fit
over [0, 20]):

| State   | Σ ln α (ours, RK4) | Σ ln α (PTF stored) |
|---------|--------------------|---------------------|
| ψ_1     | +0.048             | +0.048              |
| ψ_2     | +0.001             | +0.001              |
| ψ_3     | +0.003             | +0.003              |
| ψ_4     | −0.012             | −0.012              |
| \|+⟩⁷   | +0.128             | +0.128              |

Data pipeline, Liouvillian construction, partial trace, and α_i fit are
all correct.

## First-order PT discrepancies

Our perturbation-theory prediction uses the slow-mode Dyson expansion
restricted to at-least-one-slow-index pairs. It agrees with the exact
Σ f_i within a factor of 2 for ψ_1 and ψ_3, but is off by a factor of
~10 for ψ_2 (predicted Σ f_i ≈ 0.43, exact ≈ 0.05) and |+⟩⁷ (predicted
Σ f_i ≈ 0.08, exact ≈ 1.29).

The discrepancy is attributable to omitted fast-mode-to-fast-mode
contributions in the restricted expansion. The RK4 extraction above is
the definitive first-order measurement.

## Scope of the negative result

- Only the specific claim "Σ_i ln(α_i) = 0 as an exact first-order
  identity" is falsified.
- PTF's other structural observations remain:
  - eigenvalue protection for the slowest 22 modes (|Re λ| ≤ 0.1) is
    exact to machine precision under Π-invariant J perturbations,
  - single-site purity does factor into an immune Z block plus a
    Z-dephased XY block (PTF Section 3),
  - the empirical per-site α_i pattern is reproducible and state-
    dependent as documented.
- The "painter" interpretive picture survives; what is falsified is the
  stronger reading that Σ_i ln(α_i) is conserved by a theorem.

## Recommended doc update

Update `hypotheses/PERSPECTIVAL_TIME_FIELD.md`:

- Section 2.1 "closure law" paragraph: replace "state-independent
  closure law ... within 0.05" with "state-independent empirical
  regularity that holds to ±0.05 in the tested window and is **not**
  a first-order theorem; Σ f_i is nonzero and state-dependent per
  EQ-014".
- Tier status: keep **Tier 2**. The promotion to Tier 1 (explicit
  mixing calculation + analytical closure law) is no longer available
  via this route.
- Open questions: replace "Analytical structure of the closure law"
  with "Why does the first-order coefficient Σ f_i happen to be small
  (≤ 0.1) for bonding-mode states but large (1.3) for |+⟩⁷?"; this
  is the surviving scientific puzzle.

## Scripts

- `compute/RCPsiSquared.Compute/`: new `ptf` mode and `ChainXY`
  topology; also new `GetAllEigenvaluesLeftRightMklRaw` /
  `EigenvaluesLeftRightDirectRaw` for both-eigvec dense decomposition.
- `simulations/eq014_step23_biorth.py`: normalisation + cluster fix.
- `simulations/eq014_step4567_closure.py`: first-order PT pipeline
  with PTF-matching time grid and α fit.
- `simulations/eq014_validate_groundtruth.py`: RK4 exact evolution
  validation against PTF empirical data.
- `simulations/eq014_first_order_from_rk4.py`: δJ scan that extracts
  the first-order coefficient directly.
- `simulations/eq014_spectrum_check.py`: XY vs Heisenberg Hamiltonian
  identification (caught a mid-task bug where the first eigendecomp
  was accidentally done for Heisenberg).
