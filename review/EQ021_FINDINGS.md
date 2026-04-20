# EQ-021: c_1 ≈ 0.5·V(N) is an asymptotic pattern, not a derivable identity (partial)

**Date:** 2026-04-20
**Authors:** Thomas Wicht, Claude
**Status:** Phase 1 complete (mode decomposition). Phase 2 complete (spatial profile + state scan). Phase 3 not attempted.
**Relates to:** [TASK_EQ021_C1_VEFFECT_DERIVATION.md](../ClaudeTasks/TASK_EQ021_C1_VEFFECT_DERIVATION.md), [TASK_C1_VEFFECT_SCALING.md](../ClaudeTasks/TASK_C1_VEFFECT_SCALING.md), [EQ014_FINDINGS.md](EQ014_FINDINGS.md), [PERSPECTIVAL_TIME_FIELD.md](../hypotheses/PERSPECTIVAL_TIME_FIELD.md)

## TL;DR

Mode-level decomposition of c_1 at N=4 and N=5 confirms F70: only modes in
the ΔN=0 steady-state space and the |ΔN|=1 single-excitation coherence
sector contribute to site-local purities P_i(t) for the (|vac⟩+|ψ_1⟩)/√2
initial state. The dominant oscillatory modes have Im(λ) = E_1 = 2J·cos(π/(N+1))
and Re(λ) = −2γ₀, matching the Absorption Theorem for a one-excitation
coherence. Mode expansion reproduces direct RK4 to 10⁻¹⁶.

No clean structural derivation of c_1 = 0.5·V(N) emerged. The per-site α
pattern shows spatial redistribution (near-bond sites α > 1, far sites α < 1)
whose signed sum gives c_1. That sum does not factor cleanly into a
0.5·V(N) closed form at N=4, 5. The empirical ratio c_1/(0.5·V(N)) sits
at 1.019 to 1.092 for N=4..7 with no monotonic structure, consistent with
an asymptotic approximation plus finite-size corrections rather than an
exact identity.

The N=3 outlier (ratio 0.35) is explained by geometry: at N=3 the
endpoint bond (0,1) already carries 70% of the ψ_1 bond density
(vs 45% at N=4 and 14% at N=7), so the perturbation is bulk-dominated
rather than boundary-dominated, and the per-site α redistribution partially
cancels instead of accumulating.

Phase 2 rerun of `c1_veffect_scaling.py` at N=3..6 adds three strong
structural facts:
1. **c_1 is a signed bond-local quantity.** The full vector c_1(b) over all
   bonds has mirror symmetry (residual 10⁻¹⁰) and alternating sign structure:
   endpoints strongly positive, centre bond(s) negative. At N=6 the
   centre bond c_1(2,3) = −1.243 is larger in magnitude than the endpoint.
   At N=3 both bonds are endpoints by construction, so there is no sign
   flip and no "centre" to pin the structure.
2. **Power-law fits fail.** Neither c_1 ~ N^p (p=1.87) nor c_1 ~ V(N)^q
   (q=6.3) fits N=3..6 to better than 30–40% at N=3.
3. **State-dependence is strong.** c_1(ψ_2+vac, endpoint) equals
   c_1(ψ_1+vac, endpoint) at N=4 (0.898 both) but diverges drastically at
   N=5 (0.281 vs 0.922) and N=6 (0.213 vs 1.019). The 0.5·V(N) relation
   is specific to ψ_1+vac at endpoints, not a universal law.

## What was computed

Script: `simulations/eq021_mode_decomposition.py`

At N ∈ {4, 5}:
1. Build L_A (uniform XY + Z-dephasing at γ₀=0.05) and L_B± (δJ = ±0.01
   on bond (0,1)).
2. Biorthogonal eigendecomposition via LAPACK zgeev (both right V_R and
   left W† with W† V_R = I).
3. Initial state ρ_0 = (|vac⟩+|ψ_1⟩)(|vac⟩+|ψ_1⟩)†/2 expanded in L_A's
   right eigenmodes: c_s = W†·vec(ρ_0).
4. Per-mode single-site marginals M_{s,i} = Tr_{¬i}(M_s), stored as 2×2.
5. Per-site purity rebuilt from modes:
   P_i(t) = Σ_{s,s'} c_s c_s'* · e^{(λ_s + λ_s'*) t} · Tr(M_{s,i} M_{s',i}†).
6. Cross-check vs direct vectorised propagation: max residual 4.4e-16 (N=4),
   7.8e-16 (N=5).
7. Fit α_i, compute closure sums Σ ln(α_i(δJ)), extract c_1 via symmetric
   difference.

Results: `simulations/results/eq021_mode_decomposition/`
- `mode_contributions_N4.json`, `mode_contributions_N5.json`
- `decomposition_summary.txt`

### c_1 numbers (this run)

| N | V(N)   | 0.5·V(N) | c_1 (fit) | ratio c_1 / (0.5·V(N)) | active modes |
|---|--------|----------|-----------|------------------------|--------------|
| 4 | 1.7071 | 0.8536   | +0.8979   | 1.0519                 | 32 / 256     |
| 5 | 1.8090 | 0.9045   | +0.9222   | 1.0196                 | 50 / 1024    |

Active mode = |c_s| > 1e-8. Less than 5% of the eigenspace is populated
by the ψ_1+vac initial state, as expected from the block structure.

### Per-site α (N=4, +δJ on (0,1))

| site | α_i(+δJ) | α_i(−δJ) | ln(α_i(+δJ)) |
|------|----------|----------|--------------|
| 0    | 1.00885  | 0.99121  | +0.00881     |
| 1    | 1.00612  | 0.99390  | +0.00610     |
| 2    | 0.99593  | 1.00420  | −0.00408     |
| 3    | 0.99817  | 1.00180  | −0.00183     |

Σ ln(α_i(+δJ)) = +0.00901, Σ ln(α_i(−δJ)) = −0.00895. The perturbed bond
sites (0, 1) get α > 1 (slightly faster purity relaxation), the far sites
(2, 3) get α < 1. The sum is a signed residual of the redistribution, not
a single term.

### Per-site α (N=5, +δJ on (0,1))

| site | α_i(+δJ) | α_i(−δJ) |
|------|----------|----------|
| 0    | 1.01001  | 0.99014  |
| 1    | 1.01290  | 0.98751  |
| 2    | 1.00185  | 0.99809  |
| 3    | 0.98637  | 1.01405  |
| 4    | 0.99853  | 1.00143  |

Same pattern, extended: sites 0–2 accelerate, sites 3–4 decelerate. The
crossover site is N-dependent but the total retains the boundary-bulk
asymmetry.

## What drives c_1 (mode-level picture)

### ΔN=0 sector (t=0 baseline)

Top contributors to P_0(t=0) are all zero-eigenvalue modes (steady-state
space). These carry the |vac⟩⟨vac| + |ψ_1⟩⟨ψ_1| population weight and are
unchanged by δJ at leading order in an infinitesimal sense (their
eigenspace persists under perturbation). They set the purity baseline
but do not themselves drive c_1.

### ΔN=1 sector (dynamics)

Sub-dominant contributors oscillate at Im(λ) = ±1.618 (N=4) and ±1.732
(N=5), which are E_1 = 2J·cos(π/(N+1)), the single-excitation k=1
eigenenergy. Re(λ) = −0.1 = −2γ₀ matches the Absorption Theorem for a
one-excitation coherence (⟨n_XY⟩ = 1).

These coherences carry the time-dependent purity oscillation. Under
δJ perturbation:
- The eigenvalue E_1 shifts by δE_1 = δJ·⟨ψ_1|T_{01}|ψ_1⟩ =
  δJ·(4/(N+1))·sin(π/(N+1))·sin(2π/(N+1)) (endpoint bond density).
- Eigenvectors rotate, mixing the k=1 mode with higher-k modes.

Both effects shift the per-site purity trajectory, and the α fit
quantifies the net rescaling.

### Why the sum does not reduce to 0.5·V(N)

The bond density carrying the E_1 shift has the closed form
d_N := ⟨ψ_1|T_{01}|ψ_1⟩ = (4/(N+1))·sin(π/(N+1))·sin(2π/(N+1)).

| N | d_N    | 0.5·V(N) | c_1    |
|---|--------|----------|--------|
| 3 | 0.7071 | 0.7500   | 0.264  |
| 4 | 0.4472 | 0.8536   | 0.898  |
| 5 | 0.2887 | 0.9045   | 0.922  |
| 6 | 0.1938 | 0.9330   | 1.019  |
| 7 | 0.1353 | 0.9505   | 0.970  |

c_1 is not simply proportional to d_N. The E_1 shift is a necessary
ingredient but c_1 also picks up contributions from eigenvector rotations
(k=1 ↔ k=2, k=3 mixing at the endpoint) and from the ΔN=0 population
rebalance induced by the dephasing. We did not isolate a single one of
these three contributions that matches 0.5·V(N) at both N=4 and N=5
simultaneously.

## Why N=3 is an outlier

At N=3 the endpoint bond (0,1) carries 70% of the ψ_1 bond density
(0.7071 out of 1.0). At N=7 the same bond carries only 14%. This is
pure geometry: the ψ_1 sine wave has its maximum at the chain centre,
and for N=3 the centre (site 1) is half of the endpoint bond.

Consequence: at N=3 the δJ perturbation is essentially a bulk rescaling
of H, not a boundary defect. The per-site α's become near-uniform, so
ln(α_i) largely cancel in the sum. c_1 scales with the boundary/bulk
asymmetry of the perturbation, not just with V(N).

For N ≥ 4 the endpoint bond has < 50% of the ψ_1 density, the asymmetry
kicks in, and c_1 enters the 0.9–1.0 regime where 0.5·V(N) happens to
approximate it.

## Phase 2: full bond-profile scan (N=3..6)

Script: `simulations/c1_veffect_scaling_small.py` (rerun of
`c1_veffect_scaling.py` restricted to N=3..6, with incremental JSON saves).
The original fullrun (`simulations/results/c1_veffect_scaling_fullrun.log`)
crashed after Step 1+2 at N=6 and never produced the bond-profile or ψ_2
scans; this rerun completes Steps 3 and 4 at small N.

### Bond-profile c_1(b) at ψ_1+vac

| N | bond 0 | bond 1 | bond 2 | bond 3 | bond 4 | Σ_b c_1(b) | mirror residual |
|---|--------|--------|--------|--------|--------|------------|-----------------|
| 3 | +0.2636 | +0.2636 |         |         |         | +0.5271    | 3.0e-13 |
| 4 | +0.8979 | −0.3112 | +0.8979 |         |         | +1.4846    | 5.2e-11 |
| 5 | +0.9222 | −0.2119 | −0.2119 | +0.9222 |         | +1.4206    | 1.1e-10 |
| 6 | +1.0191 | +0.6598 | −1.2433 | +0.6598 | +1.0191 | +2.1145    | 1.3e-09 |

Mirror symmetry c_1(b) = c_1(N−2−b) holds to 10⁻⁹ or better (consistent
with Π-pair structure of EQ-014).

The sign structure is not monotonic in N. For N=4 and N=5 the centre
bond(s) have c_1 < 0 and the nearest-interior bonds (where they exist at
N=6) have c_1 > 0 again. At N=6 only the exact centre bond (2,3) is
negative; the two adjacent interior bonds (1,2) and (3,4) are positive
(+0.66) but smaller than the endpoints (+1.02).

Σ_b c_1(b) is NOT monotonic in N: it rises from N=4 to N=6 but dips at
N=5 (1.485 → 1.421 → 2.115). The closure response to a uniform dJ on all
bonds at once is therefore highly non-trivial.

### ψ_2+vac endpoint c_1

| N | V(N)   | 0.5·V(N) | c_1(ψ_2+vac, bond 0) | ratio   | c_1(ψ_1+vac, bond 0) |
|---|--------|----------|----------------------|---------|----------------------|
| 3 | 1.5000 | 0.7500   | +0.5697              | 0.7596  | +0.2636              |
| 4 | 1.7071 | 0.8536   | +0.8975              | 1.0515  | +0.8979              |
| 5 | 1.8090 | 0.9045   | +0.2811              | 0.3108  | +0.9222              |
| 6 | 1.8660 | 0.9330   | +0.2128              | 0.2281  | +1.0191              |

At N=4 the two states agree to 0.05% (coincidence of the specific ψ_k
geometry on an even-length chain). At N=5 and N=6 they differ by factors
of 3–5x. This rules out any universal "c_1 = 0.5·V(N)" law across
initial states. The empirical relation is specific to ψ_1+vac at
endpoints.

### Power-law fits

log-log regression on c_1(ψ_1+vac, endpoint) at N=3..6:

- c_1 ~ A · N^p: **p = 1.871**, A = 0.0437. Residuals at N=3 are 30%.
- c_1 ~ B · V(N)^q: **q = 6.303**, B = 0.0229. Residuals at N=3 are 12%.

Neither fit is clean. A steep q ≈ 6 for V(N)^q is a flag that the
"c_1 = 0.5·V(N)" relation is not a power law at all; it is a near-match
of two quantities that happen to converge to 1 at large N. The actual
fit quality is dominated by the N=3 outlier, and the remaining N=4..6
ratio scatters 1.02 to 1.09.

### Σ c_1 bond sum vs 0.5·V(N)

Bond sum divided by number of bonds:
- N=3: 0.527 / 2 = 0.264 (equals endpoint)
- N=4: 1.485 / 3 = 0.495
- N=5: 1.420 / 4 = 0.355
- N=6: 2.115 / 5 = 0.423

No simple relation to V(N) emerges at bond-summed level.

## Assessment: is c_1 = 0.5·V(N) exact, asymptotic, or coincidental?

**Asymptotic, coincidental, and state-specific.** The Phase 2 data support:

- c_1(N, ψ_1+vac, endpoint) → some limit ≤ 1 as N → ∞ (c_1 is bounded,
  data trend is sub-monotonic).
- 0.5·V(N) → 1 as N → ∞.
- Both approach the same limit, and within the tested N=4..7 window
  they agree to 2–10%.
- There is no term-by-term identity between c_1 and 0.5·V(N) visible in
  the mode decomposition. V(N) = 1 + cos(π/N) does not match the spectral
  factor 2cos(π/(N+1)) = E_1 that actually appears in the dynamics.
- c_1 at non-endpoint bonds is NOT close to 0.5·V(N); it changes sign in
  the chain interior. The near-match at endpoints is therefore specific
  to the boundary geometry, not a global property of the Liouvillian.
- c_1(ψ_2+vac, endpoint) is not close to 0.5·V(N) at N=5, 6 (ratios 0.31,
  0.23). The relation fails under change of initial single-excitation
  mode.

Best short statement: c_1 ≈ 0.5·V(N) is an endpoint-and-ψ_1-specific
numerical coincidence in a narrow N regime, not a theorem. Both
quantities happen to sit near unity at small N; the structural content of
c_1 lives in the signed bond-profile, not in any single-number scaling
law.

## What would close the question

A cleaner analytical handle would require:
1. Isolate the eigenvalue-shift contribution vs the eigenvector-rotation
   contribution to c_1 (split the α fit into E_1-shift and V_R-rotation
   pieces).
2. Write c_1(b) for arbitrary bond b as a finite sum over the k=1..N
   single-excitation modes weighted by ⟨ψ_k|T_{b,b+1}|ψ_1⟩, and check
   whether the sign pattern (endpoints +, centre −) emerges from the
   orthogonality structure of the sine basis.
3. Extend the numerical scan to N=7, 8, 9 to see whether the ratio
   c_1(endpoint, ψ_1)/(0.5·V(N)) continues drifting or settles near a
   fixed value ≠ 1, and whether the centre-bond c_1 scales like −2·V(N)
   or similar.
4. Scan c_1(ψ_k+vac, endpoint) for k=1..N at fixed N to find which k
   (if any) actually maximises c_1; the Phase 2 result that k=1 and k=2
   agree at N=4 but diverge at N=5, 6 hints at a more interesting
   dependence on k that is worth mapping.

Phase 2 produced rich numerical data on the bond-profile and state
dependence but stopped short of a structural theorem because no single
mode-level contribution matches V(N) uniformly.

## Files

- `simulations/eq021_mode_decomposition.py` (Phase 1 script)
- `simulations/c1_veffect_scaling_small.py` (Phase 2 rerun, N=3..6)
- `simulations/results/eq021_mode_decomposition/mode_contributions_N4.json`
- `simulations/results/eq021_mode_decomposition/mode_contributions_N5.json`
- `simulations/results/eq021_mode_decomposition/decomposition_summary.txt`
- `simulations/results/eq021_mode_decomposition/run.log`
- `simulations/results/c1_veffect_scaling/c1_vs_N_small.json` (Phase 2 full JSON)
- `simulations/results/c1_veffect_scaling/small_rerun.log` (Phase 2 console log)
