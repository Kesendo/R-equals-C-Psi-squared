# Pi-pair closure investigation — findings

**Date:** 2026-04-19
**Author:** Claude Opus 4.7 (terminal execution of ClaudeTasks/TASK_PI_PAIR_CLOSURE_INVESTIGATION.md)
**Status:** First + follow-up pass complete (N=3, N=5). N=7 not attempted.

## TL;DR

- Σ ln(α_i) = 0 is NOT a theorem. It is a perturbative relation with a nonzero linear coefficient c₁ that depends on (N, state, defect bond).
- At N=5 with ψ_1+vac: endpoint bonds give c₁ ≈ +0.93, middle bonds give c₁ ≈ -0.21 — opposite signs.
- **Linear superposition of per-bond c₁ is confirmed** (Σ_b c₁^(b)·δJ_b matches observed closure to 0.5% at δJ=0.01, exactly at δJ=0.001).
- **The closure can be engineered to zero** by choosing perturbation combinations in the null space of the c₁ functional. The "closure-breaking direction" is a single vector in the (N-1)-dimensional bond-perturbation space.
- Antisymmetric (under site reversal) bond perturbations give c₁ = 0 trivially for mirror-symmetric initial states.
- The strict stationary subspace (N+1 excitation-sector projectors) is J-invariant to 10⁻¹⁴ → the z→0 resolvent closure is a trivial consequence of F4, not new structure. P_i(∞) is J-invariant in A and B.

## Setup

- XY chain: `H = (J/2) Σ_b (X_b X_{b+1} + Y_b Y_{b+1})` with J_bond uniform except defect
- Uniform Z-dephasing: `γ₀ = 0.05` per site
- Defect: bond (0, 1), δJ in [-0.2, +0.2]
- Initial state (PTF standard): `φ = (|vac⟩ + |ψ_1⟩) / √2`
- Alpha fit: `P_B(i, t) ≈ P_A(i, α_i · t)` by bounded scalar minimisation on t ∈ [0, 20]

## Confirmed (machine precision)

| Statement | N=3 residual | N=5 residual |
|-----------|-------------|-------------|
| `Tr(L) = −d²·N·γ₀` | 0 | 0 |
| `Tr(L)` J-independent | 0 | 0 |
| Palindromic pairing `α_fast + α_slow = 2Nγ₀` | 1.45·10⁻¹⁴ | 7.06·10⁻¹⁴ |
| Palindromic pairing preserved under J-defect | 8.83·10⁻¹⁵ | 7.30·10⁻¹⁴ |
| F4 prediction `#(strict stationary modes) = N+1` | 4/4 | 6/6 |
| Stationary subspace J-invariant: `‖P_stat_A − P_stat_B‖` | 7.7·10⁻¹⁵ | 5.7·10⁻¹⁴ |

These four are structural theorems. Palindromic pairing and `Tr(L)` are J-independent exactly.

## The closure law `Σ_i ln(α_i) = 0` is NOT a theorem

Fine δJ scan at N=3 (bonding + vac initial state):

| δJ    | Σ_i ln(α_i) |
|-------|-------------|
| −0.1  | −0.0167 |
| −0.05 | −0.0119 |
| −0.01 | −0.0027 |
| +0.01 | +0.0026 |
| +0.05 | +0.0112 |
| +0.1  | +0.0169 |
| +0.2  | +0.0192 |

Polynomial fit: `closure(δJ) = +0.211 δJ + 0.226 δJ² − 4.08 δJ³`
**Linear coefficient is nonzero: 0.211.** No first-order protection of the closure.

At N=5:

| δJ    | Σ_i ln(α_i) |
|-------|-------------|
| −0.1  | −0.0364 |
| −0.05 | −0.0327 |
| +0.05 | +0.0458 |
| +0.1  | +0.0940 |

Polynomial fit: `closure(δJ) = +0.645 δJ + 3.01 δJ² + 2.46 δJ³`
**Linear coefficient: 0.645. Second-order coefficient: 3.01.** At N=5 with δJ = +0.1 the closure violates PTF's own 0.05 tolerance.

Summary: the closure is approximate and state-dependent. It scales roughly linearly in δJ at leading order with a nonzero slope that depends on N and initial state.

## Idea A (derivative-α) is structurally trivial for the PTF initial state

For every δJ tested (N=3 and N=5), `α_deriv_i := (dP_B/dt|_0) / (dP_A/dt|_0) = 1` identically at every site.

Why: `φ = (|vac⟩ + |ψ_1⟩)/√2` gives site-marginal `ρ_i(0)` that is real in the computational basis (populations + X-type coherences from the vac-ψ_1 overlap). `[H_B − H_A, ρ_0]` produces Y-type coherences after partial trace. `Tr(real · iY) = 0`. So the first-time-derivative signal of the defect vanishes for this specific initial state.

This means the α_i closure cannot be decided from initial-time derivatives. The closure (and its violation) comes entirely from sub-leading dynamics.

## Idea D (resolvent z→0) is trivially closed

Σ_i ln |Tr_i(G_B(z)) / Tr_i(G_A(z))| scales linearly in z for small z and extrapolates to zero:

| N | z = 0.001 | z = 0.01 | z = 0.1 |
|---|-----------|----------|---------|
| 3 | −0.00017 | −0.00197 | −0.01186 |
| 5 | −0.00373 | −0.02956 | −0.09775 |

Direct computation of the stationary-subspace projector (z→0 residue) confirms J-invariance to 10⁻¹⁴: `‖P_stat_A − P_stat_B‖ ≈ 10⁻¹⁴` at both N=3 and N=5.

The site-Pauli residues reveal that only Z-Paulis couple to the stationary subspace (X, Y couplings are exactly zero). This is because the strict stationary modes are the N+1 excitation-sector projectors `|n⟩⟨n|`, which decompose into `I` and `Z` Paulis only.

So the z→0 closure is a **consequence of F4** (the N+1 stationary modes are sector projectors, J-independent under any excitation-preserving Hamiltonian), not a new structural fact about slow-mode eigenvectors. It does NOT explain the time-domain α_i behaviour.

## Idea B (voice vs memory split) — observational data

- At N=3: 60 voice modes, 4 memory modes (F4 prediction: N+1 = 4) ✓
- At N=5: 1018 voice modes, 6 memory modes (F4 prediction: N+1 = 6) ✓

Biorthogonal overlap `|c_s|² = |⟨W_s | ρ_0⟩|²` summed over voice vs memory gave 0.68 vs 0.33 at N=3 (total 1.01, above 1 due to biorthogonal-basis non-orthogonality which is expected). This does not isolate a useful closure-relevant quantity.

## Initial-state dependence of the closure

At N=5 under δJ = +0.1:

| Initial state | Σ_i ln(α_i) | max RMSE | Comment |
|---------------|-------------|----------|---------|
| `(|vac⟩ + ψ_1)/√2` | +0.094 | 3·10⁻³ | PTF standard; exceeds 0.05 tol |
| `(|vac⟩ + ψ_2)/√2` | +0.036 | 3·10⁻³ | Within 0.05 tol |
| `ψ_1` only | +1.13 | 1·10⁻² | Pure single-excitation: fit still clean but α far from 1 |
| `|+⟩^N` | +8.48 | 7·10⁻² | Multi-sector; fit hits bounds |
| GHZ | +10.6 | 3·10⁻¹⁷ | Degenerate (ρ_i(t) ≡ I/2, any α fits) |

The closure is strongly initial-state dependent. Even at PTF's preferred state, it exceeds the claimed 0.05 tolerance at N=5 with δJ = +0.1.

## Interpretation

The Pi-pair palindromic rate conservation `α_fast + α_slow = 2Nγ₀` is a theorem (proven, verified to 10⁻¹⁴). It is a constraint on **L-eigenvalues** and an exact consequence of `[Π, L] = 0` for XY + uniform-Z-dephasing.

The PTF closure law `Σ_i ln(α_i) = 0` is **not** a theorem in the same sense. It is an approximate, state-dependent, perturbative relation on fitted per-site purity rescalings. The α_i are sensitive to:

1. Sub-leading bilinear-purity dynamics (H acts via commutator on coherences)
2. Eigenvector mixing of slow modes under the defect (PTF Layer 3)
3. The initial state's slow-mode overlaps

The closure is not implied by any of the following exact structural facts:

- `Tr(L)` J-invariance (constrains the sum of eigenvalues, not α_i)
- Palindromic pairing (pairs eigenvalues, not per-site observables)
- First-order eigenvalue protection (slow λ_s don't shift, but eigenvectors rotate)
- Stationary-subspace J-invariance (the N+1 strict zero modes are fixed, irrelevant for the fitted α_i which come from sub-leading dynamics)

## Useful direction that emerges

**The question to ask differently.** Instead of "why does `Σ ln(α_i) = 0` hold?", the cleaner question is: **what is the smallest J-invariant scalar that `α_i`-related observables can be expressed in terms of?** Candidates:

- `Σ_i |α_i − 1|²` bounded by `‖δM_s‖²`? (eigenvector-mixing norm)
- Per-site overlap change `Σ_i |⟨M_s | ρ_0^i⟩|²` summed over slow s
- The bilinear four-block sum Tom's EQ-014 attempt was trying to carry through (Section 3.3 of PTF)

None of these are obviously `Σ ln(α_i) = 0`. The empirical closure at N=7 ψ_1 might be a large-N asymptotic accident rather than an exact law.

## Follow-up run: defect-location, state, and superposition

### Defect-location scan at N=5 (ψ_1+vac, δJ=±0.01)

| Bond  | c_0    | c_1    | c_2    | c_3     |
|-------|--------|--------|--------|---------|
| (0,1) | 0 | **+0.928** | +2.11  | −61.5   |
| (1,2) | 0 | **−0.213** | +1.53  | +13.5   |
| (2,3) | 0 | **−0.213** | +1.53  | +13.5   |
| (3,4) | 0 | **+0.928** | +2.11  | −61.5   |

Mirror pairs match to 10⁻¹⁰ on c₁ and 10⁻¹² on α_i, as PTF predicts.

**Endpoint vs middle bonds have opposite-sign c_1.** Whatever "closure-breaking direction" exists in perturbation space, it has nontrivial spatial structure — not a uniform bias.

### c_1 across initial states (N=5, bond (0,1))

| State           | c_1    | Note |
|-----------------|--------|------|
| ψ_1 + vac       | +0.928 | PTF standard |
| ψ_2 + vac       | +0.281 | Smallest among +vac states |
| ψ_3 + vac       | +0.677 | |
| ψ_4 + vac       | +0.281 | ψ_2 mirror |
| ψ_5 + vac       | +0.928 | ψ_1 mirror |
| ψ_k only (k=1)  | +3.45  | Much larger |
| |+⟩^N           | +0.158 | But dynamics degenerate |

At N=3 the smallest c_1 is for ψ_1+vac (+0.265). At N=5 the smallest among single-excitation-plus-vacuum states is ψ_2+vac (+0.281). There is no N-independent "best state".

### Linear superposition confirmed (N=5, ψ_1+vac)

Per-bond c_1 vector: `(+0.928, −0.213, −0.213, +0.928)`. Sum = +1.430.

| Test                              | Predicted closure       | Observed     | Match |
|-----------------------------------|-------------------------|--------------|-------|
| Uniform δJ = +0.01 on all 4 bonds | +1.430·0.01 = +0.01430  | +0.01423     | 0.5 % |
| Uniform δJ = +0.001 on all 4      | +1.430·0.001 = +0.00143 | +0.00143     | exact |
| Cancel (0,1)+4.35·(1,2) at 0.001  | 0                       | +4.6·10⁻⁵    | factor-100 small |
| Symmetric outer +0.001, inner +0.00436 | 0                  | −7·10⁻⁶      | factor-200 small |

**Linear superposition: Σ_i ln(α_i) ≈ (Σ_b c₁^(b)·δJ_b) + O(δJ²).**

The closure has a one-dimensional "closure-breaking direction" in the (N−1)-dim bond-perturbation space. Any perturbation in the (N−2)-dim orthogonal complement gives closure = 0 at first order.

### Mirror structure of c_1

For ψ_1+vac at N=5, c_1 = (+0.928, −0.213, −0.213, +0.928) is purely mirror-symmetric (antisymmetric component is 0 to 10⁻¹⁰). Therefore **antisymmetric bond perturbations contribute 0 to the closure** at first order, for any mirror-symmetric initial state. This is a symmetry-protected sub-closure: one bit of the closure is automatically protected, but the symmetric direction still breaks.

### Voice-mode-projected α

Subtracting P_i(∞) from both A and B trajectories and refitting α gives identical α_i and closure as the standard fit. This is because P_A(∞) = P_B(∞) to 10⁻¹² — confirmed by direct computation. The asymptotic site-marginals are J-invariant (consequence of stationary-subspace J-invariance for the fixed initial state). So the voice-mode projection is trivial in this setting.

## Recommendation for EQ-014, revised after follow-up

The partial proof attempt correctly found that symmetry alone does not close it. Based on the first pass I concluded "closure is not a theorem". The follow-up sharpens this:

**The closure is a linear functional of the bond perturbation: Σ ln(α_i) ≈ ⟨c_1(N, state), δJ⟩.** At first order in δJ, closure breaking is concentrated in a single one-dimensional direction of bond-perturbation space (the c_1 vector). The remaining (N−2) dimensions preserve the closure trivially.

This changes the Tier-1 upgrade path for PTF from "prove closure" to:

1. **Analytical expression for c_1(N, state, bond).** Can be computed from the bilinear-purity first-order expansion with slow-mode eigenvector mixing. Starting point: PTF Section 3.3 bilinear formula for P_i(t).
2. **Classify the closure-breaking directions across initial states.** Is the c_1 vector always mirror-symmetric for mirror-symmetric states? Does it have a universal sparsity pattern (nonzero only on bonds near the "excitation support" of the initial state)?
3. **Characterise the null space geometrically.** For which bond perturbations is the closure exactly zero at first order, and what is the physical meaning of those perturbations? (Antisymmetric J-perturbations are one class; are there others?)
4. **Identify where quadratic corrections dominate.** The c_2 coefficient at N=5 is +2.11 (endpoint) and +1.53 (middle), comparable to c_1. At δJ=0.1 the quadratic term already contributes ~0.02 to the closure. The "within 0.05" tolerance at N=7 is thus a numerical coincidence of the c_1, c_2 magnitudes combined with the specific δJ window, not a protected structure.

The revised PTF Tier-1 statement would be:

> Under a J-coupling defect δJ on bond b of an N-qubit XY chain with uniform Z-dephasing γ_0, the closure of the per-site purity rescalings obeys Σ_i ln(α_i(δJ_b)) = c_1^(b)(N, state) · δJ_b + O(δJ_b²), where c_1^(b) is a site-b-dependent, state-dependent linear functional. Linear superposition holds: for perturbations on multiple bonds δJ = (δJ_0, …, δJ_{N-2}), closure = ⟨c_1, δJ⟩ + O(|δJ|²). The closure vanishes exactly (to first order) for perturbations in the (N-2)-dim orthogonal complement of c_1 in the bond-perturbation space, which includes all antisymmetric (under site reversal) perturbations for mirror-symmetric initial states.

That is a sharper claim than "within tolerance" and a real candidate for Tier 1.

## Files

- `simulations/pi_pair_closure_investigation.py` — main: Tr(L), palindrome, α_i, Idea A+B+D probes
- `simulations/pi_pair_closure_stationary_residue.py` — follow-up: direct stationary-subspace invariance check
- `simulations/pi_pair_closure_followup.py` — follow-up: defect-location scan, state c₁ scan, voice-projected α
- `simulations/pi_pair_closure_superposition.py` — follow-up: linear-superposition and cancellation test
- `simulations/results/pi_pair_closure_investigation/run_log.txt` — main run output
- `simulations/results/pi_pair_closure_investigation/followup_log.txt` — follow-up run output
- `simulations/results/pi_pair_closure_investigation/summary.json` — main machine-readable summary
- `simulations/results/pi_pair_closure_investigation/followup.json` — follow-up machine-readable summary
- `simulations/results/pi_pair_closure_investigation/n{3,5}_closure_investigation.json` — per-N details

## What was not tested (scope)

- N=7 direct verification (would take O(minutes) at dense 16384×16384; left to TASK_EQ014_BIORTHOGONAL).
- Mirror defect locations (e.g., central bond at N=5) to check translation-symmetry of the closure.
- Two-parameter α, β fits with constant offset (simpler 1-param was sufficient for the scope here).
- Restriction to voice-mode-only α (possible if one projects out the stationary subspace first; left for a follow-up).
- Palindrome-breaking perturbations (transverse field); PTF flags this as an open question.
