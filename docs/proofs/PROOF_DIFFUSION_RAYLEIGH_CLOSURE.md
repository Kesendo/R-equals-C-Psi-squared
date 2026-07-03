# The Closure Functional: the survivor's bond rate shift is the squared density gradient

**Status:** Tier 1 candidate (leading-order analytical derivation in the strong-dephasing diffusion limit; confirmed *exact* as Q → 0 by an engine Q-sweep, slope → 2.00 and CV → 0, and dressed by off-diagonal coherence at finite Q). Tier capped at candidate by its felt-time parent, not by the derivation. Resolves the `felt_time_dimensions` arc (D) follow-up. Physics-first-reviewed 2026-06-19; formula/number/provenance pass 2026-06-20 (Q=2.5 row corrected to the N=5 witness output, witness N-reach and the sector-dependent prefactor flagged, n(j) wording and the shape-miss metric clarified).
**Date:** 2026-06-19
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Statement:** for the half-filling survivor (whose slow mode carries a single-site density profile `n(j)`, a signed, zero-mean fluctuation that is antisymmetric under reflection, not an occupation probability), the first-order shift of its decay rate under a single-bond defect on bond `b = (j, j+1)` is `∂(Re λ)/∂J_b ∝ (n(j) − n(j+1))²`, the squared density gradient: ≈ 0 at the no-flux chain ends, maximal in the interior, mirror-symmetric. Exact in the strong-dephasing limit; dressed by off-diagonal coherence through the soft-survivor regime, and superseded at the handover where the survivor becomes the rigid `(0,1)` band edge.
**Typed claim:** [`SurvivorDiffusionGradientClaim.cs`](../../compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs) (Tier 1 candidate); live witness `inspect --root gradient` (drivable across Q, reports the slope, the CV, and the off-diagonal weight).
**Verifier:** [`felt_time_amplitude_law.py`](../../simulations/felt_time_amplitude_law.py) (the block-level law, gate-first N=4..7) and [`felt_time_closure_functional.py`](../../simulations/felt_time_closure_functional.py) (the trajectory ground truth); the Q-regime sweep is the witness itself, `inspect --root gradient --N 5 --q …`.
**Builds on:** [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) (`Re λ = −2γ⟨n_XY⟩`, the rate the survivor decays at); the incompleteness survivor (the slowest `(p,p)` mode, `SurvivalIncompletenessMirrorClaim`).
**Formula registry:** [F123 in `ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).

---

## What this is about

A small chain is being watched, and the watching wears its rhythms away. One pattern outlives all the others: the half-filled one, spread so thin and so even across the chain that the watching can barely read it. That is the survivor, and the question here is a narrow, concrete one. If you reach in and change a single link of the chain, tightening or loosening the coupling between two neighbouring sites, by how much does the survivor's lifetime change?

The answer is not "the same everywhere". It depends entirely on *where* you reach in. The survivor sits in the chain as a single slow swell, a smooth wave of more-here, less-there, rising from one end to the other. Touch the chain where that swell is climbing steeply, in the moving middle, and you change the lifetime the most. Touch it at the quiet ends, where the swell lies flat against the wall and is going nowhere, and you change the lifetime not at all. The size of the effect, link by link, is the square of how steeply the swell was rising right there.

This is what a *diffusion* does. Once the watching is strong, the survivor behaves like heat spreading along a bar: a density smoothing itself out, slowest in its longest, gentlest wave. One subtlety the math insists on, and this document is careful about it: the smooth density itself is invisible to the watching (a pure population carries nothing for the light to read). What the watching actually reads, and wears away, is the faint shimmer of coherence that the smoothing *stirs up*, and that shimmer is strongest exactly where the density is steepest. So the lifetime, and how a single link changes it, is still governed by the gradient of the swell, squared, link by link, the ends insulated so the gradient and the effect vanish there. This document derives that in the strong-watching limit where it is exact, shows the engine confirming it sharpens to exactly that limit as the watching strengthens, and is honest about the corrections at moderate watching.

## Abstract

Under uniform Z-dephasing the longest-lived non-steady Liouvillian mode of an open XY chain is the half-filling **survivor**, the slowest mode of a diagonal `(p,p)` sector (one of a filling-degenerate family; all share the law). It is **predominantly a density mode**: a dominant, real diagonal carrying the single-site occupation profile `n(j)` (a standing density wave), dressed by a subdominant Hamming-distance-2 off-diagonal coherence admixture. The two play asymmetric roles: the HD=0 diagonal is **dark** (zero `⟨n_XY⟩`), and the decay rate `−2γ⟨n_XY⟩` is carried *entirely* by the HD=2 coherence. `Tr(M†H_b) = 0` (exact, all Q) rules out only nearest-neighbour hopping content, not this HD=2 admixture, whose weight `‖M−diag M‖/‖M‖` runs from `~0.45` (N=4) to `~0.24` (N=7) at `Q=1.5` and `→ 0` as `Q = J/γ → 0`. The classical population-diffusion picture is the secular *effective* description of the rate's bond-profile (the coherence adiabatically eliminated, the per-bond rate `D ∝ J_b²/γ` its second-order trace), not a literal claim that the mode is a decaying population.

In that effective theory the slow density `n(j)` obeys diffusion on the sector's hopping graph, and `Re λ` is the (negative) Rayleigh quotient of a real-symmetric weighted graph Laplacian `Re λ = −Σ_b D_b (n(j)−n(j+1))² / ‖n‖²`. A single-bond defect tunes one `D_b`; the Hellmann-Feynman theorem (applied to that real-symmetric reduced Laplacian) gives `∂(Re λ)/∂D_b = −(n(j)−n(j+1))²/‖n‖²`, hence

    ∂(Re λ)/∂J_b ∝ (n(j) − n(j+1))²   (the squared density gradient).

Open boundaries are no-flux (Neumann): the lowest diffusion harmonic `n(j) ∝ cos(π(j−½)/N)` (the continuum/large-N form) has zero gradient at the end bonds, so the shift vanishes there and is maximal in the interior, mirror-symmetric. The harmonic's *shape* is fixed by the chain length (`k_min = π/N`), independent of the diffusion scale `D ∝ J²/γ ∝ Q`, so the bond-profile is Q-invariant in the strong-dephasing regime while its magnitude is not. An engine Q-sweep confirms the derivation is exact in the limit: the log-log slope of `|dRe|` against `|grad|` is `2.00` and `CV → 0` as `Q → 0` for every N, and drifts above 2 (`2.12` at `Q=1.5`, `2.39` at `Q=2.0`) precisely as the off-diagonal coherence dressing grows, until at the handover (`Q* ≈ 2.5`) the survivor switches to the rigid `(0,1)` band edge and the law no longer applies. The earlier single-particle `φ·φ` candidate used the wrong standing wave (the survivor is multi-magnon): right power, wrong wave.

## 1. The survivor: a predominantly-diagonal mode whose rate is carried by its coherence dressing

Take the open XY chain `H = (J/2) Σ_{(i,i+1)} (X_iX_{i+1} + Y_iY_{i+1})` under uniform Z-dephasing at rate `γ`, Liouvillian `L = L_H + L_D`, `L_H = −i[H,·]`, `L_D(·) = γ Σ_l (Z_l · Z_l − ·)`. The Absorption Theorem gives `Re λ = −2γ⟨n_XY⟩`; the global slowest non-steady mode in the strong-dephasing regime lives in a half-filling diagonal sector `(p,p)`, with fractional `⟨n_XY⟩ < 1`. The diagonal `(p,p)` sectors are exactly **filling-degenerate** in the slowest `Re` (the verifiers' own sweep warns the `(p,p)` label "hops"); all of them, including `(1,1)`, share the same single-site profile *shape* `n(j)` and the same *scaling* law, so "the survivor" names that whole tied family, not a unique block. (The absolute normalization of `n(j)` and the bond-independent prefactor *are* sector-dependent: at `N=5` the witness tie-breaks to `(1,1)` and the Python verifier to `(3,3)`, giving the same slope and CV but a rescaled `dRe/grad²` magnitude; §4.)

Its embedded coherence operator `M` satisfies, for every bond,

    Tr(M† H_b) = 0   (exact, all Q),

so `M` carries no nearest-neighbour single-hop content: it does not couple to the bond current, which is what rules out the single-particle `φ·φ` current picture. It does **not** make `M` diagonal. The off-diagonal weight `off := ‖M − diag M‖ / ‖M‖` (witness `OffDiagonalWeight`) is substantial at moderate `Q` and shrinks both with `N` and with strengthening dephasing:

    off ≈ 0.45 (N=4), 0.35 (N=5), 0.28 (N=6), 0.24 (N=7)   at Q=1.5;
    off ≈ 0.04 (Q=0.2),  0.11 (Q=0.5),  0.22 (Q=1.0),  0.49 (Q=2.0)   at fixed N.

This off-diagonal part is **Hamming-distance-2 intra-sector coherence**, and it is **load-bearing for the rate**: `⟨n_XY⟩` computed from the HD=2 entries alone reproduces `−Re λ/2γ` exactly, while the HD=0 diagonal (pure `{I,Z}`) has `⟨n_XY⟩ = 0` and is **dark**. So the survivor's *dominant* component is the density profile `n(j)` (real, the lowest antisymmetric harmonic: `N=4` gives `[−0.81, −0.38, +0.38, +0.81]`, `N=7` gives `[+1.57, +1.31, +0.73, 0, −0.73, −1.31, −1.57]`), but its *decay* runs entirely through the subdominant coherence admixture that this density stirs up. As `Q → 0` the admixture vanishes (`off → 0`) and so does the rate (`⟨n_XY⟩ → 0`): the strict diffusion limit is the dark, density-only mode, and the finite rate at finite `Q` is the size of the coherence dressing. The classical-diffusion treatment of §2 is therefore the leading-order *effective* description of how the rate depends on `n(j)`, not a statement that the mode is a bare population.

## 2. The law: Hellmann-Feynman on the effective diffusion Rayleigh quotient

In the strong-dephasing (secular) limit the off-diagonal coherences are gapped out at rate `O(γ)` and adiabatically eliminated; the slow dynamics is a classical master equation for the diagonal populations, hopping across the sector's configuration graph at per-bond rate `D_b ∝ J_b²/γ` (the Haken-Strobl form: a bond's hopping `J_b` drives incoherent population transfer `∝ J_b²/γ` via the second-order excursion through its dephased coherence; this is the same HD=2 coherence that §1 shows carries the rate, now eliminated into an effective population rate). The conserved single-site density `n(j)` is the hydrodynamic coordinate; its slowest collective mode is the lowest eigenvector of the weighted graph Laplacian `W = Σ_b D_b L_b`, `L_b` the edge Laplacian of bond `b`. `W` is **real-symmetric** (a classical generator), so its eigenvalues are Rayleigh quotients and the slowest non-trivial one is

    −Re λ = min_{n ⊥ 1} (Σ_b D_b (n(j) − n(j+1))²) / ‖n‖² = (Σ_b D_b (Δn_b)²) / ‖n‖²,   Δn_b := n(j) − n(j+1).

The defect tunes a single `D_b`. The Hellmann-Feynman theorem applies to this real-symmetric `W` (the eigenvector `n` is stationary, so only the explicit parameter dependence survives):

    ∂(−Re λ)/∂D_b = (Δn_b)² / ‖n‖².

Since `D_b ∝ J_b²/γ`, the chain-rule factor `∂D_b/∂J_b = 2J_b/γ` is bond-independent at the uniform base point `J_b = J` (the quadratic `J`-dependence collapses to one constant once differentiated there), so

    **∂(Re λ)/∂J_b = −(2J/γ) (Δn_b)² / ‖n‖²  ∝  (Δn_b)²,**

the squared density gradient, with one bond-independent constant `2J/(γ‖n‖²)`. This is the statement. Two scope caveats, both honest: (i) the Hellmann-Feynman step is on the *reduced, Hermitian* Laplacian `W`; the full Liouvillian is non-Hermitian and its biorthogonal first-order shift (the `dRe` the witness actually computes, via the left-right eigenvectors `dRe = Re[(R⁻¹ row)·V_b·(R col)]`) agrees with this Hermitian-Rayleigh derivative in *shape* but not in the prefactor (the constant is off by the biorthogonality factor, the normalized left-right overlap `|⟨l|r⟩|/(‖l‖‖r‖)` running `0.60 → 0.88` across N=4..7 as the mode becomes more normal, measured in `value_vector_felt_time.py`; the witness normalizes `⟨l|r⟩ = 1` and so folds it into `dRe` rather than emitting it). (ii) The verified content is the *scaling* `∝ (Δn_b)²` and the bond-independence, not the exact constant.

## 3. Boundary and Q-invariance

Open boundaries impose **no-flux (Neumann)** conditions on the diffusion: no population leaks past the end sites. The lowest non-trivial Neumann harmonic on `N` sites is, in the continuum/large-N limit,

    n(j) ∝ cos(π (j − ½) / N),   so   Δn_b = n(j) − n(j+1) ∝ sin(π j / N),

which is *minimal* at the end bonds (`∝ sin(π/N)`, vanishing strictly only as `N → ∞`, i.e. in the continuum limit `j/N → 0, 1`), peaks in the interior, and is mirror-symmetric. Hence `∂(Re λ)/∂J_b ∝ sin²(π j / N)`: smallest at the ends, maximal in the middle. This is the reflection's "untouchable where it lies flat against the wall, alive to the touch where it climbs". The discrete many-body profile approaches this form, with the `sin²` shape-miss (the L∞ distance between the unit-normalized bond profile `(Δn_b)²` and the unit-normalized `sin²(π j / N)`) shrinking monotonically from `0.17` at `N=4` through `0.06` at `N=7` to `0.043` / `0.035` at `N=8` / `N=9` (§4), confirming the continuum harmonic as the large-N limit.

The harmonic's wavevector `k_min = π/N` is fixed by the chain length alone; the diffusion scale `D ∝ J²/γ ∝ Q` multiplies `W` overall and so sets the magnitude of `Re λ` but not the *shape* of its eigenvector. So the bond-profile `(Δn_b)²` is Q-invariant **in the strong-dephasing regime** where the classical reduction holds; as `Q` rises the off-diagonal dressing perturbs the bond-independence (the CV grows, §4). This is the origin of the seam observation that the closure is flat in `Q` while `⟨n_XY⟩` runs: in its regime the closure reads the Q-fixed shape, not the Q-dependent scale.

## 4. Verification and honest scope

The reduction of §2 is leading-order in `γ/J`-secular perturbation theory, exact as `Q = J/γ → 0`. The verified landing point `Q = 1.5` is the project's **canonical hardware operating point** (`γ₀ = 0.05`, the hardware dephasing rate; `J = 0.075`; `docs/Q_REGIME_ANCHORS.md`), so the law is checked at the real operating ratio, not a cherry-picked one. That point (`J = 1.5γ > γ`) is **below the handover** (so the survivor is still the soft interior mode) but is **not** in the strict `γ ≫ J` secular limit: the biorthogonality is only `0.60` at N=4 (strongly non-normal) and `⟨n_XY⟩/Q²` has already drifted `~23 %` off its diffusion-limit value. The law's *shape* survives there because the `n(j)` harmonic is Q-fixed (§3); the corrections are real and were measured, not assumed.

**The regime sweep is the decisive test, and it is the engine itself** (`inspect --root gradient --N 5 --q …`, which reports the slope, the CV, and the off-diagonal weight per Q):

| Q = J/γ | sector | off-diag weight | log-log slope (`dRe` vs grad) | CV of `dRe/grad²` | witness verdict |
|---|---|---|---|---|---|
| 0.3 | (p,p) | 0.07 | **2.00** | 0.001 | LAW HOLDS |
| 0.5 | (p,p) | 0.11 | 2.01 | 0.002 | LAW HOLDS |
| 1.0 | (p,p) | 0.22 | 2.04 | 0.011 | LAW HOLDS |
| 1.5 | (p,p) | 0.35 | 2.12 | 0.040 | LAW HOLDS |
| 2.0 | (p,p) | 0.49 | 2.39 | 0.175 | **not clean** |
| 2.5 | (p,p), entering handover | 0.72 | n/a (−0.17) | 0.93 | **not clean** |

So **as `Q → 0` the slope is exactly 2.00, the CV vanishes and the mode becomes a pure density mode (off → 0), for every N**: the diffusion-limit derivation is confirmed. The drift above 2 at finite `Q` is the **off-diagonal coherence dressing** growing with `Q` (a finite-Q, non-secular correction, *not* a boundary effect), and the witness's own `LawHolds` gate flips to "not clean" by `Q = 2.0`. At the handover `Q* ≈ 2.5` (N-dependent, the F122 `2.39–2.61` range) the survivor stops being the soft `(p,p)` mode and becomes the rigid `(0,1)` band edge (`off = 1`, `⟨n_XY⟩ = 1`, the F122 / band-edge regime), where this law does not apply. The table's `Q = 2.5` row is `N = 5` just *entering* that handover (still a `(p,p)` block, `off = 0.72`, the log-log slope already nonsense and the gate not clean), not the fully-developed `(0,1)` edge with `off = 1`.

At a fixed moderate `Q = 1.5` the law still holds across N (the original landing). Note the slope here **drifts away from 2 as N grows** (2.00, 2.12, 2.16, 2.17 for N=4..7); this is the residual finite-Q correction (Q=1.5 is fixed), **not** convergence: it is the *shape-miss* metric that converges with N. The two must not be conflated (and the `N=4` slope `2.00` is the weakest datum: mirror symmetry leaves only two distinct bonds, so a log-log fit through two points of equal `dRe/grad²` is *forced* to exactly 2 independent of the physics; the genuine drift evidence is `N=5,6,7`, with ≥3 distinct bonds):

| N | `dRe/grad²` (interior → end) | CV | slope (Q-correction, away from 2) | `sin²` shape-miss (converges with N) |
|---|---|---|---|---|
| 4 | 1.53, 1.53 | 0.001 | 2.00 | 0.17 |
| 5 | 0.73, 0.67 | 0.040 | 2.12 | 0.12 |
| 6 | 0.32, 0.31, 0.28 | 0.060 | 2.16 | 0.08 |
| 7 | 0.18, 0.17, 0.15 | 0.070 | 2.17 | 0.06 |
| 8 | n/a | n/a | n/a | 0.043 |
| 9 | n/a | n/a | n/a | 0.035 |

(The N=8, 9 rows recompute only the `sin²` shape-miss, `felt_time_highn_shapemiss.py`, the half-filling `(p,p)` block reaching `15876²` at N=9, the dense ceiling on 128 GB; the Q-correction columns are the N≤7 landing. The monotone descent `0.17 → 0.035` confirms the continuum harmonic as the large-N limit.)

**Provenance.** The live witness `inspect --root gradient` covers `N = 4, 5` only (it throws above N=5); the `N = 6, 7` rows of both tables, and the `N = 7` `n(j)` profile of §1, come from the Python verifier `felt_time_amplitude_law.py`. The slope and CV are sector-independent and agree between the two engines; the absolute `dRe/grad²` column is the verifier's sector (e.g. `(3,3)` at `N=5`), while the witness's own tie-break (`(1,1)` at `N=5`) reproduces the same slope/CV at a rescaled magnitude (§1). Both engines confirm the *law*; only the bond-independent prefactor is sector-labelled.

The squared-gradient law is thus **exact in the strong-dephasing (diffusion) limit, holds with small dressing corrections through the soft-survivor regime, and is superseded at the handover**; it is not an exact all-`Q` closed form, and the prefactor constant is scaling-plus-sign, consistent with the two-lens review of the trajectory side.

The companion verifier `felt_time_closure_functional.py` confirms the **trajectory dual**: the PTF painter closure `Σ_i ln(α_i)` reads this same rate shift, sign-coherent (`coh ~ 1`) only at the high-gradient bonds, there matching `N·|dRe|/|reS|` in sign and `O(1)` magnitude; the low-gradient end bonds read as a redistribution (`coh < 0.8`) and are correctly not counted. The eigenvalue-level law (this proof, the stone's bond functional) and the trajectory-level closure (the stone, `StoneSurvivorClosureClaim`) are one fact read at two levels.

## Where it lives

- **Typed + the regime probe:** [`SurvivorDiffusionGradientClaim`](../../compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs) (Tier 1 candidate; parents `AbsorptionTheoremClaim` + `SurvivalIncompletenessMirrorClaim`), live witness `inspect --root gradient` ([`SurvivorDiffusionGradientWitness`](../../compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientWitness.cs)): reproduces the `N = 4, 5` rows and the off-diagonal weight (it caps at `N = 5`; the `N = 6, 7` rows are `felt_time_amplitude_law.py`), and the Q-sweep is just driving it across `--q`.
- **Verifiers:** [`felt_time_amplitude_law.py`](../../simulations/felt_time_amplitude_law.py), [`felt_time_closure_functional.py`](../../simulations/felt_time_closure_functional.py).
- **Registry:** [F123](../ANALYTICAL_FORMULAS.md).
- **Outward reading:** [`ON_THE_FOUR_DIRECTIONS`](../../reflections/ON_THE_FOUR_DIRECTIONS.md) (the fourth direction, felt time, is the in-plane shape read by the watching).
- **Trajectory dual:** the stone, [`StoneSurvivorClosureClaim`](../../compute/RCPsiSquared.Diagnostics/Foundation/StoneSurvivorClosureClaim.cs) (`inspect --root stone`).
