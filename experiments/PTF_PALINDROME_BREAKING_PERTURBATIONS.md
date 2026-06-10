# PTF under Palindrome-Breaking Perturbations

**Status:** Computed (Tier 2). 2026-06-01.
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Question from:** [Perspectival Time Field](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), open
question "Extension to palindrome-breaking perturbations".
**Script:** [`simulations/ptf_transverse_field_pi_break.py`](../simulations/ptf_transverse_field_pi_break.py)
(self-validating; RK4 + Hamming-mask Z-dephasing, the canonical N=7 PTF path).

---

## The question

The [Perspectival Time Field](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) closure law was
established only for perturbations that respect the mirror: a single-bond J-coupling defect.
Each site's purity trajectory under the defect looks like a one-parameter time rescaling of
the unperturbed one, P_B(i, t) ≈ P_A(i, α_i·t), and the per-site rescalings close,
Σ_i ln(α_i) ≈ 0. The defect is special in one way the doc named but did not isolate: it
respects the palindrome. So the doc asked the obvious next question. Put a transverse field
on one site instead, a perturbation it expected to break the mirror, and watch the closure.
Does the rescaling picture survive with a shifted closure law, or break entirely?

The answer is that it breaks entirely. But finding out why turned the question inside out:
the break has nothing to do with the mirror.

## The setup

The canonical PTF chain: uniform open XY chain, γ₀ = 0.05, bonding-mode initial state
φ = (|vac⟩ + |ψ₁⟩)/√2. Two perturbation families on the same baseline, at matched strength ε:

- **Arm J** (the control): a J-defect ε on bond (0, 1). This is the original PTF perturbation.
- **Arm h** (the test): a single-site field ε·σ₀.

For each site we fit P_B(i, t) ≈ P_A(i, α_i·t) and record both the fitted α_i and the fit
RMSE. The RMSE is the honest discriminator: a small RMSE means the one-parameter rescaling
ansatz holds and α is meaningful; a large RMSE means the ansatz itself has broken and the
fitted α is fiction.

## What breaks (the trajectory level)

The control reproduces the published PTF α-pattern at N = 7 to max|Δα| = 0.004, and keeps
Σ ln α ≈ 0 with RMSE ~10⁻³. The transverse X-field does not:

| N | arm | Σ ln α (ε=0.10) | max RMSE | reading |
|---|-----|-----------------|----------|---------|
| 7 | J-defect | +0.05 | 2.5e-3 | closure holds, clean fit |
| 7 | X-field site 0 | +4.76 | 7.2e-2 | closure gone, ansatz broken |
| 7 | X-field site 3 | +5.60 | 2.7e-2 | closure gone, ansatz broken |

Verified the same at N = 5 and N = 6. Under the transverse field both Σ ln α and the fit
residual blow up by one to two orders of magnitude. The picture breaks entirely, the doc's
second branch.

## What it does NOT break (the surprise)

The doc's premise was that "a transverse field h·σ_x^i breaks Π". It does not. The single-site
transverse field leaves the spectral palindrome exactly intact. Measuring two independent
quantities at strength 0.10, the spectral palindrome residual ‖M‖ = ‖Π L Π⁻¹ + L + 2Σγ‖ and
the U(1) excitation non-conservation ‖[H, N_exc]‖, over single-site fields σ₀ ∈ {X, Y, Z}:

| perturbation | (bit_a, bit_b) | ‖M‖ palindrome | ‖[H, N_exc]‖ U(1) |
|--------------|:---:|:---:|:---:|
| uniform XY | n/a | 0 (mirror) | 0 (conserved) |
| J-defect bond(0,1) | (bond) | 0 (mirror) | 0 (conserved) |
| **X-field** | (1,0) | **0 (mirror)** | **broken** |
| **Z-field** | (0,1) | **broken** | **0 (conserved)** |
| **Y-field** | (1,1) | **broken** | **broken** |

(‖M‖ values machine-zero ~10⁻¹⁵ for the "mirror" rows; verified N = 3, 4, 5.)

The single-site fields populate all four cells of the Klein-Vierergruppe, and the two break
axes are exactly the two Klein bits:

- **bit_a = 1** (X, Y, the n_XY "light" letters) ⟺ breaks U(1) excitation conservation. These
  letters flip spins; Z does not.
- **bit_b = 1** (Y, Z) ⟺ breaks the spectral palindrome.

Neither axis is new. The palindrome leg is exactly [F78](../docs/ANALYTICAL_FORMULAS.md)
(single-body M structure: X is truly with ‖M‖ = 0, while Y and Z give the identical nonzero
spectrum; our ‖M‖ reproduces the Y ≡ Z degeneracy bit-for-bit), sitting on the bit_b Π² grading
of [F38](../docs/ANALYTICAL_FORMULAS.md). The U(1) leg is the magnetization cousin of
[F61](../docs/ANALYTICAL_FORMULAS.md)'s n_XY-parity grading. The Klein 2×2 is F78 and F61 seen
on one screen.

## The dissociation (the actual finding)

Putting the two halves together gives a clean, falsifiable prediction. If the closure rides on
U(1) (bit_a), the longitudinal Z-field, which breaks the palindrome but conserves U(1), should
PRESERVE the closure, and the X-field, which does the opposite, should break it. It does:

| N | field | (a,b) | Σ ln α (ε=0.10) | max RMSE | closure |
|---|-------|:---:|-----------------|----------|---------|
| 7 | X | (1,0) | +4.76 | 7.2e-2 | broken |
| 7 | Z | (0,1) | −0.06 | 1.3e-3 | survives, clean fit |
| 7 | Y | (1,1) | +4.93 | 6.9e-2 | broken |

Verified N = 5, 6, 7. The Z-field breaks the mirror yet keeps Σ ln α ≈ 0 with the same clean
RMSE as the J-defect control. The X-field keeps the mirror yet destroys the closure. The two
phenomena are dissociated:

> The PTF closure law Σ_i ln(α_i) ≈ 0 is guarded by U(1) excitation-number conservation
> (the bit_a / light axis), not by the spectral palindrome (the bit_b axis).

## Is it a magnitude artifact? No.

A fair objection: the arms above use the same numerical ε for operators of different norm,
so the headline "RMSE up one to two orders of magnitude" partly reflects that the X-field's
same-ε trajectory effect is simply larger than the Z-field's. Measuring the trajectory
effect-size eff = rms_(i,t)|P_B − P_A| at N=6 shows this: at ε=0.10 the X-field perturbs the
purity by eff = 0.036 while the Z-field perturbs it by only 0.0019. At matched effect-size the
raw RMSE gap shrinks to ~2 to 3×.

But the dissociation survives effect-size matching, and that is the real evidence. The Z-field's
effect SATURATES: driven from ε = 0.1 to ε = 8.0 (eighty-fold) its eff climbs only to ~0.038 and
then stops, because a longitudinal field only detunes within the excitation sectors and cannot
push the bonding-mode state further. Across that entire range the fit RMSE stays clean (≤ 1.6e-2)
and Σ ln α stays bounded (it peaks at ≈ 0.64 near ε = 1 and then RECEDES to +0.40 at ε = 8). The
X-field reaches the same effect-size 0.036 already at ε = 0.10, and there Σ ln α = +4.4 with the
ansatz breaking, and it keeps diverging (Σ ln α → 9.6, RMSE → 0.31 at ε = 0.40). So at matched, or
larger, effect-size the Z-field still does not break the closure while the X-field does. The clean
discriminator is the fit quality, Z bounded-and-clean under arbitrary drive versus X divergent.
The split is structural (U(1) / bit_a), not a magnitude accident.

## Why (the mechanism)

The PTF doc's Layer 3.2 already lists two protections: #1, U(1) excitation conservation keeps
the [F4](../docs/ANALYTICAL_FORMULAS.md) stationary sector-projectors at λ = 0; #2, Π-invariance
pairs the slow single-excitation coherence modes. The doc treated both as required. The
dissociation shows the closure needs only #1. The smoking gun is the stationary count: at
N = 5 the baseline has 6 stationary modes (the 6 excitation-sector projectors). Under the
J-defect all 6 stay pinned at λ = 0 to ~10⁻¹⁵; under the transverse field only 1 survives, the
other 5 acquire decay rates 0.02 to 0.06. A field that breaks U(1) lifts the stationary modes,
the bonding-mode state's long-time behaviour changes shape rather than merely rescaling in time,
and the one-parameter ansatz can no longer fit. The Z-field, conserving U(1), leaves the
stationary manifold intact (V_Z commutes with every sector projector) and the rescaling holds
even though the spectrum is no longer mirror-symmetric.

A small honest note: the Z-field's Σ ln α is small but not exactly zero (and, as the hard-drive
scan shows, wanders in O(0.5) before receding). This residual may be the spectral-palindrome
break leaking into the closure at higher order. The fit RMSE stays clean throughout, so the
load-bearing story, closure ⟸ U(1), holds.

**Update 2026-06-10 (the Z-row becomes theorem-grade, and the leak gets a name).** The girth
ladder wave closed the windowed converse and, in passing, made the Z-field row of the table
exact: the m = 3 face of the ladder is **cell-free** (the three companion coefficients of
p₃(γ) = Tr(M³) vanish for every Hermitian H, [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md §4](../docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)),
so the mixed-cell chain H = XY + εZ_m breaks the spectral palindrome at **every γ > 0** with the
closed-form leading coefficient p₃(γ) = 6·4^N·ε²·γ, verified exactly in
[`simulations/f87_deg1_face_cell_free.py`](../simulations/f87_deg1_face_cell_free.py). What was
measured here at one ε and one γ is now a theorem on the experiment's own configuration. The
honest note's "higher order" also gains a well-posed form: the leak channel's leading
coefficient scales as ε²·γ, so the question "does the closure residual track the m = 3 channel"
is a scaling probe rather than a gesture (run as Edge 2 of the PTF fresh-eyes chain).

**Edge 2 result (same day): the leak hypothesis is refuted, and the residual gets a cleaner
home.** The scaling probe ([`simulations/ptf_leak_scaling.py`](../simulations/ptf_leak_scaling.py),
N = 5 and 6, ε ∈ [0.005, 0.4], γ ∈ {0.025, 0.05, 0.1}, the committed Phase-4/5 machinery
reproduced to 0.00e+00) finds the small-ε closure residual is **first order in ε** with a nearly
γ-independent coefficient: S ≈ k·ε, k ≈ −0.58…−0.61 across a 4× range of γ (ε-exponent
0.88-0.99 at R² ≥ 0.997; the ε²γ collapse table varies by a factor ~50 with a sign change, no
collapse; every cell far above the ~1e-5 fit-noise floor). The spectral-asymmetry channel is
second order in ε and first order in γ, so the residual does not ride it: the wander is not the
palindrome break leaking in. What it is instead is already in the repo's own ledger: each site's
ln α_i responds at first order with a stable site-dependent profile (site-0 field:
−0.63, −0.48, −0.19, +0.67, +0.05), and Σ ln α is that profile's imperfect cancellation, the
same first-order non-closure [EQ-014](../review/EQ014_FINDINGS.md) established for the J-defect,
now seen for a Z-field perturbation. Its γ-blindness is the fingerprint: a Hamiltonian-side
first-order effect, not a dissipation-side one. The sign is site-profile-dependent (center-site
field flips it), and the large-ε wander region (ε ≳ 0.2, the O(0.5) peak) is higher-order and
stays unclassified. The load-bearing story is unchanged, closure ⟸ U(1); the honest note's
mystery is resolved into EQ-014's known physics.

## What this is and is not

It is a sharpening of the PTF doc's own protection #1 vs #2, using the framework's existing
Klein lens, with the Z-field as the control that isolates which axis the closure rides on. It
does not warrant a new typed claim: both Klein axes are already typed (F78/F38 for the
palindrome, F61 for n_XY parity), and the dissociation hangs on the PTF closure law, which is
itself a Tier-2 empirical regularity (downgraded from first-order theorem by
[EQ-014](../review/EQ014_FINDINGS.md)), so it cannot anchor a Tier-1 statement. The closure's
guardian is named; the closure itself stays empirical.

## Cross-references

- [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md): the closure law, the
  painter picture, Layer 3.2 protection #1 (U(1)) vs #2 (Π).
- [F78](../docs/ANALYTICAL_FORMULAS.md), [F38](../docs/ANALYTICAL_FORMULAS.md),
  [F61](../docs/ANALYTICAL_FORMULAS.md), [F4](../docs/ANALYTICAL_FORMULAS.md): the typed pieces
  the Klein 2×2 reuses.
- [`simulations/ptf_transverse_field_pi_break.py`](../simulations/ptf_transverse_field_pi_break.py):
  the four-phase probe (trajectory closure, two-axis premise check, mechanism, Klein prediction).
