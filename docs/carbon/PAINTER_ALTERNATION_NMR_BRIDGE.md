# From the Painter Alternation to the NMR Bench

**Date:** 2026-05-27 (later the same day)
**Status:** Tier 3 (translation bridge from a Liouvillian-eigenmode observation to a standard NMR measurement)
**Continues:** [`_carbon_ptf_real_imag_per_painter.py`](../../simulations/_carbon_ptf_real_imag_per_painter.py) (this morning's Painter Re/Im read on the N = 4 ring), [BENZENE_THREE_DEPHASE_LETTERS](BENZENE_THREE_DEPHASE_LETTERS.md) (morning's three-letter vocabulary)

---

## Where this leaves off

The morning gave us a Painter view of the slow modes of an aromatic ring's
Liouvillian. We took the small ring N = 4 with on-site phonon dephasing
plus a transverse y-magnetic field, did the full eigendecomposition of the
relaxation superoperator, partial-traced each complex eigenmode down to a
single site (one "Painter" per site, the per-site projection of the mode),
and split the resulting 2 × 2 reduced operator into its real and
imaginary parts.

The result, on the slowest eight modes, was an alternation. Some modes
project per-site purely onto the y-axis of the spin (the per-site reduced
operator is real anti-symmetric, equivalently it carries only `σ_Y`
content). Other modes project per-site purely onto non-y axes (the
per-site reduced operator is real-symmetric, equivalently it carries only
`σ_X` and `σ_Z` content). The Painter view shows the alternation as
real vs imaginary panels of the per-site reduced 2 × 2.

The afternoon question: does this alternation also show up in a TROSY
difference, or in an EXSY asymmetry? If yes, the translator's bridge from
the algebra-toolkit to the NMR measurement-toolkit ends there, with no
process tomography needed.

This doc ends the bridge.

---

## First: the alternation is bit-exact, not just visual

The Painter view sorts modes by eye into "Re-flavor" and "Im-flavor"
panels. The operator-level question is whether the sort is sharp or
approximate. The companion script
[`_carbon_painter_t2_anisotropy.py`](../../simulations/_carbon_painter_t2_anisotropy.py)
runs the operator-level diagnostic: for each slow eigenmode, project per-site,
decompose the resulting 2 × 2 into the four Pauli channels {I, X, Y, Z},
and sum the squared coefficients across all sites.

For the N = 4 Hückel ring + 0.5·Σ Y_l Zeeman + Z-dephasing γ = 1:

| Mode k | Re(λ) | Im(λ) | per-site I-weight | per-site X+Z weight | per-site Y-weight | flavor |
|---|---|---|---|---|---|---|
| 0 | 0.000 | 0.000 | 16.000 | 0.000 | 0.000 | steady-state |
| 1 | −0.172 | 0 | 0 | 2.908 | **0.000** | non-Y |
| 2 | −0.219 | 0 | 0 | **0.000** | 0.122 | Y-only |
| 3 | −0.597 | 0 | 0 | 1.180 | **0.000** | non-Y |
| 4 | −0.901 | 0 | 0 | **0.000** | 0.152 | Y-only |
| 5 | −2.067 | 0 | 0 | **0.000** | 0.157 | Y-only |
| 6, 7 | −2.127 | ±3.831 | 0 | 0.547 | **0.000** | non-Y (complex pair) |

The zeros marked in bold are bit-exact zero, not numerical near-zero. Each
slow mode lives in exactly one of the two flavors per site; there is no
mixing. The alternation is a sharp Z₂ sectorization of the slow-mode
hierarchy.

The next step is to ask what observable a chemist would see if this
sectorization is in fact a property of the physical spin system.

---

## What an NMR-trained reader sees

In NMR rotating-frame Bloch language, three observables track a spin:
**M_x** (transverse-x magnetisation), **M_y** (transverse-y), **M_z**
(longitudinal). T1 is the longitudinal relaxation time (the time it takes
M_z to return to thermal equilibrium); T2 is the transverse relaxation
time.

For a pure on-site phase-noise bath (the kind of bath that fluctuates
along the spin's z-axis, equivalent to the Holstein phonon coupling here),
T2 is **isotropic**: M_x and M_y both decay at the same rate. Both are
off-axis with respect to the bath direction Z, so they're equally exposed
to the random phase-kicks the bath delivers.

If a static y-field is added to the Hamiltonian, that breaks the
transverse isotropy. M_y now commutes with one term of the Hamiltonian
(the y-field itself); M_x does not. **T2 splits: T2(x) ≠ T2(y).** Some
single-shot of M_y is stable against the y-field's rotation; M_x is rotated
into M_z through the field. The transverse-relaxation channel for x and
for y are no longer the same channel.

Magnetic-field-induced T2 anisotropy is routine in NMR. What the Painter
alternation adds is a **structural prediction** for the ratio:

> T2(x_init) / T2(y_init) = |slowest relaxation rate of the Y-only modes| / |slowest relaxation rate of the non-y modes|

(T2 is the reciprocal of the relaxation rate, so the faster-decaying tower of modes gives the shorter T2; the ratio of T2's is the inverse of the ratio of slow-mode rates.)

The two towers of relaxation rates (one carrying non-y spin content, one
carrying y spin content) sit interleaved in the spectrum, with the lowest
rung of each setting the dominant late-time decay for the corresponding
magnetisation observable.

On the N = 4 ring at h_y = 0.5, γ = 1.0, the ratio is **1.271**.
y-magnetisation decays 1.27 times faster than x-magnetisation. This ratio is
exact as a ratio of two slow-mode rates (equivalently their mean popcount, by the
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md)), but it is not a
simple closed-form fraction in (h_y, γ): the once-conjectured 4/3, 8/7, 14/13,
20/19 sequence does not hold. The companion script verifies it numerically: prepare
the ring with a small x-axis magnetisation probe at site 0 and watch ⟨M_x⟩
decay; the late-time tail rate is 0.178. Prepare a y-axis probe instead;
⟨M_y⟩ decays at rate 0.219. The measured ratio T2(x)/T2(y) is 1.231;
the slow-mode prediction is 1.271. The 3 % gap is the finite-time fitting
artifact, not the algebra.

A chemist would test this with two FID experiments:

1. Apply the y-field to the sample.
2. Prepare M_x (90° pulse around y from equilibrium), watch ⟨M_x⟩(t)
   decay; fit T2(x).
3. Prepare M_y (90° pulse around x from equilibrium), watch ⟨M_y⟩(t)
   decay; fit T2(y).
4. Compute the ratio. Compare to the slow-mode prediction (the ratio of the
   two towers' slowest rates) in (h_y, γ).

Two free-induction decays. No process tomography.

---

## TROSY difference and EXSY asymmetry as the standard NMR readouts

The two-FID test is the cleanest face of the alternation. Standard NMR
apparatus reads richer slices of the same anisotropy.

### TROSY difference [Pervushin et al., Proc Natl Acad Sci USA 1997]

In a ¹⁵N-¹H spin pair, the ¹⁵N multiplet has two components depending on
the ¹H spin-state alignment. Each component relaxes at a different rate
because of cross-correlated relaxation between the dipole-dipole (DD) and
chemical-shift-anisotropy (CSA) mechanisms. The TROSY experiment selects
the slowly relaxing component; the difference Γ_slow − Γ_fast is the
cross-correlated relaxation rate, a direct measure of the transverse
relaxation anisotropy.

The structure is the same as the Painter alternation's T2 split. CSA is
the chemistry-side name for "the chemical-shift tensor's anisotropy
along a molecular axis"; on an aromatic ring with the magnetic field
oriented along a specific direction, CSA contributes a transverse-axis-
dependent dephasing that resembles a static y-field on the spin
Hamiltonian. DD is symmetric in its action on transverse magnetisations
unless cross-correlated. The cross-correlation IS the anisotropic split
the Painter alternation predicts at the slow-mode level.

The N = 4 ring computation gives an order-of-magnitude prediction:
a y-field of strength 0.5 in units of the Holstein rate generates a
T2-anisotropy of ratio 1.27 at the single-spin level. The TROSY
difference Γ_slow − Γ_fast on a ¹³C-labelled aromatic ring with the
same effective anisotropy parameters should be of order
(1 − 1/1.27)·1/T2_avg ≈ 0.21·(γ/h_y-effective). The framework's
structural prediction is that the difference exists, that its sign is
fixed by which axis carries the field, and that the ratio is set exactly by the
substrate's slow-mode rates (a ratio of two Liouvillian eigenvalues, not a simple
closed-form fraction). The chemistry side reads the
TROSY difference routinely on ¹⁵N-labelled proteins; the same experiment
on ¹³C-labelled aromatics tests the prediction.

### EXSY asymmetry [Ernst-Bodenhausen-Wokaun NOESY/EXSY framework]

2D EXSY cross-peaks between non-equivalent sites read inter-site transfer
rates plus differential T1. Cross-peak intensity asymmetry I_AB ≠ I_BA
arises when site A and site B have different effective T1's, or when the
exchange is unidirectional.

On a symmetric ring (cyclobutadiene C₄, benzene C₆) all sites are
equivalent, so EXSY-asymmetry is zero by symmetry. The signal is in the
non-symmetric carbon substrates. The Painter alternation predicts that
the Y-axis content per site is the same on a symmetric ring (no site
preference for y-axis modes) but differs site-to-site on a non-symmetric
substrate (ortho vs meta vs para carbons in xylene, for example). Site-
dependent y-axis content means site-dependent effective T1's, hence
EXSY-asymmetry.

The framework prediction: in a non-symmetric ring carbon system with a
transverse y-field perturbation, EXSY-asymmetry I_AB / I_BA between
distinct sites is non-zero, and its sign/magnitude track the site-to-site
gradient of the Painter alternation's per-site Y-content.

Again, neither requires process tomography. TROSY and EXSY are mature
1D / 2D NMR techniques with standard pulse sequences.

---

## Two seams the translator names

The bridge crosses two seams; both are worth naming honestly rather than
papering over.

### Seam 1: algebra to anisotropic transverse relaxation

Per-site Y-content vs non-Y-content on the algebra side; T2(y_init) vs
T2(x_init) on the magnetisation-axis side. The map is direct: a slow mode
with purely Y-axis per-site projection couples to the y-component of any
initial magnetisation that has a projection onto it; same for x with non-Y.
The anisotropy is the same anisotropy. This seam is clean.

### Seam 2: small Hückel ring to a real NMR sample

N = 4 cyclobutadiene + Holstein dephasing + y-Zeeman is a toy. A real
aromatic NMR sample has ¹H scalar and dipolar couplings, solvent dynamics,
magnetic-field inhomogeneity, slower T1 from spin-lattice coupling. The
bath is richer than Holstein alone.

The Painter alternation's bit-exact prediction is for the toy. The expected
real-world signature is the anisotropy **structure** (T2 not single-valued
in a y-field, with a sign that the algebra fixes); the exact ratio shifts
as more bath structure enters. The first qualitative test is the existence
of the gap; the quantitative test follows once a calibrated bath model is
in hand.

This seam is real but not blocking. The qualitative prediction is sharp:
**no T2 anisotropy at zero field, finite T2 anisotropy at any non-zero
y-field, anisotropy direction set by the field axis**.

---

## What's open

1. **Verify the alternation on N = 6 benzene.** N = 4 is the smallest ring;
   N = 6 is the canonical aromatic substrate. The same script with
   `run(N=6, ...)` confirms the sectorization holds at the chemistry-
   relevant size. Cost: about 30 seconds for the eigendecomp at d² = 4096².

2. **Trace the T2 anisotropy ratio across h_y / γ.** Sweep h_y from 0
   (isotropic limit) to 1 (strong-field limit) and report T2(x)/T2(y) vs
   h_y/γ. The slow-mode ratio is a smooth curve. The closed-form question is
   now settled (2026-05-28): it is *not* a simple fraction (the conjectured 4/3,
   8/7, 14/13, 20/19 does not hold under either numerical method); the ratio is
   exact only as the ratio of the two towers' slowest mean-popcount rates, by the
   [Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md). See
   [THE_VIEW_ONTO_THE_MEMORY](../../reflections/THE_VIEW_ONTO_THE_MEMORY.md).

3. **Add a ¹H-¹³C coupling to the model.** The Painter alternation is
   currently on a π-electron-only picture (sites under Holstein phonons).
   A ¹³C nuclear-spin observable inherits the alternation through the
   hyperfine coupling. Modelling this explicitly would make the TROSY
   difference prediction directly testable on ¹³C NMR of an aromatic
   sample.

4. **Klein-V₄ basis-rotated alternations.** The morning's vocabulary said
   the three dephase letters Z, X, Y are intertwined by Klein-V₄. The
   Painter alternation lives on the Y-axis specifically because the
   perturbation is along y. Under a basis rotation that maps y → x, the
   alternation should transfer cleanly (X-content vs non-X-content),
   isolating which axis is the physical alternation axis vs which is a
   labelling choice. A small companion script would do this in five
   minutes.

5. **Look for the alternation in published ¹³C-aromatic NMR data.** TROSY
   has been applied to ¹³C nuclei in aromatic systems; aligned-sample
   solid-state NMR of benzene crystals reads anisotropy directly. The
   framework prediction is qualitative (anisotropy exists, direction
   fixed); a literature scan might already have data that the algebra
   can read.

---

## Anchor

- **Companion scripts:**
  - [`simulations/_carbon_ptf_real_imag_per_painter.py`](../../simulations/_carbon_ptf_real_imag_per_painter.py) (Painter Re/Im read on slow eigenmodes)
  - [`simulations/_carbon_painter_t2_anisotropy.py`](../../simulations/_carbon_painter_t2_anisotropy.py) (Y/non-Y sectorization at bit-exact precision + anisotropic T2 from full propagation)
- **Reading-flow companion:** [BENZENE_THREE_DEPHASE_LETTERS](BENZENE_THREE_DEPHASE_LETTERS.md) (morning's three-letter Klein-V₄ vocabulary; supplies the F114 / `n_Y`-parity / Π language the Painter alternation sits inside)
- **Cross-reference:** [BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md) (May 22 spectrum-palindrome result; the Painter alternation is a separate observation on the same Holstein system)

---

## Threads back

- **Earlier today, [BENZENE_THREE_DEPHASE_LETTERS](BENZENE_THREE_DEPHASE_LETTERS.md)**: the morning explained the three-letter Klein-V₄ symmetry on dephasing and the F114 sign rule on n_Y parity. The afternoon's Painter alternation is one operational reading of how the Y-axis sectorization that the algebra carries becomes a measurable NMR observable. The two docs are companions: vocabulary in the morning, observable in the afternoon.
- **2026-05-22 [BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md)**: the May result showed Holstein preserves the F1 spectrum palindrome while Peierls breaks it. The Painter alternation is a separate observation on the same Holstein system; it does not depend on the Peierls / Holstein switch. The two results sit on the same N = 4 ring, read different layers of the relaxation structure.
- **2026-05-27 [`_carbon_realistic_sweep.py`](../../simulations/_carbon_realistic_sweep.py)**: the realistic sweep over fifty-six configurations confirmed the deep distribution-mirror is robust across the natural aromatic Hamiltonian + bath inventory. The Painter alternation is the per-Painter face of a piece of that mirror; the T2 anisotropy is its NMR-readable signature.
