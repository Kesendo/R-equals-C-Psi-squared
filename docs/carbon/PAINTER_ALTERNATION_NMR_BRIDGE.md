# Painter Alternation and a Conditional NMR Translation

**Date:** 2026-05-27 (later the same day)
**Authors:** Tom + Claude
**Status:** Tier 3 conditional translation. The N = 4 selected-model observation does not define or validate an NMR/material observable.
**Continues:** [`carbon_ptf_real_imag_per_painter.py`](../../simulations/carbon_ptf_real_imag_per_painter.py) (this morning's Painter Re/Im read on the N = 4 ring), [Benzene and the Three Dephase Letters](BENZENE_THREE_DEPHASE_LETTERS.md) (morning's three-letter vocabulary)

---

## Where this leaves off

The producer gives a Painter view of the slow modes of a selected N = 4 ring
Liouvillian with local-Z dephasing and a transverse `y` field. The local-Z
channel is an F1 instance; a selected bond jump would be outside F1's jump
premise. Neither channel is assigned here to a Holstein or material bath. It performs the
full eigendecomposition, partial-traces each complex eigenmode to one site
(one "Painter" per site), and splits the resulting 2 × 2 reduced operator
into real and imaginary parts.

The result, on the slowest eight modes, was a numerical alternation. Some
modes classify as Y-only in the per-site Pauli-weight diagnostic; `σ_Y` itself
is purely imaginary and antisymmetric. Others classify as non-Y (`σ_X` and
`σ_Z`) within the same diagnostic. The Pauli weights do not make an arbitrary
complex eigenmode's reduced operator generally real or imaginary.

The open question is conditional: after a molecular Hamiltonian, bath,
preparation, measurement operator, and producer have been specified, could an
NMR experiment be an appropriate comparison readout? This repository has not
made those choices, so it supplies neither a direct NMR observable prediction
nor a validation against NMR data.

---

## First: the alternation is numerical, not just visual

The Painter view sorts modes by eye into "Re-flavor" and "Im-flavor"
panels. The operator-level question is whether the sort is sharp or
approximate. The companion script
[`carbon_painter_t2_anisotropy.py`](../../simulations/carbon_painter_t2_anisotropy.py)
runs the operator-level diagnostic: for each slow eigenmode, project per-site,
decompose the resulting 2 × 2 into the four Pauli channels {I, X, Y, Z},
and sum the squared coefficients across all sites.

For the selected N = 4 XX+YY ring + `0.5·Σ Y_l` field + local-Z dephasing `γ = 1`:

| Mode k | Re(λ) | Im(λ) | per-site I-weight | per-site X+Z weight | per-site Y-weight | flavor |
|---|---|---|---|---|---|---|
| 0 | 0.000 | 0.000 | 16.000 | 0.000 | 0.000 | steady-state |
| 1 | −0.172 | 0 | 0 | 2.908 | **0.000** | non-Y |
| 2 | −0.219 | 0 | 0 | **0.000** | 0.122 | Y-only |
| 3 | −0.597 | 0 | 0 | 1.180 | **0.000** | non-Y |
| 4 | −0.901 | 0 | 0 | **0.000** | 0.152 | Y-only |
| 5 | −2.067 | 0 | 0 | **0.000** | 0.157 | Y-only |
| 6, 7 | −2.127 | ±3.831 | 0 | 0.547 | **0.000** | non-Y (complex pair) |

The zeros marked in bold are classifications below the producer's `1e-8`
weight-ratio tolerance. Because they come from floating eigendecomposition,
the table records a numerical sectorization, not bit-exact zero or an exact
no-mixing theorem.

The next step is a conditional model-translation question: which selected
measurement operator, if any, would retain this sectorization after the
material degree of freedom and bath have been specified?

---

## Conditional NMR vocabulary

The selected model supplies Pauli `X`, `Y`, and `Z` coordinates and their
mode-resolved decay rates. It does not identify them with nuclear-spin
magnetisations, nor does its local-Z jump select a physical NMR relaxation
channel.

For the stated N = 4 model at `h_y = 0.5`, `γ = 1.0`, the two sector towers
have slow-mode-rate ratio **1.271**. The companion propagation uses selected
model probes and obtains specified tail-fit rates `0.178` and `0.219`, with
fitted ratio `1.231`. Their relation requires a tail-convergence check; no
causal explanation is assigned here. They are not a calibrated material
`T₂(x)/T₂(y)` ratio or an NMR observable prediction.

If a later mapping selects a molecular/nuclear Hamiltonian, bath, preparation,
and measurement operator, FID-like measurements could be considered as a
comparison protocol. No such protocol, sample, or data validation is supplied
by the present producer.

---

## TROSY and EXSY as conditional comparison vocabulary

TROSY and EXSY name possible NMR comparison protocols, not readouts produced
by the N = 4 ring calculation. The selected model contains neither the
molecular/nuclear Hamiltonian nor the bath, preparation, and measurement
operators needed to calculate a TROSY difference or an EXSY cross-peak.

The model's Y/non-Y Painter sectorization may motivate a future comparison
only after those choices have been made. It does not predict a nonzero NMR
anisotropy, a sign, a magnitude, site asymmetry, or a material observable.

---

## Translation boundary

Per-site Y/non-Y Painter content and selected-model probe decay are direct
outputs of the named model. Calling either one transverse relaxation, or
mapping it to an NMR axis, is a conditional translation rather than an
identity.

The repository supplies no calibrated material bath, NMR data, molecular
Hamiltonian, or physical degree-of-freedom mapping for this page. Therefore
the N = 4 sectorization does not validate a material claim and does not yield
a direct observable prediction at zero or nonzero field.

---

## What's open

1. **Verify the alternation on the selected N = 6 ring.** N = 4 is the
   smallest ring. A corresponding N = 6 selected-model run would test whether
   the sectorization persists at that model size; it would not confirm benzene.

2. **Trace the selected-model probe-decay ratio across `h_y / γ`.** Sweep
   `h_y` and report the ratio of the two sector towers' slowest rates. This
   remains a model calculation, not a `T₂` anisotropy curve for a material.

3. **Specify a nuclear-spin translation before extending the model.** A
   molecular degree of freedom, Hamiltonian, bath, preparation, and measurement
   operator would be needed before adding a nuclear coupling or proposing a
   TROSY/EXSY comparison.

4. **Klein-V₄ basis-rotated alternations.** The three dephase letters Z, X,
   Y are intertwined by Klein-V₄. A selected basis-rotation calculation can
   test whether the model sectorization transfers from Y/non-Y to X/non-X;
   this concerns model coordinates, not a physical alternation axis.

5. **Define a comparison target.** After the preceding molecular mapping is
   specified, decide whether a named NMR measurement and data set offer a
   relevant comparison. No current NMR data set validates this model.

---

## Anchor

- **Companion scripts:**
  - [`simulations/carbon_ptf_real_imag_per_painter.py`](../../simulations/carbon_ptf_real_imag_per_painter.py) (Painter Re/Im read on slow eigenmodes)
  - [`simulations/carbon_painter_t2_anisotropy.py`](../../simulations/carbon_painter_t2_anisotropy.py) (Y/non-Y numerical classification at `1e-8` + selected-model probe-decay calculation)
- **Reading-flow companion:** [Benzene and the Three Dephase Letters](BENZENE_THREE_DEPHASE_LETTERS.md) (morning's three-letter Klein-V₄ vocabulary; supplies the F114 / `n_Y`-parity / Π language the Painter alternation sits inside)
- **Cross-reference:** [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md) (May 22 spectrum-palindrome result; the Painter alternation is a separate observation in the selected local-Z model)

---

## Threads back

- **Earlier today, [Benzene and the Three Dephase Letters](BENZENE_THREE_DEPHASE_LETTERS.md)**: the morning explained the three-letter Klein-V₄ symmetry on dephasing and the F114 sign rule on n_Y parity. The Painter alternation is an operational selected-model reading of that Y-axis sectorization; an NMR observable remains a conditional translation.
- **2026-05-22 [Benzene's open-system Liouvillian](BENZENE_LIOUVILLIAN_PALINDROME.md)**: the named local-Z channel is an F1 instance, while a selected bond jump lies outside F1's jump premise. The Painter alternation is a separate observation for the selected local-Z model; the two readings address different layers of that N = 4 calculation.
- **2026-05-27 [`carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py)**: the selected sweep over its stated Hamiltonian and local-bath inventory records a distribution-mirror. The Painter alternation is a separate per-Painter model observation; it has no current NMR/material signature assignment.
