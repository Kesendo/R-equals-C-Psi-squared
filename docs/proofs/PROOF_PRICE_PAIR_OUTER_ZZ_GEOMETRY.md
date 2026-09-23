# Price Pair run 4: the outer ZZ arm reads half in the ideal circuit

**Object.** The conditional-Ramsey frequency difference archived by the July 2026
[Price Pair run 4](../../experiments/PRICE_PAIR_HARDWARE_PREDICTION.md), on the
Marrakesh line [93, 94, 95]. This is a statement about its preparation and
estimator in the no-relaxation model, not an inversion of the hardware result.
The F89d registry is structural context; the experiment and the
two Confirmations registries own the flight. A search of the formula registry,
proofs, earlier hardware records, glossary, OpenArcs and error ledger found no
tracked proof of this circuit's outer-arm gain. A local ignored reanalysis
also tests the geometry against synthetic counts; this page retains the
short algebra and the archived-data boundary. The
[ledger](../CAUGHT_ERRORS.md) records how the reading was caught.

## Circuit and estimator input

The external pipeline's `run_price_pair.py`, inspected at
`AIEvolution.UI/experiments/ibm_quantum_tomography/run_price_pair.py`
(`build_zz_circuits`, `analyze_zz`; 2026-09-23 workspace snapshot, SHA-256
`71c68724343c1434708dbb5b158110c55a33709ab21e6a6b5b0ea87eebce46a3`), has these
load-bearing lines:

```python
preps = {'n0': [], 'n1': [1], 'f1': [2]}
for q in (0, 2) if prep != 'f1' else (0,):
    qc.h(q)
for q in excited:
    qc.x(q)
```

After the delay the runner reads X/Y Ramsey quadratures on (0,2), or on
q0 alone for `f1`. Its estimator uses these selected lines:

```python
f_n0_q0, e1 = freq_of('n0', 0)
f_n1_q0, e2 = freq_of('n1', 0)
f_n0_q2, e3 = freq_of('n0', 2)
f_n1_q2, e4 = freq_of('n1', 2)
f_f1_q0, e5 = freq_of('f1', 0)
out = {
    'zeta01_khz': (f_n1_q0 - f_n0_q0) * 1e3, 'zeta01_err': np.hypot(e1, e2) * 1e3,
    'zeta12_khz': (f_n1_q2 - f_n0_q2) * 1e3, 'zeta12_err': np.hypot(e3, e4) * 1e3,
    'zeta02_khz': (f_f1_q0 - f_n0_q0) * 1e3, 'zeta02_err': np.hypot(e1, e5) * 1e3,
}
```

Thus `n0` puts q0 and q2 in |+⟩, with q1 in |0⟩. The `f1` arm puts q0 in
|+⟩ and q2 in |1⟩, with q1 still in |0⟩. The run-4
[JSON](../../data/ibm_price_pair_july2026/price_pair_zztest_ibm_marrakesh_20260704_083938.json)
archives the six fitted values and their job ID, but no counts or X/Y fringe
points. Consequently the hardware fit, its uncertainty and its arm envelopes
cannot be recomputed from that file.

## Ideal-arm calculation

The ζ symbols in this section are **operational frequency shifts**; their
sign relative to a bare `Z_i Z_j` Hamiltonian coefficient depends on the
phase and quadrature convention.

Let ζ₀₂ be the difference between q0's phase frequencies with q2 fixed in
|1⟩ and |0⟩, and let f₀ be the latter frequency. With a stationary |+⟩
neighbour, the q0 complex Ramsey fringe in the `n0` arm is proportional to

\[
C_{n0}(t) = \tfrac12(e^{i2\pi f_0t}+e^{i2\pi(f_0+\zeta_{02})t})
          = e^{i2\pi(f_0+\zeta_{02}/2)t}\cos(\pi\zeta_{02}t).
\]

The `f1` arm has frequency f₀+ζ₀₂. On a delay range without a cosine zero
or a phase-unwrap branch change, the free-intercept phase slopes therefore
obey `f_f1_q0 - f_n0_q0 = ζ₀₂/2`. The H on q2 in `n0` is essential: preparing
q2 instead in |0⟩ gives a reference slope f₀ and full gain ζ₀₂. For a chain
bond, both `n0` and `n1` carry the same superposed outer end. Exciting q1
multiplies the corresponding target fringe by `exp(i2πζ₀₁t)` for q0 and by
`exp(i2πζ₁₂t)` for q2; the common midpoint factor cancels. Their ideal
frequency differences have unit gain.

The archived outer-arm fit is `+0.016348127877902074 ±
0.11981306183811406 kHz`. The pre-registration used the runner's fitted
`ζ02_khz` for its `0.5 kHz` screen. That stored value is the outer-arm
frequency difference above: its point estimate, and its nominal two-sided
normal 95 % upper endpoint `0.25118172908060565 kHz`, are below the screen.
The pre-registration did not specify a confidence decision rule for that
screen. Its stronger physical reading, a bound on the fixed-neighbour
conditional shift or bare `Z₀Z₂` coupling, is a separate inference.

Twice the fitted numbers gives the **ideal-geometry equivalent**
`+0.03269625575580415 ± 0.23962612367622813 kHz`. Under the same nominal
95 % reading, its upper end is `0.5023634581612113 kHz`, just beyond
`0.5 kHz`. The observed fit is zero-consistent and much smaller than the
chain-bond readings; it does not establish a strict physical bound.

## Hardware boundary

Relaxation of the excited neighbour changes the branch weights during the
delay, so the phase slope need not have the ideal gain one half. Without the run-4
counts, a run-specific gain, covariance and physical ζ₀₂ estimate cannot be
extracted from the archived result. This page derives the circuit's ideal
geometry and the raw-versus-derived labels only.
