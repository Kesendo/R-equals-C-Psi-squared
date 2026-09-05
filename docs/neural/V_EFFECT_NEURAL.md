# The neural V-effect: coupling and drive frequency censuses

Last refreshed: 2026-09-05 (change history lives in git).

The V-effect names the nonmonotone frequency-bin census in the synthetic
coupling and drive sweeps below. Its mechanism is open. Exact conditional
palindrome pairs full complex eigenvalues; it implies neither realness,
stability, silence, optimality nor biology. P is an external drive parameter,
not temperature, heat or metabolic input. These experiments establish no
mechanism of life or Yerkes-Dodson law.

The named stores checked are [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)
(F36/F37 supply conditional complex pairing; F137 separates trace from
pairing), [docs/proofs](../proofs/MIRROR_SYMMETRY_PROOF.md) (quantum
conjugation), [the neural proof store](proofs/PROOF_PALINDROME_NEURAL.md)
(scalar entry conditions), and [experiments](../../experiments/NEURAL_GAMMA_CAVITY.md)
(connectome support null and matcher controls). The
[neural clock record](../../experiments/NEURAL_CLOCK_TWO_HANDS.md) supplies
the trace identity, not invariance of individual rates. The hardware-flight
search, including [the IBM synthesis](../../experiments/IBM_HARDWARE_SYNTHESIS.md),
and [fw.Confirmations](../../simulations/framework/confirmations.py) supply
no neural hardware confirmation. [GLOSSARY](../GLOSSARY.md) distinguishes
evidence grades; [OpenArcs](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
records conditional F36/F37 carriage with substrate questions open;
[CAUGHT_ERRORS](../CAUGHT_ERRORS.md) records diagonal, normalization and
matcher hazards. The current [README](README.md),
[algebra page](ALGEBRAIC_PALINDROME_NEURAL.md), producers and
[tests](../../simulations/neural/tests/test_neural_palindrome.py) distinguish
constructed algebra from synthetic censuses. The committed
[connectome controls](../../simulations/results/celegans_pairing_controls.txt)
own a separate support/normalization audit; the
[cavity output](../../simulations/results/neural_gamma_cavity.txt) explicitly
rejects its headline pairing interpretation. Neither supplies the counts
below, which come from the current stdout producers.

## Observable and counting rules

For eigenvalues λ_j of a real Jacobian J, activity frequencies are |Im λ_j|.
The correlation census forms |Im(λ_i+λ_j)| for every ordered pair, including
i=j: the spectrum of the formal Kronecker sum J⊗I+I⊗Jᵀ. It computes these
sums directly, without measuring correlations or simulating co-activation.

At absolute resolution ε, the shared
[`frequency_counts`](../../simulations/neural/veffect_and_heat.py) discards
frequencies ≤ ε, applies NumPy `round(frequency/ε)`, and counts distinct bins.
K_act and K_corr use the same ε. The separate n_osc counts eigenvalues with
|Im λ|>ε, retaining multiplicity and both conjugate partners. Zero means
no resolved frequency under this predicate, not exact silence. These are
not counts of sustained nonlinear oscillations. Frequencies have inverse
model-time units; no biological Hz calibration is supplied. There is no
Fourier/time grid: ε is the frequency-bin grid and the parameter grids are
listed below.

The scalar diagnostic is
`r_s = ||QJQ+J+2sI||_F / ||J||_F`, with `ei_swap` supplying Q and
s=(1/5+1/10)/2=0.15. It retains the diagonal and one scalar centre.
Small coupling can give a small residual without establishing F36.

## Coupling two constructed networks

[`veffect_exact.py`](../../simulations/neural/veffect_exact.py) uses
`build_exact_weights` from [the shared primitives](../../simulations/neural/neural_palindrome.py).
Each constituent has N/2 E and N/2 I seats, τ_E=5, τ_I=10, linear gain
α=0.5. A uses `RandomState(42)`, B `RandomState(99)`. The constructor uses
a Bernoulli mask parameter 0.3, exponential base scale 0.3, paired-weight
assignments and division by max|W|. The input mask density is not the final
support density after partner overwrites. Both balanced constituents satisfy
F36 by construction. At ε=10⁻⁶ both return K_act=K_corr=0 for N=10 and 20;
the seed-42 N=30 single-network check also returns zero.

The coupled matrix has N_c=2N+1 seats: A, B and excitatory mediator M=2N.
For each offset o∈{0,N}, both directions between M and seats o and o+N−1
receive +c. The Jacobian has J_ii=−1/τ_i and J_ij=αW_ij/τ_i; no fixed-point
solver is involved. The mediator is fixed by Q and has diagonal −0.2
rather than −s=−0.15. Thus this matrix fails F36 even at c=0. At c>0 the
positive bridge also violates Dale source signs on one inhibitory column
at N=10 and two at N=20. The exact constituents do not make this bridge
an exact-palindromic or Dale-preserving perturbation.

| c | N=10: K_act | K_corr | r_s | N=20: K_act | K_corr | r_s |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0 | 0 | 0.132127 | 0 | 0 | 0.095731 |
| 0.01 | 2 | 6 | 0.132222 | 6 | 48 | 0.095793 |
| 0.05 | 3 | 12 | 0.134480 | 7 | 62 | 0.097267 |
| 0.10 | 3 | 12 | 0.141293 | 6 | 47 | 0.101734 |
| 0.30 | 1 | 2 | 0.199466 | 5 | 34 | 0.140645 |
| 0.50 | 1 | 2 | 0.279759 | 5 | 33 | 0.195902 |
| 1.00 | 0 | 0 | 0.490095 | 5 | 31 | 0.347980 |

“0+0=48” means zero constituent bin counts and 48 correlation bins for
N_c=41, c=0.01, ε=10⁻⁶. It does not count modes released by a symmetry.
N=10 returns to zero at c=1; N=20 still returns 31.

The N=10 table survives ε/4, ε/16 and transposition. N=20 is numerically
sensitive: at c=0.05 its (K_act,K_corr) reads (7,62), (8,72), (8,73) at
ε=10⁻⁶, 2.5×10⁻⁷, 6.25×10⁻⁸. At the original ε, Jᵀ returns (7,64).
Even c=0 reads (0,0) for J and (2,7) for Jᵀ. Transposition preserves the
exact spectrum: these changes expose numerical sensitivity, not different
physical modes. Neither 48 nor 62 is an established resolution-independent
count.

## External drive on one random Dale network

[`veffect_and_heat.py`](../../simulations/neural/veffect_and_heat.py) uses one
synthetic network: N=50, 25E+25I, `RandomState(42)`, mask density 0.3,
exponential weight scale 0.3, zero self-coupling, source-column Dale signs
and division by max|W|. No F36 magnitude relation is imposed. With
τ_E=5, τ_I=10 and α=0.3, its model is

```
τ_i dx_i/dt = −x_i + S_i(α(Wx)_i + P),
S_i(z) = 1 / (1 + exp(−a_i(z−θ_i))),
(a_E,θ_E) = (1.3,4.0),  (a_I,θ_I) = (2.0,3.7).
```

Each P starts independently at x_i=0.3. The synchronous solver allows
5000 iterations, requiring max|x_next−x|<10⁻¹² and the fresh equation
residual max|x−S(αWx+P)|≤10⁻¹². Exhaustion raises an error. At the accepted
point, g_i=a_i S_i(1−S_i) and
J=diag(1/τ_i)[−I+α diag(g_i)W]. Drive changes row gains; this dependency
alone does not explain the particular frequency census.

The full P grid below uses ε=10⁻⁴. Every endpoint converges; the largest
fresh equation residual is 5.17×10⁻¹⁴. Tightening the solver tolerance to
10⁻¹⁴ leaves all listed K_act/K_corr unchanged.

| P | n_osc | K_act | K_corr |
|---:|---:|---:|---:|
| 0.0 | 10 | 2 | 4 |
| 0.5 | 14 | 3 | 7 |
| 1.0 | 18 | 6 | 13 |
| 1.5 | 32 | 8 | 22 |
| 2.0 | 38 | 13 | 39 |
| 2.5 | 40 | 17 | 60 |
| 3.0 | 40 | 19 | 90 |
| 3.5 | 40 | 15 | 107 |
| 4.0 | 40 | 19 | 124 |
| 5.0 | 38 | 15 | 81 |
| 6.0 | 30 | 9 | 34 |
| 8.0 | 10 | 2 | 3 |
| 10.0 | 0 | 0 | 0 |

124 is the largest K_corr on this declared P grid for this N=50, seed-42
network, solver and ε. It is not a crossing value, universal constant,
continuous-P optimum or biological frequency count. At P=0 the count is
four, not zero. Resolution controls give:

| P | K_corr, ε=10⁻⁴ | ε=2.5×10⁻⁵ | ε=6.25×10⁻⁶ |
|---:|---:|---:|---:|
| 0 | 4 | 13 | 48 |
| 4 | 124 | 286 | 403 |
| 10 | 0 | 0 | 4 |

The sampled rise and fall survives; the counts and endpoint silence do
not. This drive network also fails F36: r_s at P=4 is 0.102517, and at
P=10 it remains 0.000159 rather than zero.

## Companion measurements

The random-network producer also couples seed-42 and seed-99 constituents
at P=1.5, α=c=0.3, with the same bridge, density, τ_E/τ_I=5/10 and solver.
At ε=10⁻⁴ it reports:

| N per constituent | A: K_act | A: K_corr | Coupled: K_act | Coupled: K_corr |
|---:|---:|---:|---:|---:|
| 10 | 1 | 2 | 2 | 6 |
| 20 | 3 | 12 | 5 | 16 |
| 50 | 8 | 22 | 11 | 23 |
| 100 | 15 | 32 | 16 | 33 |
| 200 | 17 | 34 | 18 | 37 |

Its printed ratios divide by two times A's count, although B uses another
seed; they are not ratios to a measured sum of A and B.

The constructed N=20, seed-42 network in `veffect_exact.py` separately
receives sigmoid drive at α=0.5. At P=0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5,
4, 5, 6 and 8 all endpoints converge and n_osc=K_act=K_corr=0 at ε=10⁻⁶.
Its linear common-gain matrix satisfies F36; its sigmoid row gains do not:
r_s=0.00138 at P=0 and 0.0443 at P=3.5. Failure of F36 is insufficient
to produce a positive census.

The decay table in `veffect_and_heat.py` classifies by |Im λ|>10⁻⁶, not by
palindrome pairing. For seed-42 random constituents at P=1.5, α=0.3 and
the same density/time constants, mean(−Re λ) in the real-at-resolution
group divided by the oscillatory group is 1.626, 0.875, 1.000 and 0.836
for N=10, 20, 50 and 100. This establishes no 2× decay law.

The exact 200-seed F36 ensemble is separate: N=10, a linear gain scan,
no mediator and no sigmoid drive. Its counterexamples and the next gates
belong to the [mechanism constraint page](proofs/PROOF_VEFFECT_MECHANISM.md).
No replacement mechanism is established.

## Reproduce the census

From the repository root with NumPy, SciPy and pytest. These commands print
to stdout and do not rewrite tracked results.

The recorded tables use Python 3.12.9, NumPy 2.4.2 and SciPy 1.16.3 on
64-bit Windows. NumPy reports scipy-openblas/OpenBLAS 0.3.31.dev for BLAS
and LAPACK, with `USE64BITINT` (ILP64). This numerical backend is part of
the recorded protocol because the N=20 bins show eigensolver/roundoff
sensitivity; these versions are not universal dependency requirements.

POSIX:

```bash
export PYTHONIOENCODING=utf-8
python simulations/neural/veffect_exact.py
python simulations/neural/veffect_and_heat.py
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
```

PowerShell:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/neural/veffect_exact.py
python simulations/neural/veffect_and_heat.py
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
```

For the drive resolution table, run the Python block below using POSIX
`python - <<'PY'` with a closing `PY`, or a PowerShell single-quoted
here-string (`@'` and `'@`) piped to `python -`:

```python
import sys
sys.path.insert(0, 'simulations/neural')
import numpy as np
from veffect_and_heat import (
    make_balanced_dale_network, build_jacobian_with_sigmoid, frequency_counts,
)
W, signs = make_balanced_dale_network(50, 25, 0.3, 42)
for P in [0, 4, 10]:
    J, x, residual = build_jacobian_with_sigmoid(W, signs, 5, 10, .3, P)
    values = np.linalg.eigvals(J)
    print(P, residual,
          [frequency_counts(values, e) for e in [1e-4, 2.5e-5, 6.25e-6]])
```

For coupling controls, import `veffect_exact`, replace its module binding
`frequency_counts` with a wrapper calling the shared function at ε/4 or
ε/16, then call `main()`; for the transpose control replace its `np.linalg.eigvals`
binding with a wrapper applying a saved original function to `J.T`. Run each
in a fresh Python process. For solver refinement, call `solve_fixed_point`
with `tol=1e-14`, then rebuild J using the displayed row-gain equation.
