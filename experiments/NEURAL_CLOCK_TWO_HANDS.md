# The Neural Clock: a trace invariant and its negative controls

**Current scope:** The spectral mean is fixed by the diagonal leak rates;
individual decay rates and frequencies can both move. The tables below
record the clock probe's numerical outputs, not a palindrome, biological
or V-effect mechanism. Use the [canonical neural account](../docs/neural/README.md)
and [mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
for current gates. The historical producer still prints interpretations
that these controls do not support.

**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Script:** [neural_clock_two_hands.py](../simulations/neural/neural_clock_two_hands.py)

## The question and the matrix

The probe asks what coupling and external input change in a neural
Jacobian spectrum. It uses zero self-coupling and τ_E=5, τ_I=10 in model
time units. For the linear sweeps,

```
J_ii = −1/τ_i,  J_ij = α W_ij/τ_i  (i ≠ j),
mean Re λ = trace(J)/n = −(n_E/τ_E+n_I/τ_I)/n.
```

For equal E/I populations this is −s=−0.15, where
s=(1/τ_E+1/τ_I)/2. The invariant is the mean real part, the negative
of the mean decay rate for a stable spectrum. It is not a statement
that every real part is unchanged.

The graph comparison uses N=10, density parameter 0.3, balanced random
seed 3, α=0.5, and `RandomState(7)` for its rewired and random controls:

```
        graph              mean Re(λ)     −s target
  balanced random           −0.150000     −0.1500
  degree-preserved          −0.150000     −0.1500
  Erdos-Renyi (Dale)         −0.150000     −0.1500
```

The trace identity explains these rows for any zero-diagonal W.
It does not establish F36: that requires one involutive Q with
d_i+d_Q(i)+2s=0 and W_eff[Q(i),Q(j)]+W_eff[i,j]=0 entrywise.
See [F36/F37](../docs/ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra).

## Coupling and drive readings

The coupling sweep uses two N=20 balanced paired-weight constructions,
seeds 42 and 99, density parameter 0.3, α=0.5, and one excitatory mediator,
so n=41. Each constituent is normalized by its own max|W|. For each
constituent, both directions between the mediator and local seats 0 and
19 receive the displayed coupling. The mean real part is −6.2/41,
which rounds to −0.1512. The combined E/I swap fixes the mediator:
its leak −0.2 fails the required −s=−0.15 even at zero coupling.
The positive bridge also need not preserve Dale source signs.

The clock counts eigenvalues, with multiplicity and both conjugate
partners, satisfying |Im λ|>10⁻⁶. It computes
θ=atan2(|Im λ|,|Re λ|) in degrees and reports the maximum over the full
spectrum if any mode passes the count threshold; otherwise it reports zero.
Using |Re λ| hides its sign, so this angle does not test stability.

```
COUPLING SWEEP  (two N=20 constructions + 1 mediator)
  coupling   n_rotating   θ_max     mean Re(λ)
    0.00          0        0.0°      −0.1512
    0.01         12        1.0°      −0.1512
    0.05         14        1.9°      −0.1512
    0.10         12        2.6°      −0.1512
    0.30         10        4.9°      −0.1512
    1.00         10        3.4°      −0.1512
```

The drive sweep is a separate random Dale network: N=50, 25E+25I,
seed 42, density parameter 0.3, max|W| normalization and α=0.3.
It uses S_i(z)=1/(1+exp(−a_i(z−θ_i))), with
(a_E,θ_E)=(1.3,4.0), (a_I,θ_I)=(2.0,3.7).
At each P it starts x_i=0.3, performs 500 synchronous sigmoid updates,
and builds J from the resulting row slopes. The historical script does
not test or report the fixed-point equation residual. Its endpoint
spectra alone therefore certify no equilibrium stability or Hopf crossing.
P is external input, without a temperature or metabolic calibration.

This table uses |Im λ|>10⁻⁵, a different count threshold from the coupling
table. A displayed 0.0° can be rounding, as the nonzero counts show.

```
EXTERNAL-DRIVE SWEEP
    P     n_rotating   θ_max     mean Re(λ)
   0.0        24        0.1°      −0.1500
   1.0        40        0.2°      −0.1500
   2.0        40        0.6°      −0.1500
   4.0        40        2.6°      −0.1500
   6.0        38        0.6°      −0.1500
  10.0         4        0.0°      −0.1500
```

These are eigenvalue counts and angles on the stated grids, not distinct
correlation-frequency counts or measurements of sustained neural rhythms.
The [current V-effect report](../docs/neural/V_EFFECT_NEURAL.md) specifies
a converged drive solver and a shared frequency-bin instrument. Its
48/62/124 bins belong to named protocols; the drive peak's 124 becomes
286 and 403 under refinement. Neither that census nor these small angles
establishes a symmetry-release mechanism, a 2× law, or biological behavior.

## The fitted-residual negative control

The connectome comparison samples 200 balanced 5E+5I blocks using
`RandomState(trial+100)`, τ_E/τ_I=5/10 and α=0.3. Its source-row chemical
matrix is transposed into source-column Dale weights before building J.
The following numerical outputs are preserved as the record of this
instrument, not evidence of a wiring advantage:

```
                      diag (fitted)     off-diag
  C. elegans            0.0e+00          0.0128
  degree-preserved      identical        0.0129      ratio 1.00
  Erdos-Renyi (Dale)    identical        0.1119      ratio 0.11
```

Here the ratios are connectome/control mean residuals, rounded by the
producer. Its `residual_split` forms

```
S_diag[i] = −((QJQ)_ii + J_ii)/2,
R = QJQ + J + 2 diag(S_diag),
r_diag = ||diag(R)||₂/||J||_F,
r_off = ||R−diag(diag(R))||_F/||J||_F.
```

Thus r_diag=0 by construction for every J and Q. It is no scalar-centre
test. The separate algebraic check asks whether every fitted S_diag[i]
equals one prescribed s; the type swap and balanced leaks ensure that
here, not the measured zero.

If the nonzero supports of W_eff and QW_effQ are disjoint, then
r_off=√2||W_eff||_F/||J||_F. Sparsity alone does not ensure disjointness.
On that set the instrument reads coupling magnitude. The comparison
also divides the connectome by its full-matrix maximum before cutting
out a block, whereas the random control uses its own maximum. The two
arms therefore have different scales. Row-preserving weight rewiring
preserves the row-scaled coupling norm, making it insensitive on the
disjoint-support set; off that set the residual can change.

The [canonical matched audit](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md#3-empirical-c-elegans-null)
uses matched normalization and checks the support condition explicitly.
The full committed chemical matrix has 253 nonempty excitatory outputs
and 18 inhibitory outputs, so no sign-reversing support bijection exists
under its stored Dale labels and nonzero positive scales. No biological
network in this repository is known to satisfy F36. Selecting balanced
subcircuits does not by itself supply the missing identity.

## The quantum comparison and what remains open

The Hamiltonian commutator has zero superoperator trace:
trace(H⊗I−I⊗Hᵀ)=0. Zero-diagonal neural coupling also contributes no trace.
This shared trace fact fixes a spectral mean on each side; pairing on the
quantum side needs [F1's conjugation](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
and pairing on the neural side needs F36. The trace is not the reason
either full palindrome holds. Coupling can move individual real and
imaginary parts while preserving their mean.

The [crown-switch companion](../simulations/neural/neural_crown_switch.py)
is a separate parameter-scan probe of which mode has the largest real
part. Its character changes are not F37 partner transport, and a
leading-eigenvalue crossover is not by itself a Hopf bifurcation.
A biological frequency interpretation needs a specified circuit,
operating point, converged dynamics and time calibration. Those remain
open; neither the clock angle nor the fixed trace supplies them.

## Rerun

From the repository root:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/neural/neural_clock_two_hands.py
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q
```

The first command reproduces the historical probe and its stale printed
interpretations. It runs on import, so use it as a script. The latter two
run the current scalar-identity, complex-pairing and subspace-transport
controls. For the current coupling/drive producers and their declared
resolution, use [the neural operator's manual](../simulations/neural/README.md).
