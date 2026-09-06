# Finite-Occupation Amplitude Channels: Spectral Census and Scope

<!-- Keywords: thermal amplitude damping sigma minus sigma plus finite
occupation spectral census F137 palindrome centre frequency bins -->

**Status:** finite numerical census plus one-qubit direction gate; `n_bar` is an external channel parameter
**Producers:** [v_effect_thermal.py](../simulations/v_effect_thermal.py), [self_heating_fixpoint.py](../simulations/self_heating_fixpoint.py)
**Outputs:** [spectral census](../simulations/results/v_effect_thermal.txt), [channel audit](../simulations/results/self_heating_fixpoint.txt)

## The channel that was run

On each qubit the producer uses

```text
L_down = sqrt(gamma_amp * (n_bar + 1)) * sigma_minus
L_up   = sqrt(gamma_amp * n_bar)       * sigma_plus
sigma_minus = |0><1|,  sigma_plus = |1><0|.
```

Thus spontaneous emission remains at `n_bar=0`; absorption vanishes there.
The one-qubit endpoint gate checks directly that `|0><0|` is fixed at
`n_bar=0`, that `|1><1|` loses excited population at `gamma_amp`, and that a
positive `n_bar` drives population upward from `|0><0|`. Swapping the two
operators fails the same gate with residual `0.2`, while the implemented
direction has residual `0.0` in the spectral producer and `2.78e-17` in the
independent channel audit.

The parameter `n_bar` is supplied from outside. Neither producer contains a
law by which decay raises `n_bar`, so neither computes self-heating.

## What survives from the spectral run

The finite run uses an open isotropic-Heisenberg chain with `J=1`. A
"frequency" is one distinct rounded value of `|Im lambda|` after rounding to
four decimals; `Q_max=max |Im lambda|/(-Re lambda)` is likewise a property of
the sampled eigenvalues. These are protocol readings, not resolution-free
mode counts.

Selected rows from the correctly directed channel are:

| Channels | `n_bar` | `Q_max(N=2)` | `Q_max(N=5)` | ratio | frequency bins, N=5 |
|---|---:|---:|---:|---:|---:|
| Z dephasing `gamma_z=0.1` | 0 | 20.0 | 36.2 | 1.81 | 111 |
| Z dephasing `0.1` + amplitude `0.05` | 0 | 17.8 | 32.2 | 1.81 | 111 |
| same combined channel | 0.5 | 14.5 | 21.0 | 1.44 | 403 |
| same combined channel | 2.0 | 9.4 | 12.5 | 1.33 | 423 |
| same combined channel | 5.0 | 5.5 | 7.1 | 1.29 | 445 |

The values say only that changing the external channel parameter changes this
finite eigenspectrum and this binning protocol. They do not show energy being
converted into frequencies, and they do not define temperature, metabolic
cost, life, or complexity.

For the non-uniform Z profile `[0.5, 0.01, 0.01, 0.01, 0.01]`, compared with
the equal-total uniform profile `[0.108]*5`, the sampled `Q_max` ratio falls
from `2.97` at `n_bar=0` to `1.02` at `n_bar=10`. This is a parameter sweep of
two specified generators. It neither identifies `n_bar` with a physical
temperature nor proves a general profile theorem.

## Palindrome reading

The pure amplitude channel remains spectrally palindromic at the eigensolver
floor in this run. Its centre is

```text
-Sum(gamma_down + gamma_up) / 2,
```

the F137 centre. The code also reads the only possible candidate centre from
`trace(L)/dim` and evaluates the canonical multiset distance there. Co-axial
Z-dephasing plus the amplitude channel has distances of order `1e14` to `1e15`
in units of machine epsilon times spectral radius in the N=5 rows, so the
combined channel is not F1-paired. A centre alone is not a pairing test.

See [F137](../docs/ANALYTICAL_FORMULAS.md#f137) and the exact channel-scope
discussion in [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

## What the second producer establishes

For the selected local rates, one qubit has stationary excited-state
population

```text
p_excited = n_bar / (2*n_bar + 1).
```

The tensor product of that local state is stationary for the tested uniform
Heisenberg chains because it commutes with the Hamiltonian and each local
channel fixes it. The audit measures `||L rho_target||` from `0.0` to
`5.56e-17` for N=3 and N=5, including one non-uniform Z profile. A wrong
one-site population is the negative control and gives residual `0.141`.

This local-channel target is not, in general, the Gibbs state of the
interacting Heisenberg Hamiltonian. Comparing their energies does not close a
feedback equation. The generators establish no runaway self-heating, required
external cooling, or biological metabolism analogue.

## Reproduction

```bash
python simulations/v_effect_thermal.py
python simulations/self_heating_fixpoint.py
```
