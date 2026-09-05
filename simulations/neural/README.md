# Neural simulations: conditional algebra and synthetic censuses

The current entry point is
[neural_translation_gate.py](neural_translation_gate.py), backed by
[neural_palindrome.py](neural_palindrome.py). For the theorem, evidence
grades and biological scope, read [the canonical neural account](../../docs/neural/README.md).
Constructed synthetic matrices satisfy F36; no biological network in
this repository is known to satisfy it. The full committed C. elegans
chemical matrix fails its necessary support condition.

## What the gate tests

For J=D+W_eff, with D diagonal and W_eff zero-diagonal, one involutive
permutation Q and one scalar s must satisfy

```
d_i+d_Q(i)+2s = 0,
W_eff[Q(i),Q(j)]+W_eff[i,j] = 0  (i ≠ j).
```

Then QJQ+J+2sI=0 and the full complex eigenvalue multiset pairs under
λ↦−λ−2s, preserving multiplicity. With paired leaks −1/τ_E and −1/τ_I,
s=(1/τ_E+1/τ_I)/2. Dale signs alone provide neither support symmetry nor
the scaled magnitudes. A fitted per-seat centre is not this scalar test.
Exact pairing implies neither a real spectrum, silence, stability nor
biological behavior. See [F36/F37's proof](../../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md).

The gate includes the exact nonreal two-seat example, fixed-seat
rejection, full complex assignment and multiplicity controls, a
constructed ensemble, and Q transport into partner invariant subspaces
with a wrong-Q negative control. The historical 96% eigenvector
character-match reading is not a current result.

## Current producers and controls

| File | Purpose and scope |
|---|---|
| [neural_palindrome.py](neural_palindrome.py) | Shared constructed-network builders, scalar residual, full complex assignment, subspace transport, frequency-independent ensemble census and fixed-point solver. |
| [neural_translation_gate.py](neural_translation_gate.py) | Canonical executable gate for the conditional translation and its negative controls. |
| [tests/test_neural_palindrome.py](tests/test_neural_palindrome.py) | Focused algebra, defective/degenerate transport, input-validation, producer and convergence tests. |
| [veffect_exact.py](veffect_exact.py) | Linear coupling and sigmoid-drive censuses on constructed constituents. The added fixed excitatory mediator fails F36 even at zero coupling; the bridge need not preserve Dale signs. |
| [veffect_and_heat.py](veffect_and_heat.py) | Synthetic random Dale-network coupling/drive censuses. P is external input. The sigmoid builder returns (J, x, fresh equation residual) and raises if its fixed-point solve does not converge. |
| [find_quarter.py](find_quarter.py) | Candidate-quarter probes. Its bounded-iteration endpoints include large equation residuals and support no equilibrium-stability or Hopf verdict. |
| [cpsi_two_perspectives.py](cpsi_two_perspectives.py) | Normalized E/I amplitude readings. A product of 1/4 at equal fractions is arithmetic, not quantum CΨ or a neural threshold; its greedy real-part matcher is not the canonical complex/transport gate. |
| [celegans_pairing_controls.py](celegans_pairing_controls.py) | Full-connectome support gate G0b, matcher/tolerance and normalization controls, dynamics probes and exact-rank checks. Its run rewrites [the result file](../results/celegans_pairing_controls.txt). |
| [celegans_connectome.json](celegans_connectome.json), [neuron IDs](celegans_neuron_ids.txt) | Stored chemical/electrical matrices, Dale labels and names. Chemical sources are rows; respect transposition and rate-scaling conventions. |

The two V-effect producers count |Im λ| activity bins and
|Im(λ_i+λ_j)| correlation bins at a declared absolute resolution.
The latter is a formal eigenvalue-sum census, not a measured correlation
signal. Zero means no resolved frequency under that predicate.
The [protocol report](../../docs/neural/V_EFFECT_NEURAL.md) gives every
seed, parameter grid, solver setting and backend. Its 48/62/124 counts
are not invariants: the drive row at P=4 changes 124→286→403 under
frequency-grid refinement. These scripts establish no V-effect
mechanism, 2× decay law, thermal, metabolic or life interpretation.

## Historical exploration scripts

These files are useful for locating how a question was probed. Their
printed headlines and local matchers are not substitutes for the current
gate. In particular, several execute their entire experiment on import.

| Files | Reading constraint |
|---|---|
| `celegans_palindrome.py`, `algebraic_palindrome.py` | Early eigenvalue matching and fitted off-diagonal residuals; neither alone tests the full conditional identity. |
| `celegans_balanced.py`, `celegans_inhibitory_position.py`, `balance_vs_size.py` | Population and subnetwork probes. Balance alone does not establish F36. |
| `random_network_controls.py`, `dense_balanced_test.py`, `validation_checks.py` | Early random and sensitivity comparisons; check matching, normalization and support against the current audit. |
| `wilson_cowan_palindrome.py`, `classical_oscillator_palindrome.py` | Model-specific Jacobian/oscillator explorations; apply the scalar conditions to each proposed matrix. |
| `neural_heartbeat.py`, `fragile_bridge_neural.py` | Particular time-domain and bridge probes. An example with real eigenvalues supplies no silence or biological-safety theorem. |
| `cpsi_candidates.py`, `cpsi_deep_dive.py`, `cpsi_interference.py` | Candidate neural quantities; none establishes a quantum CΨ observable or universal neural quarter boundary. |
| `hopf_threshold.py`, `complexity_threshold.py`, `exact_pairing_test.py` | Threshold and matcher explorations; any bifurcation claim needs a converged equilibrium branch and crossing/nondegeneracy checks. |
| `neural_gamma_cavity.py`, `neural_gamma_cavity_unpaired.py` | Absolute-tolerance/order-sensitive pairing instruments. Read the [cavity event record](../../experiments/NEURAL_GAMMA_CAVITY.md) and current support controls before using their output. |
| `neural_clock_two_hands.py`, `neural_crown_switch.py` | Trace/angle and leading-mode probes. The [clock record](../../experiments/NEURAL_CLOCK_TWO_HANDS.md) states their scope: a fixed trace does not fix individual rates, and its fitted diagonal residual is vacuous. |
| `celegans_trichotomy.py`, `neural_flavor_rule.py` | Exploratory classifications; classification labels do not establish F36 or biological pairing. |

The canonical gate, shared primitives, `veffect_exact.py`,
`veffect_and_heat.py`, `find_quarter.py` and
`cpsi_two_perspectives.py` can be imported without running their experiment.
For direct imports of the producers, place this directory on `sys.path`;
their sibling imports use that path. Do not assume the same import safety
for other scripts.

## Run from the repository root

NumPy, SciPy and pytest are needed for the canonical gate and tests.
Some historical plots/probes additionally use matplotlib or networkx.

PowerShell:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/ -q

# Current stdout-only frequency censuses
python simulations/neural/veffect_exact.py
python simulations/neural/veffect_and_heat.py

# Connectome audit: rewrites simulations/results/celegans_pairing_controls.txt
python simulations/neural/celegans_pairing_controls.py
```

On POSIX shells use `export PYTHONIOENCODING=utf-8` before the same
`python` commands. The [typed counterpart](../../compute/MirrorWorld/NeuralPalindrome.cs)
can be checked with .NET 10:

```bash
dotnet test compute/MirrorWorld.Tests --filter NeuralPalindromeTests
dotnet run --project compute/MirrorWorld -- neural
```

Quantum hardware results belong to their quantum operators and protocols.
They do not confirm this neural translation on an animal.
