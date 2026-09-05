# Neural systems: the conditional palindrome and its tests

Last refreshed: 2026-09-05 (change history lives in git).

This folder translates an operator identity owned by the quantum side into
conditions on a classical neural Jacobian. It holds exactly on networks built
to meet those conditions. No biological neural network in this repository is
known to satisfy F36; the committed C. elegans chemical-connectome model fails
its support condition.

The stores checked for this page are [the F-registry](../ANALYTICAL_FORMULAS.md)
(F36/F37 own the conditional theorem; F137 separates a trace from pairing),
[docs/proofs](../proofs/MIRROR_SYMMETRY_PROOF.md) (the quantum conjugation),
[the neural proof store](proofs/PROOF_PALINDROME_NEURAL.md) (the scalar conditions
and mode transport), and [experiments](../../experiments/NEURAL_GAMMA_CAVITY.md)
(the connectome support null and matcher controls). Searches of the hardware
flight records and [fw.Confirmations](../../simulations/framework/confirmations.py)
returned no neural hardware confirmation. [The glossary](../GLOSSARY.md) fences
the evidence grades; [OpenArcs](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
records F36/F37's typed carriage and remaining substrate questions;
[CAUGHT_ERRORS](../CAUGHT_ERRORS.md) records the diagonal-condition,
normalisation and matcher failures. The current
[Python primitives](../../simulations/neural/neural_palindrome.py),
[translation gate](../../simulations/neural/neural_translation_gate.py),
[tests](../../simulations/neural/tests/test_neural_palindrome.py),
[MirrorWorld owner](../../compute/MirrorWorld/NeuralPalindrome.cs) and
[neural run description](../../compute/MirrorWorld/README.md) supply constructed
checks. The [V-effect page](V_EFFECT_NEURAL.md), its
[mechanism analysis](proofs/PROOF_VEFFECT_MECHANISM.md), and the
[coupling](../../simulations/neural/veffect_exact.py) and
[drive](../../simulations/neural/veffect_and_heat.py) producers concern synthetic
frequency censuses, with the mechanism open.

## What transfers exactly

The [quantum owner, F1](../proofs/MIRROR_SYMMETRY_PROOF.md), is the conjugation
ΠLΠ⁻¹ = −L − 2Σγ I under its stated Hamiltonian and dephasing hypotheses.
The neural transfer is finite-dimensional algebra:

```
J = D + W_eff,  D = diag(d_i),  W_eff[i,i] = 0,
Q = an involutive permutation, Q² = I,  s = one scalar.

Q J Q⁻¹ = −J − 2s I
    iff d_i + d_Q(i) + 2s = 0                    for every i
    and W_eff[Q(i),Q(j)] + W_eff[i,j] = 0        for every i ≠ j.
```

This is [F36](../ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra).
It pairs the full complex eigenvalue multiset by λ ↦ −λ − 2s, preserving
algebraic multiplicity: the centre is −s and the partner shift is −2s.
[F37](../ANALYTICAL_FORMULAS.md#f37-neural-eigenvalue-pairing-tier-1-from-f36)
gives λ + λ′ = −(1/τ_E + 1/τ_I) when d_i = −1/τ_i and Q exchanges the two
time-constant populations, so s = (1/τ_E + 1/τ_I)/2.

For the common-gain model W_eff = α diag(1/τ_i)W, with column j the source,
the W-only condition is W[Q(i),Q(j)] = −(τ_Q(i)/τ_i)W[i,j] **only when
α ≠ 0**. At α = 0 the effective coupling vanishes whatever W is. Dale signs
supply neither the paired support nor the scaled magnitudes. Equal time
constants make the diagonal condition automatic for any involutive Q; they
do not supply the coupling condition.

The [algebra page](ALGEBRAIC_PALINDROME_NEURAL.md) gives the construction,
fixed-seat condition and invariant-subspace transport. The
[proof](proofs/PROOF_PALINDROME_NEURAL.md) derives each step. Exact palindrome
alone implies neither real eigenvalues, stability, silence, nor a biological
mechanism.

## What has been tested

The canonical [Python gate](../../simulations/neural/neural_translation_gate.py)
and [MirrorWorld tests](../../compute/MirrorWorld.Tests/NeuralPalindromeTests.cs)
include the exactly representable example
J = [[−0.5, −0.25], [0.25, −0.25]], Q = (0 1), s = 0.375.
Its entrywise residual is zero and its eigenvalues are
−0.375 ± (√3/8)i: an exact palindrome with oscillation.

The Python construction's current census uses N = 10, five E and five I seats,
τ_E = 5, τ_I = 10, density 0.3, and seeds 0…199 at each α:

| α | Scalar identity passes | Complex spectrum | Unstable |
|---|---:|---:|---:|
| 0.5 | 200/200 | 24 | 0 |
| 1.5 | 200/200 | 110 | 1 |
| 3.0 | 200/200 | 149 | 15 |
| 5.0 | 200/200 | 159 | 24 |
| 10.0 | 200/200 | 167 | 45 |

Here a scalar pass means relative Frobenius residual < 10⁻¹³; complex means
some |Im λ| > 10⁻⁸; unstable means max Re λ > 0. These are finite synthetic
ensemble counts, not population estimates or measurements of brains. The gate
also rejects an incorrect fixed-seat leak, mismatched multiplicity, discarded
imaginary parts, and a wrong Q even when the spectrum pairs.

The drive window is a separate measurement on **one synthetic random Dale
network**: N = 50, 25E + 25I, seed 42, density 0.3, τ_E = 5, τ_I = 10,
α = 0.3. In [the drive producer](../../simulations/neural/veffect_and_heat.py),
external input P changes the sigmoid operating point and hence the row gains
of the effective coupling. Its frequency census has a window in P. The name
“thermal window” refers here to this external-drive sweep; P has no calibrated
meaning as temperature, metabolism or life. The producer declares its frequency
resolution and reports the fixed-point equation residual. See
[the coupling/drive account](V_EFFECT_NEURAL.md).

## What failed

The C. elegans support gate is null for the proposed biological landing.
The committed chemical matrix has **271 nonempty output rows: 253 excitatory
and 18 inhibitory** under the stored Dale labels. Sign reversal by Q requires
a bijection between those two sets, so no qualifying E/I swap exists for the
full matrix with nonzero gain and positive rate scales. The
[G0b code](../../simulations/neural/celegans_pairing_controls.py),
[stored output](../../simulations/results/celegans_pairing_controls.txt), and
[event record](../../experiments/NEURAL_GAMMA_CAVITY.md) own this result.
It does not decide every selected subnetwork or every neural model.

A fitted off-diagonal residual can be zero while the scalar diagonal condition
fails. On disjoint Q-partner supports it reduces to a coupling-norm reading;
unequal normalisations therefore cannot establish a wiring advantage.
The [algebra page](ALGEBRAIC_PALINDROME_NEURAL.md#3-empirical-c-elegans-null)
states the precise collapse condition and matched measurements.

The proposed V-effect mechanism fails: an exact palindrome can oscillate,
and the coupled construction's fixed mediator violates the scalar condition
even at zero coupling. Coupling and drive change frequency counts, but those
censuses establish no symmetry-release mechanism and no 2× decay law. A fixed
spectral mean is a trace identity; it does not keep individual decay rates fixed.

## What remains open

The next translation must identify the quantum owner, the neural candidate,
the added hypotheses, and a gate that can reject it. The
[translation roads](ALGEBRAIC_PALINDROME_NEURAL.md#translation-roads)
retain concrete candidates from the repo's eleven-axis rulebook, including
local blindness, topology, dynamics after spectral transport, and the
two-mirror construction. The rulebook generates directions; an axis is not
evidence. Each untested landing stays open.

The biological task is to find a specified circuit, operating point and
involution satisfying both F36 conditions. The dynamics task is to determine
what the coupling/drive census tracks with consistent resolution, converged
operating points, and controls that separate changes of support, sign and
magnitude. Neither question is settled by a good pair-sum average.

## How to rerun

From the repository root, with NumPy, SciPy, pytest and .NET 10 available:

```bash
PYTHONIOENCODING=utf-8 python simulations/neural/neural_translation_gate.py
PYTHONIOENCODING=utf-8 python -m pytest simulations/neural/tests/test_neural_palindrome.py -q
dotnet test compute/MirrorWorld.Tests --filter NeuralPalindromeTests
dotnet run --project compute/MirrorWorld -- neural
```

PowerShell:

```powershell
$env:PYTHONIOENCODING = 'utf-8'
python simulations/neural/neural_translation_gate.py
python -m pytest simulations/neural/tests/test_neural_palindrome.py -q
dotnet test compute/MirrorWorld.Tests --filter NeuralPalindromeTests
dotnet run --project compute/MirrorWorld -- neural
```

For the synthetic censuses, use the same UTF-8 setting and run
`python simulations/neural/veffect_exact.py` and
`python simulations/neural/veffect_and_heat.py`; both print to stdout.
`python simulations/neural/celegans_pairing_controls.py` reruns the larger
connectome audit and rewrites its linked result file.

Related readings: [the neural clock](../../experiments/NEURAL_CLOCK_TWO_HANDS.md)
for the trace invariant; [the connectome event](../../experiments/NEURAL_GAMMA_CAVITY.md)
for its separate cavity and zero-multiplicity questions; and
[the initiating hypothesis](../../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)
for the cross-domain question.
